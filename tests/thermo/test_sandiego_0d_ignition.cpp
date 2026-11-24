/* 0D constant-volume autoignition of C3H8/O2 with the San Diego 59-species
   mechanism, integrated by CVODE, cross-checked against Cantera 3.1
   IdealGasReactor ignition delays (generated offline in WSL):

     T0=1000 K, p0=4 MPa, phi=1.0 (C3H8/O2): IDT(T0+400 K) = 4.515e-4 s
     T0=1200 K, p0=4 MPa, phi=1.0 (C3H8/O2): IDT(T0+400 K) = 3.50e-5 s

   The reactor is adiabatic constant-volume: dY_i/dt = omega_i/rho and
   dT/dt = -sum_i omega_i e_sens_i / (rho cv).  CVODE BDF/dense handles the
   ignition stiffness. */

#include "SandiegoChemistry.H"
#include "SandiegoNasa7Data.H"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_types.h>

using namespace Foam;

namespace
{
constexpr double Runiv = 8.314462618;   // J/(mol K)
constexpr double Tref = 298.15;
int failures = 0;

void check(bool condition, const std::string& label)
{
    std::printf("  [%s] %s\n", condition ? "PASS" : "FAIL", label.c_str());
    if (!condition) ++failures;
}

struct ReactorContext
{
    const SandiegoChemistry* chem;
    double rho;                  // kg/m3, constant
    std::vector<double> omega;   // scratch
    double qdot;
};

int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* userData)
{
    (void)t;
    auto* ctx = static_cast<ReactorContext*>(userData);
    double Y[59];
    double T = NV_Ith_S(y, 59);
    for (int i = 0; i < 59; ++i) Y[i] = NV_Ith_S(y, i);
    std::vector<double> yVec(Y, Y + 59);
    ctx->chem->computeSourceTerms(T, 0.0, ctx->rho, yVec, ctx->omega, ctx->qdot);

    double energyRate = 0.0;
    double cv = 0.0;
    for (int i = 0; i < 59; ++i)
    {
        const double mw = SandiegoNasa7::species[i].mw;  // kg/mol
        const double eSens = (SandiegoNasa7::hMole(i, T)
                - SandiegoNasa7::hMole(i, Tref)) / mw
            - Runiv * T / mw;
        energyRate += ctx->omega[i] * eSens;
        cv += Y[i] * (SandiegoNasa7::cpMole(i, T) - Runiv) / mw;
        NV_Ith_S(ydot, i) = ctx->omega[i] / ctx->rho;
    }
    // dT/dt = (-sum omega_i e_i + qdot) / (rho cv): sensible-energy
    // convention, explicit chemical heat release.
    NV_Ith_S(ydot, 59) = (-energyRate + ctx->qdot) / (ctx->rho * cv);
    return 0;
}

double ignitionDelay(double T0, double rho, const std::vector<double>& Y0)
{
    SUNContext sunctx = nullptr;
    SUNContext_Create(SUN_COMM_NULL, &sunctx);
    N_Vector y = N_VNew_Serial(60, sunctx);
    for (int i = 0; i < 59; ++i) NV_Ith_S(y, i) = Y0[i];
    NV_Ith_S(y, 59) = T0;

    SandiegoChemistry chem;
    ReactorContext ctx{&chem, rho, std::vector<double>(59, 0.0), 0.0};

    void* mem = CVodeCreate(CV_BDF, sunctx);
    CVodeInit(mem, rhs, 0.0, y);
    CVodeSetUserData(mem, &ctx);
    // Component tolerances: mass fractions may drop to machine zero (minor
    // species) while the temperature must be tracked to ~0.1 K for the
    // ignition-threshold detection.
    N_Vector atol = N_VNew_Serial(60, sunctx);
    for (int i = 0; i < 59; ++i) NV_Ith_S(atol, i) = 1.0e-20;
    NV_Ith_S(atol, 59) = 0.5;
    CVodeSVtolerances(mem, 1.0e-9, atol);
    N_VDestroy(atol);
    SUNMatrix matrix = SUNDenseMatrix(60, 60, sunctx);
    SUNLinearSolver solver = SUNLinSol_Dense(y, matrix, sunctx);
    CVodeSetLinearSolver(mem, solver, matrix);
    CVodeSetMaxNumSteps(mem, 1000000);

    const double targetT = T0 + 400.0;
    double idt = 0.0;
    double t = 0.0, tStop = 0.5;
    double tPrev = 0.0, tPrevT = T0;
    bool refined = false;
    while (t < tStop)
    {
        // Coarse exponential stepping through the induction period; refine to
        // fixed 1e-7 s steps once the temperature approaches the threshold.
        const double dtCoarse =
            std::min(2.0 * std::max(t, 1.0e-7), 2.0e-6);
        const double dt = (NV_Ith_S(y, 59) > T0 + 100.0 && !refined)
            ? 1.0e-7 : dtCoarse;
        if (NV_Ith_S(y, 59) > T0 + 100.0) refined = true;
        double tNext = std::min(t + dt, tStop);
        double tActual = 0.0;
        int flag = CVode(mem, tNext, y, &tActual, CV_NORMAL);
        if (flag < 0)
        {
            std::printf("    CVODE fail flag=%d at t=%g T=%g\n",
                flag, tActual, NV_Ith_S(y, 59));
            break;
        }
        t = tActual;
        if (NV_Ith_S(y, 59) >= targetT)
        {
            // Linear interpolation between the bracketing samples.
            const double slope = (NV_Ith_S(y, 59) - tPrevT)
                / std::max(t - tPrev, 1.0e-300);
            idt = tPrev + (targetT - tPrevT) / slope;
            break;
        }
        tPrev = t;
        tPrevT = NV_Ith_S(y, 59);
    }
    if (NV_Ith_S(y, 59) < targetT) idt = t;
    N_VDestroy(y);
    SUNLinSolFree(solver);
    SUNMatDestroy(matrix);
    CVodeFree(&mem);
    SUNContext_Free(&sunctx);
    return idt;
}

double initialDensity(double T0, const std::vector<double>& X0)
{
    double meanMw = 0.0;
    for (int i = 0; i < 59; ++i)
        meanMw += X0[i] * SandiegoNasa7::species[i].mw;   // kg/mol
    return 4.0e6 / (Runiv / meanMw * T0);                 // kg/m3
}
} // namespace

int main()
{
    std::printf("=== San Diego 0D autoignition vs Cantera ===\n");

    // phi=1.0 C3H8/O2; mole fractions from stoichiometry: C3H8 + 5 O2.
    std::vector<double> X(59, 0.0);
    const int iC3 = SandiegoNasa7::index("C3H8");
    const int iO2 = SandiegoNasa7::index("O2");
    X[iC3] = 1.0 / 6.0;
    X[iO2] = 5.0 / 6.0;
    std::vector<double> Y(59, 0.0);
    double sum = 0.0;
    for (int i = 0; i < 59; ++i)
        Y[i] = X[i] * SandiegoNasa7::species[i].mw;
    sum = 0.0;
    for (double y : Y) sum += y;
    for (double& y : Y) y /= sum;

    // T0 = 1200 K
    {
        const double T0 = 1200.0;
        const double rho = initialDensity(T0, X);
        const double idt = ignitionDelay(T0, rho, Y);
        std::printf("  T0=1200K rho=%.6g IDT=%.6g s (Cantera 3.5e-5 s, ratio=%.4f)\n",
            rho, idt, idt / 3.5e-5);
        check(std::abs(idt / 3.5e-5 - 1.0) < 0.15,
            "1200 K ignition delay matches Cantera within 15%");
    }
    // T0 = 1000 K
    {
        const double T0 = 1000.0;
        const double rho = initialDensity(T0, X);
        const double idt = ignitionDelay(T0, rho, Y);
        std::printf("  T0=1000K rho=%.6g IDT=%.6g s (Cantera 4.515e-4 s, ratio=%.4f)\n",
            rho, idt, idt / 4.515e-4);
        check(std::abs(idt / 4.515e-4 - 1.0) < 0.15,
            "1000 K ignition delay matches Cantera within 15%");
    }

    std::printf("%s\n", failures == 0
        ? "SANDIEGO_0D_IGNITION_PASS"
        : "SANDIEGO_0D_IGNITION_FAIL");
    return failures == 0 ? 0 : 1;
}
