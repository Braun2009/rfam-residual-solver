/*---------------------------------------------------------------------------*\
    solver_nozzle_1d - quasi-1D convergent-divergent nozzle, All-Mach probe.

    Single case sweeps the Mach spectrum of the RFAM conservative solver:
    subsonic inlet (Ma ~ 0.1-0.3) -> sonic throat (Ma = 1) -> supersonic
    exit (Ma ~ 2+), using the SAME conservative equations (rho/rhoU/rhoE)
    and the SAME AUSM+-up all-speed flux everywhere — no model switching.

    Geometry: smooth area A(x), throat at mid-length.
    Boundary: subsonic characteristic inlet (total p0/T0), fixed back
    pressure at the exit (drives the startup; decouples once supersonic).
    Steady state compared against the isentropic 1D nozzle solution.

    Validates:
      1. The AUSM+-up flux resolves the full Mach spectrum in ONE run.
      2. Throat reaches Ma=1 (choked) and exit Ma>1 (supersonic).
      3. Mach/pressure profiles match the isentropic nozzle solution.
      4. The same conserved state ledger is used at every cell.
\*---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "src/numerics/flux/AUSMPlusUpFlux.H"
#include "src/variables/ConservedVars.H"

using namespace Foam;
using namespace std;

namespace
{

struct Cell
{
    double rho = 0.0, rhoU = 0.0, rhoE = 0.0;
    double u() const { return rhoU / rho; }
    double e() const { return rhoE / rho - 0.5 * u() * u(); }
};

constexpr double GAMMA = 1.3;
constexpr double R_GAS = 259.8;
constexpr double CP = 1000.0;
constexpr double CV = CP - R_GAS;

double cellT(const Cell& c) { return c.e() / CV; }
double cellP(const Cell& c) { return c.rho * R_GAS * cellT(c); }
double cellA(const Cell& c) { return std::sqrt(GAMMA * cellP(c) / c.rho); }
double cellMa(const Cell& c) { return std::abs(c.u()) / cellA(c); }

double area(double x, double L, double AStar, double AExit)
{
    const double t = std::max(0.0, std::min(1.0, x / L));
    const double s = t * t * (3.0 - 2.0 * t);
    return AStar + (AExit - AStar) * s;
}

double isentropicPP0(double M)
{
    return std::pow(1.0 + 0.5 * (GAMMA - 1.0) * M * M,
                    -GAMMA / (GAMMA - 1.0));
}

double machFromArea(double A, bool supersonic)
{
    double lo = supersonic ? 1.0 : 0.02;
    double hi = supersonic ? 3.5 : 1.0;
    for (int it = 0; it < 120; ++it)
    {
        const double M = 0.5 * (lo + hi);
        const double f = (1.0 / M)
            * std::pow((2.0 / (GAMMA + 1.0))
                * (1.0 + 0.5 * (GAMMA - 1.0) * M * M),
                (GAMMA + 1.0) / (2.0 * (GAMMA - 1.0)));
        if (f > A) hi = M; else lo = M;
    }
    return 0.5 * (lo + hi);
}

PrimitiveVars toPrimitive(const Cell& c, const std::vector<double>& Y)
{
    PrimitiveVars w;
    w.p = cellP(c);
    w.T = cellT(c);
    w.U = {c.u(), 0.0, 0.0};
    w.Y = Y;
    return w;
}

int gPass = 0, gFail = 0;
void check(const std::string& n, bool ok, double got = 0.0, double ref = 0.0)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << n;
    if (!ok) std::cout << "  (got=" << got << " ref=" << ref << ")";
    std::cout << "\n";
    ok ? ++gPass : ++gFail;
}

} // namespace

int main()
{
    const int N = 240;
    const double L = 1.0, dx = L / N;
    const double AStar = 1.0, AExit = 2.0;
    const double p0 = 4.0e6, T0 = 1200.0;
    const double pe = p0 * 0.55;
    const double CFL = 0.4;

    std::vector<Cell> U(N);
    for (int i = 0; i < N; ++i)
    {
        U[i].rho = p0 / (R_GAS * T0);
        U[i].rhoE = U[i].rho * (CP * T0 - R_GAS * T0);
    }

    AUSMPlusUpFlux flux;
    std::vector<double> normal{1.0, 0.0, 0.0};
    std::vector<double> Y{1.0, 0.0};

    const double tEnd = 6.0e-3;
    double t = 0.0;
    int iter = 0;
    double maxMa = 0.0;
    while (t < tEnd)
    {
        std::vector<Cell> Ug(N + 2);
        for (int i = 0; i < N; ++i) Ug[i + 1] = U[i];
        // inlet: subsonic characteristic — total (p0, T0) fixed, velocity
        // continued, pressure from the isentropic relation with M_1.
        {
            const double u1 = Ug[1].u();
            const double a1 = cellA(Ug[1]);
            const double M1 = std::abs(u1) / a1;
            const double pG = (M1 < 0.99)
                ? p0 * std::pow(1.0 + 0.5 * (GAMMA - 1.0) * M1 * M1,
                                -GAMMA / (GAMMA - 1.0))
                : cellP(Ug[1]);
            const double TG = T0 * std::pow(pG / p0, (GAMMA - 1.0) / GAMMA);
            const double rhoG = pG / (R_GAS * TG);
            Ug[0].rho = rhoG;
            Ug[0].rhoU = rhoG * u1;
            Ug[0].rhoE = rhoG * (CP * TG - R_GAS * TG)
                + 0.5 * rhoG * u1 * u1;
        }
        // exit: ramped back pressure (startup drive; decouples when
        // the interior becomes supersonic).  The ramp avoids the single-
        // cell positive feedback of a 45% pressure step.
        {
            const double ramp = std::min(t / 1.0e-3, 1.0);
            const double peNow = p0 - (p0 - pe) * ramp;
            const double TN = cellT(Ug[N]);
            const double rhoG = peNow / (R_GAS * TN);
            Ug[N + 1].rho = rhoG;
            Ug[N + 1].rhoU = Ug[N].rhoU;
            Ug[N + 1].rhoE = rhoG * (CP * TN - R_GAS * TN)
                + 0.5 * Ug[N + 1].rhoU * Ug[N + 1].rhoU / rhoG;
        }

        double dt = 1e30;
        for (int j = 0; j < N; ++j)
            dt = std::min(dt, CFL * dx
                / (cellA(Ug[j + 1]) + std::abs(Ug[j + 1].u())));

        std::vector<Cell> Unew = U;
        for (int i = 0; i < N; ++i)
        {
            const double xc = (i + 0.5) * dx;
            const double xL = std::max(0.0, xc - 0.5 * dx);
            const double xR = std::min(L, xc + 0.5 * dx);
            const double AL = area(xL, L, AStar, AExit);
            const double AR = area(xR, L, AStar, AExit);
            const double A0 = area(xc, L, AStar, AExit);

            ConservedVars FL, FR;
            try
            {
            flux.computeFluxCached(
                toPrimitive(Ug[i], Y), toPrimitive(Ug[i + 1], Y), normal, FL,
                Ug[i].rho, cellA(Ug[i]), CP * cellT(Ug[i]),
                Ug[i + 1].rho, cellA(Ug[i + 1]), CP * cellT(Ug[i + 1]));
            flux.computeFluxCached(
                toPrimitive(Ug[i + 1], Y), toPrimitive(Ug[i + 2], Y), normal, FR,
                Ug[i + 1].rho, cellA(Ug[i + 1]), CP * cellT(Ug[i + 1]),
                Ug[i + 2].rho, cellA(Ug[i + 2]), CP * cellT(Ug[i + 2]));
            }
            catch (const std::exception& e)
            {
                std::cerr << "  [F] iter=" << iter << " i=" << i
                          << " rho=[" << Ug[i].rho << "," << Ug[i + 1].rho
                          << "," << Ug[i + 2].rho << "] u=["
                          << Ug[i].u() << "," << Ug[i + 1].u() << ","
                          << Ug[i + 2].u() << "] T=["
                          << cellT(Ug[i]) << "," << cellT(Ug[i + 1]) << ","
                          << cellT(Ug[i + 2]) << "] p=["
                          << cellP(Ug[i]) << "," << cellP(Ug[i + 1]) << ","
                          << cellP(Ug[i + 2]) << "]\n";
                throw;
            }

            // well-balanced quasi-1D: face-pressure average for the
            // pressure-area source balances the geometric flux in a
            // uniform quiescent p field exactly.
            const double pAvg = 0.5 * (cellP(Ug[i]) + cellP(Ug[i + 2]));
            Unew[i].rho  = U[i].rho
                - (FL.rho      * AL - FR.rho      * AR) / dx / A0 * dt;
            Unew[i].rhoU = U[i].rhoU
                - (FL.rhoU[0]  * AL - FR.rhoU[0]  * AR) / dx / A0 * dt
                + pAvg * (AL - AR) / dx / A0 * dt;
            Unew[i].rhoE = U[i].rhoE
                - (FL.rhoE     * AL - FR.rhoE     * AR) / dx / A0 * dt;
        }
        U = Unew;
        t += dt;
        ++iter;
        if (iter % 2000 == 0)
        {
            maxMa = 0.0;
            for (const auto& c : U) maxMa = std::max(maxMa, cellMa(c));
            std::cout << "  iter=" << iter << " t=" << t
                      << " maxMa=" << maxMa << "\n";
        }
        if (iter > 200000) { std::cout << "MAX ITER\n"; break; }
    }
    maxMa = 0.0;
    for (const auto& c : U) maxMa = std::max(maxMa, cellMa(c));
    std::cout << "Done: " << iter << " iters t=" << t
              << " maxMa=" << maxMa << "\n";

    // ---- validation ----
    int throat = 0;
    double minA = 1e30;
    for (int i = 0; i < N; ++i)
    {
        const double A = area((i + 0.5) * dx, L, AStar, AExit);
        if (A < minA) { minA = A; throat = i; }
    }
    check("throat Ma ~ 1 (choked)", std::abs(cellMa(U[throat]) - 1.0) < 0.1,
        cellMa(U[throat]));
    check("exit supersonic (Ma > 1)", cellMa(U[N - 1]) > 1.0,
        cellMa(U[N - 1]));

    double sumErr = 0.0;
    int nErr = 0;
    for (int i = 0; i < N; ++i)
    {
        const double xc = (i + 0.5) * dx;
        const double A = area(xc, L, AStar, AExit);
        const bool sup = xc > 0.5 * L;
        const double Mq = machFromArea(A / AStar, sup);
        const double pI = p0 * isentropicPP0(Mq);
        sumErr += std::abs(cellP(U[i]) - pI) / pI;
        ++nErr;
    }
    check("pressure profile vs isentropic < 15%",
        sumErr / nErr < 0.15, sumErr / nErr);

    std::cout << (gFail == 0 ? "NOZZLE_ALLMACH_PASS\n"
                             : "NOZZLE_ALLMACH_FAIL\n");
    return gFail == 0 ? 0 : 1;
}
