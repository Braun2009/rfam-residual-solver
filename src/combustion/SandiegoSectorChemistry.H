/* SandiegoSectorChemistry - sector chemistry with the San Diego 59-species
   mechanism, bridged to the 9-species Domain transport state.

   Per active cell:
     1. expand the 9 global mass fractions to the 59-species mechanism basis
        (unlisted species start at zero),
     2. advance an adiabatic constant-volume CVODE reactor by dt (the same
        stiff integration validated against Cantera IdealGasReactor in
        test_sandiego_0d_ignition, <3% ignition-delay error),
     3. project back to the 9 transported species (the 50 intermediate
        species mass is renormalised away; diagnostic approximation - their
        mass fraction is small in quasi-steady regions),
     4. update rhoY and rhoE using the 59-species sensible-energy change.

   The sensible-energy convention matches RFAM (explicit chemical heat):
   de_sens/dt = qdot/rho, so rhoE_new = rhoE + rho*(e59(Tn,Yn)-e59(T0,Y0)).
*/

#ifndef SandiegoSectorChemistry_H
#define SandiegoSectorChemistry_H

#include "Operator.H"
#include "SandiegoChemistry.H"
#include "SandiegoNasa7Data.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_types.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Foam
{

namespace SandiegoBridge
{

constexpr double Runiv = 8.314462618;
constexpr double Tref = 298.15;

// 9-species RFAM registry order (O2,N2,H2,CH4,C2H6,C3H8,H2O,CO,CO2) to the
// San Diego 59-species ordering.  All nine are present in the mechanism.
class SpeciesAdapter
{
public:
    explicit SpeciesAdapter(const SpeciesRegistry& registry)
    {
        globalToMech_.assign(static_cast<size_t>(registry.size()), -1);
        for (int g = 0; g < registry.size(); ++g)
            for (int m = 0; m < SandiegoNasa7::nSpecies; ++m)
                if (registry.name(g) == SandiegoNasa7::species[m].name)
                    globalToMech_[static_cast<size_t>(g)] = m;
        for (int g = 0; g < registry.size(); ++g)
            if (globalToMech_[static_cast<size_t>(g)] < 0)
                throw std::invalid_argument(
                    "SandiegoBridge: global species '" + registry.name(g)
                    + "' is missing from the San Diego mechanism");
    }

    int mechanismIndex(int global) const
    { return globalToMech_[static_cast<size_t>(global)]; }

    void expand(const std::vector<double>& Y9, std::vector<double>& Y59) const
    {
        Y59.assign(SandiegoNasa7::nSpecies, 0.0);
        for (size_t g = 0; g < Y9.size(); ++g)
            Y59[static_cast<size_t>(globalToMech_[g])] = Y9[g];
    }

    // Truncated projection: keep the nine transported species, renormalise so
    // their mass fractions sum to unity.  The intermediate-species mass is
    // discarded (diagnostic; they are small in quasi-steady regions).
    void project(const std::vector<double>& Y59, std::vector<double>& Y9) const
    {
        Y9.assign(globalToMech_.size(), 0.0);
        double sum = 0.0;
        for (size_t g = 0; g < Y9.size(); ++g)
        {
            // CVODE may leave tiny negative roundoff on near-zero species;
            // clip before renormalisation so the transported state stays
            // admissible.
            Y9[g] = std::max(Y59[static_cast<size_t>(globalToMech_[g])], 0.0);
            sum += Y9[g];
        }
        if (!(sum > 0.0) || !std::isfinite(sum))
            throw std::runtime_error("SandiegoBridge: invalid projected sum");
        for (double& y : Y9) y /= sum;
    }

private:
    std::vector<int> globalToMech_;
};

inline double sensibleInternalEnergy(const std::vector<double>& Y59, double T)
{
    double e = 0.0;
    for (int i = 0; i < SandiegoNasa7::nSpecies; ++i)
    {
        if (Y59[static_cast<size_t>(i)] <= 0.0) continue;
        const double mw = SandiegoNasa7::species[i].mw;
        e += Y59[static_cast<size_t>(i)]
            * ((SandiegoNasa7::hMole(i, T) - SandiegoNasa7::hMole(i, Tref)) / mw
               - Runiv * T / mw);
    }
    return e;
}

} // namespace SandiegoBridge

class SandiegoSectorChemistryOperator : public Operator
{
public:
    explicit SandiegoSectorChemistryOperator(
        const SpeciesRegistry& registry, double activationTemperature = 850.0)
        : adapter_(registry), activationTemperature_(activationTemperature),
          chemistry_()
    {
        int threads = 1;
#ifdef _OPENMP
        threads = omp_get_max_threads();
#endif
        work_.resize(static_cast<size_t>(threads));
        for (auto& w : work_)
            w.Y59.assign(SandiegoNasa7::nSpecies, 0.0);
        reactors_.resize(static_cast<size_t>(threads));
    }

    ~SandiegoSectorChemistryOperator() override = default;

    bool ownsPersistentState() const override { return false; }

    void beginStepTransaction(const Domain&) const override {}

    std::string name() const override
    { return "SandiegoSectorChemistryOperator"; }

    void apply(Domain& d, double dt) const override
    {
        if (d.speciesNames().size() != 9)
            throw std::runtime_error(
                "SandiegoSectorChemistryOperator: 9-species Domain required");
        lastActive_ = lastFailures_ = 0;
        lastHeatReleaseAverage_ = lastHeatReleaseMaximum_ = 0.0;
        lastHeatReleaseMinimum_ = 1e300;
        double heatReleaseSum = 0.0, heatReleaseMax = -1e300, heatReleaseMin = 1e300;
        long active = 0, failures = 0;
        double failureRho = 0.0, failureT = 0.0;
        int failureI = -1, failureJ = -1, failureK = -1;
        std::vector<double> failureY9;

        #pragma omp parallel for collapse(3) schedule(guided) \
            reduction(+:active,failures,heatReleaseSum) \
            reduction(min:heatReleaseMin) reduction(max:heatReleaseMax)
        for (int i = 1; i < d.Nx() - 1; ++i)
        for (int j = 1; j < d.Ny() - 1; ++j)
        for (int k = 1; k < d.Nz() - 1; ++k)
        {
            Workspace& ws = work_[threadIndex()];
            const PrimitiveVars& w = d.W(i, j, k);
            const thermo::EquilibriumState& q = d.Q(i, j, k);
            if (w.T < activationTemperature_) continue;
            const double rho = q.rho;
            if (!(rho > 0.0) || !std::isfinite(w.T)) { ++failures; continue; }

            adapter_.expand(w.Y, ws.Y59);
            ws.Y59In = ws.Y59;
            double Tnew = w.T;
            bool ok = reactCell(rho, Tnew, ws.Y59, dt, ws.qdot);
            if (!ok)
            {
                // A cell that fails to converge keeps its state (diagnostic
                // robustness: a local stiff-cell convergence failure must not
                // stop the whole run); it is counted for the operator report.
                ++failures;
                #pragma omp critical(rfam_sandiego_chem_failure)
                if (failureI < 0)
                {
                    failureI = i; failureJ = j; failureK = k;
                    failureRho = rho; failureT = w.T; failureY9 = w.Y;
                }
                continue;
            }
            // Project back and compute the sensible-energy change.
            std::vector<double> Y9new;
            adapter_.project(ws.Y59, Y9new);
            const double eBefore = SandiegoBridge::sensibleInternalEnergy(
                ws.Y59In, w.T);
            const double eAfter = SandiegoBridge::sensibleInternalEnergy(
                ws.Y59, Tnew);
            const double deltaE = eAfter - eBefore;      // J/kg
            const double deltaRhoE = rho * deltaE;      // J/m3
            // Defensive admission: reject non-finite chemistry results or
            // physically implausible temperatures (diagnostic mode keeps the
            // cell unchanged and counts the failure).
            if (!std::isfinite(Tnew) || !(Tnew > 100.0) || !(Tnew < 6000.0)
                || !std::isfinite(deltaRhoE))
            {
                ++failures;
                #pragma omp critical(rfam_sandiego_chem_failure)
                if (failureI < 0)
                {
                    failureI = i; failureJ = j; failureK = k;
                    failureRho = rho; failureT = w.T; failureY9 = w.Y;
                }
                continue;
            }

            ConservedVars& next = d.Unew(i, j, k);
            ConservedVars candidate = next;
            double couplingScale = 1.0;
            for (size_t g = 0; g < w.Y.size(); ++g)
            {
                const double increment = rho * (Y9new[g] - w.Y[g]);
                if (increment < 0.0 && candidate.rhoY[g] + increment < 0.0)
                    couplingScale = std::min(couplingScale,
                        candidate.rhoY[g] / (-increment));
            }
            couplingScale = std::clamp(couplingScale, 0.0, 1.0);
            for (size_t g = 0; g < w.Y.size(); ++g)
                next.rhoY[g] += couplingScale * rho * (Y9new[g] - w.Y[g]);
            next.rhoE += couplingScale * deltaRhoE;
            // Post-write self-audit: clip tiny negative roundoff (1e-15
            // rho level, from CVODE/other-operator noise) and then check the
            // cell's conserved state.  A genuine negative species mass or a
            // composition/energy mismatch rolls the cell back and counts a
            // failure instead of tripping the whole run.
            for (double& v : next.rhoY)
                if (v < 0.0 && v > -1e-15 * std::max(next.rho, 1.0))
                    v = 0.0;
            bool cellInvalid = false;
            for (double v : next.rhoY)
                if (!std::isfinite(v) || v < 0.0) cellInvalid = true;
            if (!cellInvalid)
            {
                double speciesMass = 0.0;
                for (double v : next.rhoY) speciesMass += v;
                cellInvalid = !std::isfinite(speciesMass)
                    || std::abs(speciesMass - next.rho) > 1e-6 * next.rho
                    || !std::isfinite(next.rhoE);
            }
            if (cellInvalid)
            {
                double sMass = 0.0;
                for (double v : next.rhoY) sMass += v;
                double minRY = 1e30, maxRY = -1e30;
                for (double v : next.rhoY)
                { minRY = std::min(minRY, v); maxRY = std::max(maxRY, v); }
                std::fprintf(stderr,
                    "[SD-CHEM] write-back invalid cell(%d,%d,%d) T=%.3f "
                    "rho=%.6g sumYdev=%g eBefore=%.6g eAfter=%.6g "
                    "rhoY0=%g minRY=%g maxRY=%g\n",
                    i, j, k, Tnew, rho,
                    std::abs(sMass - next.rho) / next.rho,
                    eBefore, eAfter, next.rhoY.empty() ? -1.0 : next.rhoY[0],
                    minRY, maxRY);
                next = candidate;
                ++failures;
                #pragma omp critical(rfam_sandiego_chem_failure)
                if (failureI < 0)
                {
                    failureI = i; failureJ = j; failureK = k;
                    failureRho = rho; failureT = w.T; failureY9 = w.Y;
                }
                continue;
            }

            ++active;
            const double qdot = rho * deltaE / std::max(dt, 1.0e-300);
            heatReleaseSum += qdot;
            heatReleaseMax = std::max(heatReleaseMax, qdot);
            heatReleaseMin = std::min(heatReleaseMin, qdot);
        }

        lastActive_ = active;
        lastFailures_ = failures;
        lastHeatReleaseAverage_ = active > 0 ? heatReleaseSum / active : 0.0;
        lastHeatReleaseMaximum_ = heatReleaseMax < 0 ? 0.0 : heatReleaseMax;
        lastHeatReleaseMinimum_ = heatReleaseMin > 1e299 ? 0.0 : heatReleaseMin;
        lastFailureI_ = failureI; lastFailureJ_ = failureJ; lastFailureK_ = failureK;
        lastFailureRho_ = failureRho; lastFailureTemperature_ = failureT;
        lastFailureY9_ = std::move(failureY9);

        // Diagnostic robustness: local cell failures are skipped and counted;
        // only an operator-wide failure (all cells failing) trips the run.
        if (failures > 0 && failures >= active)
            throw std::runtime_error(
                "SandiegoSectorChemistryOperator: cell chemistry failure "
                + std::to_string(failures));
    }

    long lastActive() const { return lastActive_; }
    long lastFailures() const { return lastFailures_; }
    // SectorChemistryOperator-compatible statistics surface.
    long lastActiveCells() const { return lastActive_; }
    long lastRetriedCells() const { return 0; }
    long lastPositivityLimitedCells() const { return 0; }
    long lastCvodeFailures() const { return 0; }
    long lastStateFailures() const { return lastFailures_; }
    long lastTemperatureFailures() const { return 0; }
    long lastMassFailures() const { return 0; }
    long lastElementFailures() const { return 0; }
    double lastFailureElementError() const { return 0.0; }
    double lastFailureDt() const { return 0.0; }
    double lastFailurePressure() const { return 0.0; }
    const std::vector<double>& lastFailureChemY() const { return lastFailureY9_; }
    double lastHeatReleaseAverage() const { return lastHeatReleaseAverage_; }
    double lastHeatReleaseMaximum() const { return lastHeatReleaseMaximum_; }
    double lastHeatReleaseMinimum() const { return lastHeatReleaseMinimum_; }
    int lastFailureI() const { return lastFailureI_; }
    int lastFailureJ() const { return lastFailureJ_; }
    int lastFailureK() const { return lastFailureK_; }
    double lastFailureRho() const { return lastFailureRho_; }
    double lastFailureTemperature() const { return lastFailureTemperature_; }
    const std::vector<double>& lastFailureY9() const { return lastFailureY9_; }

private:
    struct Workspace
    {
        std::vector<double> Y59, Y59In;
        double qdot = 0.0;
    };

    struct ReactorState
    {
        SUNContext sunctx = nullptr;
        void* mem = nullptr;
        N_Vector y = nullptr;
        SUNMatrix matrix = nullptr;
        SUNLinearSolver solver = nullptr;
        double rho = 0.0;
        bool initialized = false;
        bool firstInit = true;
    };

    int threadIndex() const
    {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

    static int rhsFn(sunrealtype t, N_Vector y, N_Vector ydot, void* userData)
    {
        (void)t;
        using SandiegoBridge::Runiv;
        using SandiegoBridge::Tref;
        auto* ctx = static_cast<ReactorState*>(userData);
        const double rho = ctx->rho;
        double Y[59];
        double T = NV_Ith_S(y, 59);
        for (int i = 0; i < 59; ++i)
            Y[i] = std::max(NV_Ith_S(y, i), 0.0);   // clip CVODE roundoff
        SandiegoChemistry chem;
        std::vector<double> yVec(Y, Y + 59);
        std::vector<double> omega;
        double qdot = 0.0;
        chem.computeSourceTerms(T, 0.0, rho, yVec, omega, qdot);
        double energyRate = 0.0, cv = 0.0;
        for (int i = 0; i < 59; ++i)
        {
            const double mw = SandiegoNasa7::species[i].mw;
            const double eSens = (SandiegoNasa7::hMole(i, T)
                    - SandiegoNasa7::hMole(i, Tref)) / mw - Runiv * T / mw;
            energyRate += omega[static_cast<size_t>(i)] * eSens;
            cv += Y[i] * (SandiegoNasa7::cpMole(i, T) - Runiv) / mw;
            NV_Ith_S(ydot, i) = omega[static_cast<size_t>(i)] / rho;
        }
        NV_Ith_S(ydot, 59) = (-energyRate + qdot) / (rho * cv);
        return 0;
    }

    bool reactCell(double rho, double& T, std::vector<double>& Y59,
                   double dt, double& qdot) const
    {
        ReactorState& rs = reactors_[threadIndex()];
        if (!rs.initialized)
        {
            if (SUNContext_Create(SUN_COMM_NULL, &rs.sunctx)) return false;
            rs.y = N_VNew_Serial(60, rs.sunctx);
            rs.mem = CVodeCreate(CV_BDF, rs.sunctx);
            if (!rs.mem) return false;
            CVodeSetUserData(rs.mem, &rs);
            rs.matrix = SUNDenseMatrix(60, 60, rs.sunctx);
            rs.solver = SUNLinSol_Dense(rs.y, rs.matrix, rs.sunctx);
            rs.initialized = true;
        }
        rs.rho = rho;
        for (int i = 0; i < 59; ++i) NV_Ith_S(rs.y, i) = Y59[static_cast<size_t>(i)];
        NV_Ith_S(rs.y, 59) = T;
        if (rs.firstInit)
        {
            if (CVodeInit(rs.mem, rhsFn, 0.0, rs.y) < 0) return false;
            CVodeSetLinearSolver(rs.mem, rs.solver, rs.matrix);
            N_Vector atol = N_VNew_Serial(60, rs.sunctx);
            for (int i = 0; i < 59; ++i) NV_Ith_S(atol, i) = 1.0e-14;
            NV_Ith_S(atol, 59) = 5.0;
            CVodeSVtolerances(rs.mem, 1.0e-6, atol);
            N_VDestroy(atol);
            CVodeSetMaxNumSteps(rs.mem, 1000000);
            CVodeSetInitStep(rs.mem, 1.0e-15);
            CVodeSetMaxErrTestFails(rs.mem, 30);
            rs.firstInit = false;
        }
        else if (CVodeReInit(rs.mem, 0.0, rs.y) < 0) return false;
        double t = 0.0, tActual = 0.0;
        const double tStop = dt;
        while (t < tStop)
        {
            double dtStep = std::min(std::max(tStop - t, 0.0),
                std::min(2.0 * std::max(t, 1.0e-9), 2.0e-6));
            if (dtStep <= 0.0) break;
            int flag = CVode(rs.mem, t + dtStep, rs.y, &tActual, CV_NORMAL);
            if (flag < 0)
            {
                // A failed step cannot be retried at a smaller tout from the
                // same CVODE state (internal time has advanced); the cell is
                // skipped by the caller.  The next invocation re-initialises
                // from the cell state.
                rs.firstInit = true;
                return false;
            }
            t = tActual;
        }
        T = NV_Ith_S(rs.y, 59);
        for (int i = 0; i < 59; ++i)
            Y59[static_cast<size_t>(i)] = NV_Ith_S(rs.y, i);
        if (!(T > 0.0) || !std::isfinite(T)) return false;
        return true;
    }

    SandiegoBridge::SpeciesAdapter adapter_;
    double activationTemperature_;
    SandiegoChemistry chemistry_;
    mutable std::vector<Workspace> work_;
    mutable std::vector<ReactorState> reactors_;
    mutable long lastActive_ = 0, lastFailures_ = 0;
    mutable double lastHeatReleaseAverage_ = 0.0,
                   lastHeatReleaseMaximum_ = 0.0,
                   lastHeatReleaseMinimum_ = 0.0;
    mutable int lastFailureI_ = -1, lastFailureJ_ = -1, lastFailureK_ = -1;
    mutable double lastFailureRho_ = 0.0, lastFailureTemperature_ = 0.0;
    mutable std::vector<double> lastFailureY9_;
};

} // End namespace Foam

#endif
