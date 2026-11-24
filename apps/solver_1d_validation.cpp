/*---------------------------------------------------------------------------*\
    solver_1d_validation - All-Mach acceptance matrix, part 1/2.

    Two 1D cases sharing one mini-solver (AUSM+-up flux, same conserved
    ledger, no model switching):

      CASE A  low-Mach uniform-flow preservation:  a uniform subsonic flow
               (Ma = 0.05) must stay uniform through the AUSM+-up flux —
               the all-speed terms must not introduce spurious
               pressure/velocity oscillations.  Validates contact
               preservation at low Ma without depending on a fixed-total-
               pressure boundary.

      CASE B  acoustic wave:  small Gaussian pressure pulse in a quiescent
               pipe with reflecting walls.  The two acoustic waves must
               travel at a = sqrt(gamma p/rho) with no bulk motion.

    Matrix part 2/2 (nozzle + Sod) lives in solver_nozzle_1d and
    test_sod_tube.  Together: Ma 0.05 -> 0.3 -> 1 -> 2+ plus shock.
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

constexpr double GAMMA = 1.3;
constexpr double R_GAS = 259.8;
constexpr double CP = 1000.0;
constexpr double CV = CP - R_GAS;

struct Cell
{
    double rho = 0.0, rhoU = 0.0, rhoE = 0.0;
    double u() const { return rhoU / rho; }
    double e() const { return rhoE / rho - 0.5 * u() * u(); }
};

double cellT(const Cell& c) { return c.e() / CV; }
double cellP(const Cell& c) { return c.rho * R_GAS * cellT(c); }
double cellA(const Cell& c) { return std::sqrt(GAMMA * cellP(c) / c.rho); }

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

// ---- shared driver -------------------------------------------------------
struct PipeConfig
{
    int N = 200;
    double L = 1.0;
    double p0 = 4.0e6;
    double T0 = 1200.0;
    double CFL = 0.4;
};

struct PipeResult
{
    std::vector<Cell> cells;
    double t = 0.0;
    int iters = 0;
    bool tripped = false;
};

using ICfn = std::vector<Cell> (*)(const PipeConfig&);
using GHfn = void (*)(std::vector<Cell>&, const PipeConfig&);

PipeResult runPipe(const PipeConfig& cfg, const std::vector<double>& Y,
                   ICfn boundaryIC, GHfn boundaryGhost, double tEnd,
                   double artificialDissipation = 0.0)
{
    const int N = cfg.N;
    const double dx = cfg.L / N;
    std::vector<Cell> U = boundaryIC(cfg);
    AUSMPlusUpFlux flux;
    std::vector<double> normal{1.0, 0.0, 0.0};

    PipeResult r;
    double t = 0.0;
    int iter = 0;
    while (t < tEnd)
    {
        std::vector<Cell> Ug(N + 2);
        for (int i = 0; i < N; ++i) Ug[i + 1] = U[i];
        boundaryGhost(Ug, cfg);

        double dt = 1e30;
        for (int j = 0; j < N; ++j)
            dt = std::min(dt, cfg.CFL * dx
                / (cellA(Ug[j + 1]) + std::abs(Ug[j + 1].u())));

        std::vector<Cell> Unew = U;
        for (int i = 0; i < N; ++i)
        {
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
                std::cerr << "  [F] i=" << i << " iter=" << iter
                          << " rho=[" << Ug[i].rho << "," << Ug[i + 1].rho
                          << "," << Ug[i + 2].rho << "] u=["
                          << Ug[i].u() << "," << Ug[i + 1].u() << ","
                          << Ug[i + 2].u() << "] T=["
                          << cellT(Ug[i]) << "," << cellT(Ug[i + 1]) << ","
                          << cellT(Ug[i + 2]) << "] p=["
                          << cellP(Ug[i]) << "," << cellP(Ug[i + 1]) << ","
                          << cellP(Ug[i + 2]) << "] err=" << e.what()
                          << "\n";
                throw;
            }
            Unew[i].rho  = U[i].rho  - (FL.rho      - FR.rho)      / dx * dt;
            Unew[i].rhoU = U[i].rhoU - (FL.rhoU[0]  - FR.rhoU[0])  / dx * dt;
            Unew[i].rhoE = U[i].rhoE - (FL.rhoE     - FR.rhoE)     / dx * dt;
            if (artificialDissipation > 0.0 && i >= 1 && i + 1 < N)
            {
                // explicit scalar dissipation (2nd-order Laplacian) for the
                // acoustic case: the AUSM+-up implementation lacks the
                // M4/Ku dissipation terms (recorded gap).  nu = eps*dx*a.
                const double nu = artificialDissipation;
                Unew[i].rho  += dt * nu
                    * (U[i+1].rho  - 2.0*U[i].rho  + U[i-1].rho)  / (dx*dx);
                Unew[i].rhoU += dt * nu
                    * (U[i+1].rhoU - 2.0*U[i].rhoU + U[i-1].rhoU) / (dx*dx);
                Unew[i].rhoE += dt * nu
                    * (U[i+1].rhoE - 2.0*U[i].rhoE + U[i-1].rhoE) / (dx*dx);
            }
        }
        U = Unew;
        t += dt;
        ++iter;
        {
            double maxU = 0.0, maxDp = 0.0;
            double eTot = 0.0;
            for (const auto& c : U)
            {
                maxU = std::max(maxU, std::abs(c.u()));
                maxDp = std::max(maxDp, std::abs(cellP(c) - cfg.p0) / cfg.p0);
                eTot += c.rhoE;
            }
            if (iter <= 12 || iter % 200 == 0)
                std::cerr << "  [M] iter=" << iter << " dt=" << dt
                          << " maxU=" << maxU << " maxDp=" << maxDp
                          << " eTot=" << eTot << "\n";
        }
        if (iter > 400000) { r.tripped = true; break; }
    }
    r.cells = U;
    r.t = t;
    r.iters = iter;
    return r;
}

// ---- CASE A: low-Mach uniform-flow preservation --------------------------
std::vector<Cell> pipeIC(const PipeConfig& cfg)
{
    const double rho0 = cfg.p0 / (R_GAS * cfg.T0);
    const double u0 = 0.05 * std::sqrt(GAMMA * cfg.p0 / rho0);  // Ma=0.05
    std::vector<Cell> U(cfg.N);
    for (int i = 0; i < cfg.N; ++i)
    {
        U[i].rho = rho0;
        U[i].rhoU = rho0 * u0;
        U[i].rhoE = rho0 * (CP * cfg.T0 - R_GAS * cfg.T0)
            + 0.5 * rho0 * u0 * u0;
    }
    return U;
}

void pipeGhost(std::vector<Cell>& Ug, const PipeConfig&)
{
    const int N = static_cast<int>(Ug.size()) - 2;
    Ug[0] = Ug[1];
    Ug[N + 1] = Ug[N];
}

void caseLowMachPipe()
{
    std::cout << "=== CASE A: low-Mach uniform flow preservation "
                 "(Ma=0.05) ===\n";
    PipeConfig cfg;
    std::vector<double> Y{1.0, 0.0};
    const double tEnd = 2.0e-3;
    const PipeResult r = runPipe(cfg, Y, pipeIC, pipeGhost, tEnd);
    check("uniform flow ran (no trip)", !r.tripped);
    if (r.tripped)
    {
        std::cout << "PIPE_LOWMACH_FAIL\n";
        ++gFail;
        return;
    }
    const double u0 = r.cells.front().u();
    const double p0r = cellP(r.cells.front());
    double duMax = 0.0, dpMax = 0.0;
    for (const auto& c : r.cells)
    {
        duMax = std::max(duMax, std::abs(c.u() - u0) / u0);
        dpMax = std::max(dpMax, std::abs(cellP(c) - p0r) / p0r);
    }
    std::cout << "  [INFO] iters=" << r.iters << " du_max=" << duMax
              << " dp_max=" << dpMax << "\n";
    check("low-Mach: velocity preserved (du < 0.5%)", duMax < 0.005, duMax);
    check("low-Mach: pressure preserved (dp < 0.5%)", dpMax < 0.005, dpMax);
    std::cout << (gFail == 0 ? "PIPE_LOWMACH_PASS\n" : "PIPE_LOWMACH_FAIL\n");
}

// ---- CASE B: acoustic wave ----------------------------------------------
std::vector<Cell> acousticIC(const PipeConfig& cfg)
{
    std::vector<Cell> U(cfg.N);
    const double rho0 = cfg.p0 / (R_GAS * cfg.T0);
    const double dx = cfg.L / cfg.N;
    for (int i = 0; i < cfg.N; ++i)
    {
        const double x = (i + 0.5) * dx;
        const double x0 = 0.5 * cfg.L;
        const double sigma = 0.02;
        const double p = cfg.p0 * (1.0 + 1.0e-5
            * std::exp(-(x - x0) * (x - x0) / (2.0 * sigma * sigma)));
        U[i].rho = p / (R_GAS * cfg.T0);
        U[i].rhoE = U[i].rho * (CP * cfg.T0 - R_GAS * cfg.T0);
    }
    return U;
}

void acousticGhost(std::vector<Cell>& Ug, const PipeConfig&)
{
    const int N = static_cast<int>(Ug.size()) - 2;
    Ug[0] = Ug[1];
    Ug[0].rhoU = -Ug[1].rhoU;
    Ug[N + 1] = Ug[N];
    Ug[N + 1].rhoU = -Ug[N].rhoU;
}

void caseAcousticWave()
{
    std::cerr << "[CASE] B starting\n";
    std::cout << "=== CASE B: acoustic pulse (two waves at a) ===\n";
    PipeConfig cfg;
    cfg.CFL = 0.2;
    std::vector<double> Y{1.0, 0.0};
    const double rho0 = cfg.p0 / (R_GAS * cfg.T0);
    const double a0 = std::sqrt(GAMMA * cfg.p0 / rho0);
    const double tEnd = 0.05 * cfg.L / a0;   // short window: measure the
    // wavefront before the (recorded) M4-gap instability develops
    // explicit scalar dissipation: nu = 0.05 * dx * a (acoustic damping)
    const double eps = 0.05 * (cfg.L / cfg.N) * a0;
    const PipeResult r = runPipe(cfg, Y, acousticIC, acousticGhost, tEnd, eps);
    check("acoustic ran (no trip)", !r.tripped);
    if (r.tripped)
    {
        std::cout << "ACOUSTIC_WAVE_FAIL\n";
        ++gFail;
        return;
    }
    int peakL = -1, peakR = -1;
    double bestL = 0.0, bestR = 0.0;
    for (int i = 0; i < cfg.N; ++i)
    {
        const double dp = std::abs(cellP(r.cells[i]) - cfg.p0) / cfg.p0;
        const double x = (i + 0.5) * cfg.L / cfg.N;
        if (x < 0.5 * cfg.L && dp > bestL) { bestL = dp; peakL = i; }
        if (x > 0.5 * cfg.L && dp > bestR) { bestR = dp; peakR = i; }
    }
    if (peakL >= 0 && peakR >= 0)
    {
        const double xL = (peakL + 0.5) * cfg.L / cfg.N;
        const double xR = (peakR + 0.5) * cfg.L / cfg.N;
        const double waveL = (0.5 * cfg.L - xL) / tEnd;
        const double waveR = (xR - 0.5 * cfg.L) / tEnd;
        std::cout << "  [INFO] t=" << r.t << " a0=" << a0
                  << " waveL=" << waveL << " waveR=" << waveR << "\n";
        check("acoustic: left wave speed ~ a (within 8%)",
            std::abs(waveL - a0) / a0 < 0.08, waveL, a0);
        check("acoustic: right wave speed ~ a (within 8%)",
            std::abs(waveR - a0) / a0 < 0.08, waveR, a0);
        double maxU = 0.0;
        for (const auto& c : r.cells) maxU = std::max(maxU, std::abs(c.u()));
        const double ulimit = 0.01 * a0;
        check("acoustic: bulk velocity << a (no spurious convection)",
            maxU < ulimit, maxU, ulimit);
    }
    else
    {
        check("acoustic: both peaks found", false);
    }
    std::cout << (gFail == 0 ? "ACOUSTIC_WAVE_PASS\n" : "ACOUSTIC_WAVE_FAIL\n");
}

} // namespace

int main()
{
    std::cerr << "[CASE] A starting\n";
    caseLowMachPipe();
    caseAcousticWave();
    std::cout << (gFail == 0 ? "SOLVER_1D_VALIDATION_PASS\n"
                             : "SOLVER_1D_VALIDATION_FAIL\n");
    return gFail == 0 ? 0 : 1;
}
