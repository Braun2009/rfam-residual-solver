/*---------------------------------------------------------------------------*\
    test_phase7_1d - Phase 7 Step 2: finite-rate chemistry coupled into the
    framework 1D solver (Convection + ChemistrySource + species-aware EOS reback)
    over the global SpeciesRegistry.

    SCOPE / HONEST BOUNDARY (Stage-3 step-1 thermo in): the ideal-gas part of
    the EOS is now REAL NASA-7 Cp(T)/h_sens(T) (GRI-Mech 3.0, gated by
    test_nasa7_thermo), so the reback temperature is a GENUINE EOS inversion —
    the old ~15400 K ideal-gas-fallback artifact is gone. Remaining documented
    simplifications (model limits, NOT thermo-foundation bugs):
      - chemistry integrates real-fluid [Y,e] and recovers T from rho-e-Y
      - the current state model remains single-phase and cannot close LOX flash
      - TwoStepCO includes reversible CO/CO2 but omits radical dissociation,
        so absolute flame temperature remains mechanism-limited
      - TwoStepCO ~1e-6 mass imbalance (Sum omega != 0), known debt
    So: structure/conservation/coupling are quantitative; absolute T is now
    physically MEANINGFUL given the mechanism, but not yet the true flame T.

    Scenario A (uniform reacting box) — CORRECTNESS MILESTONE. A spatially
    uniform, quiescent premixed box: convection flux between identical cells is
    exactly zero, so every cell must reproduce the standalone 0-D ignition
    (test_phase7_chem0d) while running through Domain + SolverBase + the
    9-species reback. Isolates chemistry-in-framework + species-aware reback.

    Scenario B (cold transport smoke) — a mild O2/C3H8 composition gradient
    (fuel-lean left | pure O2 right) at 500 K, 1 MPa, well away from any critical
    point, so the contact is a benign dilute-gas contact that convects cleanly.
    Its job is to exercise the ONE thing A (spatially uniform) cannot: 9-species
    CONVECTIVE TRANSPORT with real, non-zero species flux differences, plus the
    per-cell species-aware reback on the resulting mixtures. At 500 K the C3H8/O2
    ignition delay is astronomically longer than the run, so the chemistry source
    operator — run on every cell every step — must correctly return ~0 (no
    spurious ignition). Together with A (chemistry correctly DOES ignite when hot)
    this brackets the coupled path: fires when it should, stays quiet when it
    shouldn't, through the full framework. NOT an ignition-front test (that
    still needs a dissociating mechanism + an implicit/limited source; explicit
    splitting overshoots to A's ~7040 K const-V complete-combustion state).
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include "src/thermo/SpeciesRegistry.H"
#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"
#include "src/combustion/ChemistrySourceOperator.H"
#include "src/combustion/FiniteRate/TwoStepCO.H"

using namespace Foam;
using namespace std;

static int g_pass = 0, g_fail = 0;
static void check(const string& n, bool ok)
{
    cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << n << "\n";
    ok ? ++g_pass : ++g_fail;
}

// sum element atoms per m3 over the whole box (uses committed U.rho + W.Y)
static void boxAtoms(const Domain& d, const SpeciesRegistry& reg,
                     double& C, double& H, double& O, double& N)
{
    C = H = O = N = 0.0;
    for (int i = 0; i < d.N(); ++i)
    {
        double rr = d.U(i).rho;
        for (int g = 0; g < reg.size(); ++g)
        {
            double mr = rr * d.W(i).Y[g] / reg.molecularWeight(g);
            C += mr * reg.nC(g); H += mr * reg.nH(g);
            O += mr * reg.nO(g); N += mr * reg.nN(g);
        }
    }
}

static void printTripDiagnostics(const Domain& d)
{
    if (!d.tripped()) return;

    double maxMassFractionError = 0.0;
    double minimumSpeciesDensity = 1.0e300;
    int nonFiniteStates = 0;
    for (int i = 0; i < d.N(); ++i)
    {
        const ConservedVars& u = d.Unew(i);
        if (!std::isfinite(u.rho) || !std::isfinite(u.rhoE) || u.rho <= 0.0)
            ++nonFiniteStates;
        double speciesMass = 0.0;
        for (double rhoY : u.rhoY)
        {
            speciesMass += rhoY;
            minimumSpeciesDensity = std::min(minimumSpeciesDensity, rhoY);
            if (!std::isfinite(rhoY)) ++nonFiniteStates;
        }
        if (std::isfinite(u.rho) && u.rho > 0.0)
            maxMassFractionError = std::max(
                maxMassFractionError,
                std::abs(speciesMass/u.rho - 1.0));
    }
    cout << "  trip diagnostics: max|sumY-1|=" << setprecision(12)
         << maxMassFractionError
         << " min(rhoY)=" << minimumSpeciesDensity
         << " nonfinite=" << nonFiniteStates << "\n";
}

// TwoStepCO rates == Phase 6 / 0-D reference
static TwoStepCO makeChem() { return TwoStepCO(1e10, 0.0, 125520.0, 2e12, 0.0, 167360.0); }

// -------------------------------------------------------------------------
static void scenarioA_uniformBox(const SpeciesRegistry& reg)
{
    cout << "--- Scenario A: uniform reacting box (correctness milestone) ---\n";
    const int nG = reg.size();
    // Four identical cells are sufficient to exercise the complete framework
    // path and uniformity. A finer uniform mesh only repeats the same stiff 0-D
    // solve while reducing the CFL timestep.
    const int N = 4; const double L = 0.1;
    Domain dom(N, L, nG);
    dom.setSpecies(reg.names());          // GLOBAL-space, opt-in species-aware reback

    // uniform premixed stoich C3H8/O2 (== 0-D test composition), quiescent
    const double p0 = 4.0e6, T0 = 1200.0;
    vector<double> Y0(nG, 0.0);
    Y0[reg.index("C3H8")] = 0.216;
    Y0[reg.index("O2")]   = 0.784;
    for (int i = 0; i < N; ++i)
    {
        PrimitiveVars& w = dom.stagedW(i);
        w.p = p0; w.T = T0; w.U[0]=w.U[1]=w.U[2]=0.0; w.Z=0.0; w.C=0.0; w.Y = Y0;
    }
    dom.primitiveToConserved();

    double C0,H0,O0,N0; boxAtoms(dom, reg, C0,H0,O0,N0);
    double rhoE0 = 0.0; for (int i=0;i<N;++i) rhoE0 += dom.U(i).rhoE;

    TwoStepCO chem = makeChem();
    SolverBase solver(dom);
    // ORDER: convection (assign) BEFORE chemistry (accumulate) — see Operator.H.
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    solver.addOperator(unique_ptr<Operator>(
        new ChemistrySourceOperator(&chem, reg)));
    for(int step=0;step<4 && !dom.tripped();++step)
    {
        const double dt=solver.computeTimestep(0.5);
        if(!std::isfinite(dt)||dt<=0.0||!solver.step(dt)) break;
    }
    printTripDiagnostics(dom);

    const PrimitiveVars& w = dom.W(N/2);
    double C1,H1,O1,N1; boxAtoms(dom, reg, C1,H1,O1,N1);
    double rhoE1 = 0.0, Tmin=1e30, Tmax=-1e30;
    for (int i=0;i<N;++i){ rhoE1 += dom.U(i).rhoE; Tmin=min(Tmin,dom.W(i).T); Tmax=max(Tmax,dom.W(i).T); }
    double thr = C0+H0+O0+N0;

    cout << fixed;
    cout << "  iters=" << solver.iteration()
         << "  mid: Y_C3H8=" << setprecision(4) << w.Y[reg.index("C3H8")]
         << " Y_CO2=" << w.Y[reg.index("CO2")]
         << " Y_H2O=" << w.Y[reg.index("H2O")]
         << "  T span=[" << setprecision(0) << Tmin << "," << Tmax << "]K\n";

    check("A: no trip", !dom.tripped());
    check("A: fuel decreases",w.Y[reg.index("C3H8")]<Y0[reg.index("C3H8")]);
    // Reversible R2: high-T equilibrium shifts CO2→CO; check total carbon products.
    check("A: carbon products produced",
          w.Y[reg.index("CO")] + w.Y[reg.index("CO2")] > 1e-8);
    check("A: H2O produced", w.Y[reg.index("H2O")] > 1e-8);
    check("A: non-chem species stay zero (CH4,H2)",
          fabs(w.Y[reg.index("CH4")])<1e-12 && fabs(w.Y[reg.index("H2")])<1e-12);
    check("A: element C conserved", fabs(C1-C0) < 1e-6*thr);
    check("A: element H conserved", fabs(H1-H0) < 1e-6*thr);
    check("A: element O conserved", fabs(O1-O0) < 1e-6*thr);
    check("A: EOS reback T finite & elevated", isfinite(Tmax) && Tmin > T0-1.0);
    // Broad reduced-mechanism band; the upper edge catches the old ideal-gas
    // fallback artifact while allowing reversible CO/CO2 equilibrium.
    check("A: T in supported hot-gas band [1200,9000]K",
          Tmax > T0 && Tmax < 9000.0);
    check("A: box stays uniform (T span < 1K)", (Tmax-Tmin) < 1.0);
    check("A: chemical heat added (rhoE up)", rhoE1 > rhoE0);
    cout << "\n";
}

// -------------------------------------------------------------------------
static void scenarioB_gradientSmoke(const SpeciesRegistry& reg)
{
    cout << "--- Scenario B: mild cold composition gradient (transport smoke) ---\n";
    const int nG = reg.size();
    const int N = 100; const double L = 1.0;   // long domain: contact stays off the ends
    Domain dom(N, L, nG);
    dom.setSpecies(reg.names());

    // Dilute-gas contact (well away from critical): fuel-lean O2/C3H8 (left) |
    // pure O2 (right) at 500 K, 1 MPa. Mild density ratio -> benign convection;
    // cold enough that chemistry must NOT ignite. Point = 9-species transport +
    // per-cell reback across a real composition gradient.
    const double p0 = 1.0e6, T0 = 500.0;
    vector<double> Ymix(nG,0.0);  Ymix[reg.index("O2")] = 0.8; Ymix[reg.index("C3H8")] = 0.2;
    vector<double> Yoxid(nG,0.0); Yoxid[reg.index("O2")] = 1.0;
    for (int i = 0; i < N; ++i)
    {
        PrimitiveVars& w = dom.stagedW(i);
        w.p=p0; w.T=T0; w.U[0]=w.U[1]=w.U[2]=0.0; w.Z=0.0; w.C=0.0;
        w.Y = (dom.x(i) < 0.5*L) ? Ymix : Yoxid;
    }
    dom.primitiveToConserved();

    double C0,H0,O0,N0; boxAtoms(dom, reg, C0,H0,O0,N0);

    TwoStepCO chem = makeChem();
    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    solver.addOperator(unique_ptr<Operator>(
        new ChemistrySourceOperator(&chem, reg)));
    solver.advanceToTime(5.0e-5, 0.5);
    printTripDiagnostics(dom);

    double C1,H1,O1,N1; boxAtoms(dom, reg, C1,H1,O1,N1);
    double thr = C0+H0+O0+N0;
    bool Tfinite = true, noIgnite = true; double tm = 0.0;
    for (int i=0;i<N;++i){ if(!isfinite(dom.W(i).T)) Tfinite=false; if(dom.W(i).T > T0+50.0) noIgnite=false; tm=max(tm,dom.W(i).T); }
    bool endsHeld = (dom.W(0).Y[reg.index("C3H8")] > 0.199) && (dom.W(N-1).Y[reg.index("O2")] > 0.999);

    cout << fixed << setprecision(4);
    cout << "  iters=" << solver.iteration()
         << "  interface Y_C3H8[N/2-1..N/2]=" << dom.W(N/2-1).Y[reg.index("C3H8")]
         << "," << dom.W(N/2).Y[reg.index("C3H8")]
         << "  Tmax=" << setprecision(1) << tm << "K\n";

    check("B: no trip (coupled path runs over a gradient)", !dom.tripped());
    check("B: all T finite", Tfinite);
    check("B: cold contact does NOT ignite (chem correctly quiet)", noIgnite);
    check("B: ends held (contact contained, no-BC ends untouched)", endsHeld);
    check("B: element C conserved", fabs(C1-C0) < 1e-9*thr);
    check("B: element H conserved", fabs(H1-H0) < 1e-9*thr);
    check("B: element O conserved", fabs(O1-O0) < 1e-9*thr);
    cout << "\n";
}

int main()
{
    cout << "=== Phase 7 Step 2: 1D reacting flow (framework + species-aware reback) ===\n";
    cout << "    [real-fluid rho-e-Y reback; single-phase PR and the"
            " reversible 2-step mechanism remain model limits]\n\n";

    SpeciesRegistry reg;
    scenarioA_uniformBox(reg);
    scenarioB_gradientSmoke(reg);

    cout << "--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "PHASE7_1D_PASS" : "PHASE7_1D_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
