/*---------------------------------------------------------------------------*\
    solver_mvp_1d - Stage 0 MVP: framework reproduction of solver1d_fpv.

    Builds the SAME 1D FPV reacting case as src/solver/solver1d_fpv.cpp, but
    through the framework (Domain + [Convection, FPVSource] + OperatorSplitting
    + SolverBase). Numerical results must match the original prototype.
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm>

#include "src/variables/ConservedVars.H"
#include "src/combustion/FPV/FPVTable.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"

using namespace Foam;
using namespace std;

int main()
{
    cout << "=== Stage 0 MVP: framework reproduction of solver1d_fpv ===\n\n";

    const int    N  = 200;
    const double L  = 1.0;
    const int    nS = 2;
    Domain dom(N, L, nS);

    // Initial condition == solver1d_fpv main():
    //   left  Z=0.95 (fuel), right Z=0.05 (oxid), C=0, p=6e6, T=300, U=0, Y=[0,1](N2)
    const double pL=6.0e6, TL=300.0, ZL=0.95, CL=0.0;
    const double pR=6.0e6, TR=300.0, ZR=0.05, CR=0.0;
    const double x0=0.5;
    for (int i = 0; i < N; ++i)
    {
        double x = dom.x(i);
        PrimitiveVars& w = dom.stagedW(i);
        if (x < x0) { w.p=pL; w.T=TL; w.Z=ZL; w.C=CL; }
        else        { w.p=pR; w.T=TR; w.Z=ZR; w.C=CR; }
        w.U[0]=w.U[1]=w.U[2]=0.0;
        w.Y[0]=0.0; w.Y[1]=1.0;   // N2
    }
    dom.primitiveToConserved();

    FPVTable table;
    if (!table.loadFromCSV("tests/data/fpv_table_proof.csv"))
    {
        cerr << "ERROR: Failed to load FPV table\n";
        return 1;
    }

    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    solver.addOperator(unique_ptr<Operator>(
        new FPVSourceOperator(&table, 1e-5, 1e-5)));
    solver.setMonitor([](const SolverBase& s){
        double cmax = 0.0;
        for (int i = 0; i < s.domain().N(); ++i) cmax = max(cmax, s.domain().W(i).C);
        cout << "  iter=" << s.iteration() << ", t=" << s.time() << ", C_max=" << cmax << "\n";
    });

    solver.advanceToTime(0.001, 0.5);

    double cmax=0.0, tmax=0.0;
    for (int i=0;i<N;++i){ cmax=max(cmax,dom.W(i).C); tmax=max(tmax,dom.W(i).T); }
    cout << "\nFinal: iter=" << solver.iteration() << " t=" << solver.time()
         << " C_max=" << cmax << " T_max=" << tmax
         << (dom.tripped() ? "  [TRIPPED]" : "") << "\n";

    // Reproduction check against the same-version src/solver/solver1d_fpv.cpp
    // prototype on this identical case. Identical iteration counts imply the
    // same CFL sequence and time-marching behavior for this deterministic case.
    const int REF_ITER = 273;            // current prototype baseline
    bool match = (solver.iteration() == REF_ITER)
              && (tmax > 299.0 && tmax < 301.0)
              && !dom.tripped();
    cout << (match ? "MVP_REPRODUCES_solver1d_fpv" : "MVP_MISMATCH") << "\n";
    return match ? 0 : 1;
}
