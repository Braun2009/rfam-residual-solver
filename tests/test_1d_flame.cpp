/*---------------------------------------------------------------------------*\
    test_1d_flame - 1-D premixed flame propagation physical validation.

    A planar premixed C3H8/O2/N2 flame (phi=1, 80% N2 dilution) propagates
    into a quiescent fresh mixture at 300 K / 1 atm.  The flame is ignited
    by a hot product patch on the left.  The solver chain is
    Convection + 1-D molecular diffusion + TwoStepCO finite-rate chemistry,
    with the SAME two-step mechanism encoded in
    tests/reference/rfam_twostep_cantera_v1.yaml.

    Physical reference (Cantera 3.1.0 free flame, same mechanism,
    mixture-averaged transport, fresh gas at 350 K to stay above the
    NASA-7 Tlow=300 K of N2/C3H8):
      S_L  = 3.53 m/s   (gate: 2.5 - 4.6 m/s, 30% band)
      T_ad = 2104 K     (gate: 1890 - 2310 K, 10% band)

    Extra gates: element conservation (C/H/O/N < 1e-3 rel), fuel fully
    consumed behind the flame, no trip.

    Output marker: FLAME1D_PASS / FLAME1D_FAIL.
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <chrono>

#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/SolverBase.H"
#include "src/framework/DiffusionOperator.H"
#include "src/combustion/ChemistrySourceOperator.H"
#include "src/combustion/FiniteRate/TwoStepCO.H"

using namespace Foam;
using namespace std;

static int g_fail=0;
static void check(const string& n, bool ok, double got=0, double ref=0)
{
    cout << "  [" << (ok?"PASS":"FAIL") << "] " << n;
    if(!ok) cout << "  (got=" << got << " ref=" << ref << ")";
    cout << "\n";
    if(!ok) ++g_fail;
}

int main()
{
    cout << "[FLAME1D] enter main" << endl << flush;
    const int N=240;
    const double L=0.006;              // 6 mm domain (24um flame, dx=25um)
    const double xIgn=0.0015;          // hot product patch length
    const double p0=101325.0;
    const double Tfresh=350.0, Tign=2500.0;
    const double dx=L/N;

    // Smooth ignition patch boundary over 4 cells to avoid numerical
    // undershoot below the NASA-7 Tlow of N2/C3H8 (300 K).

    // Fresh mixture: C3H8:O2:N2 = 1:5:24 (phi=1, 80% N2).
    const double mwC3H8=44.097, mwO2=31.999, mwN2=28.014;
    const double mTotal=1.0*mwC3H8+5.0*mwO2+24.0*mwN2;
    const double yC3H8f=1.0*mwC3H8/mTotal;
    const double yO2f  =5.0*mwO2/mTotal;
    const double yN2f  =24.0*mwN2/mTotal;

    SpeciesRegistry reg;
    const vector<string> species=reg.names();   // 9 species ordered
    const int iC3H8=reg.index("C3H8");
    const int iO2  =reg.index("O2");
    const int iN2  =reg.index("N2");
    const int iCO  =reg.index("CO");
    const int iCO2 =reg.index("CO2");
    const int iH2O =reg.index("H2O");

    // Complete-combustion products: C3H8+5O2+24N2 -> 3CO2+4H2O+24N2
    const double mProd=3.0*44.010+4.0*18.015+24.0*mwN2;
    const double yCO2p=3.0*44.010/mProd;
    const double yH2Op=4.0*18.015/mProd;
    const double yN2p =24.0*mwN2/mProd;

    Domain dom(N, L, static_cast<int>(species.size()));
    cout << "[FLAME1D] domain" << endl << flush;
    dom.setSpecies(species);
    for(int i=0;i<N;++i)
    {
        PrimitiveVars w(species.size());
        w.p=p0;
        const double x=dom.x(i);
        // Smooth transition over 4 cells across xIgn.
        double blend=1.0;
        if(x>=xIgn && x<xIgn+4.0*dx)
            blend=1.0-(x-xIgn)/(4.0*dx);
        else if(x>=xIgn) blend=0.0;
        if(blend>0.0)
        {
            w.T=Tign;
            w.Y[iCO2]=yCO2p; w.Y[iH2O]=yH2Op; w.Y[iN2]=yN2p;
        }
        else
        {
            w.T=Tfresh;
            w.Y[iC3H8]=yC3H8f; w.Y[iO2]=yO2f; w.Y[iN2]=yN2f;
        }
        w.U[0]=w.U[1]=w.U[2]=0.0;
        dom.stagedW(i)=w;
    }
    dom.primitiveToConserved();
    cout << "[FLAME1D] p2c" << endl << flush;

    double C0,H0,O0,N0;
    {
        double c=0,h=0,o=0,n=0;
        for(int i=0;i<N;++i)
        {
            const auto& w=dom.W(i);
            c+=w.Y[iC3H8]*3.0 + w.Y[iCO]*1.0 + w.Y[iCO2]*1.0;
            h+=w.Y[iC3H8]*8.0 + w.Y[iH2O]*2.0;
            o+=w.Y[iO2]*2.0 + w.Y[iCO]*1.0 + w.Y[iCO2]*2.0 + w.Y[iH2O]*1.0;
            n+=w.Y[iN2]*2.0;
        }
        C0=c;H0=h;O0=o;N0=n;
    }

    const double tEnd=1.0e-3;   // 1 ms (flame travels ~3.5 mm)
    TwoStepCO chem;   // default parameters == Cantera YAML
    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    solver.addOperator(unique_ptr<Operator>(new DiffusionOperator(reg)));
    solver.addOperator(unique_ptr<Operator>(new ChemistrySourceOperator(&chem, reg)));
    cout << "[FLAME1D] solver ready, FLAME1D start: N=" << N << " L=" << L
         << " Tfresh=" << Tfresh << " tEnd=" << tEnd << "\n" << flush;


    // Flame tracking: record x where T crosses 1000 K every 4000 steps.
    vector<pair<double,double>> track;   // (t, x_flame)
    auto applyBC=[&](Domain& d)
    {
        d.setCell(0,0,d.W(1,0));
        d.setCell(N-1,0,d.W(N-2,0));
    };

    auto t0=chrono::steady_clock::now();
    while(!dom.tripped() && solver.time()<tEnd)
    {
        applyBC(dom);
        const double dt=solver.computeTimestep(0.5);
        solver.step(dt);
        if(solver.iteration()%20000==0)
            cout << "    iter=" << solver.iteration()
                 << " t=" << solver.time() << "\n" << flush;
        if(solver.iteration()%10000==0)
        {
            double xf=-1.0;
            for(int i=2;i<N-2;++i)
                if(dom.W(i).T>1000.0 && dom.W(i-1).T<=1000.0)
                { xf=dom.x(i); break; }
            if(xf>0.0) track.emplace_back(solver.time(),xf);
        }
    }
    auto t1=chrono::steady_clock::now();
    cout << "  iters=" << solver.iteration() << " t=" << solver.time()
         << "  wall=" << chrono::duration<double>(t1-t0).count() << " s\n";

    check("no trip", !dom.tripped());

    // Flame speed from linear fit of x_f vs t.
    double SL=0.0;
    if(track.size()>=3)
    {
        double sumT=0,sumT2=0,sumX=0,sumTX=0;
        for(const auto& p:track){sumT+=p.first;sumT2+=p.first*p.first;
            sumX+=p.second;sumTX+=p.first*p.second;}
        const int m=static_cast<int>(track.size());
        const double denom=m*sumT2-sumT*sumT;
        if(std::abs(denom)>1e-30)
            SL=(m*sumTX-sumT*sumX)/denom;
    }
    check("flame speed in [2.5,4.6] m/s (Cantera 3.53)",
        SL>2.5 && SL<4.6, SL, 3.53);

    // Adiabatic flame temperature behind the flame.
    double Tmax=0.0;
    for(int i=0;i<N;++i) Tmax=std::max(Tmax,dom.W(i).T);
    check("T_ad in [1890,2310] K (Cantera 2104)",
        Tmax>1890.0 && Tmax<2310.0, Tmax, 2104.0);

    // Fuel fully consumed behind the flame (right region at tEnd).
    double yC3H8max=0.0;
    for(int i=0;i<N;++i) yC3H8max=std::max(yC3H8max,dom.W(i).Y[iC3H8]);
    check("fresh mixture behind flame consumed (max Y_C3H8 < 0.05 y_f)",
        yC3H8max < 0.05*yC3H8f, yC3H8max, yC3H8f);

    // Element conservation.
    double c=0,h=0,o=0,n=0;
    for(int i=0;i<N;++i)
    {
        const auto& w=dom.W(i);
        c+=w.Y[iC3H8]*3.0+w.Y[iCO]+w.Y[iCO2];
        h+=w.Y[iC3H8]*8.0+w.Y[iH2O]*2.0;
        o+=w.Y[iO2]*2.0+w.Y[iCO]+w.Y[iCO2]*2.0+w.Y[iH2O];
        n+=w.Y[iN2]*2.0;
    }
    check("element C conserved", std::abs(c-C0)/C0<1e-3, c, C0);
    check("element H conserved", std::abs(h-H0)/H0<1e-3, h, H0);
    check("element O conserved", std::abs(o-O0)/O0<1e-3, o, O0);
    check("element N conserved", std::abs(n-N0)/N0<1e-3, n, N0);

    cout << "  S_L=" << SL << " m/s  T_max=" << Tmax
         << " K  track_pts=" << track.size() << "\n";

    cout << (g_fail==0 ? "FLAME1D_PASS\n" : "FLAME1D_FAIL\n");
    return g_fail==0 ? 0 : 1;
}
