/*---------------------------------------------------------------------------*\
    test_phase7_chem0d - Phase 7 core: finite-rate chemistry in global-space
    via SpeciesMapper + CVODE, validated as a 0-D reactor.

    Goal: prove the chem<->global bridge + explicit-qdot energy bookkeeping is
    correct, i.e. the SAME ignition as Phase 6 (0-D CVODE) but now expressed in
    SpeciesRegistry global-space, with element + energy conservation.

    Placeholder thermo: chem-space temperature uses ideal Cv; CFD-side T (EOS
    reback) not exercised here. Asserts are structural/conservation only.
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "src/thermo/SpeciesRegistry.H"
#include "src/combustion/FiniteRateReactor.H"
#include "src/combustion/FiniteRate/TwoStepCO.H"

using namespace Foam;
using namespace std;

int main()
{
    cout << "=== Phase 7 core: finite-rate chemistry 0-D (global-space + mapper) ===\n\n";
    int pass=0, fail=0;
    auto check=[&](const string& n, bool ok){ cout<<"  ["<<(ok?"PASS":"FAIL")<<"] "<<n<<"\n"; ok?++pass:++fail; };

    SpeciesRegistry reg;
    TwoStepCO chem(1e10, 0.0, 125520.0, 2e12, 0.0, 167360.0);   // == Phase 6 rates
    FiniteRateReactor reactor(
        &chem,reg,FiniteRateReactor::LegacyFixedCvDiagnostic(1500.0));

    const int nG = reg.size();
    const double rho=15.0, p=4e6;
    double T=1200.0;
    // global Y0: stoichiometric C3H8/O2 (Phase 6: chem Y0={0.216 C3H8, 0.784 O2})
    vector<double> Yg(nG, 0.0);
    Yg[reg.index("C3H8")] = 0.216;
    Yg[reg.index("O2")]   = 0.784;

    // element atoms (per m3) helper
    auto atoms=[&](const vector<double>& Y,double rr,double&C,double&H,double&O,double&N){
        C=H=O=N=0; for(int g=0;g<reg.size();++g){double mr=rr*Y[g]/reg.molecularWeight(g);
            C+=mr*reg.nC(g);H+=mr*reg.nH(g);O+=mr*reg.nO(g);N+=mr*reg.nN(g);} };
    double C0,H0,O0,N0; atoms(Yg,rho,C0,H0,O0,N0);

    cout << "  global species (" << nG << "): "; for(int g=0;g<nG;++g) cout<<reg.name(g)<<" "; cout<<"\n";
    cout << "  init: Y[C3H8]="<<Yg[reg.index("C3H8")]<<" Y[O2]="<<Yg[reg.index("O2")]<<" T="<<T<<"K\n\n";

    // integrate 0..2e-4 s (== Phase 6), explicit qdot bookkeeping
    const double dt=2e-6, t_end=2e-4;
    double rhoE_chem = 0.0;   // accumulated chemical heat [J/m3]
    bool dE_nonneg = true;
    int nsteps = (int)(t_end/dt + 0.5);
    cout << "  t(us)   T_chem(K)  Y_C3H8   Y_CO2    Y_H2O\n";
    for (int s=0; s<nsteps; ++s)
    {
        double dE;
        reactor.reactLegacyFixedCv(rho, p, T, Yg, dt, dE);
        rhoE_chem += dE;
        if (dE < -1e-3) dE_nonneg = false;     // exothermic: dE >= 0
        if (s % 20 == 0 || s==nsteps-1)
            cout << "  " << fixed << setprecision(1) << setw(5) << (s+1)*dt*1e6
                 << "   " << setprecision(0) << setw(6) << T
                 << "   " << setprecision(4) << setw(6) << Yg[reg.index("C3H8")]
                 << "   " << setw(6) << Yg[reg.index("CO2")]
                 << "   " << setw(6) << Yg[reg.index("H2O")] << "\n";
    }

    double C1,H1,O1,N1; atoms(Yg,rho,C1,H1,O1,N1);
    double thr = (C0+H0+O0+N0);
    cout << "\n";

    // ---- assertions ----
    check("fuel consumed (Y_C3H8 -> ~0)", Yg[reg.index("C3H8")] < 1e-3);
    // With reversible R2 (CO2 ⇌ CO + 0.5O2), high-T equilibrium shifts CO2→CO.
    // Check total carbon products (CO+CO2) rather than CO2 alone.
    double Y_CO_CO2 = Yg[reg.index("CO")] + Yg[reg.index("CO2")];
    check("carbon products produced (CO+CO2 > 0.3)", Y_CO_CO2 > 0.3);
    check("H2O produced", Yg[reg.index("H2O")] > 0.1);
    check("non-chem species stay zero (e.g. CH4,H2)",
          fabs(Yg[reg.index("CH4")])<1e-15 && fabs(Yg[reg.index("H2")])<1e-15);
    cout << scientific << setprecision(2);
    check("element C conserved", fabs(C1-C0) < 1e-9*thr);
    check("element H conserved", fabs(H1-H0) < 1e-9*thr);
    check("element O conserved", fabs(O1-O0) < 1e-9*thr);
    check("exothermic (all dE>=0)", dE_nonneg);
    check("chemical heat released (rhoE_chem>0)", rhoE_chem > 0.0);
    cout << "  accumulated chemical heat rhoE_chem = " << rhoE_chem << " J/m3\n";

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << pass << " passed, " << fail << " failed\n";
    cout << (fail==0 ? "PHASE7_CHEM0D_PASS" : "PHASE7_CHEM0D_FAIL") << "\n";
    return fail==0 ? 0 : 1;
}
