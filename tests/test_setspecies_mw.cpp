/*---------------------------------------------------------------------------*\
    test_setspecies_mw - species-mapping smoke test (Phase 7 Step 2).

    SMALL & FATAL. Before wiring chemistry into the 1D framework, prove that
    RealFluidState::setSpecies routes the caller's EXACT species order down to
    the EOS mixing rule:
      - a pure-N2 mixture in a 6-species vector must report MW(N2)=0.0280,
        NOT the legacy default {O2,C3H8} answer (0.044) that the old code gave
        when it silently looped only the first two slots.
      - an unknown species must THROW (no silent fallback).
    Far easier to localise a mixing/ordering bug here than inside a 1D run.

    MW values are the SpeciesCriticalDB (createDefault) values the mixing rule
    actually uses. Placeholder thermo elsewhere is irrelevant: MW is a pure
    mole-weighted sum over the resolved species, independent of a_mix/b_mix.
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include "RealFluidState.H"

using namespace Foam;
using namespace std;

int main()
{
    cout << "=== Phase 7 Step 2: setSpecies MW-routing smoke test ===\n\n";
    int pass=0, fail=0;
    auto check=[&](const string& n, bool ok){ cout<<"  ["<<(ok?"PASS":"FAIL")<<"] "<<n<<"\n"; ok?++pass:++fail; };
    auto rel =[&](double a,double b){ return fabs(a-b) < 1e-4*fabs(b); };

    // chem-space order — deliberately NOT the default {O2,C3H8}.
    const vector<string> sp = {"C3H8","O2","CO","CO2","H2O","N2"};
    const int n = (int)sp.size();

    auto mwOfPure=[&](const string& name)->double{
        RealFluidState s; s.setSpecies(sp);
        vector<double> Y(n,0.0);
        int idx=-1; for(int i=0;i<n;++i) if(sp[i]==name) idx=i;
        Y[idx]=1.0;
        s.updatePTY(4.0e6, 1200.0, Y);
        return s.MW();
    };

    double mwN2   = mwOfPure("N2");
    double mwC3H8 = mwOfPure("C3H8");
    double mwO2   = mwOfPure("O2");
    cout << fixed << setprecision(6);
    cout << "  pure N2   -> MW = " << mwN2   << "  (expect 0.028014)\n";
    cout << "  pure C3H8 -> MW = " << mwC3H8 << "  (expect 0.044096)\n";
    cout << "  pure O2   -> MW = " << mwO2   << "  (expect 0.031999)\n\n";

    check("pure N2 MW = 0.028014 (NOT default C3H8 0.044)", rel(mwN2, 0.028014));
    check("pure C3H8 MW = 0.044096", rel(mwC3H8, 0.044096));
    check("pure O2 MW = 0.031999", rel(mwO2, 0.031999));

    // stoichiometric C3H8/O2 (== 0-D test): mole-weighted MW must sit strictly
    // between O2 and C3H8 (O2 dominates by moles -> closer to O2).
    {
        RealFluidState s; s.setSpecies(sp);
        vector<double> Y(n,0.0); Y[0]=0.216; Y[1]=0.784;   // C3H8, O2
        s.updatePTY(4.0e6, 1200.0, Y);
        double mw = s.MW();
        cout << "  stoich C3H8/O2 -> MW = " << mw << "  (expect O2<MW<C3H8, ~0.0340)\n\n";
        check("stoich mixture MW between O2 and C3H8", mw > mwO2 && mw < mwC3H8);
    }

    // Explicit failure: unknown species must throw, never silently fall back.
    bool threw=false;
    try { RealFluidState s; s.setSpecies({"O2","NOSUCH"}); }
    catch (const std::exception&) { threw=true; }
    check("unknown species throws (no silent fallback)", threw);

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << pass << " passed, " << fail << " failed\n";
    cout << (fail==0 ? "SETSPECIES_MW_PASS" : "SETSPECIES_MW_FAIL") << "\n";
    return fail==0 ? 0 : 1;
}
