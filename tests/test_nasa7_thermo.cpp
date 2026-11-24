/*---------------------------------------------------------------------------*\
    test_nasa7_thermo - consistency test for NASA-7 ideal-gas thermodynamics
    (Stage 3, step 1). This is the acceptance gate that guards against the
    N2-JANAF class of coefficient / range-boundary discontinuity bugs
    [[n2-thermo-discontinuity]] before NASA-7 is trusted in the reacting flow.

    Checks:
      1. Cp(298.15 K) per species vs textbook values (+/- 2%).
      2. hsens(298.15 K) == 0 per species (sensible reference self-consistent).
      3. Range-boundary (Tcommon = 1000 K) continuity of Cp and h (< 0.1%) —
         this is the direct guard against a JANAF-style enthalpy jump.
      4. Cp > 0 and enthalpy strictly monotone increasing over 200..4000 K.
      5. Mixture weighting: pure species == that species; 50/50 == mean.

    Pure std C++. Build: g++ -O2 -std=c++17 -Isrc -Isrc/thermo test_nasa7_thermo.cpp
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include "src/thermo/NASA7Thermo.H"
#include "src/thermo/SpeciesRegistry.H"

using namespace Foam;
using namespace std;

static int g_pass = 0, g_fail = 0;
static void check(const string& n, bool ok, double got = 0, double exp = 0)
{
    cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << n;
    if (!ok) cout << "   (got " << got << ", expected " << exp << ")";
    cout << "\n";
    ok ? ++g_pass : ++g_fail;
}

int main()
{
    cout << "=== NASA-7 ideal-gas thermodynamics consistency test ===\n\n";
    const NASA7Thermo& th = nasa7Thermo();
    SpeciesRegistry reg;

    // ---- 1. Cp(298.15) mass-specific [J/(kg K)] vs textbook (+/- 2%) --------
    cout << "-- Cp(298.15 K) vs textbook [J/(kg K)] --\n";
    struct CpRef { const char* name; double cp; };
    // standard 298 K ideal-gas Cp (mass basis)
    const vector<CpRef> refs = {
        {"O2", 918.0}, {"N2", 1040.0}, {"H2", 14310.0}, {"CH4", 2220.0},
        {"C2H6", 1750.0}, {"C3H8", 1670.0}, {"H2O", 1864.0}, {"CO", 1040.0}, {"CO2", 844.0}
    };
    for (const auto& r : refs)
    {
        double cp = th.cp_mole(r.name, 298.15) / th.MW(r.name);   // -> J/(kg K)
        double relerr = fabs(cp - r.cp) / r.cp;
        cout << "    " << setw(5) << left << r.name << " Cp=" << fixed << setprecision(1)
             << setw(9) << right << cp << "  ref=" << setw(8) << r.cp
             << "  err=" << setprecision(2) << 100*relerr << "%\n";
        check(string("Cp(298) ") + r.name + " within 2%", relerr < 0.02, cp, r.cp);
    }

    // ---- 2. hsens(298.15) == 0 per species ---------------------------------
    cout << "\n-- hsens(298.15 K) == 0 --\n";
    for (const auto& r : refs)
    {
        double h = th.hsens_mole(r.name, 298.15);
        check(string("hsens(298) ") + r.name + " ~ 0", fabs(h) < 1e-6, h, 0.0);
    }

    // ---- 3. Range-boundary continuity at Tcommon (Cp and h) ----------------
    cout << "\n-- Continuity at range boundary (JANAF-jump guard) --\n";
    for (const auto& r : refs)
    {
        // approach 1000 K from both sides
        double Tm = 999.999, Tp = 1000.001;
        double cpm = th.cp_mole(r.name, Tm), cpp = th.cp_mole(r.name, Tp);
        double hm  = th.h_mole(r.name, Tm),  hp  = th.h_mole(r.name, Tp);
        double dcp = fabs(cpp - cpm) / fabs(cpm);
        double dh  = fabs(hp - hm) / (fabs(hm) + 1.0);
        check(string("Cp continuous @1000K ") + r.name, dcp < 1e-3, dcp, 0.0);
        check(string("h  continuous @1000K ") + r.name, dh  < 1e-3, dh, 0.0);
    }

    // ---- 4. Cp > 0 and enthalpy monotone over 200..4000 K ------------------
    cout << "\n-- Cp>0 and h monotone (200..4000 K) --\n";
    for (const auto& r : refs)
    {
        bool cpPos = true, hMono = true;
        double hPrev = th.h_mole(r.name, 200.0);
        for (double T = 210.0; T <= 4000.0; T += 10.0)
        {
            if (th.cp_mole(r.name, T) <= 0.0) cpPos = false;
            double h = th.h_mole(r.name, T);
            if (h <= hPrev) hMono = false;
            hPrev = h;
        }
        check(string("Cp>0 ") + r.name, cpPos);
        check(string("h monotone ") + r.name, hMono);
    }

    // ---- 5. Mixture weighting ----------------------------------------------
    cout << "\n-- Mixture weighting --\n";
    vector<string> names = reg.names();
    int nS = reg.size();
    // pure CO2
    {
        vector<double> Y(nS, 0.0); Y[reg.index("CO2")] = 1.0;
        double cpMix = th.cp_mass_mix(names, Y, 1500.0);
        double cpPure = th.cp_mole("CO2", 1500.0) / th.MW("CO2");
        check("pure-CO2 mixture == CO2", fabs(cpMix - cpPure) < 1e-6*cpPure, cpMix, cpPure);
    }
    // 50/50 O2/N2 by mass
    {
        vector<double> Y(nS, 0.0); Y[reg.index("O2")] = 0.5; Y[reg.index("N2")] = 0.5;
        double cpMix = th.cp_mass_mix(names, Y, 800.0);
        double cpO2 = th.cp_mole("O2", 800.0) / th.MW("O2");
        double cpN2 = th.cp_mole("N2", 800.0) / th.MW("N2");
        double cpAvg = 0.5*(cpO2 + cpN2);
        check("50/50 O2/N2 == mass-mean", fabs(cpMix - cpAvg) < 1e-6*cpAvg, cpMix, cpAvg);
    }
    // hsens(298) of any mixture == 0
    {
        vector<double> Y(nS, 0.0); Y[reg.index("CO2")] = 0.4; Y[reg.index("H2O")] = 0.6;
        double h = th.hsens_mass_mix(names, Y, 298.15);
        check("mixture hsens(298) ~ 0", fabs(h) < 1e-6, h, 0.0);
    }

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail == 0 ? "NASA7_THERMO_PASS" : "NASA7_THERMO_FAIL") << "\n";
    return g_fail == 0 ? 0 : 1;
}
