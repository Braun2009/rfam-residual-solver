/*---------------------------------------------------------------------------*\
    test_pr_mixing - Stage 3 step 3: true quadratic VdW mixing validation.

    Checks (PR_MIXING_PASS sentinel):
      1. Pure species: a_quad(Y_pure_i) == a_i(T),  b_lin(Y_pure_i) == b_i
      2. Mass->mole conversion sanity (pure and binary)
      3. O2/C3H8 binary: rho > 0, Z in (0,2], Cp/Cv/a^2 > 0  (subcritical gas)
      4. 9-species combustion products: rho/Z/Cp/Cv/a^2 all positive, physical
      5. Legacy opt-out: RealFluidState without setSpecies still gives old path
      6. Analytic da_mix/dT and d2a_mix/dT2 match five-point differences
      7. Departure functions remain finite and sign-consistent

    Build (standalone, no SUNDIALS needed):
      g++ -O2 -std=c++17 -Isrc -Isrc/thermo -Isrc/thermo/EquationOfState \
          -Isrc/thermo/MixingRules -Isrc/thermo/RootSelection \
          test_pr_mixing.cpp \
          src/thermo/RealFluidState.C \
          src/thermo/EquationOfState/PengRobinsonEOS.C \
          src/thermo/MixingRules/VdWOneFluidMixing.C \
          src/thermo/MixingRules/SpeciesCriticalDB.C \
          src/thermo/MixingRules/BinaryInteractionDB.C \
          src/thermo/RootSelection/CubicRoots.C
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

#include "src/thermo/SpeciesRegistry.H"
#include "src/thermo/RealFluidState.H"
#include "src/thermo/MixingRules/VdWOneFluidMixing.H"
#include "src/thermo/MixingRules/SpeciesCriticalDB.H"
#include "src/thermo/MixingRules/BinaryInteractionDB.H"
#include "src/thermo/EquationOfState/PengRobinsonEOS.H"

using namespace Foam;
using namespace std;

static int g_pass = 0, g_fail = 0;
static void check(const string& name, bool ok, double got = 0, double ref = 0)
{
    cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
    if (!ok) cout << "   (got=" << got << " ref=" << ref << ")";
    cout << "\n";
    ok ? ++g_pass : ++g_fail;
}

// Compute per-species a_i(T) and b_i the same way VdWOneFluidMixing helper
// does so we have a reference for the "pure species" regression check.
static double pureA(const CriticalProperties& p, double T)
{
    constexpr double R = 8.314462618;
    double kappa = 0.37464 + 1.54226*p.omega - 0.26992*p.omega*p.omega;
    double alpha = std::pow(1.0 + kappa*(1.0 - std::sqrt(T/p.Tc)), 2);
    return 0.45724 * R*R * p.Tc*p.Tc / p.Pc * alpha;
}
static double pureB(const CriticalProperties& p)
{
    constexpr double R = 8.314462618;
    return 0.07780 * R * p.Tc / p.Pc;
}

int main()
{
    cout << "=== PR true quadratic VdW mixing test (Stage 3 step 3) ===\n\n";

    SpeciesCriticalDB  critDB  = SpeciesCriticalDB::createDefault();
    BinaryInteractionDB kijDB  = BinaryInteractionDB::createDefault();
    SpeciesRegistry     reg;

    const vector<string>& allNames = reg.names();
    const int nS = reg.size();

    // ------------------------------------------------------------------ 1/2
    cout << "-- 1. Pure-species a_quad == a_i(T), b_lin == b_i --\n";
    const double T_test = 1000.0;
    const vector<string> testSpecies = {"O2","N2","H2","CH4","C3H8","H2O","CO","CO2"};
    for (const auto& sp : testSpecies)
    {
        // Build a 1-component EOS
        VdWOneFluidMixing mix(critDB, kijDB);
        mix.setSpeciesNames({sp});
        mix.enableMultiSpeciesMixing();

        vector<double> Y1 = {1.0};
        double aGot = mix.a_quad(Y1, T_test);
        double bGot = mix.b_lin(Y1);

        const CriticalProperties& cp = critDB.get(sp);
        double aRef = pureA(cp, T_test);
        double bRef = pureB(cp);

        double aErr = std::fabs(aGot - aRef) / (std::fabs(aRef) + 1e-30);
        double bErr = std::fabs(bGot - bRef) / (std::fabs(bRef) + 1e-30);

        check("a_quad pure " + sp + " (rel<1e-10)", aErr < 1e-10, aGot, aRef);
        check("b_lin  pure " + sp + " (rel<1e-10)", bErr < 1e-10, bGot, bRef);
    }

    // ------------------------------------------------------------------ 3
    cout << "\n-- 2. Mass->mole roundtrip (50/50 O2/N2) --\n";
    {
        VdWOneFluidMixing mix(critDB, kijDB);
        mix.setSpeciesNames({"O2","N2"});
        mix.enableMultiSpeciesMixing();

        // Y mass = {0.5, 0.5} -> x_O2 = (0.5/0.032)/(0.5/0.032 + 0.5/0.028)
        const double MWO2 = critDB.get("O2").MW, MWN2 = critDB.get("N2").MW;
        double nO2 = 0.5/MWO2, nN2 = 0.5/MWN2, nTot = nO2+nN2;
        double xO2 = nO2/nTot, xN2 = nN2/nTot;

        double aRef = xO2*xO2*pureA(critDB.get("O2"),T_test)
                    + 2.0*xO2*xN2*std::sqrt(pureA(critDB.get("O2"),T_test)*pureA(critDB.get("N2"),T_test))*1.0 // kij=0 for O2/N2
                    + xN2*xN2*pureA(critDB.get("N2"),T_test);
        double bRef = xO2*pureB(critDB.get("O2")) + xN2*pureB(critDB.get("N2"));

        vector<double> Y2 = {0.5, 0.5};
        double aGot = mix.a_quad(Y2, T_test);
        double bGot = mix.b_lin(Y2);
        check("a_quad 50/50 O2/N2 vs hand (rel<1e-8)",
              std::fabs(aGot-aRef)/(std::fabs(aRef)+1e-30) < 1e-8, aGot, aRef);
        check("b_lin  50/50 O2/N2 vs hand (rel<1e-8)",
              std::fabs(bGot-bRef)/(std::fabs(bRef)+1e-30) < 1e-8, bGot, bRef);
    }

    // ------------------------------------------------------------------ 4
    cout << "\n-- 3. O2/C3H8 subcritical gas: rho/Z/Cp/Cv/a^2 positive (opt-in) --\n";
    {
        const double p0 = 4.0e6, T0 = 400.0;
        vector<string> sp2 = {"O2","C3H8"};
        vector<double> Y2  = {0.784, 0.216};

        RealFluidState s;
        s.setSpecies(sp2);
        auto status = s.updatePTY(p0, T0, Y2);

        check("O2/C3H8 updatePTY success",
              status == RealFluidState::Status::SUCCESS);
        check("O2/C3H8 rho > 0",  s.rho() > 0.0,  s.rho());
        check("O2/C3H8 Z in (0,2]", s.Z() > 0.0 && s.Z() <= 2.0, s.Z());
        check("O2/C3H8 Cp > 0",   s.Cp() > 0.0,   s.Cp());
        check("O2/C3H8 Cv > 0",   s.Cv() > 0.0,   s.Cv());
        check("O2/C3H8 a > 0",    s.a() > 0.0,    s.a());
        cout << "    p=" << p0 << " T=" << T0
             << " rho=" << fixed << setprecision(3) << s.rho()
             << " Z=" << s.Z() << " Cp=" << s.Cp() << " a=" << s.a() << "\n";
    }

    // ------------------------------------------------------------------ 5
    cout << "\n-- 4. 9-species combustion products: all physical (opt-in) --\n";
    {
        // Approximate post-combustion products from stoich C3H8/O2 reaction
        const double p0 = 4.0e6, T0 = 3000.0;
        vector<double> Yprod(nS, 0.0);
        Yprod[reg.index("CO2")] = 0.35;
        Yprod[reg.index("H2O")] = 0.20;
        Yprod[reg.index("CO")]  = 0.15;
        Yprod[reg.index("O2")]  = 0.10;
        Yprod[reg.index("N2")]  = 0.20;

        RealFluidState s;
        s.setSpecies(allNames);
        auto status = s.updatePTY(p0, T0, Yprod);

        check("combustion products updatePTY success",
              status == RealFluidState::Status::SUCCESS);
        check("combustion products rho > 0",   s.rho() > 0.0,  s.rho());
        check("combustion products Z in (0,3]", s.Z() > 0.0 && s.Z() <= 3.0, s.Z());
        check("combustion products Cp > 0",     s.Cp() > 0.0,  s.Cp());
        check("combustion products Cv > 0",     s.Cv() > 0.0,  s.Cv());
        check("combustion products a > 0",      s.a()  > 0.0,  s.a());
        cout << "    p=" << p0 << " T=" << T0
             << " rho=" << fixed << setprecision(3) << s.rho()
             << " Z=" << s.Z() << " Cp=" << setprecision(1) << s.Cp()
             << " a=" << s.a() << "\n";
    }

    // ------------------------------------------------------------------ 6
    cout << "\n-- 5. Legacy path unchanged: RealFluidState default (no setSpecies) --\n";
    {
        // The default 2-species {O2,C3H8} pseudo-critical path must still work.
        RealFluidState sLeg;
        vector<double> Y2 = {0.784, 0.216};
        auto status = sLeg.updatePTY(4.0e6, 400.0, Y2);
        check("legacy default updatePTY success",
              status == RealFluidState::Status::SUCCESS);
        check("legacy rho > 0", sLeg.rho() > 0.0, sLeg.rho());
        check("legacy a > 0",   sLeg.a()   > 0.0, sLeg.a());
    }

    // ------------------------------------------------------------------ 7
    cout << "\n-- 6. a_mix temperature derivatives match finite differences --\n";
    {
        VdWOneFluidMixing mix(critDB, kijDB);
        mix.setSpeciesNames(allNames);
        mix.enableMultiSpeciesMixing();

        vector<double> Yprod(nS, 0.0);
        Yprod[reg.index("CO2")] = 0.35;
        Yprod[reg.index("H2O")] = 0.20;
        Yprod[reg.index("N2")]  = 0.45;

        for (double T : {300.0, 1000.0, 3000.0})
        {
            const double dT = 1.0e-3*T;
            const double am2 = mix.a_quad(Yprod, T - 2.0*dT);
            const double am1 = mix.a_quad(Yprod, T - dT);
            const double a0 = mix.a_quad(Yprod, T);
            const double ap1 = mix.a_quad(Yprod, T + dT);
            const double ap2 = mix.a_quad(Yprod, T + 2.0*dT);
            const double daNumeric =
                (am2 - 8.0*am1 + 8.0*ap1 - ap2)/(12.0*dT);
            const double d2aNumeric =
                (-ap2 + 16.0*ap1 - 30.0*a0 + 16.0*am1 - am2)
                /(12.0*dT*dT);
            const double da = mix.dadT_quad(Yprod, T);
            const double d2a = mix.d2adT2_quad(Yprod, T);
            const double daRel = std::abs(da-daNumeric)
                / std::max(std::abs(daNumeric), 1.0e-12);
            const double d2aRel = std::abs(d2a-d2aNumeric)
                / std::max(std::abs(d2aNumeric), 1.0e-12);
            check("da_mix/dT finite-difference at T=" + to_string((int)T) + "K",
                  daRel < 1.0e-5, da, daNumeric);
            check("d2a_mix/dT2 finite-difference at T=" + to_string((int)T) + "K",
                  d2aRel < 2.0e-3, d2a, d2aNumeric);
        }
    }

    // ------------------------------------------------------------------ 8
    cout << "\n-- 7. updateRhoEY reback with true mixing: T finite and positive --\n";
    {
        const double p0 = 4.0e6, T0 = 2000.0;
        vector<double> Yprod(nS, 0.0);
        Yprod[reg.index("CO2")] = 0.35;
        Yprod[reg.index("H2O")] = 0.20;
        Yprod[reg.index("CO")]  = 0.15;
        Yprod[reg.index("O2")]  = 0.10;
        Yprod[reg.index("N2")]  = 0.20;

        RealFluidState s;
        s.setSpecies(allNames);
        s.updatePTY(p0, T0, Yprod);
        double rho0 = s.rho(), e0 = s.e();

        RealFluidState s2;
        s2.setSpecies(allNames);
        auto status = s2.updateRhoEY(rho0, e0, Yprod, T0);
        double Tback = s2.T(), pback = s2.p();

        check("rhoEY reback status SUCCESS",
              status == RealFluidState::Status::SUCCESS);
        check("rhoEY T roundtrip within 1K",
              std::fabs(Tback - T0) < 1.0, Tback, T0);
        check("rhoEY p roundtrip within 0.1%",
              std::fabs(pback - p0) / p0 < 1e-3, pback, p0);
        cout << "    T_in=" << T0 << " T_back=" << fixed << setprecision(2)
             << Tback << "  p_in=" << p0 << " p_back=" << setprecision(0) << pback << "\n";
    }

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail == 0 ? "PR_MIXING_PASS" : "PR_MIXING_FAIL") << "\n";
    return g_fail == 0 ? 0 : 1;
}
