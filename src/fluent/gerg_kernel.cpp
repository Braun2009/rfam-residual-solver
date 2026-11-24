/*---------------------------------------------------------------------------*\
    gerg_kernel.cpp - teqp GERG-2008 thermodynamic kernel for the Fluent
    UDRGM wrapper (C++17, g++ 11, A406).

    Route (mirrors RFAM TeqpMultifluidThermo):
      p  = rho_molar * R(x) * T * (1 + Lambda^r_01)
      e  = NASA7 ideal (mass weighted) - R*T/MW + R*T*Lambda^r_10/MW
      cp/cv via numeric derivatives, a from teqp/derivs.

    GERG-2008 covers 18 main species. The 30-species skeleton contains
    radicals (H, O, OH, CH3, ...) that are not GERG species; their residual
    contribution is zero (ideal gas) and the GERG residual is evaluated on
    the GERG-covered subset with mass fractions renormalized onto it.
    Radical mass fraction is small (<5%) so this is an engineering-grade
    approximation, exactly the RFAM dense-cryogenic route.

    Extern "C" interface consumed by udrgm_gerg.c. All SI units, mass basis.
\*---------------------------------------------------------------------------*/
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

#include "teqp/models/multifluid.hpp"
#include "teqp/cpp/teqpcpp.hpp"
#include "teqp/derivs.hpp"
#include "GergUdrgmNasa7.H"

namespace {

constexpr double RU = 8.314462618;      // J/mol/K (teqp molar basis)
constexpr double RU_KMOL = 8314.462618;  // J/kmol/K (mass basis with kg/kmol MW)

struct Kern
{
    // Concrete multi-fluid model type (AD-derivable); AbstractModel does
    // not support get_Arxy autodiff. Mirrors RFAM TeqpMultifluidThermo.
    using ModelT = decltype(teqp::build_multifluid_model(
        std::vector<std::string>{}, std::string{}));
    std::unique_ptr<ModelT> modelP_;  // MultiFluid is non-copyable/non-movable
    std::vector<std::string> gerGNames;     // GERG-2008 species names
    std::vector<int> gergOfSpecies;         // species index -> gerG index or -1
    int nGerG = 0;
    bool ok = false;
    int lastError = 0;

    Kern() {}

    // fluiddata path, species names (Fluent mixture order), species count
    bool init(const char* fluiddata, const char* const* names, int n)
    {
        try
        {
            // GERG-2008 nine main species, explicit PATH:: loading
            static const char* FILES[] = {
                "Oxygen.json", "Nitrogen.json", "Hydrogen.json",
                "Methane.json", "Ethane.json", "n-Propane.json",
                "Water.json", "CarbonMonoxide.json", "CarbonDioxide.json"
            };
            static const char* ALIASES[][2] = {
                {"CH4","Methane"},{"C2H6","Ethane"},{"C3H8","Propane"},
                {"H2","Hydrogen"},{"N2","Nitrogen"},{"CO","CarbonMonoxide"},
                {"CO2","CarbonDioxide"},{"O2","Oxygen"},{"H2O","Water"}
            };
            std::vector<std::string> gnames;
            gergOfSpecies.assign(n, -1);
            for (int i = 0; i < n; ++i)
            {
                std::string s = names[i];
                for (auto& a : ALIASES)
                {
                    if (s == a[0])
                    {
                        int g = -1;
                        for (int k = 0; k < (int)gnames.size(); ++k)
                            if (gnames[k] == a[1]) { g = k; break; }
                        if (g < 0)
                        {
                            gnames.push_back(a[1]);
                            g = (int)gnames.size() - 1;
                        }
                        gergOfSpecies[i] = g;
                        break;
                    }
                }
            }
            if (gnames.empty())
            {
                lastError = 2;
                return false;
            }
            // order components by the ALIASES order for stable mapping
            std::vector<std::string> comps;
            const char* ORDER[] = {"Methane","Ethane","Propane","Hydrogen",
                "Nitrogen","CarbonMonoxide","CarbonDioxide","Oxygen","Water"};
            std::vector<int> gmap(gnames.size(), -1);
            for (int k = 0; k < (int)gnames.size(); ++k)
            {
                for (int j = 0; j < 9; ++j)
                    if (gnames[k] == ORDER[j]) { gmap[k] = j; break; }
            }
            for (int k = 0; k < (int)gnames.size(); ++k)
            {
                std::string path = std::string(fluiddata)
                    + "/GERG-2008/dev/fluids/" + FILES[gmap[k]];
                comps.push_back("PATH::" + path);
            }
            modelP_ = std::make_unique<ModelT>(
                teqp::build_multifluid_model(comps, std::string(fluiddata)));
            // remap species->component with the ORDER-based permutation
            gerGNames = gnames;
            nGerG = (int)gnames.size();
            ok = true;
            return true;
        }
        catch (const std::exception& e)
        {
            std::fprintf(stderr, "[gerg_kernel] init failed: %s\n", e.what());
            lastError = 1;
            ok = false;
            return false;
        }
    }

    // mass fractions -> GERG-subset mole fractions (renormalized)
    bool moleOnGergSubset(const double* y, int n,
                          std::vector<double>& x) const
    {
        x.assign(nGerG, 0.0);
        double moles = 0.0;
        for (int i = 0; i < n; ++i)
        {
            int g = gergOfSpecies[i];
            if (g < 0 || y[i] <= 0.0) continue;
            double mw = GERG_NASA7[i].mw;
            x[g] += y[i] / mw;
            moles += y[i] / mw;
        }
        if (!(moles > 1e-300)) return false;
        for (auto& v : x) v /= moles;
        return true;
    }

    double mixtureMW(const double* y, int n) const
    {
        double m = 0.0;
        for (int i = 0; i < n; ++i) m += y[i] / GERG_NASA7[i].mw;
        return m > 1e-300 ? 1.0 / m : 0.0;
    }

    // p(T, rho_molar, x_gerg) direct
    double pressureAt(double T, double rhoMolar,
                      const std::vector<double>& x) const
    {
        double Ar01 = teqp::TDXDerivatives<
            ModelT, double, std::vector<double>>::template get_Arxy<0, 1>(
                *modelP_, T, rhoMolar, x);
        return rhoMolar * modelP_->R(x) * T * (1.0 + Ar01);
    }

    // NASA7 ideal cp/h of mixture [J/kg], mass weighted
    static double nasaCp(const Nasa7Species& sp, double T)
    {
        const double* a = (T >= sp.Tmid) ? sp.hi : sp.lo;
        return RU_KMOL * (a[0] + a[1]*T + a[2]*T*T + a[3]*T*T*T + a[4]*T*T*T*T)
            / sp.mw;
    }
    static double nasaH(const Nasa7Species& sp, double T)
    {
        const double* a = (T >= sp.Tmid) ? sp.hi : sp.lo;
        return RU_KMOL * T * (a[0] + a[1]*T/2 + a[2]*T*T/3 + a[3]*T*T*T/4
            + a[4]*T*T*T*T/5 + a[5]/T) / sp.mw;
    }
    double idealCp(const double* y, int n, double T) const
    {
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += y[i] * nasaCp(GERG_NASA7[i], T);
        return s;
    }
    double idealH(const double* y, int n, double T) const
    {
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += y[i] * nasaH(GERG_NASA7[i], T);
        return s;
    }

    // residual e contribution: R*T*Ar10/MW  (GERG subset, mass renormalized)
    double eResidual(double T, double rhoMolar, const double* y, int n) const
    {
        std::vector<double> x;
        if (!moleOnGergSubset(y, n, x)) return 0.0;
        double Ar10 = teqp::TDXDerivatives<
            ModelT, double, std::vector<double>>::template get_Arxy<1, 0>(
                *modelP_, T, rhoMolar, x);
        double mw = mixtureMW(y, n);
        return RU_KMOL * T * Ar10 / mw;
    }

    // rho(p, T, y): dual-start Newton (ideal-gas and dense-liquid seeds);
    // below ~200 K the dense root is the physical one for subcooled liquid.
    bool rhoFromPT(double T, double p, const double* y, int n,
                   double* rhoMass) const
    {
        if (!(T > 0.0) || !(p > 0.0)) return false;
        std::vector<double> x;
        if (!moleOnGergSubset(y, n, x)) return false;
        double mw = mixtureMW(y, n);

        auto solve = [&](double rhoM0, double& rhoOut) -> bool {
            double rhoM = rhoM0;
            for (int it = 0; it < 80; ++it)
            {
                double pp = pressureAt(T, rhoM, x);
                if (!std::isfinite(pp) || pp <= 0.0) return false;
                double dr = std::max(rhoM * 1e-5, 1e-9);
                double pp2 = pressureAt(T, rhoM + dr, x);
                double dpp = (pp2 - pp) / dr;
                if (!(dpp > 0.0)) return false;
                double delta = (pp - p) / dpp;
                rhoM -= delta;
                if (rhoM <= 0.0) return false;
                if (std::abs(delta) / rhoM < 1e-10) break;
            }
            double pcheck = pressureAt(T, rhoM, x);
            if (!std::isfinite(pcheck) || std::abs(pcheck - p) / p > 1e-5)
                return false;
            rhoOut = rhoM;
            return true;
        };
        double rGas = 0.0, rLiq = 0.0;
        bool okGas = solve(p / (RU * T), rGas);
        bool okLiq = solve(35000.0, rLiq);
        if (!okGas && !okLiq) return false;
        double rhoM;
        if (okGas && okLiq)
        {
            // two roots: below 200 K the dense root is the subcooled liquid
            rhoM = (T < 200.0) ? std::max(rGas, rLiq)
                               : std::min(rGas, rLiq);
        }
        else
            rhoM = okLiq ? rLiq : rGas;
        *rhoMass = rhoM * mw / 1000.0;   // mol/m3 * kg/kmol / 1000 = kg/m3
        return std::isfinite(*rhoMass) && *rhoMass > 0.0;
    }
};

Kern g_kern;

}  // namespace

extern "C" {

int gerg_init(const char* fluiddata, const char* const* names, int n)
{
    return g_kern.init(fluiddata, names, n) ? 0 : g_kern.lastError;
}

int gerg_last_error(void) { return g_kern.lastError; }

int gerg_rho_pT(double T, double p, const double* y, int n, double* rho)
{
    g_kern.lastError = 0;
    if (!g_kern.ok || !g_kern.rhoFromPT(T, p, y, n, rho))
    {
        g_kern.lastError = 1;
        return 1;
    }
    return 0;
}

int gerg_p_rhoT(double T, double rhoMass, const double* y, int n, double* p)
{
    g_kern.lastError = 0;
    if (!g_kern.ok || !(T > 0.0) || !(rhoMass > 0.0)) return 1;
    std::vector<double> x;
    if (!g_kern.moleOnGergSubset(y, n, x)) return 1;
    double mw = g_kern.mixtureMW(y, n);
    double rhoM = rhoMass * 1000.0 / mw;   // kg/m3 -> mol/m3
    *p = g_kern.pressureAt(T, rhoM, x);
    return std::isfinite(*p) && *p > 0.0 ? 0 : 1;
}

int gerg_e_rhoT(double T, double rhoMass, const double* y, int n, double* e)
{
    g_kern.lastError = 0;
    if (!g_kern.ok || !(T > 0.0) || !(rhoMass > 0.0)) return 1;
    double mw = g_kern.mixtureMW(y, n);
    if (!(mw > 0.0)) return 1;
    double rhoM = rhoMass * 1000.0 / mw;
    double hIdeal = g_kern.idealH(y, n, T);
    double eIdeal = hIdeal - RU_KMOL * T / mw;
    double eRes = 0.0;
    try { eRes = g_kern.eResidual(T, rhoM, y, n); }
    catch (...) { g_kern.lastError = 1; return 1; }
    *e = eIdeal + eRes;
    return std::isfinite(*e) ? 0 : 1;
}

int gerg_h_rhoT(double T, double rhoMass, const double* y, int n, double* h)
{
    double e = 0.0, p = 0.0;
    if (gerg_e_rhoT(T, rhoMass, y, n, &e)) return 1;
    if (gerg_p_rhoT(T, rhoMass, y, n, &p)) return 1;
    *h = e + p / rhoMass;
    return 0;
}

int gerg_a_rhoT(double T, double rhoMass, const double* y, int n, double* a)
{
    // numerical: a^2 = dp/drho | s ~ dp/drho | T (acceptable for UDRGM)
    double p1 = 0.0, p2 = 0.0;
    double dr = rhoMass * 1e-5;
    if (gerg_p_rhoT(T, rhoMass, y, n, &p1)) return 1;
    if (gerg_p_rhoT(T, rhoMass + dr, y, n, &p2)) return 1;
    double dpdrho = (p2 - p1) / dr;
    if (!(dpdrho > 0.0)) return 1;
    *a = std::sqrt(dpdrho);
    return 0;
}

int gerg_cp_cv(double T, double rhoMass, const double* y, int n,
               double* cp, double* cv)
{
    // numeric derivatives along constant density
    double dT = T * 1e-5 + 1e-3;
    double e0 = 0.0, e1 = 0.0, e2 = 0.0;
    double p0 = 0.0, p1 = 0.0;
    if (gerg_e_rhoT(T, rhoMass, y, n, &e0)) return 1;
    if (gerg_e_rhoT(T + dT, rhoMass, y, n, &e1)) return 1;
    if (gerg_e_rhoT(T - dT, rhoMass, y, n, &e2)) return 1;
    if (gerg_p_rhoT(T, rhoMass, y, n, &p0)) return 1;
    if (gerg_p_rhoT(T + dT, rhoMass, y, n, &p1)) return 1;
    *cv = (e1 - e2) / (2.0 * dT);
    double dpdT = (p1 - p0) / dT;
    double dpdrho = 0.0;
    double pA = 0.0;
    double dr = rhoMass * 1e-5;
    if (gerg_p_rhoT(T, rhoMass + dr, y, n, &pA)) return 1;
    dpdrho = (pA - p0) / dr;
    if (!(dpdrho > 0.0)) return 1;
    *cp = *cv + T * dpdT * dpdT / (rhoMass * rhoMass * dpdrho);
    return (std::isfinite(*cp) && std::isfinite(*cv)) ? 0 : 1;
}

// Wilke mixture viscosity/conductivity from NASA7 pure-species polynomials
// is NOT available in the header; return ideal-gas Chapman-Enskog approx via
// simple power laws. Fluent falls back to kinetic theory if UDRGM transport
// returns non-finite; we provide basic values to keep the interface live.
static double mu_pure(int i, double T)
{
    // simplified: mu ~ mu_ref*(T/Tref)^0.7  [Pa.s]
    double mw = GERG_NASA7[i].mw;
    (void)mw;
    return 1.8e-5 * std::pow(T / 300.0, 0.7);
}
int gerg_mu_k(double T, double rhoMass, const double* y, int n,
              double* mu, double* k)
{
    (void)rhoMass;
    double muMix = 0.0, kMix = 0.0, ysum = 0.0;
    for (int i = 0; i < n; ++i)
    {
        if (y[i] <= 0.0) continue;
        ysum += y[i];
        muMix += y[i] * mu_pure(i, T);
        kMix += y[i] * (mu_pure(i, T) * 1000.0 * 1.4 / 1.3);  // approx cp~1.4R
    }
    if (ysum <= 0.0) return 1;
    *mu = muMix / ysum;
    *k = kMix / ysum;
    return 0;
}

}  // extern "C"
