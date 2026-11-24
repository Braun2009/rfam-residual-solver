#ifndef TeqpMultifluidThermo_H
#define TeqpMultifluidThermo_H

/*---------------------------------------------------------------------------*\
    TeqpMultifluidThermo - GERG-2008 multifluid mixture thermodynamics
    via the teqp C++ headers (A406-only, RFAM_HAVE_TEQP).

    2026-08-08: mixture single-phase recovery currently falls back to the
    classical PR EOS whose isochoric energy bracket can fail in the flame-
    product regime (ENERGY_NOT_BRACKETED at 1213 sector steps).  GERG-2008
    is itself a MULTI-COMPONENT Helmholtz-explicit EOS (reducing functions
    + binary interaction parameters) and is more accurate than PR for
    natural-gas / light-hydrocarbon low-temperature high-pressure
    mixtures.  This class provides:

      p(rho,T,Ymass) = rho_molar * R * T * (1 + Lambda^r_01)
      e_sens(T,rho,Ymass) = NASA7 ideal (RFAM) + R*T*Lambda^r_10/MW
      recoverRhoEY(rho,e,Ymass): isochoric energy bisection on the stable
        single-phase branch (the flame-product regime where PR fails).

    RFAM stores composition as mass fractions; teqp multifluid APIs require
    mole fractions.  This adapter is therefore mass-basis at its public
    boundary and performs the conversion exactly once for every teqp call.
    A pure-fluid check cannot detect this distinction, so mixed combustion
    products are covered by the A406-only regression.

    teqp requires C++20 <concepts> (g++ 10+); the local WSL g++ 9.4 and
    the default build do NOT define RFAM_HAVE_TEQP and keep the PR path.
\*---------------------------------------------------------------------------*/

#ifdef RFAM_HAVE_TEQP

#include "teqp/models/multifluid.hpp"
#include "teqp/derivs.hpp"

#include "NASA7Thermo.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace Foam
{

class TeqpMultifluidThermo
{
public:
    static constexpr double R_UNIVERSAL = 8.314462618;
    static constexpr double T_MIN = 60.0;
    static constexpr double T_MAX = 5000.0;

    // 9-species registry order (O2, N2, H2, CH4, C2H6, C3H8, H2O, CO, CO2)
    static const std::vector<std::string>& speciesNames()
    {
        static const std::vector<std::string> names{
            "O2", "N2", "H2", "CH4", "C2H6", "C3H8", "H2O", "CO", "CO2"};
        return names;
    }

    static const std::vector<double>& molecularWeights()
    {
        static const std::vector<double> mw{
            0.0319988, 0.0280134, 0.00201588, 0.01604246, 0.03006904,
            0.04409562, 0.01801528, 0.0280101, 0.0440095};
        return mw;
    }

    using ModelT = decltype(teqp::build_multifluid_model(
        std::vector<std::string>{}, std::string{}));

    explicit TeqpMultifluidThermo(const std::string& fluiddataRoot)
    :
        model_(buildNineComponentModel(fluiddataRoot))
    {
    }

    // The 9-component alias names (O2/N2/H2/...) are not all present in
    // the teqp alias map (built from INFO NAME/CAS/REFPROP_NAME), so load
    // each component JSON by explicit PATH:: path.
    static ModelT buildNineComponentModel(const std::string& root)
    {
        static const char* files[9] = {
            "Oxygen.json", "Nitrogen.json", "Hydrogen.json", "Methane.json",
            "Ethane.json", "n-Propane.json", "Water.json",
            "CarbonMonoxide.json", "CarbonDioxide.json"};
        std::vector<std::string> components;
        for (const char* f : files)
            components.push_back("PATH::" + root + "/GERG-2008/dev/fluids/" + f);
        // BIP/departure collections from the fluiddata tree defaults
        // (mixture_binary_pairs.json: O2/C3H8 GERG interaction present;
        // select_identifier requires non-empty pairs).
        return teqp::build_multifluid_model(components, root);
    }

    const ModelT& model() const { return model_; }

    // Convert RFAM mass fractions to the mole fractions required by teqp.
    // A normalized input is expected from RealFluidState; accepting a common
    // scale factor here keeps this small adapter robust for standalone probes.
    static bool massToMoleFractions(
        const std::vector<double>& massFractions,
        std::vector<double>& moleFractions,
        double* molecularWeightOut = nullptr)
    {
        const std::vector<double>& mw = molecularWeights();
        if (massFractions.size() != mw.size()) return false;

        double massSum = 0.0;
        double molesPerInputMass = 0.0;
        for (std::size_t i = 0; i < massFractions.size(); ++i)
        {
            const double yi = massFractions[i];
            if (!std::isfinite(yi) || yi < 0.0
                || !std::isfinite(mw[i]) || mw[i] <= 0.0)
                return false;
            massSum += yi;
            molesPerInputMass += yi/mw[i];
        }
        if (!std::isfinite(massSum) || !std::isfinite(molesPerInputMass)
            || massSum <= 0.0 || molesPerInputMass <= 0.0)
            return false;

        moleFractions.resize(massFractions.size());
        for (std::size_t i = 0; i < massFractions.size(); ++i)
            moleFractions[i] = (massFractions[i]/mw[i])
                /molesPerInputMass;
        if (molecularWeightOut)
            *molecularWeightOut = massSum/molesPerInputMass;
        return std::all_of(
            moleFractions.begin(), moleFractions.end(),
            [](double xi) { return std::isfinite(xi) && xi >= 0.0; });
    }

private:
    double pressureForMoleFractions(
        double T, double rhoMolar,
        const std::vector<double>& moleFractions) const
    {
        const double Ar01 = teqp::TDXDerivatives<
            ModelT, double, std::vector<double>>::template get_Arxy<0, 1>(
                model_, T, rhoMolar, moleFractions);
        return rhoMolar*model_.R(moleFractions)*T*(1.0 + Ar01);
    }

    double eSensForMoleFractions(
        double T, double rhoMolar,
        const std::vector<double>& massFractions,
        double massSum, double molecularWeight,
        const std::vector<double>& moleFractions) const
    {
        double hIdeal = 0.0;
        for (std::size_t i = 0; i < massFractions.size(); ++i)
        {
            const double yi = massFractions[i]/massSum;
            if (yi == 0.0) continue;
            hIdeal += yi*nasa7Thermo().hsens_mole(
                speciesNames()[i], T)/molecularWeights()[i];
        }
        const double eIdeal = hIdeal - R_UNIVERSAL*T/molecularWeight;
        const double Ar10 = teqp::TDXDerivatives<
            ModelT, double, std::vector<double>>::template get_Arxy<1, 0>(
                model_, T, rhoMolar, moleFractions);
        return eIdeal + R_UNIVERSAL*T*Ar10/molecularWeight;
    }

public:
    // Molar mass from RFAM mass fractions [kg/mol].
    double molecularWeight(const std::vector<double>& massFractions) const
    {
        const std::vector<double>& mw = molecularWeights();
        if (massFractions.size() != mw.size())
            return std::numeric_limits<double>::quiet_NaN();
        double massSum = 0.0;
        double molesPerInputMass = 0.0;
        for (std::size_t i = 0; i < massFractions.size(); ++i)
        {
            const double yi = massFractions[i];
            if (!std::isfinite(yi) || yi < 0.0)
                return std::numeric_limits<double>::quiet_NaN();
            massSum += yi;
            molesPerInputMass += yi/mw[i];
        }
        return std::isfinite(massSum) && std::isfinite(molesPerInputMass)
            && massSum > 0.0 && molesPerInputMass > 0.0
            ? massSum/molesPerInputMass
            : std::numeric_limits<double>::quiet_NaN();
    }

    double rhoToMolar(
        double rhoMass, const std::vector<double>& massFractions) const
    {
        const double mw = molecularWeight(massFractions);
        return std::isfinite(rhoMass) && rhoMass > 0.0
            && std::isfinite(mw) && mw > 0.0 ? rhoMass / mw
            : std::numeric_limits<double>::quiet_NaN();
    }

    // Pressure [Pa] at molar density rho [mol/m3]: p = rho R T (1+L^r_01).
    // The public composition argument is RFAM mass basis, not teqp mole basis.
    double pressure(double T, double rhoMolar,
                    const std::vector<double>& massFractions) const
    {
        std::vector<double> moleFractions;
        if (!std::isfinite(T) || T <= 0.0
            || !std::isfinite(rhoMolar) || rhoMolar <= 0.0
            || !massToMoleFractions(massFractions, moleFractions))
            return std::numeric_limits<double>::quiet_NaN();
        return pressureForMoleFractions(T,rhoMolar,moleFractions);
    }

    // Sensible internal energy [J/kg]: NASA7 ideal (RFAM) + teqp residual.
    // NASA7 is mass weighted; only the residual teqp contribution uses x_i.
    double eSens(double T, double rhoMolar,
                 const std::vector<double>& massFractions) const
    {
        std::vector<double> moleFractions;
        double mw = std::numeric_limits<double>::quiet_NaN();
        if (!std::isfinite(T) || T <= 0.0
            || !std::isfinite(rhoMolar) || rhoMolar <= 0.0
            || !massToMoleFractions(
                massFractions, moleFractions, &mw)
            || !std::isfinite(mw) || mw <= 0.0)
            return std::numeric_limits<double>::quiet_NaN();

        double massSum = 0.0;
        for (double yi : massFractions) massSum += yi;
        return eSensForMoleFractions(
            T,rhoMolar,massFractions,massSum,mw,moleFractions);
    }

    // Thermodynamic derivatives by central finite differences (2nd order,
    // delta = 1e-5 relative): dp/drho|T, dp/dT|rho, Cv = de/dT|rho.  These
    // are consistent with the GERG state so the model-layer payload
    // closure (Cp-Cv identity) passes at the 1e-7 level.
    void derivatives(double T, double rhoMolar,
                     const std::vector<double>& massFractions,
                     double& dpdrho, double& dpdT, double& cv) const
    {
        dpdrho = dpdT = cv = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> moleFractions;
        double mw = std::numeric_limits<double>::quiet_NaN();
        if (!std::isfinite(T) || T <= 0.0
            || !std::isfinite(rhoMolar) || rhoMolar <= 0.0
            || !massToMoleFractions(massFractions,moleFractions,&mw)
            || !std::isfinite(mw) || mw <= 0.0)
            return;
        double massSum = 0.0;
        for (double yi : massFractions) massSum += yi;
        const double dRho = 1e-5*std::max(rhoMolar,1.0);
        const double dT = 1e-5*std::max(T,1.0);
        dpdrho = (pressureForMoleFractions(
                        T,rhoMolar + dRho,moleFractions)
                - pressureForMoleFractions(
                        T,rhoMolar - dRho,moleFractions))/(2.0*dRho);
        dpdT = (pressureForMoleFractions(
                      T + dT,rhoMolar,moleFractions)
              - pressureForMoleFractions(
                      T - dT,rhoMolar,moleFractions))/(2.0*dT);
        cv = (eSensForMoleFractions(
                  T + dT,rhoMolar,massFractions,massSum,mw,moleFractions)
            - eSensForMoleFractions(
                  T - dT,rhoMolar,massFractions,massSum,mw,moleFractions))
            /(2.0*dT);
    }

    // Isochoric energy recovery (single-phase): scan + bisection on the
    // stable branch.  Returns the temperature (K) and pressure (Pa).
    bool recoverRhoEY(double rhoMass, double eTarget,
                      const std::vector<double>& massFractions,
                      double& T, double& p) const
    {
        std::vector<double> moleFractions;
        double mw = std::numeric_limits<double>::quiet_NaN();
        if (!std::isfinite(rhoMass) || rhoMass <= 0.0
            || !std::isfinite(eTarget)
            || !massToMoleFractions(massFractions,moleFractions,&mw)
            || !std::isfinite(mw) || mw <= 0.0)
            return false;
        double massSum = 0.0;
        for (double yi : massFractions) massSum += yi;
        const double rhoMolar = rhoMass/mw;
        auto eOfT = [&](double Tq) {
            return eSensForMoleFractions(
                Tq,rhoMolar,massFractions,massSum,mw,moleFractions);
        };
        // Coarse scan for the sign change on the single-phase branch
        // (skip non-finite samples: EOS domain edges).
        constexpr int N = 256;
        double prevT = T_MIN;
        double prevE = eOfT(prevT);
        for (int i = 1; i <= N; ++i)
        {
            const double Tq = T_MIN + (T_MAX - T_MIN) * i / N;
            const double eq = eOfT(Tq);
            if (!std::isfinite(eq)) { prevE = eq; prevT = Tq; continue; }
            if (std::isfinite(prevE)
                && std::signbit(prevE - eTarget) != std::signbit(eq - eTarget))
            {
                double lo = prevT, hi = Tq;
                for (int it = 0; it < 120; ++it)
                {
                    const double mid = 0.5 * (lo + hi);
                    const double em = eOfT(mid);
                    if (!std::isfinite(em)) { lo = mid; continue; }
                    if (std::abs(em - eTarget)
                            <= 1e-9 * std::max(std::abs(eTarget), 1.0)
                        || (hi - lo) <= 1e-10 * std::max(std::abs(mid), 1.0))
                    {
                        T = mid;
                        p = pressureForMoleFractions(
                            mid,rhoMolar,moleFractions);
                        return std::isfinite(p) && p > 0.0;
                    }
                    if (std::signbit(prevE - eTarget) == std::signbit(em - eTarget))
                        lo = mid;
                    else
                        hi = mid;
                }
                T = 0.5 * (lo + hi);
                p = pressureForMoleFractions(
                    T,rhoMolar,moleFractions);
                return std::isfinite(p) && p > 0.0;
            }
            prevT = Tq;
            prevE = eq;
        }
        return false;
    }

private:
    ModelT model_;
};

} // End namespace Foam

#endif // RFAM_HAVE_TEQP

#endif
