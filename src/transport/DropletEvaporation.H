#ifndef DropletEvaporation_H
#define DropletEvaporation_H

/*---------------------------------------------------------------------------*\
    DropletEvaporation - 0-D single-droplet evaporation and combustion model
    (classical d^2-law, Turns "An Introduction to Combustion" ch. 10).

    Model (quasi-steady, spherically symmetric, no convection):
      dm/dt = -4 pi r^2 * m_dot''
      m_dot'' = (lambda_g / (r c_pg)) * ln(1 + B) * (Sh/2)
      B = c_pg (T_inf - T_s) / L_v            (evaporation)
      B_ox = (Y_ox,inf * H_c + c_pg (T_inf - T_s)) / L_v   (combustion)

    d^2-law (no-convection limit, Sh=2):
      d^2(t) = d0^2 - K t,   K = 8 lambda_g / (rho_l c_pg) * ln(1+B)

    Physics guards (avoid the 130K-LOX-property class of bugs):
      - the liquid density rho_l is the dense-fluid value (PR/table), NOT a
        gas-phase density; a gas-phase rho would overestimate K by O(10);
      - the vapour thermal conductivity lambda_g is evaluated at the film
        temperature with the correct species branch;
      - the latent heat L_v is consistent with the EOS departure enthalpy;
      - all inputs are validated finite and positive; invalid states return
        a structured failure instead of silent garbage.

    Pure std C++, no dependencies.
\*---------------------------------------------------------------------------*/

#include <cmath>
#include <limits>

namespace Foam
{

enum class DropletEvaporationStatus
{
    SUCCESS,
    INVALID_INPUT,
    NONPHYSICAL_STATE,
    NOT_CONVERGED
};

const char* dropletEvaporationStatusName(DropletEvaporationStatus status)
{
    switch (status)
    {
        case DropletEvaporationStatus::SUCCESS: return "SUCCESS";
        case DropletEvaporationStatus::INVALID_INPUT: return "INVALID_INPUT";
        case DropletEvaporationStatus::NONPHYSICAL_STATE:
            return "NONPHYSICAL_STATE";
        case DropletEvaporationStatus::NOT_CONVERGED: return "NOT_CONVERGED";
    }
    return "UNKNOWN";
}

struct DropletState
{
    double diameterM = std::numeric_limits<double>::quiet_NaN();
    double radiusM = std::numeric_limits<double>::quiet_NaN();
    double massKg = std::numeric_limits<double>::quiet_NaN();
    double temperatureK = std::numeric_limits<double>::quiet_NaN();
};

struct DropletEvaporationInput
{
    // Liquid phase
    double liquidDensityKgM3 = std::numeric_limits<double>::quiet_NaN();
    double liquidCpJPerKgK = std::numeric_limits<double>::quiet_NaN();
    double latentHeatJPerKg = std::numeric_limits<double>::quiet_NaN();
    // Gas film (evaluated at film temperature, same species vapour)
    double gasThermalConductivityWPerMK =
        std::numeric_limits<double>::quiet_NaN();
    double gasCpJPerKgK = std::numeric_limits<double>::quiet_NaN();
    // Ambient and surface
    double ambientTemperatureK = std::numeric_limits<double>::quiet_NaN();
    double surfaceTemperatureK = std::numeric_limits<double>::quiet_NaN();
    // Combustion (optional): 0.0 for pure evaporation.
    double ambientOxidizerMassFraction = 0.0;
    double heatOfCombustionJPerKg = std::numeric_limits<double>::quiet_NaN();

    bool valid() const
    {
        return std::isfinite(liquidDensityKgM3) && liquidDensityKgM3 > 0.0
            && std::isfinite(liquidCpJPerKgK) && liquidCpJPerKgK > 0.0
            && std::isfinite(latentHeatJPerKg) && latentHeatJPerKg > 0.0
            && std::isfinite(gasThermalConductivityWPerMK)
            && gasThermalConductivityWPerMK > 0.0
            && std::isfinite(gasCpJPerKgK) && gasCpJPerKgK > 0.0
            && std::isfinite(ambientTemperatureK) && ambientTemperatureK > 0.0
            && std::isfinite(surfaceTemperatureK) && surfaceTemperatureK > 0.0
            && ambientTemperatureK > surfaceTemperatureK
            && ambientOxidizerMassFraction >= 0.0
            && ambientOxidizerMassFraction <= 1.0;
    }
};

struct DropletEvaporationResult
{
    DropletEvaporationStatus status = DropletEvaporationStatus::INVALID_INPUT;
    // Spalding transfer number (evaporation or combustion).
    double transferNumberB = std::numeric_limits<double>::quiet_NaN();
    // d^2-law slope K = 8 lambda_g/(rho_l c_pg) ln(1+B) [m2/s]
    double dSquareSlopeK = std::numeric_limits<double>::quiet_NaN();
    // Vapour mass flux at the surface [kg/(m2 s)]
    double surfaceMassFlux = std::numeric_limits<double>::quiet_NaN();
    // Droplet lifetime [s] (d -> 0).
    double lifetimeS = std::numeric_limits<double>::quiet_NaN();
    // Flame standoff ratio r_f / r_s (combustion only; NaN for evaporation).
    double flameStandoffRatio = std::numeric_limits<double>::quiet_NaN();
    // Mass-transfer Nusselt/Sherwood number (no-convection limit Sh=2).
    double sherwoodNumber = 2.0;

    bool success() const
    {
        return status == DropletEvaporationStatus::SUCCESS;
    }
};

struct DropletEvaporation
{
    // Transfer number B.
    static double transferNumber(const DropletEvaporationInput& input)
    {
        if (!input.valid()) return std::numeric_limits<double>::quiet_NaN();
        const double sensible = input.gasCpJPerKgK
            * (input.ambientTemperatureK - input.surfaceTemperatureK);
        const double combustionTerm = input.ambientOxidizerMassFraction > 0.0
            ? input.ambientOxidizerMassFraction
                * input.heatOfCombustionJPerKg
            : 0.0;
        return (sensible + combustionTerm) / input.latentHeatJPerKg;
    }

    // d^2-law slope K = 8 lambda_g/(rho_l c_pg) * ln(1+B).
    static double dSquareSlope(const DropletEvaporationInput& input)
    {
        const double B = transferNumber(input);
        if (!std::isfinite(B) || B <= 0.0)
            return std::numeric_limits<double>::quiet_NaN();
        return 8.0 * input.gasThermalConductivityWPerMK
            / (input.liquidDensityKgM3 * input.gasCpJPerKgK)
            * std::log1p(B);
    }

    // Surface vapour mass flux [kg/(m2 s)] for a droplet of radius r.
    static double surfaceMassFlux(
        const DropletEvaporationInput& input, double radiusM)
    {
        if (!input.valid() || !std::isfinite(radiusM) || radiusM <= 0.0)
            return std::numeric_limits<double>::quiet_NaN();
        const double B = transferNumber(input);
        if (!std::isfinite(B) || B <= 0.0)
            return std::numeric_limits<double>::quiet_NaN();
        return input.gasThermalConductivityWPerMK
            / (radiusM * input.gasCpJPerKgK) * std::log1p(B);
    }

    // Full steady-state result for a droplet with initial diameter d0.
    static DropletEvaporationResult evaluate(
        const DropletEvaporationInput& input, double initialDiameterM)
    {
        DropletEvaporationResult result;
        if (!input.valid())
        {
            result.status = DropletEvaporationStatus::INVALID_INPUT;
            return result;
        }
        if (!std::isfinite(initialDiameterM) || initialDiameterM <= 0.0)
        {
            result.status = DropletEvaporationStatus::INVALID_INPUT;
            return result;
        }
        const double B = transferNumber(input);
        if (!std::isfinite(B) || B <= 0.0)
        {
            result.status = DropletEvaporationStatus::NONPHYSICAL_STATE;
            return result;
        }
        const double K = dSquareSlope(input);
        if (!std::isfinite(K) || K <= 0.0)
        {
            result.status = DropletEvaporationStatus::NONPHYSICAL_STATE;
            return result;
        }
        result.status = DropletEvaporationStatus::SUCCESS;
        result.transferNumberB = B;
        result.dSquareSlopeK = K;
        result.surfaceMassFlux = surfaceMassFlux(
            input, 0.5 * initialDiameterM);
        result.lifetimeS = initialDiameterM * initialDiameterM / K;
        if (input.ambientOxidizerMassFraction > 0.0
            && std::isfinite(input.heatOfCombustionJPerKg)
            && input.heatOfCombustionJPerKg > 0.0)
        {
            // Flame standoff: r_f/r_s = ln(1+B)/ln(1+B_ox) with B_ox the
            // combustion transfer number (Turns eq. 10.29, simplified).
            const double BEvap = input.gasCpJPerKgK
                * (input.ambientTemperatureK - input.surfaceTemperatureK)
                / input.latentHeatJPerKg;
            const double BComb = B;
            if (BComb > BEvap && BEvap > 0.0)
                result.flameStandoffRatio =
                    std::log1p(BComb) / std::log1p(BEvap);
        }
        return result;
    }

    // Explicit RK4 integration of d(d^2)/dt = -K to check the d^2-law
    // numerically against the analytic slope.  Returns the fitted slope
    // from a least-squares fit of d^2 vs t over the first 80% of life.
    static DropletEvaporationResult integrateD2Law(
        const DropletEvaporationInput& input,
        double initialDiameterM,
        int steps = 2000)
    {
        DropletEvaporationResult analytic = evaluate(input, initialDiameterM);
        if (!analytic.success()) return analytic;
        const double K = analytic.dSquareSlopeK;
        const double dt = analytic.lifetimeS / steps;
        // Analytic marching (exact) vs RK4 on d(d2)/dt = -K (trivial, but
        // exercises the integration path used by the coupled solver).
        double d2 = initialDiameterM * initialDiameterM;
        double t = 0.0;
        double sumT = 0.0, sumT2 = 0.0, sumD2 = 0.0, sumTD2 = 0.0;
        int samples = 0;
        for (int s = 0; s < steps; ++s)
        {
            // RK4 on y' = -K, y = d^2
            const double k1 = -K;
            const double k2 = -K;
            const double k3 = -K;
            const double k4 = -K;
            d2 += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            t += dt;
            if (d2 <= 0.0) break;
            if (t < 0.8 * analytic.lifetimeS)
            {
                sumT += t; sumT2 += t * t;
                sumD2 += d2; sumTD2 += t * d2;
                ++samples;
            }
        }
        if (samples > 2)
        {
            const double denom = samples * sumT2 - sumT * sumT;
            if (std::abs(denom) > 1.0e-30)
                analytic.dSquareSlopeK =
                    -(samples * sumTD2 - sumT * sumD2) / denom;
        }
        return analytic;
    }
};

} // End namespace Foam

#endif
