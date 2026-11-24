#ifndef MassActionThermo_H
#define MassActionThermo_H

#include "NASA7Thermo.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace Foam
{

inline double regularizedConcentrationPower(
    double concentration,
    double exponent)
{
    constexpr double transition = 1e-6; // mol/m3
    if (exponent == 0.0) return 1.0;
    if (concentration <= 0.0) return 0.0;
    if (concentration >= transition || exponent >= 1.0)
        return std::pow(concentration, exponent);

    const double scaled = concentration / transition;
    return std::pow(transition, exponent)
        * ((3.0 - exponent) * scaled * scaled
          + (exponent - 2.0) * scaled * scaled * scaled);
}

inline double concentrationEquilibriumConstant(
    double T,
    const std::vector<std::string>& speciesNames,
    const std::vector<double>& stoich)
{
    if (!std::isfinite(T) || T <= 0.0)
        throw std::runtime_error(
            "concentrationEquilibriumConstant: temperature must be finite and positive");
    if (speciesNames.size() != stoich.size() || speciesNames.empty())
        throw std::runtime_error(
            "concentrationEquilibriumConstant: invalid reaction dimensions");

    static constexpr double Runiv = 8.314462618;
    static constexpr double pStd = 101325.0;
    static constexpr double minimumKc = 1e-30;
    const NASA7Thermo& thermo = nasa7Thermo();
    double deltaH = 0.0;
    double deltaS = 0.0;
    double deltaN = 0.0;
    bool hasReactionSpecies = false;

    for (std::size_t i = 0; i < stoich.size(); ++i)
    {
        const double nu = stoich[i];
        if (nu == 0.0) continue;
        if (!std::isfinite(nu) || !thermo.has(speciesNames[i]))
            throw std::runtime_error(
                "concentrationEquilibriumConstant: unsupported reaction species '"
                + speciesNames[i] + "'");
        hasReactionSpecies = true;
        deltaH += nu * thermo.h_mole(speciesNames[i], T);
        deltaS += nu * thermo.s_mole(speciesNames[i], T);
        deltaN += nu;
    }
    if (!hasReactionSpecies)
        throw std::runtime_error(
            "concentrationEquilibriumConstant: reaction has no stoichiometry");

    const double deltaG = deltaH - T * deltaS;
    const double logKc = -deltaG / (Runiv * T)
        - deltaN * std::log(Runiv * T / pStd);
    if (!std::isfinite(logKc))
        throw std::runtime_error(
            "concentrationEquilibriumConstant: non-finite equilibrium constant");

    const double logMinimum = std::log(minimumKc);
    const double logMaximum = std::log(std::numeric_limits<double>::max());
    if (logKc <= logMinimum) return minimumKc;
    if (logKc >= logMaximum) return std::numeric_limits<double>::max();
    return std::exp(logKc);
}

inline void validateReversibleMassAction(
    const std::vector<std::string>& speciesNames,
    const std::vector<double>& concentrationExponents,
    const std::vector<double>& stoich)
{
    if (speciesNames.size() != concentrationExponents.size()
        || speciesNames.size() != stoich.size()
        || speciesNames.empty())
        throw std::runtime_error(
            "reversible mass-action reaction has invalid dimensions");

    bool hasReactant = false;
    bool hasProduct = false;
    const NASA7Thermo& thermo = nasa7Thermo();
    for (std::size_t i = 0; i < stoich.size(); ++i)
    {
        const double nu = stoich[i];
        const double exponent = concentrationExponents[i];
        if (!std::isfinite(nu) || !std::isfinite(exponent) || exponent < 0.0)
            throw std::runtime_error(
                "reversible mass-action reaction has invalid numeric data");

        const double expected = nu < 0.0 ? -nu : 0.0;
        const double tolerance = 1e-12 * std::max(1.0, std::abs(expected));
        if (std::abs(exponent - expected) > tolerance)
            throw std::runtime_error(
                "reversible mass-action exponent mismatch for species '"
                + speciesNames[i] + "'");

        if (nu < 0.0) hasReactant = true;
        if (nu > 0.0) hasProduct = true;
        if (nu != 0.0 && !thermo.has(speciesNames[i]))
            throw std::runtime_error(
                "reversible mass-action thermodynamics unavailable for species '"
                + speciesNames[i] + "'");
    }
    if (!hasReactant || !hasProduct)
        throw std::runtime_error(
            "reversible mass-action reaction requires reactants and products");
}

} // End namespace Foam

#endif
