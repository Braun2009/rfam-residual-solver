#ifndef TabulatedIdealCaloric_H
#define TabulatedIdealCaloric_H

#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace Foam
{
namespace phaseEquilibrium
{

enum class TabulatedIdealCaloricStatus
{
    SUCCESS,
    UNKNOWN_SPECIES,
    INVALID_TEMPERATURE,
    OUT_OF_RANGE,
    INVALID_FROZEN_TABLE,
    NONFINITE_PROPERTY
};

const char* tabulatedIdealCaloricStatusName(
    TabulatedIdealCaloricStatus status);

struct TabulatedIdealCaloricProperties
{
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double sensibleEnthalpy = std::numeric_limits<double>::quiet_NaN();
    double cp = std::numeric_limits<double>::quiet_NaN();
    // Integral of cp(T)/T from 298.15 K on the same frozen Hermite
    // caloric curve.  This is an ideal-gas standard-state entropy increment,
    // not a dense-fluid or mixing entropy.
    double sensibleEntropy = std::numeric_limits<double>::quiet_NaN();
};

struct TabulatedIdealCaloricResult
{
    TabulatedIdealCaloricStatus status =
        TabulatedIdealCaloricStatus::UNKNOWN_SPECIES;
    std::optional<TabulatedIdealCaloricProperties> state;
    std::string species;
    double supportedTemperatureMinimum =
        std::numeric_limits<double>::quiet_NaN();
    double supportedTemperatureMaximum =
        std::numeric_limits<double>::quiet_NaN();

    bool success() const
    {
        return status == TabulatedIdealCaloricStatus::SUCCESS
            && state.has_value();
    }
};

struct TabulatedIdealCaloricNode
{
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double sensibleEnthalpy = std::numeric_limits<double>::quiet_NaN();
    double cp = std::numeric_limits<double>::quiet_NaN();
};

struct TabulatedIdealCaloricTableData
{
    std::string species;
    std::vector<TabulatedIdealCaloricNode> nodes;
};

// Frozen mass-specific ideal-gas sensible calorics. The table contains only
// HMASS_IDEALGAS(T)-HMASS_IDEALGAS(298.15 K) and CP0MASS from CoolProp 7.2.0.
// The entropy increment is derived from the frozen interpolant by integrating
// cp/T, so dh/dT=cp and ds/dT=cp/T share one curve.
// It supplements the ideal caloric part of PR; it is not a liquid/two-phase
// model, a transport model, or an HEOS residual-property backend.
class TabulatedIdealCaloric
{
public:
    static constexpr double referenceTemperature = 298.15;

    bool has(const std::string& species) const;

    bool temperatureRange(
        const std::string& species,
        double& minimum,
        double& maximum) const;

    std::vector<TabulatedIdealCaloricTableData> data() const;

    static std::string provenanceManifest(
        const std::vector<TabulatedIdealCaloricTableData>& tables);

    static std::string provenanceSha256(
        const std::vector<TabulatedIdealCaloricTableData>& tables);

    std::string provenanceManifest() const;
    std::string provenanceSha256() const;

    TabulatedIdealCaloricResult evaluate(
        const std::string& species,
        double temperature) const;
};

const TabulatedIdealCaloric& tabulatedIdealCaloric();

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
