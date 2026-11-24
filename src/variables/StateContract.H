#ifndef StateContract_H
#define StateContract_H

#include "SpeciesRegistry.H"
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace Foam
{

// The conservative energy definition is a solver-wide contract.  It is not a
// chemistry-model option and may not vary between EOS, transport, or source
// operators within one Domain.
enum class EnergyConvention
{
    SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT
};

inline const char* energyConventionName(EnergyConvention convention)
{
    switch (convention)
    {
        case EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT:
            return "sensible-internal-energy-explicit-chemical-heat";
    }
    return "unknown";
}

class ConservativeStateContract
{
    std::vector<std::string> speciesNames_;
    EnergyConvention energyConvention_ =
        EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;

public:
    EnergyConvention energyConvention() const { return energyConvention_; }

    bool hasExplicitSpeciesLayout() const { return !speciesNames_.empty(); }
    const std::vector<std::string>& speciesNames() const { return speciesNames_; }

    void setSpeciesLayout(
        const std::vector<std::string>& names,
        std::size_t expectedSpeciesCount)
    {
        if (names.size() != expectedSpeciesCount)
            throw std::invalid_argument(
                "state contract species count does not match conserved state");

        SpeciesRegistry registry;
        std::unordered_set<std::string> unique;
        unique.reserve(names.size());
        for (const std::string& name : names)
        {
            if (name.empty() || !registry.has(name))
                throw std::invalid_argument(
                    "state contract contains an unknown species: " + name);
            if (!unique.insert(name).second)
                throw std::invalid_argument(
                    "state contract contains a duplicate species: " + name);
        }
        speciesNames_ = names;
    }
};

} // End namespace Foam

#endif
