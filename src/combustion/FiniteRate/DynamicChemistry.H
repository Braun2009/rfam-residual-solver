#ifndef DynamicChemistry_H
#define DynamicChemistry_H

/*---------------------------------------------------------------------------*\
    DynamicChemistry - ChemistryModel driven by data from MechanismReader.

    Implements the same interface as TwoStepCO but loaded at runtime from
    a YAML mechanism file. Uses ArrheniusRate for each reaction.

    Source terms per reaction j:
        r_j = A_j T^n_j exp(-Ea_j/R T) * PROD_k [C_k]^nu_k
        omega_i += stoich_ij * MW_i * r_j   (net mass rate [kg/m3/s])
        qdot    += (-DH_j) * r_j             (heat release [W/m3])
\*---------------------------------------------------------------------------*/

#include "ChemistryModel.H"
#include "ArrheniusRate.H"
#include "MassActionThermo.H"
#include "MechanismReader.H"
#include "SpeciesRegistry.H"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Foam
{

class DynamicChemistry : public ChemistryModel
{
private:
    static constexpr double R_J = 8.314462618;  // [J/(mol K)]
    static constexpr double molecularWeightTolerance = 1.0e-10;
    static constexpr double formationEnthalpyTolerance = 1.0e-8;
    static constexpr double reactionEnthalpyTolerance = 1.0e-8;
    static constexpr double elementBalanceTolerance = 1.0e-12;
    static constexpr double massBalanceTolerance = 1.0e-10;

    std::vector<ArrheniusRate>  rates_;
    // per-reaction: concExp[rxn][species_local_idx] = exponent
    // per-reaction: stoich[rxn][species_local_idx]  = net stoich (positive = product)
    // per-reaction: DH[rxn] = [J/mol_limiting_reactant]
    std::vector<std::vector<double>> concExp_;
    std::vector<std::vector<double>> stoich_;
    std::vector<double>              DH_;
    std::vector<bool>                reversible_;

    void build(const MechanismData& mech)
    {
        const SpeciesRegistry registry;
        nSpecies_    = (int)mech.species.size();
        nReactions_  = (int)mech.reactions.size();
        speciesNames_.resize(nSpecies_);
        molecularWeights_.resize(nSpecies_);
        for (int i = 0; i < nSpecies_; ++i)
        {
            const SpeciesData& species = mech.species[i];
            if (std::find(
                    speciesNames_.begin(), speciesNames_.begin() + i,
                    species.name) != speciesNames_.begin() + i)
                throw std::runtime_error(
                    "DynamicChemistry: duplicate species '" + species.name
                    + "'");

            const int globalIndex = registry.index(species.name);
            if (globalIndex < 0)
                throw std::runtime_error(
                    "DynamicChemistry: species '" + species.name
                    + "' is not present in SpeciesRegistry");

            const double referenceMW = registry.molecularWeight(globalIndex);
            if (!std::isfinite(species.MW) || species.MW <= 0.0)
                throw std::runtime_error(
                    "DynamicChemistry: invalid molecular weight for species '"
                    + species.name + "'");
            if (std::abs(species.MW - referenceMW)
                > molecularWeightTolerance
                    * std::max(1.0, std::abs(referenceMW)))
                throw std::runtime_error(
                    "DynamicChemistry: molecular weight mismatch for species '"
                    + species.name + "'");

            const double referenceHf =
                registry.formationEnthalpy(globalIndex);
            if (!std::isfinite(species.hf))
                throw std::runtime_error(
                    "DynamicChemistry: invalid formation enthalpy for species '"
                    + species.name + "'");
            if (std::abs(species.hf - referenceHf)
                > formationEnthalpyTolerance
                    * std::max(1.0, std::abs(referenceHf)))
                throw std::runtime_error(
                    "DynamicChemistry: formation enthalpy mismatch for species '"
                    + species.name + "'");

            speciesNames_[i] = species.name;
            molecularWeights_[i] = referenceMW;
        }

        for (const auto& rxn : mech.reactions)
        {
            if (rxn.name.empty())
                throw std::runtime_error(
                    "DynamicChemistry: reaction name must not be empty");
            if (!std::isfinite(rxn.A) || rxn.A <= 0.0)
                throw std::runtime_error(
                    "DynamicChemistry: forward pre-exponential must be finite and positive in reaction '"
                    + rxn.name + "'");
            if (!std::isfinite(rxn.n)
                || !std::isfinite(rxn.Ea)
                || !std::isfinite(rxn.DH))
                throw std::runtime_error(
                    "DynamicChemistry: invalid numeric data in reaction '"
                    + rxn.name + "'");

            std::vector<double> ce(nSpecies_, 0.0);
            bool hasPositiveExponent = false;
            for (const auto& kv : rxn.concExp)
            {
                int idx = mech.speciesIndex(kv.first);
                if (idx < 0)
                    throw std::runtime_error(
                        "DynamicChemistry: unknown concentration species '"
                        + kv.first + "' in reaction '" + rxn.name + "'");
                if (!std::isfinite(kv.second) || kv.second < 0.0)
                    throw std::runtime_error(
                        "DynamicChemistry: invalid concentration exponent for species '"
                        + kv.first + "' in reaction '" + rxn.name + "'");
                ce[idx] = kv.second;
                hasPositiveExponent = hasPositiveExponent || kv.second > 0.0;
            }
            if (!hasPositiveExponent)
                throw std::runtime_error(
                    "DynamicChemistry: reaction '" + rxn.name
                    + "' requires at least one positive concentration exponent");

            std::vector<double> st(nSpecies_, 0.0);
            for (const auto& kv : rxn.stoich)
            {
                int idx = mech.speciesIndex(kv.first);
                if (idx < 0)
                    throw std::runtime_error(
                        "DynamicChemistry: unknown stoichiometric species '"
                        + kv.first + "' in reaction '" + rxn.name + "'");
                if (!std::isfinite(kv.second) || kv.second == 0.0)
                    throw std::runtime_error(
                        "DynamicChemistry: invalid stoichiometric coefficient for species '"
                        + kv.first + "' in reaction '" + rxn.name + "'");
                st[idx] = kv.second;
            }

            bool hasReactant = false;
            bool hasProduct = false;
            long double massResidual = 0.0L;
            long double massScale = 0.0L;
            long double elementResidual[4] = {0.0L, 0.0L, 0.0L, 0.0L};
            long double elementScale[4] = {0.0L, 0.0L, 0.0L, 0.0L};
            long double derivedDH = 0.0L;
            auto checkedAccumulate = [&, this](
                long double& aggregate,
                long double term,
                const char* quantity)
            {
                const long double next = aggregate + term;
                if (!std::isfinite(term) || !std::isfinite(next))
                    throw std::runtime_error(
                        "DynamicChemistry: non-finite reaction aggregate '"
                        + std::string(quantity) + "' in '" + rxn.name + "'");
                aggregate = next;
            };
            for (int i = 0; i < nSpecies_; ++i)
            {
                const double nu = st[i];
                if (nu == 0.0) continue;
                hasReactant = hasReactant || nu < 0.0;
                hasProduct = hasProduct || nu > 0.0;

                const int globalIndex = registry.index(speciesNames_[i]);
                const double mw = registry.molecularWeight(globalIndex);
                const long double nuExtended = static_cast<long double>(nu);
                const long double mwExtended = static_cast<long double>(mw);
                checkedAccumulate(
                    massResidual, nuExtended * mwExtended, "mass residual");
                checkedAccumulate(
                    massScale, std::abs(nuExtended) * mwExtended, "mass scale");
                const long double elements[4] = {
                    static_cast<long double>(registry.nC(globalIndex)),
                    static_cast<long double>(registry.nH(globalIndex)),
                    static_cast<long double>(registry.nO(globalIndex)),
                    static_cast<long double>(registry.nN(globalIndex))};
                for (int element = 0; element < 4; ++element)
                {
                    checkedAccumulate(
                        elementResidual[element],
                        nuExtended * elements[element],
                        "element residual");
                    checkedAccumulate(
                        elementScale[element],
                        std::abs(nuExtended) * elements[element],
                        "element scale");
                }
                checkedAccumulate(
                    derivedDH,
                    nuExtended * static_cast<long double>(
                        registry.formationEnthalpy(globalIndex)),
                    "reaction enthalpy");
            }
            if (!hasReactant || !hasProduct)
                throw std::runtime_error(
                    "DynamicChemistry: reaction '" + rxn.name
                    + "' requires reactants and products");

            const long double maximumStored =
                static_cast<long double>(std::numeric_limits<double>::max());
            auto requireSupportedAggregate = [&](
                long double value,
                const char* quantity)
            {
                if (!std::isfinite(value) || std::abs(value) > maximumStored)
                    throw std::runtime_error(
                        "DynamicChemistry: reaction aggregate '"
                        + std::string(quantity)
                        + "' exceeds the supported numeric range in '"
                        + rxn.name + "'");
            };
            requireSupportedAggregate(massResidual, "mass residual");
            requireSupportedAggregate(massScale, "mass scale");
            requireSupportedAggregate(derivedDH, "reaction enthalpy");
            for (int element = 0; element < 4; ++element)
            {
                requireSupportedAggregate(
                    elementResidual[element], "element residual");
                requireSupportedAggregate(
                    elementScale[element], "element scale");
            }

            for (int element = 0; element < 4; ++element)
            {
                if (std::abs(elementResidual[element])
                    > static_cast<long double>(elementBalanceTolerance)
                        * std::max(1.0L, elementScale[element]))
                    throw std::runtime_error(
                        "DynamicChemistry: reaction element balance mismatch in '"
                        + rxn.name + "'");
            }
            if (std::abs(massResidual)
                > static_cast<long double>(massBalanceTolerance)
                    * std::max(1.0L, massScale))
                throw std::runtime_error(
                    "DynamicChemistry: reaction mass balance mismatch in '"
                    + rxn.name + "'");

            const long double enthalpyScale = std::max(
                {1.0L, std::abs(static_cast<long double>(rxn.DH)),
                 std::abs(derivedDH)});
            const long double enthalpyDifference =
                static_cast<long double>(rxn.DH) - derivedDH;
            if (!std::isfinite(enthalpyScale)
                || !std::isfinite(enthalpyDifference))
                throw std::runtime_error(
                    "DynamicChemistry: non-finite reaction enthalpy comparison in '"
                    + rxn.name + "'");
            if (std::abs(enthalpyDifference)
                > static_cast<long double>(reactionEnthalpyTolerance)
                    * enthalpyScale)
                throw std::runtime_error(
                    "DynamicChemistry: reaction enthalpy mismatch in '"
                    + rxn.name + "'");

            if (rxn.reversible)
                validateReversibleMassAction(speciesNames_, ce, st);
            rates_.emplace_back(rxn.A, rxn.n, rxn.Ea);
            DH_.push_back(static_cast<double>(derivedDH));
            concExp_.push_back(ce);
            stoich_.push_back(st);
            reversible_.push_back(rxn.reversible);
        }
    }

public:
    explicit DynamicChemistry(const MechanismData& mech) { build(mech); }
    explicit DynamicChemistry(const std::string& filename)
    {
        build(MechanismReader::read(filename));
    }

    void computeSourceTerms(
        double T, double p, double rho,
        const std::vector<double>& Y,
        std::vector<double>& omega_i,
        double& qdot) const override
    {
        omega_i.assign(nSpecies_, 0.0);
        qdot = 0.0;

        // Molar concentrations [mol/m3]
        std::vector<double> C(nSpecies_);
        for (int i = 0; i < nSpecies_; ++i)
        {
            double MW = molecularWeights_[i];
            C[i] = (MW > 0 && Y[i] > 0) ? rho * Y[i] / MW : 0.0;
        }

        for (int j = 0; j < nReactions_; ++j)
        {
            // Rate = A T^n exp(-Ea/RT) * PROD C_k^ce_k
            // If any concentration with positive exponent is zero → rate = 0.
            const double kForward = rates_[j].rate(T);
            double rForward = kForward;
            for (int i = 0; i < nSpecies_ && rForward > 0.0; ++i)
            {
                if (concExp_[j][i] == 0.0) continue;
                if (C[i] <= 0.0) { rForward = 0.0; break; }
                rForward *= reversible_[j]
                    ? std::pow(C[i], concExp_[j][i])
                    : regularizedConcentrationPower(C[i], concExp_[j][i]);
            }
            double r = rForward;
            if (reversible_[j])
            {
                double reverseProduct = 1.0;
                for (int i = 0; i < nSpecies_; ++i)
                {
                    if (stoich_[j][i] <= 0.0) continue;
                    if (C[i] <= 0.0) { reverseProduct = 0.0; break; }
                    reverseProduct *= std::pow(C[i], stoich_[j][i]);
                }
                const double Kc = concentrationEquilibriumConstant(
                    T, speciesNames_, stoich_[j]);
                r -= (kForward / Kc) * reverseProduct;
            }

            // Mass source terms [kg/m3/s] = stoich * MW * r
            for (int i = 0; i < nSpecies_; ++i)
                omega_i[i] += stoich_[j][i] * molecularWeights_[i] * r;

            // Heat release [W/m3] = -DH * r   (DH < 0 => exothermic => qdot > 0)
            qdot += (-DH_[j]) * r;
        }
    }

    void print() const override
    {
        std::cout << "DynamicChemistry: " << nSpecies_ << " species, "
                  << nReactions_ << " reactions\n";
        for (int i = 0; i < nSpecies_; ++i)
            std::cout << "  " << speciesNames_[i]
                      << " MW=" << molecularWeights_[i] << "\n";
    }
};

} // End namespace Foam
#endif
