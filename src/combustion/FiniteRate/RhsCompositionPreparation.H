#ifndef RhsCompositionPreparation_H
#define RhsCompositionPreparation_H

#include "SpeciesMapper.H"
#include "SpeciesRegistry.H"

#include <cmath>
#include <utility>
#include <vector>

namespace Foam
{

// Shared composition preparation for every finite-rate RHS consumer. The
// helper commits its output vectors only after the complete operation passes.
class RhsCompositionPreparation
{
public:
    static constexpr double traceConcentrationMolM3 = 1.0e-6;

    static double regularizedNonnegative(
        double value, double transition)
    {
        if (value <= 0.0) return 0.0;
        if (value >= transition) return value;
        const double scaled = value/transition;
        return transition*scaled*scaled*(2.0 - scaled);
    }

    static bool prepare(
        double rho,
        const SpeciesRegistry& registry,
        const SpeciesMapper& mapper,
        const std::vector<double>& chemistryMolecularWeights,
        const std::vector<double>& globalTemplate,
        std::vector<double>& chemistryY,
        std::vector<double>& globalY,
        double* preNormalizationMassFractionSum = nullptr)
    {
        if (!std::isfinite(rho) || rho <= 0.0
            || globalTemplate.size()
                != static_cast<std::size_t>(registry.size())
            || globalY.size() != globalTemplate.size()
            || chemistryY.size()
                != static_cast<std::size_t>(mapper.nChem())
            || chemistryMolecularWeights.size() != chemistryY.size())
            return false;

        for (int global = 0; global < registry.size(); ++global)
        {
            const double molecularWeight =
                registry.molecularWeight(global);
            if (!std::isfinite(molecularWeight)
                || molecularWeight <= 0.0)
                return false;
        }
        for (double molecularWeight : chemistryMolecularWeights)
            if (!std::isfinite(molecularWeight)
                || molecularWeight <= 0.0)
                return false;

        std::vector<double> preparedGlobal = globalTemplate;
        std::vector<double> preparedChemistry = chemistryY;
        for (int global = 0; global < registry.size(); ++global)
        {
            const double transition = traceConcentrationMolM3
                *registry.molecularWeight(global)/rho;
            preparedGlobal[global] = regularizedNonnegative(
                preparedGlobal[global], transition);
        }

        for (int chemistry = 0; chemistry < mapper.nChem(); ++chemistry)
        {
            const double transition = traceConcentrationMolM3
                *chemistryMolecularWeights[chemistry]/rho;
            preparedGlobal[mapper.globalIndex(chemistry)] =
                regularizedNonnegative(
                    preparedChemistry[chemistry], transition);
        }

        double massFractionSum = 0.0;
        for (double value : preparedGlobal)
            massFractionSum += value;
        if (!std::isfinite(massFractionSum)
            || massFractionSum <= 0.0)
            return false;

        for (double& value : preparedGlobal)
            value /= massFractionSum;
        for (int chemistry = 0; chemistry < mapper.nChem(); ++chemistry)
            preparedChemistry[chemistry] =
                preparedGlobal[mapper.globalIndex(chemistry)];

        globalY = std::move(preparedGlobal);
        chemistryY = std::move(preparedChemistry);
        if (preNormalizationMassFractionSum)
            *preNormalizationMassFractionSum = massFractionSum;
        return true;
    }
};

} // End namespace Foam

#endif
