#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "FiniteRate/RhsCompositionPreparation.H"

namespace
{

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& label)
{
    ++checks;
    if (!condition)
    {
        ++failures;
        std::cerr << "FAIL: " << label << '\n';
    }
}

bool close(double lhs, double rhs, double relative = 2.0e-15)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= relative
            *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

} // namespace

int main()
{
    using namespace Foam;

    SpeciesRegistry registry;
    const std::vector<std::string> chemistryNames = {
        "C3H8", "O2", "CO", "CO2", "H2O", "N2"};
    SpeciesMapper mapper(registry, chemistryNames);
    std::vector<double> molecularWeights;
    for (int chemistry = 0; chemistry < mapper.nChem(); ++chemistry)
        molecularWeights.push_back(
            registry.molecularWeight(mapper.globalIndex(chemistry)));

    const double rho = 1.0;
    const double fuelTransition =
        RhsCompositionPreparation::traceConcentrationMolM3
        *molecularWeights[0]/rho;
    const double fuelMassFraction = 0.5*fuelTransition;
    const double remaining = 1.0 - fuelMassFraction;
    std::vector<double> chemistryY = {
        fuelMassFraction, 0.23*remaining, 0.0, 0.0, 0.0,
        0.77*remaining};
    std::vector<double> globalTemplate = {
        0.20, 0.20, 0.00, 0.00, 0.00, 0.20, 0.10, 0.10, 0.20};
    std::vector<double> globalY(globalTemplate.size(), 0.0);
    double preNormalizationSum = -1.0;

    check(RhsCompositionPreparation::prepare(
            rho, registry, mapper, molecularWeights, globalTemplate,
            chemistryY, globalY, &preNormalizationSum),
        "the production counterexample prepares successfully");
    check(close(preNormalizationSum, 0.9999999944880474),
        "mapper overwrite precedes complete-vector normalization");
    check(close(std::accumulate(globalY.begin(), globalY.end(), 0.0), 1.0),
        "prepared global composition is normalized");
    for (int chemistry = 0; chemistry < mapper.nChem(); ++chemistry)
        check(chemistryY[chemistry]
                == globalY[mapper.globalIndex(chemistry)],
            "normalized chemistry entries are gathered from global order");
    check(close(
            chemistryY[0],
            0.375*fuelTransition/preNormalizationSum),
        "sub-transition fuel smoothing is normalized after overwrite");

    std::vector<double> atomicChemistry = chemistryY;
    std::vector<double> atomicGlobal = globalY;
    std::vector<double> invalidWeights = molecularWeights;
    invalidWeights[2] = std::numeric_limits<double>::quiet_NaN();
    check(!RhsCompositionPreparation::prepare(
            rho, registry, mapper, invalidWeights, globalTemplate,
            atomicChemistry, atomicGlobal, nullptr)
            && atomicChemistry == chemistryY
            && atomicGlobal == globalY,
        "invalid molecular weight leaves both output vectors unchanged");

    check(!RhsCompositionPreparation::prepare(
            0.0, registry, mapper, molecularWeights, globalTemplate,
            atomicChemistry, atomicGlobal, nullptr),
        "nonpositive density is rejected");

    std::vector<double> zeroChemistry(
        static_cast<std::size_t>(mapper.nChem()), 0.0);
    std::vector<double> zeroTemplate(
        static_cast<std::size_t>(registry.size()), 0.0);
    std::vector<double> zeroGlobal(zeroTemplate.size(), 1.0);
    const std::vector<double> zeroChemistryBefore = zeroChemistry;
    const std::vector<double> zeroGlobalBefore = zeroGlobal;
    check(!RhsCompositionPreparation::prepare(
            rho, registry, mapper, molecularWeights, zeroTemplate,
            zeroChemistry, zeroGlobal, nullptr)
            && zeroChemistry == zeroChemistryBefore
            && zeroGlobal == zeroGlobalBefore,
        "zero complete composition fails atomically");

    if (failures != 0)
    {
        std::cerr << "RHS_COMPOSITION_PREPARATION_FAIL checks="
                  << checks << " failures=" << failures << '\n';
        return 1;
    }
    std::cout << "RHS_COMPOSITION_PREPARATION_PASS checks="
              << checks << "\n";
    return 0;
}
