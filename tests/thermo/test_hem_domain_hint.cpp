#include "Domain.H"

#include <atomic>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{

using namespace Foam;

int failures = 0;

void check(bool condition, const char* label)
{
    std::cout << "  [" << (condition ? "PASS" : "FAIL") << "] "
              << label << '\n';
    if (!condition) ++failures;
}

class RecordingSeedModel final : public thermo::EquilibriumThermoModel
{
    mutable std::atomic<int> recoveries_{0};
    mutable std::atomic<int> continuationSeeds_{0};
    std::vector<std::string> species_{"O2"};
    std::string id_{"recording-domain-hem-seed"};

    static thermo::EquilibriumState state()
    {
        constexpr double pressure = 1.0e6;
        constexpr double temperature = 100.0;
        constexpr double vaporMassFraction = 0.4;
        constexpr double rhoLiquid = 1000.0;
        constexpr double rhoVapor = 10.0;
        constexpr double molecularWeight = 0.0319988;
        constexpr double gasConstant = 8.314462618;
        constexpr double hLiquid = 1.0e5;
        constexpr double hVapor = 5.0e5;

        thermo::EquilibriumPhaseProperties liquid;
        liquid.massFractions = {1.0};
        liquid.moleFractions = {1.0};
        liquid.molecularWeight = molecularWeight;
        liquid.rho = rhoLiquid;
        liquid.h = hLiquid;
        liquid.e = hLiquid - pressure/rhoLiquid;
        liquid.compressibility = pressure*molecularWeight
            /(rhoLiquid*gasConstant*temperature);

        thermo::EquilibriumPhaseProperties vapor;
        vapor.massFractions = {1.0};
        vapor.moleFractions = {1.0};
        vapor.molecularWeight = molecularWeight;
        vapor.rho = rhoVapor;
        vapor.h = hVapor;
        vapor.e = hVapor - pressure/rhoVapor;
        vapor.compressibility = pressure*molecularWeight
            /(rhoVapor*gasConstant*temperature);

        const double specificVolume =
            (1.0 - vaporMassFraction)/rhoLiquid
            + vaporMassFraction/rhoVapor;

        thermo::EquilibriumState result;
        result.region = thermo::EquilibriumRegion::TWO_PHASE;
        result.overallMassFractions = {1.0};
        result.p = pressure;
        result.T = temperature;
        result.rho = 1.0/specificVolume;
        result.h = (1.0 - vaporMassFraction)*liquid.h
            + vaporMassFraction*vapor.h;
        result.e = (1.0 - vaporMassFraction)*liquid.e
            + vaporMassFraction*vapor.e;
        result.vaporMoleFraction = vaporMassFraction;
        result.vaporMassFraction = vaporMassFraction;
        result.vaporVolumeFraction =
            (vaporMassFraction/rhoVapor)/specificVolume;
        result.acoustic.closure =
            thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
        result.acoustic.evidence = thermo::AcousticEvidence::UNVERIFIED;
        result.acoustic.pressureDensityDerivative = 1.0e5;
        result.liquid = liquid;
        result.vapor = vapor;
        return result;
    }

public:
    thermo::EquilibriumResult recoverRhoE(
        double density,
        double internalEnergy,
        const std::vector<double>& massFractions,
        const thermo::EquilibriumSeed* seed) const override
    {
        ++recoveries_;
        if (seed
            && seed->previousRegion == thermo::EquilibriumRegion::TWO_PHASE
            && std::abs(seed->previousTemperature - 100.0) < 1.0e-12
            && std::abs(seed->previousVaporMassFraction - 0.4) < 1.0e-12)
        {
            ++continuationSeeds_;
        }

        thermo::EquilibriumResult result;
        result.diagnostics.modelId = id_;
        result.diagnostics.backendStatus = "SUCCESS";
        result.state = state();
        result.status = thermo::EquilibriumStatus::SUCCESS;
        if (massFractions != std::vector<double>{1.0}
            || std::abs(density - result.state->rho) > 1.0e-10
            || std::abs(internalEnergy - result.state->e) > 1.0e-8)
        {
            result.status = thermo::EquilibriumStatus::RESIDUAL_CONTRACT_FAILED;
            result.state.reset();
        }
        return result;
    }

    const std::string& modelId() const override { return id_; }
    const std::vector<std::string>& speciesNames() const override
    { return species_; }
    EnergyConvention energyConvention() const override
    {
        return EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
    }

    int recoveries() const { return recoveries_.load(); }
    int continuationSeeds() const { return continuationSeeds_.load(); }
};

} // End anonymous namespace

int main()
{
    using namespace Foam;

    std::cout << "=== Domain HEM continuation-seed lifecycle ===\n";
    const auto model = std::make_shared<RecordingSeedModel>();
    Domain domain(1, 1.0, 1, model);

    thermo::EquilibriumResult initial = model->recoverRhoE(
        1.0/((1.0 - 0.4)/1000.0 + 0.4/10.0),
        (1.0 - 0.4)*(1.0e5 - 1.0e6/1000.0)
            + 0.4*(5.0e5 - 1.0e6/10.0),
        {1.0}, nullptr);
    check(initial.success(), "recording model creates a valid two-phase state");

    PrimitiveVars primitive(1);
    primitive.p = initial.state->p;
    primitive.T = initial.state->T;
    primitive.U = {0.0, 0.0, 0.0};
    primitive.Y = {1.0};
    primitive.Z = 0.0;
    primitive.C = 0.0;
    check(domain.replaceCompleteField({primitive}, {*initial.state}),
        "Domain publishes the initial complete two-phase state");

    domain.beginStep();
    check(domain.commit(), "Domain recommits the unchanged conserved state");
    check(model->continuationSeeds() == 1,
        "Domain commit forwards previous region/T/quality as a seed");
    check(domain.Q(0).region == thermo::EquilibriumRegion::TWO_PHASE
            && std::abs(domain.Q(0).vaporMassFraction - 0.4) < 1.0e-12,
        "continuation seed does not change the recovered physical state");

    thermo::EquilibriumSeed restartSeed;
    restartSeed.pressure = domain.Q(0).p;
    restartSeed.temperature = domain.Q(0).T;
    restartSeed.previousRegion = domain.Q(0).region;
    restartSeed.previousTemperature = domain.Q(0).T;
    restartSeed.previousVaporMassFraction = domain.Q(0).vaporMassFraction;
    check(domain.replaceConservedField({domain.U(0)}, {restartSeed}),
        "explicit restart equilibrium seed is accepted");
    check(model->continuationSeeds() == 2,
        "restart replacement forwards the explicit phase continuation seed");

    if (failures != 0)
    {
        std::cerr << "HEM_DOMAIN_HINT_FAIL failures=" << failures << '\n';
        return 1;
    }
    std::cout << "HEM_DOMAIN_HINT_PASS\n";
    return 0;
}
