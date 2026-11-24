#include "CondensedCarbonModel.H"

#include "ThermoDataProvenance.H"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Foam
{
namespace
{

constexpr long double relativeTolerance = 1.0e-10L;
constexpr long double absoluteMassTolerance = 1.0e-12L;
constexpr long double absoluteElementTolerance = 1.0e-12L;
constexpr long double absoluteEnergyTolerance = 1.0e-6L;

bool finiteNonnegative(double value)
{
    return std::isfinite(value) && value >= 0.0;
}

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

bool closeAggregate(
    long double residual,
    long double scale,
    long double absoluteTolerance)
{
    return std::isfinite(residual) && std::isfinite(scale)
        && std::abs(residual)
            <= absoluteTolerance
                + relativeTolerance*std::max(1.0L,scale);
}

bool supportedDouble(long double value)
{
    return std::isfinite(value)
        && value <= static_cast<long double>(
            std::numeric_limits<double>::max())
        && value >= -static_cast<long double>(
            std::numeric_limits<double>::max());
}

bool addChecked(double base,long double increment,double& result)
{
    const long double candidate = static_cast<long double>(base)+increment;
    if(!supportedDouble(candidate)) return false;
    result=static_cast<double>(candidate);
    return std::isfinite(result);
}

std::string nonemptyModelId(const std::string& id,const char* kind)
{
    if(id.empty())
        throw std::invalid_argument(
            std::string("condensed-carbon ")+kind+" model id is empty");
    return id;
}

} // namespace

const char* condensedCarbonMomentModelName(CondensedCarbonMomentModel model)
{
    switch(model)
    {
        case CondensedCarbonMomentModel::ONE_MOMENT: return "one-moment";
        case CondensedCarbonMomentModel::TWO_MOMENT: return "two-moment";
    }
    return "unknown";
}

const char* condensedCarbonApplicabilityStatusName(
    CondensedCarbonApplicabilityStatus status)
{
    switch(status)
    {
        case CondensedCarbonApplicabilityStatus::SUPPORTED:
            return "SUPPORTED";
        case CondensedCarbonApplicabilityStatus::WARNING_RESEARCH_ONLY:
            return "WARNING_RESEARCH_ONLY";
        case CondensedCarbonApplicabilityStatus::OUTSIDE_DILUTE_LIMIT:
            return "OUTSIDE_DILUTE_LIMIT";
        case CondensedCarbonApplicabilityStatus::SLIP_MODEL_REQUIRED:
            return "SLIP_MODEL_REQUIRED";
        case CondensedCarbonApplicabilityStatus::SOOT_ENERGY_MODEL_REQUIRED:
            return "SOOT_ENERGY_MODEL_REQUIRED";
        case CondensedCarbonApplicabilityStatus::FAILED:
            return "FAILED";
    }
    return "UNKNOWN";
}

bool CondensedCarbonState::valid(CondensedCarbonMomentModel model) const noexcept
{
    if(!finiteNonnegative(rhoSoot)||!finiteNonnegative(numberDensity))
        return false;
    if(model==CondensedCarbonMomentModel::ONE_MOMENT)
        return numberDensity==0.0;
    return (rhoSoot==0.0)==(numberDensity==0.0);
}

bool WallCarbonInventory::valid() const noexcept
{
    return finiteNonnegative(massPerArea);
}

bool CondensedCarbonApplicabilityLimits::valid() const noexcept
{
    return finitePositive(maximumMassLoading)
        &&finitePositive(maximumVolumeFraction)
        &&finitePositive(maximumHeatCapacityRatio)
        &&finitePositive(maximumStokesNumber);
}

void CondensedCarbonOptions::validate() const
{
    if(!finitePositive(carbonMolarMass)
       ||!std::isfinite(carbonFormationEnthalpy)
       ||!finitePositive(materialDensity)
       ||!finitePositive(diagnosticSpecificHeat)
       ||!limits.valid())
        throw std::invalid_argument(
            "invalid condensed-carbon material or applicability data");
    if(momentModel==CondensedCarbonMomentModel::ONE_MOMENT
       &&!finitePositive(oneMomentParticleDiameter))
        throw std::invalid_argument(
            "one-moment condensed carbon requires a particle diameter");
    if(!initialState.valid(momentModel)||!inletState.valid(momentModel))
        throw std::invalid_argument(
            "invalid condensed-carbon initial or inlet state");
    if(molecularDiffusionEnabled||turbulentDiffusionEnabled
       ||thermophoresisEnabled)
        throw std::invalid_argument(
            "condensed-carbon v1 does not implement diffusion or thermophoresis");
}

SootApplicabilityReport evaluateCondensedCarbonApplicability(
    const CondensedCarbonOptions& options,
    const CondensedCarbonApplicabilityInputs& inputs)
{
    SootApplicabilityReport result;
    try { options.validate(); }
    catch(...) { return result; }
    if(!finitePositive(inputs.rhoFluid)
       ||!finiteNonnegative(inputs.rhoSoot)
       ||!finitePositive(inputs.fluidCp)
       ||!finitePositive(inputs.dynamicViscosity)
       ||!finitePositive(inputs.fluidTimeScale)
       ||!finitePositive(inputs.particleDiameter))
        return result;

    result.massLoading=inputs.rhoSoot/inputs.rhoFluid;
    result.volumeFraction=inputs.rhoSoot/options.materialDensity;
    result.heatCapacityRatio=
        inputs.rhoSoot*options.diagnosticSpecificHeat
        /(inputs.rhoFluid*inputs.fluidCp);
    const double particleResponseTime=
        options.materialDensity*inputs.particleDiameter
        *inputs.particleDiameter/(18.0*inputs.dynamicViscosity);
    result.stokesNumber=particleResponseTime/inputs.fluidTimeScale;
    if(!finiteNonnegative(result.massLoading)
       ||!finiteNonnegative(result.volumeFraction)
       ||!finiteNonnegative(result.heatCapacityRatio)
       ||!finiteNonnegative(result.stokesNumber))
        return result;

    if(result.massLoading>options.limits.maximumMassLoading
       ||result.volumeFraction>options.limits.maximumVolumeFraction)
        result.status=
            CondensedCarbonApplicabilityStatus::OUTSIDE_DILUTE_LIMIT;
    else if(result.stokesNumber>options.limits.maximumStokesNumber)
        result.status=
            CondensedCarbonApplicabilityStatus::SLIP_MODEL_REQUIRED;
    else if(result.heatCapacityRatio
            >options.limits.maximumHeatCapacityRatio)
        result.status=
            CondensedCarbonApplicabilityStatus::SOOT_ENERGY_MODEL_REQUIRED;
    else
        result.status=options.productionValidated
            ?CondensedCarbonApplicabilityStatus::SUPPORTED
            :CondensedCarbonApplicabilityStatus::WARNING_RESEARCH_ONLY;
    return result;
}

double CondensedCarbonExchangeContract::expectedFluidEnergyRate(
    const PrimitiveVars& fluid,
    const std::vector<std::string>& speciesNames,
    const SpeciesRegistry& registry,
    const CondensedCarbonOptions& options,
    const CondensedCarbonExchange& exchange)
{
    options.validate();
    if(exchange.fluidSpeciesMassRate.size()!=speciesNames.size()
       ||fluid.U.size()!=3)
        throw std::invalid_argument(
            "condensed-carbon energy source layout mismatch");
    long double chemicalRate=0.0L;
    long double fluidMassRate=0.0L;
    for(std::size_t species=0;species<speciesNames.size();++species)
    {
        const double rate=exchange.fluidSpeciesMassRate[species];
        const int global=registry.index(speciesNames[species]);
        if(!std::isfinite(rate)||global<0)
            throw std::invalid_argument(
                "invalid condensed-carbon fluid species source");
        const double molecularWeight=registry.molecularWeight(global);
        if(!finitePositive(molecularWeight))
            throw std::invalid_argument(
                "invalid condensed-carbon species molecular weight");
        fluidMassRate+=static_cast<long double>(rate);
        chemicalRate-=static_cast<long double>(rate)
            /static_cast<long double>(molecularWeight)
            *static_cast<long double>(registry.formationEnthalpy(global));
    }
    if(!std::isfinite(exchange.sootMassRate))
        throw std::invalid_argument("invalid condensed-carbon soot source");
    chemicalRate-=static_cast<long double>(exchange.sootMassRate)
        /static_cast<long double>(options.carbonMolarMass)
        *static_cast<long double>(options.carbonFormationEnthalpy);
    long double velocitySquared=0.0L;
    for(double velocity:fluid.U)
    {
        if(!std::isfinite(velocity))
            throw std::invalid_argument(
                "invalid condensed-carbon source velocity");
        velocitySquared+=static_cast<long double>(velocity)*velocity;
    }
    const long double total=chemicalRate
        +0.5L*velocitySquared*fluidMassRate;
    if(!supportedDouble(total))
        throw std::overflow_error(
            "condensed-carbon expected energy source overflow");
    return static_cast<double>(total);
}

CondensedCarbonExchangeResult CondensedCarbonExchangeContract::evaluate(
    const ConservedVars& fluidBase,
    const PrimitiveVars& fluid,
    const CondensedCarbonState& carbonBase,
    const std::vector<std::string>& speciesNames,
    const SpeciesRegistry& registry,
    const CondensedCarbonOptions& options,
    const CondensedCarbonExchange& exchange,
    double dt)
{
    CondensedCarbonExchangeResult result;
    result.fluidCandidate=fluidBase;
    result.carbonCandidate=carbonBase;
    auto reject=[&](const std::string& reason)
    {
        result.accepted=false;
        result.reason=reason;
        return result;
    };
    try { options.validate(); }
    catch(const std::exception& error) { return reject(error.what()); }
    if(!std::isfinite(dt)||dt<0.0)
        return reject("invalid condensed-carbon exchange timestep");
    if(!carbonBase.valid(options.momentModel))
        return reject("invalid condensed-carbon base state");
    if(exchange.fluidSpeciesMassRate.size()!=speciesNames.size()
       ||fluidBase.rhoY.size()!=speciesNames.size()
       ||fluid.U.size()!=3||fluidBase.rhoU.size()!=3)
        return reject("condensed-carbon exchange layout mismatch");
    if(!std::isfinite(exchange.sootMassRate)
       ||!std::isfinite(exchange.sootNumberRate)
       ||!std::isfinite(exchange.fluidEnergyRate))
        return reject("non-finite condensed-carbon exchange");
    if(options.momentModel==CondensedCarbonMomentModel::ONE_MOMENT
       &&exchange.sootNumberRate!=0.0)
        return reject("one-moment exchange contains a number source");

    long double fluidMassRate=0.0L;
    long double massScale=std::abs(
        static_cast<long double>(exchange.sootMassRate));
    std::array<long double,4> elementResidual={0.0L,0.0L,0.0L,0.0L};
    std::array<long double,4> elementScale={0.0L,0.0L,0.0L,0.0L};
    for(std::size_t species=0;species<speciesNames.size();++species)
    {
        const double rate=exchange.fluidSpeciesMassRate[species];
        const int global=registry.index(speciesNames[species]);
        if(!std::isfinite(rate)||global<0)
            return reject("invalid condensed-carbon fluid species source");
        const double molecularWeight=registry.molecularWeight(global);
        if(!finitePositive(molecularWeight))
            return reject("invalid condensed-carbon molecular weight");
        const long double extendedRate=rate;
        fluidMassRate+=extendedRate;
        massScale+=std::abs(extendedRate);
        const int atoms[4]={registry.nC(global),registry.nH(global),
                            registry.nO(global),registry.nN(global)};
        const long double molarRate=
            extendedRate/static_cast<long double>(molecularWeight);
        for(int element=0;element<4;++element)
        {
            const long double contribution=molarRate*atoms[element];
            elementResidual[element]+=contribution;
            elementScale[element]+=std::abs(contribution);
        }
    }
    const long double sootMassRate=exchange.sootMassRate;
    const long double massResidual=fluidMassRate+sootMassRate;
    const long double sootCarbonMolarRate=sootMassRate
        /static_cast<long double>(options.carbonMolarMass);
    elementResidual[0]+=sootCarbonMolarRate;
    elementScale[0]+=std::abs(sootCarbonMolarRate);
    if(!supportedDouble(fluidMassRate)||!supportedDouble(massResidual))
        return reject("condensed-carbon mass aggregate overflow");
    result.diagnostics.fluidMassRate=static_cast<double>(fluidMassRate);
    result.diagnostics.massResidual=static_cast<double>(massResidual);
    if(!closeAggregate(massResidual,massScale,absoluteMassTolerance))
        return reject("condensed-carbon mass exchange is not closed");
    for(int element=0;element<4;++element)
    {
        if(!supportedDouble(elementResidual[element])
           ||!supportedDouble(elementScale[element]))
            return reject("condensed-carbon element aggregate overflow");
        result.diagnostics.elementResidual[element]=
            static_cast<double>(elementResidual[element]);
        if(!closeAggregate(elementResidual[element],elementScale[element],
                           absoluteElementTolerance))
            return reject("condensed-carbon element exchange is not closed");
    }

    double expectedEnergy=0.0;
    try
    {
        expectedEnergy=expectedFluidEnergyRate(
            fluid,speciesNames,registry,options,exchange);
    }
    catch(const std::exception& error) { return reject(error.what()); }
    result.diagnostics.expectedFluidEnergyRate=expectedEnergy;
    const long double energyResidual=
        static_cast<long double>(exchange.fluidEnergyRate)-expectedEnergy;
    const long double energyScale=std::abs(
        static_cast<long double>(exchange.fluidEnergyRate))
        +std::abs(static_cast<long double>(expectedEnergy));
    if(!supportedDouble(energyResidual))
        return reject("condensed-carbon energy aggregate overflow");
    result.diagnostics.energyResidual=static_cast<double>(energyResidual);
    if(!closeAggregate(energyResidual,energyScale,absoluteEnergyTolerance))
        return reject("condensed-carbon energy exchange is not closed");

    for(std::size_t species=0;species<speciesNames.size();++species)
    {
        if(!addChecked(fluidBase.rhoY[species],
                static_cast<long double>(dt)
                    *exchange.fluidSpeciesMassRate[species],
                result.fluidCandidate.rhoY[species])
           ||result.fluidCandidate.rhoY[species]<0.0)
            return reject("condensed-carbon exchange produces negative species");
    }
    if(!addChecked(fluidBase.rho,static_cast<long double>(dt)*fluidMassRate,
                   result.fluidCandidate.rho)
       ||result.fluidCandidate.rho<=0.0)
        return reject("condensed-carbon exchange produces invalid fluid mass");
    for(int direction=0;direction<3;++direction)
        if(!addChecked(fluidBase.rhoU[direction],
                static_cast<long double>(dt)*fluid.U[direction]*fluidMassRate,
                result.fluidCandidate.rhoU[direction]))
            return reject("condensed-carbon momentum source overflow");
    if(!addChecked(fluidBase.rhoE,
            static_cast<long double>(dt)*exchange.fluidEnergyRate,
            result.fluidCandidate.rhoE))
        return reject("condensed-carbon energy source overflow");
    if(!addChecked(fluidBase.rhoZ,
            static_cast<long double>(dt)*fluid.Z*fluidMassRate,
            result.fluidCandidate.rhoZ)
       ||!addChecked(fluidBase.rhoC,
            static_cast<long double>(dt)*fluid.C*fluidMassRate,
            result.fluidCandidate.rhoC))
        return reject("condensed-carbon scalar source overflow");
    if(!addChecked(carbonBase.rhoSoot,
            static_cast<long double>(dt)*exchange.sootMassRate,
            result.carbonCandidate.rhoSoot)
       ||!addChecked(carbonBase.numberDensity,
            static_cast<long double>(dt)*exchange.sootNumberRate,
            result.carbonCandidate.numberDensity)
       ||!result.carbonCandidate.valid(options.momentModel))
        return reject("condensed-carbon exchange produces invalid soot state");

    long double speciesMass=0.0L;
    for(double value:result.fluidCandidate.rhoY)
        speciesMass+=static_cast<long double>(value);
    const long double speciesResidual=
        speciesMass-result.fluidCandidate.rho;
    if(!closeAggregate(speciesResidual,
            std::abs(speciesMass)+std::abs(
                static_cast<long double>(result.fluidCandidate.rho)),
            absoluteMassTolerance))
        return reject("condensed-carbon candidate violates fluid mass layout");

    result.accepted=true;
    result.reason="SUCCESS";
    return result;
}

std::string condensedCarbonModelManifest(
    const CondensedCarbonOptions& options,
    const SpeciesRegistry& registry,
    const std::string& kineticsModelId,
    const std::string& wallModelId)
{
    options.validate();
    std::string manifest="schema=rfam.condensed-carbon/v1;";
    thermoProvenance::appendString(manifest,"moment.model",
        condensedCarbonMomentModelName(options.momentModel));
    thermoProvenance::appendString(manifest,"energy.convention",
        "sensible-internal-energy-explicit-chemical-heat");
    thermoProvenance::appendDouble(
        manifest,"carbon.molar-mass",options.carbonMolarMass);
    thermoProvenance::appendDouble(
        manifest,"carbon.formation-enthalpy",
        options.carbonFormationEnthalpy);
    thermoProvenance::appendDouble(
        manifest,"carbon.material-density",options.materialDensity);
    thermoProvenance::appendDouble(
        manifest,"carbon.diagnostic-cp",options.diagnosticSpecificHeat);
    thermoProvenance::appendDouble(
        manifest,"one-moment.particle-diameter",
        options.oneMomentParticleDiameter);
    thermoProvenance::appendDouble(
        manifest,"initial.rho-soot",options.initialState.rhoSoot);
    thermoProvenance::appendDouble(
        manifest,"initial.number-density",options.initialState.numberDensity);
    thermoProvenance::appendDouble(
        manifest,"inlet.rho-soot",options.inletState.rhoSoot);
    thermoProvenance::appendDouble(
        manifest,"inlet.number-density",options.inletState.numberDensity);
    thermoProvenance::appendDouble(
        manifest,"limit.mass-loading",options.limits.maximumMassLoading);
    thermoProvenance::appendDouble(
        manifest,"limit.volume-fraction",options.limits.maximumVolumeFraction);
    thermoProvenance::appendDouble(
        manifest,"limit.heat-capacity-ratio",
        options.limits.maximumHeatCapacityRatio);
    thermoProvenance::appendDouble(
        manifest,"limit.stokes",options.limits.maximumStokesNumber);
    thermoProvenance::appendSigned(
        manifest,"production-validated",options.productionValidated?1:0);
    thermoProvenance::appendSigned(
        manifest,"applicability-policy",
        options.applicabilityPolicy
            ==CondensedCarbonApplicabilityPolicy::FAIL_CLOSED?1:0);
    thermoProvenance::appendSigned(
        manifest,"molecular-diffusion",options.molecularDiffusionEnabled?1:0);
    thermoProvenance::appendSigned(
        manifest,"turbulent-diffusion",options.turbulentDiffusionEnabled?1:0);
    thermoProvenance::appendSigned(
        manifest,"thermophoresis",options.thermophoresisEnabled?1:0);
    thermoProvenance::appendString(
        manifest,"species-registry-sha256",registry.provenanceSha256());
    thermoProvenance::appendString(
        manifest,"kinetics.model-id",nonemptyModelId(kineticsModelId,"kinetics"));
    thermoProvenance::appendString(
        manifest,"wall.model-id",nonemptyModelId(wallModelId,"wall"));
    return manifest;
}

std::string condensedCarbonModelFingerprint(
    const CondensedCarbonOptions& options,
    const SpeciesRegistry& registry,
    const std::string& kineticsModelId,
    const std::string& wallModelId)
{
    return thermoProvenance::sha256Hex(condensedCarbonModelManifest(
        options,registry,kineticsModelId,wallModelId));
}

} // End namespace Foam
