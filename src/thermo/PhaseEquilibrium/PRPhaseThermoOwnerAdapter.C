#include "PRPhaseThermoOwnerAdapter.H"

#include "PRPhaseAuthorityIdentity.H"
#include "ThermoDataProvenance.H"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <utility>

namespace Foam
{
namespace phaseEquilibrium
{
namespace
{

constexpr const char* energyConvention =
    "SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT";
constexpr const char* phaseCaloricSchema = "rfam.phase-caloric/v4";
constexpr double relativeAcceptance = 1.0e-8;
constexpr double compositionAcceptance = 1.0e-10;

bool close(
    double lhs,
    double rhs,
    double absoluteTolerance,
    double relativeTolerance)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= absoluteTolerance
            + relativeTolerance
                *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

double relativeResidual(double actual, double reference, double scaleFloor)
{
    if (!std::isfinite(actual) || !std::isfinite(reference)
        || !std::isfinite(scaleFloor) || scaleFloor <= 0.0)
        return std::numeric_limits<double>::infinity();
    return std::abs(actual - reference)
        /std::max({std::abs(actual), std::abs(reference), scaleFloor});
}

bool exactDouble(double lhs, double rhs)
{
    return std::memcmp(&lhs, &rhs, sizeof(double)) == 0;
}

bool exactVector(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i)
        if (!exactDouble(lhs[i], rhs[i])) return false;
    return true;
}

bool ownerIdValid(const std::string& value)
{
    if (value.empty() || value.size() > 128U) return false;
    const auto alphaNumeric = [](unsigned char character)
    {
        return (character >= 'A' && character <= 'Z')
            || (character >= 'a' && character <= 'z')
            || (character >= '0' && character <= '9');
    };
    if (!alphaNumeric(static_cast<unsigned char>(value.front()))) return false;
    for (char raw : value)
    {
        const unsigned char character = static_cast<unsigned char>(raw);
        if (!alphaNumeric(character) && character != '_'
            && character != '.' && character != ':' && character != '/'
            && character != '-')
            return false;
    }
    return true;
}

bool nonblankPrintableAscii(const std::string& value)
{
    bool nonblank = false;
    if (value.empty()) return false;
    for (char raw : value)
    {
        const unsigned char character = static_cast<unsigned char>(raw);
        if (character < 0x20U || character > 0x7eU) return false;
        if (character != ' ' && character != '\t') nonblank = true;
    }
    return nonblank;
}

bool normalized(const std::vector<double>& values)
{
    if (values.empty()) return false;
    double sum = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value) || value < 0.0 || value > 1.0)
            return false;
        sum += value;
    }
    return close(sum, 1.0, compositionAcceptance, compositionAcceptance);
}

const char* phaseIdentityName(PRPhaseIdentity identity)
{
    switch (identity)
    {
        case PRPhaseIdentity::VAPOR: return "VAPOR";
        case PRPhaseIdentity::SUPERCRITICAL: return "SUPERCRITICAL";
    }
    return nullptr;
}

const char* rootIdentityName(PRPhaseRootIdentity identity)
{
    switch (identity)
    {
        case PRPhaseRootIdentity::PR_VAPOR_ROOT: return "PR_VAPOR_ROOT";
        case PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT:
            return "PR_UNIQUE_SUPERCRITICAL_ROOT";
    }
    return nullptr;
}

std::string jsonString(const std::string& value)
{
    std::string result;
    result.reserve(value.size() + 2U);
    result.push_back('"');
    for (char character : value)
    {
        if (character == '"' || character == '\\')
            result.push_back('\\');
        result.push_back(character);
    }
    result.push_back('"');
    return result;
}

std::string f64be(double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t),
        "phase owner fingerprints require binary64");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    std::ostringstream encoded;
    encoded.imbue(std::locale::classic());
    encoded << std::hex << std::setfill('0') << std::setw(16) << bits;
    return encoded.str();
}

std::string f64beArray(const std::vector<double>& values)
{
    std::string result = "[";
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i != 0U) result.push_back(',');
        result += jsonString(f64be(values[i]));
    }
    result.push_back(']');
    return result;
}

std::string stringArray(const std::vector<std::string>& values)
{
    std::string result = "[";
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i != 0U) result.push_back(',');
        result += jsonString(values[i]);
    }
    result.push_back(']');
    return result;
}

std::string canonicalFingerprint(
    const std::string& domain,
    const std::string& canonicalJson)
{
    std::string material = domain;
    material.push_back('\0');
    material += canonicalJson;
    return thermoProvenance::sha256Hex(material);
}

std::string orderedSpeciesFingerprint(
    const std::vector<std::string>& speciesNames)
{
    if (speciesNames.empty()) return {};
    for (const std::string& name : speciesNames)
        if (!nonblankPrintableAscii(name)) return {};
    const std::string payload =
        "{\"ordered_species\":" + stringArray(speciesNames)
        + ",\"schema\":\"rfam.phase-caloric.ordered-species/v4\"}";
    return canonicalFingerprint(
        "rfam.phase-caloric.ordered-species/v4", payload);
}

std::string optionsFingerprint(const PRPhaseRhoEOptions& options)
{
    if (!std::isfinite(options.minimumTemperature)
        || !std::isfinite(options.maximumTemperature)
        || !std::isfinite(options.derivativeRelativeStep)
        || !std::isfinite(options.energyAbsoluteTolerance)
        || !std::isfinite(options.energyRelativeTolerance)
        || !std::isfinite(options.densityRelativeTolerance)
        || !std::isfinite(options.rootRelativeTolerance)
        || !std::isfinite(options.derivativeIdentityRelativeTolerance)
        || !std::isfinite(options.derivativeIdentityAbsoluteTolerance)
        || options.temperatureScanIntervals <= 0
        || options.maximumIterations <= 0)
        return {};

    const std::string payload =
        "{\"density_relative_tolerance_f64be\":"
            + jsonString(f64be(options.densityRelativeTolerance))
        + ",\"derivative_identity_absolute_tolerance_f64be\":"
            + jsonString(f64be(options.derivativeIdentityAbsoluteTolerance))
        + ",\"derivative_identity_relative_tolerance_f64be\":"
            + jsonString(f64be(options.derivativeIdentityRelativeTolerance))
        + ",\"derivative_relative_step_f64be\":"
            + jsonString(f64be(options.derivativeRelativeStep))
        + ",\"energy_absolute_tolerance_f64be\":"
            + jsonString(f64be(options.energyAbsoluteTolerance))
        + ",\"energy_relative_tolerance_f64be\":"
            + jsonString(f64be(options.energyRelativeTolerance))
        + ",\"maximum_iterations\":"
            + std::to_string(options.maximumIterations)
        + ",\"maximum_temperature_k_f64be\":"
            + jsonString(f64be(options.maximumTemperature))
        + ",\"minimum_temperature_k_f64be\":"
            + jsonString(f64be(options.minimumTemperature))
        + ",\"root_relative_tolerance_f64be\":"
            + jsonString(f64be(options.rootRelativeTolerance))
        + ",\"schema\":\"rfam.phase-caloric.options/v4\""
        + ",\"temperature_scan_intervals\":"
            + std::to_string(options.temperatureScanIntervals)
        + "}";
    return canonicalFingerprint("rfam.phase-caloric.options/v4", payload);
}

std::string committedQFingerprint(
    const CompleteCellState& committed,
    const PRPhaseThermoParentAdmission& admission,
    const std::string& speciesFingerprint)
{
    const std::string payload =
        "{\"energy_convention\":" + jsonString(energyConvention)
        + ",\"ordered_species_count\":"
            + std::to_string(committed.U.rhoY.size())
        + ",\"ordered_species_fingerprint\":"
            + jsonString(speciesFingerprint)
        + ",\"parent_owner_id\":" + jsonString(admission.parentOwnerId)
        + ",\"rho_c_kg_m3_f64be_ordered_scalars\":"
            + f64beArray({committed.U.rhoC})
        + ",\"rho_e_j_m3_f64be\":" + jsonString(f64be(committed.U.rhoE))
        + ",\"rho_kg_m3_f64be\":" + jsonString(f64be(committed.U.rho))
        + ",\"rho_u_kg_m2_s_f64be_xyz\":"
            + f64beArray(committed.U.rhoU)
        + ",\"rho_y_kg_m3_f64be_ordered_species\":"
            + f64beArray(committed.U.rhoY)
        + ",\"rho_z_kg_m3_f64be_ordered_scalars\":"
            + f64beArray({committed.U.rhoZ})
        + ",\"schema\":\"rfam.committed-parent-q/v4\"}";
    return canonicalFingerprint(
        "rfam.phase-caloric.committed-parent-q/v4", payload);
}

std::string phasePayloadFingerprint(
    const CompleteCellState& committed,
    const thermo::EquilibriumPhaseProperties& phase,
    const PRPhaseThermoParentAdmission& admission,
    const std::string& modelId,
    const std::string& speciesFingerprint,
    const std::string& phaseOptionsFingerprint)
{
    const std::string payload =
        "{\"compressibility_z_f64be\":"
            + jsonString(f64be(phase.compressibility))
        + ",\"density_kg_m3_f64be\":" + jsonString(f64be(phase.rho))
        + ",\"energy_convention\":" + jsonString(energyConvention)
        + ",\"mass_fractions_f64be_ordered_species\":"
            + f64beArray(phase.massFractions)
        + ",\"model_id\":" + jsonString(modelId)
        + ",\"ordered_species_count\":"
            + std::to_string(phase.massFractions.size())
        + ",\"ordered_species_fingerprint\":"
            + jsonString(speciesFingerprint)
        + ",\"parent_owner_id\":" + jsonString(admission.parentOwnerId)
        + ",\"parent_phase_identity\":"
            + jsonString(phaseIdentityName(admission.phaseIdentity))
        + ",\"phase_caloric_options_fingerprint\":"
            + jsonString(phaseOptionsFingerprint)
        + ",\"phase_caloric_schema_id\":"
            + jsonString(phaseCaloricSchema)
        + ",\"pressure_pa_f64be\":" + jsonString(f64be(committed.Q.p))
        + ",\"root_identity\":"
            + jsonString(rootIdentityName(admission.rootIdentity))
        + ",\"schema\":\"rfam.parent-phase-payload/v4\""
        + ",\"specific_enthalpy_j_kg_f64be\":"
            + jsonString(f64be(phase.h))
        + ",\"specific_internal_energy_j_kg_f64be\":"
            + jsonString(f64be(phase.e))
        + ",\"temperature_k_f64be\":"
            + jsonString(f64be(committed.Q.T))
        + "}";
    return canonicalFingerprint(
        "rfam.phase-caloric.parent-phase-payload/v4", payload);
}

std::string parentBindingDigest(
    const std::string& qFingerprint,
    const std::string& phaseFingerprint,
    const PRPhaseThermoParentAdmission& admission,
    const std::string& modelId,
    std::size_t speciesCount,
    const std::string& speciesFingerprint,
    const std::string& phaseOptionsFingerprint,
    const std::string& displacementAuthorityFingerprint)
{
    const std::string payload =
        "{\"committed_parent_q_fingerprint\":" + jsonString(qFingerprint)
        + ",\"energy_convention\":" + jsonString(energyConvention)
        + ",\"model_id\":" + jsonString(modelId)
        + ",\"ordered_species_count\":" + std::to_string(speciesCount)
        + ",\"ordered_species_fingerprint\":"
            + jsonString(speciesFingerprint)
        + ",\"parent_generation\":"
            + std::to_string(admission.parentGeneration)
        + ",\"parent_owner_id\":" + jsonString(admission.parentOwnerId)
        + ",\"parent_phase_identity\":"
            + jsonString(phaseIdentityName(admission.phaseIdentity))
        + ",\"parent_phase_payload_fingerprint\":"
            + jsonString(phaseFingerprint)
        + ",\"phase_caloric_options_fingerprint\":"
            + jsonString(phaseOptionsFingerprint)
        + ",\"phase_caloric_schema_id\":"
            + jsonString(phaseCaloricSchema)
        + (displacementAuthorityFingerprint.empty()
            ? std::string()
            : ",\"phase_displacement_authority_fingerprint\":"
                + jsonString(displacementAuthorityFingerprint))
        + ",\"process_runtime_epoch\":"
            + std::to_string(admission.processRuntimeEpoch)
        + ",\"root_identity\":"
            + jsonString(rootIdentityName(admission.rootIdentity))
        + "}";
    return canonicalFingerprint(
        "rfam.phase-caloric.parent-binding/v4", payload);
}

std::string transactionId(
    const std::string& bindingDigest,
    const std::string& parentOwnerId,
    std::uint64_t ordinal)
{
    const std::string payload =
        "{\"parent_binding_digest\":" + jsonString(bindingDigest)
        + ",\"parent_owner_id\":" + jsonString(parentOwnerId)
        + ",\"transaction_ordinal\":" + std::to_string(ordinal) + "}";
    return canonicalFingerprint("rfam.phase-caloric.transaction/v4", payload);
}

std::string tokenDigest(
    const std::string& bindingDigest,
    const PRPhaseThermoParentAdmission& admission,
    const std::string& transaction,
    std::uint64_t sequence)
{
    const std::string payload =
        "{\"parent_binding_digest\":" + jsonString(bindingDigest)
        + ",\"parent_generation\":"
            + std::to_string(admission.parentGeneration)
        + ",\"parent_owner_id\":" + jsonString(admission.parentOwnerId)
        + ",\"process_runtime_epoch\":"
            + std::to_string(admission.processRuntimeEpoch)
        + ",\"root_identity\":"
            + jsonString(rootIdentityName(admission.rootIdentity))
        + ",\"transaction_id\":" + jsonString(transaction)
        + ",\"transaction_sequence\":" + std::to_string(sequence) + "}";
    return canonicalFingerprint(
        "rfam.phase-caloric.live-root-token/v4", payload);
}

} // End anonymous namespace

const char* phaseThermoOwnerStatusName(PRPhaseThermoOwnerStatus status)
{
    switch (status)
    {
        case PRPhaseThermoOwnerStatus::SUCCESS: return "SUCCESS";
        case PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case PRPhaseThermoOwnerStatus::INVALID_PARENT_ADMISSION:
            return "INVALID_PARENT_ADMISSION";
        case PRPhaseThermoOwnerStatus::INVALID_PARENT_STATE_CONTRACT:
            return "INVALID_PARENT_STATE_CONTRACT";
        case PRPhaseThermoOwnerStatus::INVALID_PARENT_PHASE_ADMISSION:
            return "INVALID_PARENT_PHASE_ADMISSION";
        case PRPhaseThermoOwnerStatus::TRANSACTION_ALREADY_ACTIVE:
            return "TRANSACTION_ALREADY_ACTIVE";
        case PRPhaseThermoOwnerStatus::NO_ACTIVE_TRANSACTION:
            return "NO_ACTIVE_TRANSACTION";
        case PRPhaseThermoOwnerStatus::PARENT_BINDING_DRIFT:
            return "PARENT_BINDING_DRIFT";
        case PRPhaseThermoOwnerStatus::STALE_TOKEN: return "STALE_TOKEN";
        case PRPhaseThermoOwnerStatus::INVALID_REQUEST: return "INVALID_REQUEST";
        case PRPhaseThermoOwnerStatus::DISPLACEMENT_NOT_AUTHORIZED:
            return "DISPLACEMENT_NOT_AUTHORIZED";
        case PRPhaseThermoOwnerStatus::DISPLACEMENT_AUTHORITY_FAILED:
            return "DISPLACEMENT_AUTHORITY_FAILED";
        case PRPhaseThermoOwnerStatus::SEQUENCE_OVERFLOW:
            return "SEQUENCE_OVERFLOW";
        case PRPhaseThermoOwnerStatus::KERNEL_FAILURE: return "KERNEL_FAILURE";
        case PRPhaseThermoOwnerStatus::RECOVERY_CONTRACT_FAILED:
            return "RECOVERY_CONTRACT_FAILED";
    }
    return "INVALID_CONFIGURATION";
}

bool PRPhaseThermoOwnerResiduals::accepted() const
{
    const double* residuals[] = {
        &pressureRelative,
        &temperatureRelative,
        &compressibilityRelative,
        &densityRelative,
        &specificInternalEnergyRelative,
        &specificEnthalpyRelative,
        &massFractionL1,
        &eosRoot,
        &departureCaloric,
        &cvMinusDedTRho,
        &rootContinuity};
    for (const double* residual : residuals)
        if (!std::isfinite(*residual) || *residual < 0.0) return false;
    return pressureRelative <= relativeAcceptance
        && temperatureRelative <= relativeAcceptance
        && compressibilityRelative <= relativeAcceptance
        && densityRelative <= relativeAcceptance
        && specificInternalEnergyRelative <= relativeAcceptance
        && specificEnthalpyRelative <= relativeAcceptance
        && massFractionL1 <= compositionAcceptance
        && eosRoot <= relativeAcceptance
        && departureCaloric <= relativeAcceptance
        && cvMinusDedTRho <= relativeAcceptance
        && rootContinuity == 0.0;
}

PRPhaseThermoOwnerAdapter::PRPhaseThermoOwnerAdapter(
    const PRPhaseRhoERecovery& recovery)
:
    recovery_(recovery),
    orderedSpeciesFingerprint_(orderedSpeciesFingerprint(
        recovery.speciesNames())),
    optionsFingerprint_(optionsFingerprint(recovery.options())),
    configurationValid_(recovery.configured()
        && !orderedSpeciesFingerprint_.empty()
        && !optionsFingerprint_.empty())
{}

PRPhaseThermoOwnerAdapter::PRPhaseThermoOwnerAdapter(
    const PRPhaseRhoERecovery& recovery,
    const PRPhaseDisplacementAuthority& displacementAuthority,
    const thermo::EquilibriumThermoModel& modelAuthority)
:
    recovery_(recovery),
    displacementAuthority_(&displacementAuthority),
    modelAuthority_(&modelAuthority),
    orderedSpeciesFingerprint_(orderedSpeciesFingerprint(
        recovery.speciesNames())),
    optionsFingerprint_(optionsFingerprint(recovery.options())),
    displacementAuthorityFingerprint_(
        displacementAuthority.authorityFingerprint()),
    configurationValid_(
        recovery.configured()
        && displacementAuthority.configured()
        && !orderedSpeciesFingerprint_.empty()
        && !optionsFingerprint_.empty()
        && validPRPhaseAuthorityFingerprint(
            displacementAuthorityFingerprint_)
        && recovery.physicalModelAuthorityFingerprint()
            == displacementAuthority.physicalModelAuthorityFingerprint()
        && recovery.authorityFingerprint()
            == displacementAuthority.recoveryAuthorityFingerprint()
        && modelAuthority.speciesNames() == recovery.speciesNames()
        && modelAuthority.energyConvention()
            == EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT
        && validPRPhaseAuthorityFingerprint(
            modelAuthority.phaseLocalPhysicalModelFingerprint())
        && modelAuthority.phaseLocalPhysicalModelFingerprint()
            == recovery.physicalModelAuthorityFingerprint())
{}

PRPhaseThermoOwnerStatus PRPhaseThermoOwnerAdapter::validateParent(
    const PRPhaseThermoParentView& parent,
    ValidatedParent& validated,
    PRPhaseThermoOwnerDiagnostics& diagnostics) const
{
    validated = {};
    if (!parent.stateContract || !parent.committed || !parent.model)
        return PRPhaseThermoOwnerStatus::INVALID_PARENT_ADMISSION;
    if (modelAuthority_ && parent.model != modelAuthority_)
        return PRPhaseThermoOwnerStatus::INVALID_PARENT_ADMISSION;
    if (!ownerIdValid(parent.admission.parentOwnerId)
        || parent.admission.processRuntimeEpoch == 0U
        || parent.admission.parentGeneration == 0U
        || !phaseIdentityName(parent.admission.phaseIdentity)
        || !rootIdentityName(parent.admission.rootIdentity))
        return PRPhaseThermoOwnerStatus::INVALID_PARENT_ADMISSION;

    const ConservativeStateContract& stateContract = *parent.stateContract;
    const CompleteCellState& committed = *parent.committed;
    const thermo::EquilibriumThermoModel& model = *parent.model;
    const std::vector<std::string>& species = recovery_.speciesNames();
    if (!stateContract.hasExplicitSpeciesLayout()
        || stateContract.speciesNames() != species
        || model.speciesNames() != species
        || stateContract.energyConvention()
            != EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT
        || model.energyConvention() != stateContract.energyConvention()
        || !nonblankPrintableAscii(model.modelId())
        || committed.U.rhoY.size() != species.size()
        || !committed.valid()
        || !model.stateValidForModel(committed.Q))
        return PRPhaseThermoOwnerStatus::INVALID_PARENT_STATE_CONTRACT;
    if (modelAuthority_
        && model.phaseLocalPhysicalModelFingerprint()
            != recovery_.physicalModelAuthorityFingerprint())
        return PRPhaseThermoOwnerStatus::INVALID_PARENT_STATE_CONTRACT;

    double speciesDensity = 0.0;
    for (double rhoY : committed.U.rhoY)
    {
        if (!std::isfinite(rhoY) || rhoY < 0.0)
            return PRPhaseThermoOwnerStatus::INVALID_PARENT_STATE_CONTRACT;
        speciesDensity += rhoY;
    }
    if (!close(
            speciesDensity,
            committed.U.rho,
            compositionAcceptance,
            relativeAcceptance))
        return PRPhaseThermoOwnerStatus::INVALID_PARENT_STATE_CONTRACT;

    diagnostics.parentStateContractValidated = true;
    if (parent.admission.phaseIdentity == PRPhaseIdentity::VAPOR)
    {
        if (parent.admission.rootIdentity
                != PRPhaseRootIdentity::PR_VAPOR_ROOT)
            return PRPhaseThermoOwnerStatus::INVALID_PARENT_PHASE_ADMISSION;
        if (committed.Q.region == thermo::EquilibriumRegion::TWO_PHASE
            && committed.Q.vapor.has_value())
            validated.phase = &*committed.Q.vapor;
        else if (committed.Q.region == thermo::EquilibriumRegion::SINGLE_VAPOR
            && committed.Q.singlePhase.has_value())
            validated.phase = &*committed.Q.singlePhase;
        else
            return PRPhaseThermoOwnerStatus::INVALID_PARENT_PHASE_ADMISSION;
    }
    else if (parent.admission.phaseIdentity == PRPhaseIdentity::SUPERCRITICAL)
    {
        if (parent.admission.rootIdentity
                != PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT
            || committed.Q.region != thermo::EquilibriumRegion::SUPERCRITICAL
            || !committed.Q.singlePhase.has_value())
            return PRPhaseThermoOwnerStatus::INVALID_PARENT_PHASE_ADMISSION;
        validated.phase = &*committed.Q.singlePhase;
    }
    else
    {
        return PRPhaseThermoOwnerStatus::INVALID_PARENT_PHASE_ADMISSION;
    }

    if (!validated.phase
        || validated.phase->massFractions.size() != species.size()
        || !normalized(validated.phase->massFractions))
        return PRPhaseThermoOwnerStatus::INVALID_PARENT_PHASE_ADMISSION;
    return PRPhaseThermoOwnerStatus::SUCCESS;
}

PRPhaseThermoOwnerAdapter::ParentBinding
PRPhaseThermoOwnerAdapter::fingerprintParent(
    const PRPhaseThermoParentView& parent,
    const ValidatedParent& validated) const
{
    ParentBinding binding;
    binding.parentOwnerId = parent.admission.parentOwnerId;
    binding.modelId = parent.model->modelId();
    binding.processRuntimeEpoch = parent.admission.processRuntimeEpoch;
    binding.parentGeneration = parent.admission.parentGeneration;
    binding.phaseIdentity = parent.admission.phaseIdentity;
    binding.rootIdentity = parent.admission.rootIdentity;
    binding.displacementAuthorityFingerprint =
        displacementAuthorityFingerprint_;
    binding.committedParentQFingerprint = committedQFingerprint(
        *parent.committed, parent.admission, orderedSpeciesFingerprint_);
    binding.parentPhasePayloadFingerprint = phasePayloadFingerprint(
        *parent.committed,
        *validated.phase,
        parent.admission,
        binding.modelId,
        orderedSpeciesFingerprint_,
        optionsFingerprint_);
    binding.parentBindingDigest = parentBindingDigest(
        binding.committedParentQFingerprint,
        binding.parentPhasePayloadFingerprint,
        parent.admission,
        binding.modelId,
        recovery_.speciesNames().size(),
        orderedSpeciesFingerprint_,
        optionsFingerprint_,
        displacementAuthorityFingerprint_);
    return binding;
}

PRPhaseThermoOwnerResult PRPhaseThermoOwnerAdapter::recoverBound(
    const PRPhaseThermoParentView& parent,
    const ValidatedParent& validated,
    const ParentBinding& binding,
    const PRPhaseThermoOwnerRequest& request,
    const std::string& transaction,
    std::uint64_t nextSequence,
    PRPhaseThermoOwnerDiagnostics diagnostics) const
{
    PRPhaseThermoOwnerResult result;
    result.diagnostics = std::move(diagnostics);
    if (!std::isfinite(request.phaseDensity) || request.phaseDensity <= 0.0
        || !std::isfinite(request.phaseSpecificInternalEnergy)
        || request.phaseMassFractions.size() != recovery_.speciesNames().size()
        || !normalized(request.phaseMassFractions)
        || (request.pressureSeed.has_value()
            && (!std::isfinite(*request.pressureSeed)
                || *request.pressureSeed <= 0.0))
        || !request.temperatureSeed.has_value()
        || !std::isfinite(*request.temperatureSeed)
        || *request.temperatureSeed <= 0.0)
    {
        result.status = PRPhaseThermoOwnerStatus::INVALID_REQUEST;
        return result;
    }

    result.diagnostics.kernelInvoked = true;
    const PRPhaseRhoEResult recovered = recovery_.recover(
        request.phaseDensity,
        request.phaseSpecificInternalEnergy,
        request.phaseMassFractions,
        binding.phaseIdentity,
        binding.rootIdentity,
        *request.temperatureSeed);
    result.diagnostics.kernelStatus = recovered.status;
    if (!recovered.success())
    {
        result.status = PRPhaseThermoOwnerStatus::KERNEL_FAILURE;
        return result;
    }

    const PRPhaseCaloricProperties& state = *recovered.state;
    const PRPhaseRhoEDerivatives& derivatives = *recovered.derivatives;
    const thermo::EquilibriumPhaseProperties& reference = *validated.phase;
    result.residuals.pressureRelative = relativeResidual(
        state.pressure, parent.committed->Q.p, 1.0);
    result.residuals.temperatureRelative = relativeResidual(
        state.temperature, parent.committed->Q.T, 1.0);
    result.residuals.compressibilityRelative = relativeResidual(
        state.compressibility, reference.compressibility, 1.0);
    result.residuals.densityRelative = relativeResidual(
        state.rho, reference.rho, 1.0);
    result.residuals.specificInternalEnergyRelative = relativeResidual(
        state.e, reference.e, 1.0);
    result.residuals.specificEnthalpyRelative = relativeResidual(
        state.h, reference.h, 1.0);
    if (state.massFractions.size() == reference.massFractions.size())
    {
        result.residuals.massFractionL1 = 0.0;
        for (std::size_t i = 0; i < state.massFractions.size(); ++i)
            result.residuals.massFractionL1 += std::abs(
                state.massFractions[i] - reference.massFractions[i]);
    }
    result.residuals.eosRoot =
        recovered.diagnostics.eosRootRelativeResidual;
    const double pressureWork = state.pressure/state.rho;
    result.residuals.departureCaloric =
        std::abs(state.e - (state.h - pressureWork))
        /std::max({
            std::abs(state.e), std::abs(state.h),
            std::abs(pressureWork), 1.0});
    result.residuals.cvMinusDedTRho = relativeResidual(
        derivatives.cv, derivatives.dedT_rho, 1.0);
    result.residuals.rootContinuity =
        recovered.phaseIdentity == binding.phaseIdentity
            && recovered.rootIdentity == binding.rootIdentity
        ? 0.0 : std::numeric_limits<double>::infinity();

    if (!recovered.derivatives->valid()
        || state.root != PhaseRoot::VAPOR
        || !result.residuals.accepted())
    {
        result.status = PRPhaseThermoOwnerStatus::RECOVERY_CONTRACT_FAILED;
        return result;
    }

    result.status = PRPhaseThermoOwnerStatus::SUCCESS;
    result.phaseState = state;
    result.derivatives = derivatives;
    result.acceptedToken = issueToken(binding, transaction, nextSequence);
    result.diagnostics.resultTokenIssued = true;
    return result;
}

PRPhaseThermoLiveToken PRPhaseThermoOwnerAdapter::issueToken(
    const ParentBinding& binding,
    const std::string& transaction,
    std::uint64_t sequence) const
{
    PRPhaseThermoParentAdmission tokenAdmission;
    tokenAdmission.parentOwnerId = binding.parentOwnerId;
    tokenAdmission.processRuntimeEpoch = binding.processRuntimeEpoch;
    tokenAdmission.parentGeneration = binding.parentGeneration;
    tokenAdmission.phaseIdentity = binding.phaseIdentity;
    tokenAdmission.rootIdentity = binding.rootIdentity;

    PRPhaseThermoLiveToken token;
    token.parentOwnerId = binding.parentOwnerId;
    token.transactionId = transaction;
    token.processRuntimeEpoch = binding.processRuntimeEpoch;
    token.parentGeneration = binding.parentGeneration;
    token.transactionSequence = sequence;
    token.parentBindingDigest = binding.parentBindingDigest;
    token.rootIdentity = binding.rootIdentity;
    token.tokenDigest = tokenDigest(
        binding.parentBindingDigest,
        tokenAdmission,
        transaction,
        sequence);
    return token;
}

PRPhaseThermoOwnerResult PRPhaseThermoOwnerAdapter::recoverDisplacedBound(
    const PRPhaseThermoParentView& parent,
    const ParentBinding& binding,
    const PRPhaseThermoOwnerRequest& request,
    const std::string& transaction,
    std::uint64_t nextSequence,
    PRPhaseThermoOwnerDiagnostics diagnostics) const
{
    PRPhaseThermoOwnerResult result;
    result.diagnostics = std::move(diagnostics);
    if (!displacementAuthority_)
    {
        result.status = PRPhaseThermoOwnerStatus::DISPLACEMENT_NOT_AUTHORIZED;
        return result;
    }

    const PRPhaseRhoERequest parentRequest(
        active_->admittedDensity,
        active_->admittedSpecificInternalEnergy,
        active_->admittedMassFractions,
        active_->admittedTemperature);
    const PRPhaseRhoERequest requestedRequest(
        request.phaseDensity,
        request.phaseSpecificInternalEnergy,
        request.phaseMassFractions,
        *request.temperatureSeed);

    result.diagnostics.displacementAuthorityInvoked = true;
    const PRPhaseDisplacementResult displaced = displacementAuthority_->evaluate(
        parentRequest,
        requestedRequest,
        binding.phaseIdentity,
        binding.rootIdentity);
    result.diagnostics.displacementStatus = displaced.status;
    result.diagnostics.displacementDiagnostics = displaced.diagnostics;
    result.diagnostics.kernelInvoked = true;
    result.diagnostics.kernelStatus = displaced.diagnostics.recoveryStatus;
    if (!displaced.success() || !displaced.requested
        || !displaced.requested->state || !displaced.requested->derivatives)
    {
        result.status = PRPhaseThermoOwnerStatus::DISPLACEMENT_AUTHORITY_FAILED;
        return result;
    }

    const PRPhaseRhoEResult& recovered = *displaced.requested;
    const PRPhaseCaloricProperties& state = *recovered.state;
    if (!recovered.derivatives->valid()
        || recovered.phaseIdentity != binding.phaseIdentity
        || recovered.rootIdentity != binding.rootIdentity
        || state.root != PhaseRoot::VAPOR)
    {
        result.status = PRPhaseThermoOwnerStatus::DISPLACEMENT_AUTHORITY_FAILED;
        result.diagnostics.displacementStatus =
            PRPhaseDisplacementStatus::REQUEST_CLOSURE_FAILED;
        return result;
    }

    result.status = PRPhaseThermoOwnerStatus::SUCCESS;
    result.phaseState = state;
    result.derivatives = *recovered.derivatives;
    result.acceptedToken = issueToken(binding, transaction, nextSequence);
    result.diagnostics.resultTokenIssued = true;
    (void)parent;
    return result;
}

PRPhaseThermoOwnerResult PRPhaseThermoOwnerAdapter::reject(
    PRPhaseThermoOwnerStatus status,
    PRPhaseThermoOwnerDiagnostics diagnostics,
    bool closeTransaction)
{
    if (closeTransaction) active_.reset();
    PRPhaseThermoOwnerResult result;
    result.status = status;
    result.diagnostics = std::move(diagnostics);
    return result;
}

bool PRPhaseThermoOwnerAdapter::sameBinding(
    const ParentBinding& lhs,
    const ParentBinding& rhs)
{
    return lhs.committedParentQFingerprint
            == rhs.committedParentQFingerprint
        && lhs.parentPhasePayloadFingerprint
            == rhs.parentPhasePayloadFingerprint
        && lhs.parentBindingDigest == rhs.parentBindingDigest
        && lhs.displacementAuthorityFingerprint
            == rhs.displacementAuthorityFingerprint
        && lhs.parentOwnerId == rhs.parentOwnerId
        && lhs.modelId == rhs.modelId
        && lhs.processRuntimeEpoch == rhs.processRuntimeEpoch
        && lhs.parentGeneration == rhs.parentGeneration
        && lhs.phaseIdentity == rhs.phaseIdentity
        && lhs.rootIdentity == rhs.rootIdentity;
}

bool PRPhaseThermoOwnerAdapter::sameToken(
    const PRPhaseThermoLiveToken& token,
    const ActiveTransaction& active)
{
    PRPhaseThermoParentAdmission admission;
    admission.parentOwnerId = active.binding.parentOwnerId;
    admission.processRuntimeEpoch = active.binding.processRuntimeEpoch;
    admission.parentGeneration = active.binding.parentGeneration;
    admission.phaseIdentity = active.binding.phaseIdentity;
    admission.rootIdentity = active.binding.rootIdentity;
    const std::string expected = tokenDigest(
        active.binding.parentBindingDigest,
        admission,
        active.transactionId,
        active.transactionSequence);
    return token.parentOwnerId == active.binding.parentOwnerId
        && token.transactionId == active.transactionId
        && token.processRuntimeEpoch == active.binding.processRuntimeEpoch
        && token.parentGeneration == active.binding.parentGeneration
        && token.transactionSequence == active.transactionSequence
        && token.parentBindingDigest == active.binding.parentBindingDigest
        && token.rootIdentity == active.binding.rootIdentity
        && token.tokenDigest == active.currentTokenDigest
        && token.tokenDigest == expected;
}

PRPhaseThermoOwnerResult PRPhaseThermoOwnerAdapter::open(
    const PRPhaseThermoParentView& parent,
    double temperatureSeed)
{
    PRPhaseThermoOwnerDiagnostics diagnostics;
    if (!configurationValid_)
        return reject(
            PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION,
            std::move(diagnostics), false);
    if (active_)
        return reject(
            PRPhaseThermoOwnerStatus::TRANSACTION_ALREADY_ACTIVE,
            std::move(diagnostics), false);

    ValidatedParent validated;
    const PRPhaseThermoOwnerStatus validation = validateParent(
        parent, validated, diagnostics);
    if (validation != PRPhaseThermoOwnerStatus::SUCCESS)
        return reject(validation, std::move(diagnostics), false);

    const ParentBinding binding = fingerprintParent(parent, validated);
    diagnostics.parentFingerprintsComputed = true;
    diagnostics.committedParentQFingerprint =
        binding.committedParentQFingerprint;
    diagnostics.parentPhasePayloadFingerprint =
        binding.parentPhasePayloadFingerprint;
    diagnostics.parentBindingDigest = binding.parentBindingDigest;
    const auto ordinalIt = transactionOrdinals_.find(binding.parentOwnerId);
    const std::uint64_t currentOrdinal =
        ordinalIt == transactionOrdinals_.end() ? 0U : ordinalIt->second;
    if (currentOrdinal == std::numeric_limits<std::uint64_t>::max())
        return reject(
            PRPhaseThermoOwnerStatus::SEQUENCE_OVERFLOW,
            std::move(diagnostics), false);

    const std::uint64_t ordinal = currentOrdinal + 1U;
    const std::string transaction = transactionId(
        binding.parentBindingDigest, binding.parentOwnerId, ordinal);
    transactionOrdinals_[binding.parentOwnerId] = ordinal;

    PRPhaseThermoOwnerRequest request;
    request.phaseDensity = validated.phase->rho;
    request.phaseSpecificInternalEnergy = validated.phase->e;
    request.phaseMassFractions = validated.phase->massFractions;
    request.pressureSeed = parent.committed->Q.p;
    request.temperatureSeed = temperatureSeed;
    PRPhaseThermoOwnerResult result = recoverBound(
        parent,
        validated,
        binding,
        request,
        transaction,
        1U,
        std::move(diagnostics));
    if (!result.success()) return result;

    ActiveTransaction active;
    active.binding = binding;
    active.transactionId = transaction;
    active.transactionSequence = 1U;
    active.currentTokenDigest = result.acceptedToken->tokenDigest;
    active.admittedDensity = validated.phase->rho;
    active.admittedSpecificInternalEnergy = validated.phase->e;
    active.admittedTemperature = result.phaseState->temperature;
    active.admittedMassFractions = validated.phase->massFractions;
    active_ = std::move(active);
    return result;
}

PRPhaseThermoProbeResult PRPhaseThermoOwnerAdapter::probe(
    const PRPhaseThermoParentView& parent,
    const PRPhaseThermoOwnerRequest& request,
    const PRPhaseThermoLiveToken& token) const
{
    PRPhaseThermoOwnerDiagnostics diagnostics;
    const auto failure = [&diagnostics](PRPhaseThermoOwnerStatus status)
    {
        PRPhaseThermoProbeResult result;
        result.status = status;
        result.diagnostics = diagnostics;
        return result;
    };

    if (!configurationValid_)
        return failure(PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION);

    ValidatedParent validated;
    const PRPhaseThermoOwnerStatus validation = validateParent(
        parent, validated, diagnostics);
    if (validation != PRPhaseThermoOwnerStatus::SUCCESS)
        return failure(validation);

    const ParentBinding binding = fingerprintParent(parent, validated);
    diagnostics.parentFingerprintsComputed = true;
    diagnostics.committedParentQFingerprint =
        binding.committedParentQFingerprint;
    diagnostics.parentPhasePayloadFingerprint =
        binding.parentPhasePayloadFingerprint;
    diagnostics.parentBindingDigest = binding.parentBindingDigest;
    if (!active_)
        return failure(PRPhaseThermoOwnerStatus::NO_ACTIVE_TRANSACTION);
    if (!sameBinding(binding, active_->binding))
        return failure(PRPhaseThermoOwnerStatus::PARENT_BINDING_DRIFT);
    if (!sameToken(token, *active_))
        return failure(PRPhaseThermoOwnerStatus::STALE_TOKEN);
    diagnostics.inputTokenAdmitted = true;

    if (!std::isfinite(request.phaseDensity) || request.phaseDensity <= 0.0
        || !std::isfinite(request.phaseSpecificInternalEnergy)
        || request.phaseMassFractions.size() != recovery_.speciesNames().size()
        || !normalized(request.phaseMassFractions)
        || (request.pressureSeed.has_value()
            && (!std::isfinite(*request.pressureSeed)
                || *request.pressureSeed <= 0.0))
        || !request.temperatureSeed.has_value()
        || !std::isfinite(*request.temperatureSeed)
        || *request.temperatureSeed <= 0.0)
        return failure(PRPhaseThermoOwnerStatus::INVALID_REQUEST);

    PRPhaseThermoOwnerRequest evaluatedRequest = request;
    constexpr double probeCanonicalizationTolerance =
        64.0*std::numeric_limits<double>::epsilon();
    const double densityRelativeDifference = std::abs(
            request.phaseDensity - active_->admittedDensity)
        /std::max({
            std::abs(request.phaseDensity),
            std::abs(active_->admittedDensity), 1.0});
    const double energyRelativeDifference = std::abs(
            request.phaseSpecificInternalEnergy
                - active_->admittedSpecificInternalEnergy)
        /std::max({
            std::abs(request.phaseSpecificInternalEnergy),
            std::abs(active_->admittedSpecificInternalEnergy), 1.0});
    double compositionL1Difference = 0.0;
    for (std::size_t i = 0; i < request.phaseMassFractions.size(); ++i)
        compositionL1Difference += std::abs(
            request.phaseMassFractions[i]
                - active_->admittedMassFractions[i]);
    const bool byteExactReplay = exactDouble(
            request.phaseDensity, active_->admittedDensity)
        && exactDouble(
            request.phaseSpecificInternalEnergy,
            active_->admittedSpecificInternalEnergy)
        && exactVector(
            request.phaseMassFractions,
            active_->admittedMassFractions);
    if (!byteExactReplay
        && densityRelativeDifference <= probeCanonicalizationTolerance
        && energyRelativeDifference <= probeCanonicalizationTolerance
        && compositionL1Difference <= probeCanonicalizationTolerance)
    {
        evaluatedRequest.phaseDensity = active_->admittedDensity;
        evaluatedRequest.phaseSpecificInternalEnergy =
            active_->admittedSpecificInternalEnergy;
        evaluatedRequest.phaseMassFractions =
            active_->admittedMassFractions;
        evaluatedRequest.temperatureSeed = active_->admittedTemperature;
        diagnostics.requestCanonicalizedToAdmittedState = true;
    }

    const bool exactReplay = exactDouble(
            evaluatedRequest.phaseDensity, active_->admittedDensity)
        && exactDouble(
            evaluatedRequest.phaseSpecificInternalEnergy,
            active_->admittedSpecificInternalEnergy)
        && exactVector(
            evaluatedRequest.phaseMassFractions,
            active_->admittedMassFractions);
    if (!exactReplay && !displacementAuthority_)
        return failure(
            PRPhaseThermoOwnerStatus::DISPLACEMENT_NOT_AUTHORIZED);
    if (active_->transactionSequence
        == std::numeric_limits<std::uint64_t>::max())
        return failure(PRPhaseThermoOwnerStatus::SEQUENCE_OVERFLOW);

    PRPhaseThermoOwnerResult evaluated = displacementAuthority_
        ? recoverDisplacedBound(
            parent,
            binding,
            evaluatedRequest,
            active_->transactionId,
            active_->transactionSequence + 1U,
            std::move(diagnostics))
        : recoverBound(
            parent,
            validated,
            binding,
            evaluatedRequest,
            active_->transactionId,
            active_->transactionSequence + 1U,
            std::move(diagnostics));

    PRPhaseThermoProbeResult result;
    result.status = evaluated.status;
    result.phaseState = std::move(evaluated.phaseState);
    result.derivatives = std::move(evaluated.derivatives);
    result.residuals = evaluated.residuals;
    result.diagnostics = std::move(evaluated.diagnostics);
    // recoverBound() constructs the token that a committing call would issue.
    // A probe never exposes or admits it.
    result.diagnostics.resultTokenIssued = false;
    return result;
}

PRPhaseThermoOwnerResult PRPhaseThermoOwnerAdapter::recover(
    const PRPhaseThermoParentView& parent,
    const PRPhaseThermoOwnerRequest& request,
    const PRPhaseThermoLiveToken& token)
{
    PRPhaseThermoOwnerDiagnostics diagnostics;
    if (!configurationValid_)
        return reject(
            PRPhaseThermoOwnerStatus::INVALID_CONFIGURATION,
            std::move(diagnostics), true);

    ValidatedParent validated;
    const PRPhaseThermoOwnerStatus validation = validateParent(
        parent, validated, diagnostics);
    if (validation != PRPhaseThermoOwnerStatus::SUCCESS)
        return reject(validation, std::move(diagnostics), true);

    const ParentBinding binding = fingerprintParent(parent, validated);
    diagnostics.parentFingerprintsComputed = true;
    diagnostics.committedParentQFingerprint =
        binding.committedParentQFingerprint;
    diagnostics.parentPhasePayloadFingerprint =
        binding.parentPhasePayloadFingerprint;
    diagnostics.parentBindingDigest = binding.parentBindingDigest;
    if (!active_)
        return reject(
            PRPhaseThermoOwnerStatus::NO_ACTIVE_TRANSACTION,
            std::move(diagnostics), false);
    if (!sameBinding(binding, active_->binding))
        return reject(
            PRPhaseThermoOwnerStatus::PARENT_BINDING_DRIFT,
            std::move(diagnostics), true);
    if (!sameToken(token, *active_))
        return reject(
            PRPhaseThermoOwnerStatus::STALE_TOKEN,
            std::move(diagnostics), true);
    diagnostics.inputTokenAdmitted = true;

    if (!std::isfinite(request.phaseDensity) || request.phaseDensity <= 0.0
        || !std::isfinite(request.phaseSpecificInternalEnergy)
        || request.phaseMassFractions.size() != recovery_.speciesNames().size()
        || !normalized(request.phaseMassFractions)
        || (request.pressureSeed.has_value()
            && (!std::isfinite(*request.pressureSeed)
                || *request.pressureSeed <= 0.0))
        || !request.temperatureSeed.has_value()
        || !std::isfinite(*request.temperatureSeed)
        || *request.temperatureSeed <= 0.0)
        return reject(
            PRPhaseThermoOwnerStatus::INVALID_REQUEST,
            std::move(diagnostics), true);

    const bool exactReplay = exactDouble(
            request.phaseDensity, active_->admittedDensity)
        && exactDouble(
            request.phaseSpecificInternalEnergy,
            active_->admittedSpecificInternalEnergy)
        && exactVector(
            request.phaseMassFractions,
            active_->admittedMassFractions);
    if (!exactReplay && !displacementAuthority_)
        return reject(
            PRPhaseThermoOwnerStatus::DISPLACEMENT_NOT_AUTHORIZED,
            std::move(diagnostics), true);

    if (active_->transactionSequence
        == std::numeric_limits<std::uint64_t>::max())
        return reject(
            PRPhaseThermoOwnerStatus::SEQUENCE_OVERFLOW,
            std::move(diagnostics), false);
    const std::uint64_t nextSequence = active_->transactionSequence + 1U;
    const std::string transaction = active_->transactionId;
    PRPhaseThermoOwnerResult result = displacementAuthority_
        ? recoverDisplacedBound(
            parent,
            binding,
            request,
            transaction,
            nextSequence,
            std::move(diagnostics))
        : recoverBound(
            parent,
            validated,
            binding,
            request,
            transaction,
            nextSequence,
            std::move(diagnostics));
    if (!result.success())
    {
        active_.reset();
        return result;
    }

    active_->currentTokenDigest = result.acceptedToken->tokenDigest;
    active_->transactionSequence = nextSequence;
    if (result.diagnostics.displacementAuthorityInvoked)
    {
        active_->admittedDensity = result.phaseState->rho;
        active_->admittedSpecificInternalEnergy = result.phaseState->e;
        active_->admittedTemperature = result.phaseState->temperature;
        active_->admittedMassFractions = result.phaseState->massFractions;
    }
    return result;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
