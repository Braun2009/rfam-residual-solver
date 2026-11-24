#ifndef SectorRestartOwnerCodec_H
#define SectorRestartOwnerCodec_H

#include "Domain.H"
#include "SectorTurbulenceOperator.H"
#include "ThermoDataProvenance.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Foam
{

enum class SectorRestartOwnerCodecStatus : std::uint8_t
{
    SUCCESS = 0,
    INVALID_SOURCE,
    INVALID_TARGET,
    INVALID_MAGIC,
    UNSUPPORTED_VERSION,
    TRUNCATED,
    TRAILING_DATA,
    STRUCTURE_MISMATCH,
    THERMO_IDENTITY_MISMATCH,
    SPECIES_MISMATCH,
    INVALID_VALUE,
    STATE_RECOVERY_FAILED,
    PHASE_ASSERTION_MISMATCH,
    TURBULENCE_REJECTED,
    INTERNAL_ERROR
};

inline const char* sectorRestartOwnerCodecStatusName(
    SectorRestartOwnerCodecStatus status) noexcept
{
    switch (status)
    {
        case SectorRestartOwnerCodecStatus::SUCCESS: return "SUCCESS";
        case SectorRestartOwnerCodecStatus::INVALID_SOURCE:
            return "INVALID_SOURCE";
        case SectorRestartOwnerCodecStatus::INVALID_TARGET:
            return "INVALID_TARGET";
        case SectorRestartOwnerCodecStatus::INVALID_MAGIC:
            return "INVALID_MAGIC";
        case SectorRestartOwnerCodecStatus::UNSUPPORTED_VERSION:
            return "UNSUPPORTED_VERSION";
        case SectorRestartOwnerCodecStatus::TRUNCATED: return "TRUNCATED";
        case SectorRestartOwnerCodecStatus::TRAILING_DATA:
            return "TRAILING_DATA";
        case SectorRestartOwnerCodecStatus::STRUCTURE_MISMATCH:
            return "STRUCTURE_MISMATCH";
        case SectorRestartOwnerCodecStatus::THERMO_IDENTITY_MISMATCH:
            return "THERMO_IDENTITY_MISMATCH";
        case SectorRestartOwnerCodecStatus::SPECIES_MISMATCH:
            return "SPECIES_MISMATCH";
        case SectorRestartOwnerCodecStatus::INVALID_VALUE:
            return "INVALID_VALUE";
        case SectorRestartOwnerCodecStatus::STATE_RECOVERY_FAILED:
            return "STATE_RECOVERY_FAILED";
        case SectorRestartOwnerCodecStatus::PHASE_ASSERTION_MISMATCH:
            return "PHASE_ASSERTION_MISMATCH";
        case SectorRestartOwnerCodecStatus::TURBULENCE_REJECTED:
            return "TURBULENCE_REJECTED";
        case SectorRestartOwnerCodecStatus::INTERNAL_ERROR:
            return "INTERNAL_ERROR";
    }
    return "UNKNOWN";
}

struct SectorRestartEncodeResult
{
    SectorRestartOwnerCodecStatus status =
        SectorRestartOwnerCodecStatus::INVALID_SOURCE;
    std::string detail;
    std::string bytes;

    bool success() const noexcept
    {
        return status == SectorRestartOwnerCodecStatus::SUCCESS;
    }
};

struct SectorRestartFluidDecodeResult
{
    SectorRestartOwnerCodecStatus status =
        SectorRestartOwnerCodecStatus::INVALID_TARGET;
    std::string detail;
    std::unique_ptr<Domain> stagedDomain;

    bool success() const noexcept
    {
        return status == SectorRestartOwnerCodecStatus::SUCCESS
            && static_cast<bool>(stagedDomain);
    }
};

struct SectorRestartFluidRecoveryAssertion
{
    std::uint64_t cell = 0;
    thermo::EquilibriumSeed seed;
    thermo::EquilibriumRegion expectedRegion =
        thermo::EquilibriumRegion::UNKNOWN;
    std::string phaseIdentitySha256;
};

struct SectorRestartTurbulenceState
{
    std::vector<double> k;
    std::vector<double> omega;
};

struct SectorRestartTurbulenceDecodeResult
{
    SectorRestartOwnerCodecStatus status =
        SectorRestartOwnerCodecStatus::INVALID_TARGET;
    std::string detail;
    std::unique_ptr<SectorRestartTurbulenceState> stagedState;

    bool success() const noexcept
    {
        return status == SectorRestartOwnerCodecStatus::SUCCESS
            && static_cast<bool>(stagedState);
    }
};

class SectorRestartOwnerCodec
{
    static constexpr std::uint32_t fluidCodecVersion = 2U;
    static constexpr std::uint32_t turbulenceCodecVersion = 1U;
    static constexpr std::size_t maximumStringBytes = 1024U * 1024U;
    static constexpr std::uint64_t maximumCells = 1000000000ULL;
    static constexpr std::size_t sha256HexBytes = 64U;

    class Writer
    {
        std::string bytes_;

    public:
        void raw(const char* data, std::size_t count)
        { bytes_.append(data,count); }

        void u8(std::uint8_t value)
        { bytes_.push_back(static_cast<char>(value)); }

        void u32(std::uint32_t value)
        {
            for (unsigned shift = 0; shift < 32; shift += 8)
                u8(static_cast<std::uint8_t>(value >> shift));
        }

        void u64(std::uint64_t value)
        {
            for (unsigned shift = 0; shift < 64; shift += 8)
                u8(static_cast<std::uint8_t>(value >> shift));
        }

        void real(double value)
        {
            std::uint64_t bits = 0;
            static_assert(sizeof(bits) == sizeof(value),
                "restart codec requires binary64 double");
            std::memcpy(&bits,&value,sizeof(bits));
            u64(bits);
        }

        void text(const std::string& value)
        {
            u32(static_cast<std::uint32_t>(value.size()));
            raw(value.data(),value.size());
        }

        std::string finish() { return std::move(bytes_); }
    };

    class Reader
    {
        const std::string& bytes_;
        std::size_t offset_ = 0;

    public:
        explicit Reader(const std::string& bytes) : bytes_(bytes) {}

        bool raw(char* target, std::size_t count)
        {
            if (count > bytes_.size() - offset_) return false;
            std::memcpy(target,bytes_.data()+offset_,count);
            offset_ += count;
            return true;
        }

        bool u8(std::uint8_t& value)
        {
            if (offset_ >= bytes_.size()) return false;
            value = static_cast<std::uint8_t>(bytes_[offset_++]);
            return true;
        }

        bool u32(std::uint32_t& value)
        {
            value = 0;
            for (unsigned shift = 0; shift < 32; shift += 8)
            {
                std::uint8_t byte = 0;
                if (!u8(byte)) return false;
                value |= static_cast<std::uint32_t>(byte) << shift;
            }
            return true;
        }

        bool u64(std::uint64_t& value)
        {
            value = 0;
            for (unsigned shift = 0; shift < 64; shift += 8)
            {
                std::uint8_t byte = 0;
                if (!u8(byte)) return false;
                value |= static_cast<std::uint64_t>(byte) << shift;
            }
            return true;
        }

        bool real(double& value)
        {
            std::uint64_t bits = 0;
            if (!u64(bits)) return false;
            std::memcpy(&value,&bits,sizeof(value));
            return true;
        }

        bool text(std::string& value)
        {
            std::uint32_t count = 0;
            if (!u32(count) || count > maximumStringBytes
                || count > bytes_.size() - offset_)
                return false;
            value.assign(bytes_.data()+offset_,count);
            offset_ += count;
            return true;
        }

        bool atEnd() const noexcept { return offset_ == bytes_.size(); }
        std::size_t remaining() const noexcept
        { return bytes_.size() - offset_; }
    };

    static bool sameScalar(double lhs, double rhs) noexcept
    {
        std::uint64_t left = 0, right = 0;
        std::memcpy(&left,&lhs,sizeof(left));
        std::memcpy(&right,&rhs,sizeof(right));
        return left == right;
    }

    static bool sameConserved(
        const ConservedVars& lhs,
        const ConservedVars& rhs) noexcept
    {
        if (!sameScalar(lhs.rho,rhs.rho)
            || lhs.rhoU.size() != rhs.rhoU.size()
            || !sameScalar(lhs.rhoE,rhs.rhoE)
            || lhs.rhoY.size() != rhs.rhoY.size()
            || !sameScalar(lhs.rhoZ,rhs.rhoZ)
            || !sameScalar(lhs.rhoC,rhs.rhoC))
            return false;
        for (std::size_t index = 0; index < lhs.rhoU.size(); ++index)
            if (!sameScalar(lhs.rhoU[index],rhs.rhoU[index])) return false;
        for (std::size_t index = 0; index < lhs.rhoY.size(); ++index)
            if (!sameScalar(lhs.rhoY[index],rhs.rhoY[index])) return false;
        return true;
    }

    static bool samePhaseScalar(double lhs, double rhs) noexcept
    {
        constexpr double absoluteTolerance = 1.0e-12;
        constexpr double relativeTolerance = 1.0e-12;
        return std::isfinite(lhs) && std::isfinite(rhs)
            && std::abs(lhs-rhs) <= absoluteTolerance
                + relativeTolerance*std::max(std::abs(lhs),std::abs(rhs));
    }

    static bool samePhaseVector(
        const std::vector<double>& lhs,
        const std::vector<double>& rhs) noexcept
    {
        if (lhs.size() != rhs.size()) return false;
        for (std::size_t index = 0; index < lhs.size(); ++index)
            if (!samePhaseScalar(lhs[index],rhs[index])) return false;
        return true;
    }

    static bool samePhaseComposition(
        const std::optional<thermo::EquilibriumPhaseProperties>& lhs,
        const std::optional<thermo::EquilibriumPhaseProperties>& rhs) noexcept
    {
        return lhs.has_value() == rhs.has_value()
            && (!lhs.has_value()
                || (samePhaseVector(lhs->massFractions,rhs->massFractions)
                    && samePhaseVector(
                        lhs->moleFractions,rhs->moleFractions)));
    }

    static bool samePhaseAssertion(
        const thermo::EquilibriumState& lhs,
        const thermo::EquilibriumState& rhs) noexcept
    {
        if (!lhs.valid() || !rhs.valid() || lhs.region != rhs.region
            || !samePhaseVector(
                lhs.overallMassFractions,rhs.overallMassFractions)
            || !samePhaseComposition(lhs.liquid,rhs.liquid)
            || !samePhaseComposition(lhs.vapor,rhs.vapor)
            || !samePhaseComposition(lhs.singlePhase,rhs.singlePhase))
            return false;
        return lhs.region != thermo::EquilibriumRegion::TWO_PHASE
            || (samePhaseScalar(
                    lhs.vaporMoleFraction,rhs.vaporMoleFraction)
                && samePhaseScalar(
                    lhs.vaporMassFraction,rhs.vaporMassFraction)
                && samePhaseScalar(
                    lhs.vaporVolumeFraction,rhs.vaporVolumeFraction));
    }

    static bool validSha256(const std::string& value) noexcept
    {
        if (value.size() != sha256HexBytes) return false;
        for (const unsigned char character : value)
            if (!((character >= '0' && character <= '9')
                || (character >= 'a' && character <= 'f')))
                return false;
        return true;
    }

    static bool appendComposition(
        std::string& material,
        const std::string& prefix,
        const std::vector<double>& composition,
        std::size_t speciesCount)
    {
        if (composition.size() != speciesCount) return false;
        thermoProvenance::appendSize(
            material,prefix + ".count",composition.size());
        for (std::size_t species = 0; species < composition.size(); ++species)
        {
            const double value = composition[species];
            if (!std::isfinite(value) || value < 0.0) return false;
            thermoProvenance::appendDouble(
                material,prefix + "[" + std::to_string(species) + "]",value);
        }
        return true;
    }

    static bool appendPhaseCaloricState(
        std::string& material,
        const std::string& prefix,
        const std::optional<thermo::EquilibriumPhaseProperties>& phase,
        std::size_t speciesCount)
    {
        thermoProvenance::appendSize(
            material,prefix + ".present",phase.has_value() ? 1U : 0U);
        if(!phase.has_value()) return true;
        if(!phase->valid()
           ||!appendComposition(
                material,prefix + ".mass",phase->massFractions,
                speciesCount)
           ||!appendComposition(
                material,prefix + ".mole",phase->moleFractions,
                speciesCount))
            return false;
        thermoProvenance::appendDouble(
            material,prefix + ".molecular-weight",phase->molecularWeight);
        thermoProvenance::appendDouble(
            material,prefix + ".rho",phase->rho);
        thermoProvenance::appendDouble(
            material,prefix + ".h",phase->h);
        thermoProvenance::appendDouble(
            material,prefix + ".e",phase->e);
        thermoProvenance::appendDouble(
            material,prefix + ".compressibility",phase->compressibility);
        thermoProvenance::appendSize(
            material,prefix + ".entropy.present",
            phase->entropyAvailable ? 1U : 0U);
        if(phase->entropyAvailable)
            thermoProvenance::appendDouble(
                material,prefix + ".s",phase->s);
        return true;
    }

    static bool phaseIdentitySha256(
        const thermo::EquilibriumState& state,
        std::size_t speciesCount,
        std::string& fingerprint)
    {
        // The conserved field remains authoritative. This digest binds the
        // complete canonical recovery that may be exposed through a live
        // restart root; it is not an alternative state owner.
        if (!state.valid()
            || state.overallMassFractions.size() != speciesCount)
            return false;

        std::string material;
        thermoProvenance::appendString(
            material,"schema","rfam.restart-phase-caloric-identity/v1");
        thermoProvenance::appendSize(
            material,"region",static_cast<std::size_t>(state.region));
        if (!appendComposition(
                material,"overall-mass",state.overallMassFractions,
                speciesCount))
            return false;
        thermoProvenance::appendDouble(material,"mixture.p",state.p);
        thermoProvenance::appendDouble(material,"mixture.T",state.T);
        thermoProvenance::appendDouble(material,"mixture.rho",state.rho);
        thermoProvenance::appendDouble(material,"mixture.h",state.h);
        thermoProvenance::appendDouble(material,"mixture.e",state.e);
        thermoProvenance::appendSize(
            material,"mixture.entropy.present",
            state.entropyAvailable ? 1U : 0U);
        if(state.entropyAvailable)
            thermoProvenance::appendDouble(material,"mixture.s",state.s);

        if(!state.acoustic.valid()) return false;
        thermoProvenance::appendSize(
            material,"acoustic.closure",
            static_cast<std::size_t>(state.acoustic.closure));
        thermoProvenance::appendSize(
            material,"acoustic.evidence",
            static_cast<std::size_t>(state.acoustic.evidence));
        thermoProvenance::appendDouble(
            material,"acoustic.dpdrho-s",
            state.acoustic.pressureDensityDerivative);

        thermoProvenance::appendSize(
            material,"thermodynamics.present",
            state.thermodynamics.has_value() ? 1U : 0U);
        if(state.thermodynamics.has_value())
        {
            const thermo::EquilibriumThermodynamicPayload& payload=
                *state.thermodynamics;
            if(!payload.valid()) return false;
            thermoProvenance::appendSize(
                material,"thermodynamics.version",payload.version);
            thermoProvenance::appendSize(
                material,"thermodynamics.closure",
                static_cast<std::size_t>(payload.closure));
            thermoProvenance::appendDouble(
                material,"thermodynamics.Cp",payload.Cp);
            thermoProvenance::appendDouble(
                material,"thermodynamics.Cv",payload.Cv);
            thermoProvenance::appendDouble(
                material,"thermodynamics.dpdrho-T",payload.dpdrho_T);
            thermoProvenance::appendDouble(
                material,"thermodynamics.dpdT-rho",payload.dpdT_rho);
            thermoProvenance::appendDouble(
                material,"thermodynamics.dedrho-T",payload.dedrho_T);
            thermoProvenance::appendDouble(
                material,"thermodynamics.dedT-rho",payload.dedT_rho);
        }

        const bool qualityPresent =
            state.region == thermo::EquilibriumRegion::TWO_PHASE;
        thermoProvenance::appendSize(
            material,"quality.present",qualityPresent ? 1U : 0U);
        if (qualityPresent)
        {
            if (!std::isfinite(state.vaporMoleFraction)
                || !std::isfinite(state.vaporMassFraction)
                || !std::isfinite(state.vaporVolumeFraction))
                return false;
            thermoProvenance::appendDouble(
                material,"quality.vapor-mole",state.vaporMoleFraction);
            thermoProvenance::appendDouble(
                material,"quality.vapor-mass",state.vaporMassFraction);
            thermoProvenance::appendDouble(
                material,"quality.vapor-volume",state.vaporVolumeFraction);
        }
        if (!appendPhaseCaloricState(
                material,"liquid",state.liquid,speciesCount)
            || !appendPhaseCaloricState(
                material,"vapor",state.vapor,speciesCount)
            || !appendPhaseCaloricState(
                material,"single",state.singlePhase,speciesCount))
            return false;

        fingerprint = thermoProvenance::sha256Hex(material);
        return validSha256(fingerprint);
    }

    static bool validDimension(int value) noexcept
    { return value > 0; }

    static bool validLength(double value, int cells) noexcept
    {
        return std::isfinite(value)
            && (cells == 1 ? value >= 0.0 : value > 0.0)
            && !std::signbit(value);
    }

    static bool writeConserved(
        Writer& writer,
        const ConservedVars& state,
        std::size_t speciesCount)
    {
        if (state.rhoU.size() != 3U
            || state.rhoY.size() != speciesCount)
            return false;
        writer.real(state.rho);
        for (double value : state.rhoU) writer.real(value);
        writer.real(state.rhoE);
        for (double value : state.rhoY) writer.real(value);
        writer.real(state.rhoZ);
        writer.real(state.rhoC);
        return true;
    }

    static bool readConserved(
        Reader& reader,
        ConservedVars& state,
        std::size_t speciesCount)
    {
        state = ConservedVars(static_cast<int>(speciesCount));
        if (!reader.real(state.rho)) return false;
        for (double& value : state.rhoU)
            if (!reader.real(value)) return false;
        if (!reader.real(state.rhoE)) return false;
        for (double& value : state.rhoY)
            if (!reader.real(value)) return false;
        return reader.real(state.rhoZ) && reader.real(state.rhoC);
    }

    template<class Result>
    static Result failure(
        SectorRestartOwnerCodecStatus status,
        const std::string& detail)
    {
        Result result;
        result.status = status;
        result.detail = detail;
        return result;
    }

public:
    static std::optional<std::string> phaseIdentityForState(
        const thermo::EquilibriumState& state,
        std::size_t speciesCount) noexcept
    {
        try
        {
            std::string fingerprint;
            if (!phaseIdentitySha256(state,speciesCount,fingerprint))
                return std::nullopt;
            return fingerprint;
        }
        catch (...)
        {
            return std::nullopt;
        }
    }

    static SectorRestartEncodeResult encodeFluid(const Domain& domain) noexcept
    {
        try
        {
            if (domain.tripped() || !domain.stateReady()
                || domain.stepTransactionActive()
                || !validDimension(domain.Nx())
                || !validDimension(domain.Ny())
                || !validDimension(domain.Nz())
                || !validLength(domain.Lx(),domain.Nx())
                || !validLength(domain.Ly(),domain.Ny())
                || !validLength(domain.Lz(),domain.Nz())
                || domain.speciesNames().size()
                    != static_cast<std::size_t>(domain.nSpecies())
                || domain.thermoModelId().empty())
                return failure<SectorRestartEncodeResult>(
                    SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                    "fluid source is not a quiescent complete Domain");

            const std::uint64_t cellCount =
                static_cast<std::uint64_t>(domain.N());
            if (cellCount == 0 || cellCount > maximumCells)
                return failure<SectorRestartEncodeResult>(
                    SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                    "fluid cell count is outside the codec range");

            Writer writer;
            const std::array<char,8> magic =
                {{'R','F','A','M','F','L','2','\0'}};
            writer.raw(magic.data(),magic.size());
            writer.u32(fluidCodecVersion);
            writer.u32(static_cast<std::uint32_t>(domain.Nx()));
            writer.u32(static_cast<std::uint32_t>(domain.Ny()));
            writer.u32(static_cast<std::uint32_t>(domain.Nz()));
            writer.u32(static_cast<std::uint32_t>(domain.nSpecies()));
            writer.real(domain.Lx());
            writer.real(domain.Ly());
            writer.real(domain.Lz());
            writer.u64(domain.acceptedSteps());
            writer.u64(domain.rejectedSteps());
            writer.text(domain.thermoModelId());
            writer.u32(static_cast<std::uint32_t>(
                domain.speciesNames().size()));
            for (const std::string& species : domain.speciesNames())
            {
                if (species.empty() || species.size() > maximumStringBytes)
                    return failure<SectorRestartEncodeResult>(
                        SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                        "fluid source has an invalid species identity");
                writer.text(species);
            }
            writer.u64(cellCount);
            for (int k = 0; k < domain.Nz(); ++k)
                for (int j = 0; j < domain.Ny(); ++j)
                    for (int i = 0; i < domain.Nx(); ++i)
                    {
                        const ConservedVars& committed = domain.U(i,j,k);
                        const ConservedVars& candidate = domain.Unew(i,j,k);
                        if (!sameConserved(committed,candidate)
                            || !writeConserved(
                                writer,committed,
                                domain.speciesNames().size()))
                            return failure<SectorRestartEncodeResult>(
                                SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                                "fluid source is not at an accepted-step boundary");
                        const double sourceTemperature = domain.W(i,j,k).T;
                        if (!std::isfinite(sourceTemperature)
                            || sourceTemperature <= 0.0)
                            return failure<SectorRestartEncodeResult>(
                                SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                                "fluid source has an invalid temperature seed");
                        const CompleteCellConversionResult canonical =
                            CompleteStateConversion::fromConserved(
                                committed,domain.W(i,j,k),domain.thermoModel());
                        // Reject a stale phase cache before freezing the
                        // canonical rho-e-Y recovery used by the decoder.
                        if (!canonical.success()
                            || !samePhaseAssertion(
                                domain.Q(i,j,k),canonical.state->Q))
                            return failure<SectorRestartEncodeResult>(
                                SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                                "fluid source phase assertion is not the "
                                "canonical rho-e-Y recovery");
                        const double temperature = canonical.state->W.T;
                        if (!std::isfinite(temperature) || temperature <= 0.0)
                            return failure<SectorRestartEncodeResult>(
                                SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                                "canonical fluid recovery has an invalid "
                                "temperature seed");
                        writer.real(temperature);
                        std::string phaseIdentity;
                        if (!phaseIdentitySha256(
                                canonical.state->Q,
                                domain.speciesNames().size(),phaseIdentity))
                            return failure<SectorRestartEncodeResult>(
                                SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                                "fluid source has an invalid phase assertion");
                        writer.text(phaseIdentity);
                    }

            SectorRestartEncodeResult result;
            result.status = SectorRestartOwnerCodecStatus::SUCCESS;
            result.bytes = writer.finish();
            return result;
        }
        catch (const std::exception& error)
        {
            return failure<SectorRestartEncodeResult>(
                SectorRestartOwnerCodecStatus::INTERNAL_ERROR,error.what());
        }
        catch (...)
        {
            return failure<SectorRestartEncodeResult>(
                SectorRestartOwnerCodecStatus::INTERNAL_ERROR,
                "unknown fluid encoding failure");
        }
    }

    static SectorRestartFluidDecodeResult decodeFluid(
        const std::string& bytes,
        const Domain& targetTemplate,
        const std::vector<SectorRestartFluidRecoveryAssertion>&
            recoveryAssertions = {}) noexcept
    {
        try
        {
            if (targetTemplate.tripped()
                || targetTemplate.stepTransactionActive()
                || targetTemplate.thermoModelId().empty()
                || targetTemplate.speciesNames().size()
                    != static_cast<std::size_t>(targetTemplate.nSpecies()))
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::INVALID_TARGET,
                    "fluid target template is not admissible");

            Reader reader(bytes);
            std::array<char,8> magic{};
            const std::array<char,8> expectedMagic =
                {{'R','F','A','M','F','L','2','\0'}};
            if (!reader.raw(magic.data(),magic.size()))
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "fluid payload is shorter than its magic");
            if (magic != expectedMagic)
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::INVALID_MAGIC,
                    "fluid payload magic mismatch");

            std::uint32_t version = 0, nx = 0, ny = 0, nz = 0, nSpecies = 0;
            double lx = 0.0, ly = 0.0, lz = 0.0;
            std::uint64_t acceptedSteps = 0, rejectedSteps = 0;
            if (!reader.u32(version))
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "fluid payload has no version");
            if (version != fluidCodecVersion)
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::UNSUPPORTED_VERSION,
                    "fluid payload version is unsupported");
            if (!reader.u32(nx) || !reader.u32(ny) || !reader.u32(nz)
                || !reader.u32(nSpecies) || !reader.real(lx)
                || !reader.real(ly) || !reader.real(lz)
                || !reader.u64(acceptedSteps)
                || !reader.u64(rejectedSteps))
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "fluid payload header is truncated");
            if (nx != static_cast<std::uint32_t>(targetTemplate.Nx())
                || ny != static_cast<std::uint32_t>(targetTemplate.Ny())
                || nz != static_cast<std::uint32_t>(targetTemplate.Nz())
                || nSpecies
                    != static_cast<std::uint32_t>(targetTemplate.nSpecies())
                || !sameScalar(lx,targetTemplate.Lx())
                || !sameScalar(ly,targetTemplate.Ly())
                || !sameScalar(lz,targetTemplate.Lz()))
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::STRUCTURE_MISMATCH,
                    "fluid payload mesh or species dimensions mismatch");

            std::string thermoModelId;
            if (!reader.text(thermoModelId))
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "fluid thermo identity is truncated");
            if (thermoModelId != targetTemplate.thermoModelId())
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::THERMO_IDENTITY_MISMATCH,
                    "fluid thermo identity mismatch");

            std::uint32_t speciesCount = 0;
            if (!reader.u32(speciesCount))
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "fluid species roster is truncated");
            if (speciesCount != nSpecies)
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::SPECIES_MISMATCH,
                    "fluid species count mismatch");
            for (std::uint32_t index = 0; index < speciesCount; ++index)
            {
                std::string species;
                if (!reader.text(species))
                    return failure<SectorRestartFluidDecodeResult>(
                        SectorRestartOwnerCodecStatus::TRUNCATED,
                        "fluid species identity is truncated");
                if (species != targetTemplate.speciesNames()[index])
                    return failure<SectorRestartFluidDecodeResult>(
                        SectorRestartOwnerCodecStatus::SPECIES_MISMATCH,
                        "fluid species order mismatch");
            }

            std::uint64_t cellCount = 0;
            if (!reader.u64(cellCount))
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "fluid cell roster is truncated");
            if (cellCount != static_cast<std::uint64_t>(targetTemplate.N())
                || cellCount == 0 || cellCount > maximumCells)
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::STRUCTURE_MISMATCH,
                    "fluid cell count mismatch");
            if (!recoveryAssertions.empty()
                && recoveryAssertions.size()
                    != static_cast<std::size_t>(cellCount))
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::STRUCTURE_MISMATCH,
                    "fluid recovery assertion count mismatch");
            const std::uint64_t doublesPerCell =
                static_cast<std::uint64_t>(nSpecies) + 8ULL;
            const std::uint64_t minimumBytesPerCell =
                doublesPerCell * sizeof(double)
                + sizeof(std::uint32_t) + sha256HexBytes;
            if (doublesPerCell == 0 || minimumBytesPerCell == 0
                || cellCount > reader.remaining() / minimumBytesPerCell)
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "fluid cell roster exceeds the remaining payload");

            std::vector<ConservedVars> conserved;
            std::vector<double> temperatureSeeds;
            std::vector<std::string> phaseIdentities;
            conserved.reserve(static_cast<std::size_t>(cellCount));
            temperatureSeeds.reserve(static_cast<std::size_t>(cellCount));
            phaseIdentities.reserve(static_cast<std::size_t>(cellCount));
            for (std::uint64_t cell = 0; cell < cellCount; ++cell)
            {
                ConservedVars state(static_cast<int>(nSpecies));
                double temperature = 0.0;
                std::string phaseIdentity;
                if (!readConserved(reader,state,nSpecies)
                    || !reader.real(temperature)
                    || !reader.text(phaseIdentity))
                    return failure<SectorRestartFluidDecodeResult>(
                        SectorRestartOwnerCodecStatus::TRUNCATED,
                        "fluid cell payload is truncated");
                if (!std::isfinite(temperature) || temperature <= 0.0)
                    return failure<SectorRestartFluidDecodeResult>(
                        SectorRestartOwnerCodecStatus::INVALID_VALUE,
                        "fluid temperature seed is invalid");
                if (!validSha256(phaseIdentity))
                    return failure<SectorRestartFluidDecodeResult>(
                        SectorRestartOwnerCodecStatus::INVALID_VALUE,
                        "fluid phase assertion is not canonical SHA-256");
                if (!recoveryAssertions.empty())
                {
                    const SectorRestartFluidRecoveryAssertion& assertion =
                        recoveryAssertions[static_cast<std::size_t>(cell)];
                    const bool validRegion =
                        assertion.expectedRegion
                            == thermo::EquilibriumRegion::SINGLE_LIQUID
                        || assertion.expectedRegion
                            == thermo::EquilibriumRegion::SINGLE_VAPOR
                        || assertion.expectedRegion
                            == thermo::EquilibriumRegion::SUPERCRITICAL
                        || assertion.expectedRegion
                            == thermo::EquilibriumRegion::TWO_PHASE;
                    if (assertion.cell != cell
                        || !std::isfinite(assertion.seed.pressure)
                        || assertion.seed.pressure <= 0.0
                        || !std::isfinite(assertion.seed.temperature)
                        || assertion.seed.temperature <= 0.0
                        || !validRegion
                        || !validSha256(assertion.phaseIdentitySha256))
                        return failure<SectorRestartFluidDecodeResult>(
                            SectorRestartOwnerCodecStatus::INVALID_VALUE,
                            "fluid recovery assertion is invalid at cell "
                                + std::to_string(cell));
                    if (!samePhaseScalar(
                            temperature,assertion.seed.temperature)
                        || phaseIdentity != assertion.phaseIdentitySha256)
                        return failure<SectorRestartFluidDecodeResult>(
                            SectorRestartOwnerCodecStatus::
                                PHASE_ASSERTION_MISMATCH,
                            "fluid payload and authoritative recovery "
                            "assertion differ at cell "
                                + std::to_string(cell));
                }
                conserved.push_back(std::move(state));
                temperatureSeeds.push_back(temperature);
                phaseIdentities.push_back(std::move(phaseIdentity));
            }
            if (!reader.atEnd())
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRAILING_DATA,
                    "fluid payload has trailing bytes");

            auto staged = std::make_unique<Domain>(
                targetTemplate.Nx(),targetTemplate.Ny(),targetTemplate.Nz(),
                targetTemplate.Lx(),targetTemplate.Ly(),targetTemplate.Lz(),
                targetTemplate.nSpecies(),targetTemplate.thermoModelOwner());
            staged->setSpecies(targetTemplate.speciesNames());
            bool recovered = false;
            if (recoveryAssertions.empty())
                recovered = staged->replaceConservedField(
                    conserved,temperatureSeeds);
            else
            {
                std::vector<thermo::EquilibriumSeed> authoritativeSeeds;
                authoritativeSeeds.reserve(recoveryAssertions.size());
                for (const SectorRestartFluidRecoveryAssertion& assertion :
                     recoveryAssertions)
                    authoritativeSeeds.push_back(assertion.seed);
                recovered = staged->replaceConservedField(
                    conserved,authoritativeSeeds);
            }
            if (!recovered)
            {
                const Domain::FailureInfo& failureInfo = staged->failure();
                return failure<SectorRestartFluidDecodeResult>(
                    SectorRestartOwnerCodecStatus::STATE_RECOVERY_FAILED,
                    "fluid rho-e-Y recovery failed at cell "
                        + std::to_string(failureInfo.cell) + ": "
                        + failureInfo.detail);
            }
            for (int cell = 0; cell < staged->N(); ++cell)
            {
                std::string recoveredPhaseIdentity;
                if (!phaseIdentitySha256(
                        staged->Q(cell),static_cast<std::size_t>(nSpecies),
                        recoveredPhaseIdentity)
                    || recoveredPhaseIdentity
                        != phaseIdentities[static_cast<std::size_t>(cell)])
                    return failure<SectorRestartFluidDecodeResult>(
                        SectorRestartOwnerCodecStatus::PHASE_ASSERTION_MISMATCH,
                        "fluid recovered phase assertion mismatch at cell "
                            + std::to_string(cell));
                if (!recoveryAssertions.empty())
                {
                    const SectorRestartFluidRecoveryAssertion& assertion =
                        recoveryAssertions[static_cast<std::size_t>(cell)];
                    const thermo::EquilibriumState& state = staged->Q(cell);
                    if (state.region != assertion.expectedRegion
                        || !samePhaseScalar(state.p,assertion.seed.pressure)
                        || !samePhaseScalar(
                            state.T,assertion.seed.temperature)
                        || recoveredPhaseIdentity
                            != assertion.phaseIdentitySha256)
                        return failure<SectorRestartFluidDecodeResult>(
                            SectorRestartOwnerCodecStatus::
                                PHASE_ASSERTION_MISMATCH,
                            "fluid recovered state violates the authoritative "
                            "p/T/phase assertion at cell "
                                + std::to_string(cell));
                }
            }
            staged->restoreDetachedRestartStepLedger(
                acceptedSteps,rejectedSteps);

            SectorRestartFluidDecodeResult result;
            result.status = SectorRestartOwnerCodecStatus::SUCCESS;
            result.stagedDomain = std::move(staged);
            return result;
        }
        catch (const std::exception& error)
        {
            return failure<SectorRestartFluidDecodeResult>(
                SectorRestartOwnerCodecStatus::INTERNAL_ERROR,error.what());
        }
        catch (...)
        {
            return failure<SectorRestartFluidDecodeResult>(
                SectorRestartOwnerCodecStatus::INTERNAL_ERROR,
                "unknown fluid decoding failure");
        }
    }

    static SectorRestartEncodeResult encodeTurbulence(
        const Domain& domain,
        const SectorTurbulenceOperator& turbulence) noexcept
    {
        try
        {
            if (domain.tripped() || !domain.stateReady()
                || domain.stepTransactionActive()
                || turbulence.stepTransactionActive()
                || !turbulence.meshMatches(domain))
                return failure<SectorRestartEncodeResult>(
                    SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                    "turbulence source is not quiescent or mesh-compatible");
            const std::vector<double>& k = turbulence.kField();
            const std::vector<double>& omega = turbulence.omegaField();
            if (!turbulence.fieldsAreRestorable(k,omega))
                return failure<SectorRestartEncodeResult>(
                    SectorRestartOwnerCodecStatus::INVALID_SOURCE,
                    "turbulence source fields are invalid");

            Writer writer;
            const std::array<char,8> magic =
                {{'R','F','A','M','T','B','1','\0'}};
            writer.raw(magic.data(),magic.size());
            writer.u32(turbulenceCodecVersion);
            writer.u32(static_cast<std::uint32_t>(domain.Nx()));
            writer.u32(static_cast<std::uint32_t>(domain.Ny()));
            writer.u32(static_cast<std::uint32_t>(domain.Nz()));
            writer.u64(static_cast<std::uint64_t>(k.size()));
            for (double value : k) writer.real(value);
            for (double value : omega) writer.real(value);

            SectorRestartEncodeResult result;
            result.status = SectorRestartOwnerCodecStatus::SUCCESS;
            result.bytes = writer.finish();
            return result;
        }
        catch (const std::exception& error)
        {
            return failure<SectorRestartEncodeResult>(
                SectorRestartOwnerCodecStatus::INTERNAL_ERROR,error.what());
        }
        catch (...)
        {
            return failure<SectorRestartEncodeResult>(
                SectorRestartOwnerCodecStatus::INTERNAL_ERROR,
                "unknown turbulence encoding failure");
        }
    }

    static SectorRestartTurbulenceDecodeResult decodeTurbulence(
        const std::string& bytes,
        const Domain& targetDomain,
        const SectorTurbulenceOperator& targetTurbulence) noexcept
    {
        try
        {
            if (targetDomain.tripped()
                || targetDomain.stepTransactionActive()
                || targetTurbulence.stepTransactionActive()
                || !targetTurbulence.meshMatches(targetDomain))
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::INVALID_TARGET,
                    "turbulence target is not quiescent or mesh-compatible");

            Reader reader(bytes);
            std::array<char,8> magic{};
            const std::array<char,8> expectedMagic =
                {{'R','F','A','M','T','B','1','\0'}};
            if (!reader.raw(magic.data(),magic.size()))
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "turbulence payload is shorter than its magic");
            if (magic != expectedMagic)
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::INVALID_MAGIC,
                    "turbulence payload magic mismatch");
            std::uint32_t version = 0, nx = 0, ny = 0, nz = 0;
            std::uint64_t count = 0;
            if (!reader.u32(version))
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "turbulence payload has no version");
            if (version != turbulenceCodecVersion)
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::UNSUPPORTED_VERSION,
                    "turbulence payload version is unsupported");
            if (!reader.u32(nx) || !reader.u32(ny) || !reader.u32(nz)
                || !reader.u64(count))
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "turbulence payload header is truncated");
            if (nx != static_cast<std::uint32_t>(targetDomain.Nx())
                || ny != static_cast<std::uint32_t>(targetDomain.Ny())
                || nz != static_cast<std::uint32_t>(targetDomain.Nz())
                || count != targetTurbulence.expectedFieldSize()
                || count == 0 || count > maximumCells)
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::STRUCTURE_MISMATCH,
                    "turbulence payload dimensions mismatch");
            if (count > reader.remaining() / (2U * sizeof(double)))
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRUNCATED,
                    "turbulence fields exceed the remaining payload");

            auto staged = std::make_unique<SectorRestartTurbulenceState>();
            staged->k.resize(static_cast<std::size_t>(count));
            staged->omega.resize(static_cast<std::size_t>(count));
            for (double& value : staged->k)
                if (!reader.real(value))
                    return failure<SectorRestartTurbulenceDecodeResult>(
                        SectorRestartOwnerCodecStatus::TRUNCATED,
                        "turbulence k payload is truncated");
            for (double& value : staged->omega)
                if (!reader.real(value))
                    return failure<SectorRestartTurbulenceDecodeResult>(
                        SectorRestartOwnerCodecStatus::TRUNCATED,
                        "turbulence omega payload is truncated");
            if (!reader.atEnd())
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::TRAILING_DATA,
                    "turbulence payload has trailing bytes");
            if (!targetTurbulence.fieldsAreRestorable(
                    staged->k,staged->omega))
                return failure<SectorRestartTurbulenceDecodeResult>(
                    SectorRestartOwnerCodecStatus::TURBULENCE_REJECTED,
                    "turbulence target rejected the detached fields");

            SectorRestartTurbulenceDecodeResult result;
            result.status = SectorRestartOwnerCodecStatus::SUCCESS;
            result.stagedState = std::move(staged);
            return result;
        }
        catch (const std::exception& error)
        {
            return failure<SectorRestartTurbulenceDecodeResult>(
                SectorRestartOwnerCodecStatus::INTERNAL_ERROR,error.what());
        }
        catch (...)
        {
            return failure<SectorRestartTurbulenceDecodeResult>(
                SectorRestartOwnerCodecStatus::INTERNAL_ERROR,
                "unknown turbulence decoding failure");
        }
    }
};

} // End namespace Foam

#endif
