#ifndef RestartClockAdmission_H
#define RestartClockAdmission_H

#include "RestartGenerationStore.H"
#include "SolverBase.H"

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>

namespace Foam
{

enum class RestartClockAdmissionStatus
{
    SUCCESS,
    INVALID_STORE_RESULT,
    INVALID_GENERATION_CLOCK,
    INVALID_CHECKPOINT_CLOCK,
    TIME_MISMATCH,
    ITERATION_OVERFLOW,
    ITERATION_MISMATCH,
    PREVIOUS_DT_PRESENCE_MISMATCH
};

struct RestartCheckpointClockObservation
{
    double physicalTime = 0.0;
    int acceptedIteration = 0;
    int threadCount = 1;
};

inline const char* restartClockAdmissionStatusName(
    RestartClockAdmissionStatus status) noexcept
{
    switch (status)
    {
        case RestartClockAdmissionStatus::SUCCESS: return "SUCCESS";
        case RestartClockAdmissionStatus::INVALID_STORE_RESULT:
            return "INVALID_STORE_RESULT";
        case RestartClockAdmissionStatus::INVALID_GENERATION_CLOCK:
            return "INVALID_GENERATION_CLOCK";
        case RestartClockAdmissionStatus::INVALID_CHECKPOINT_CLOCK:
            return "INVALID_CHECKPOINT_CLOCK";
        case RestartClockAdmissionStatus::TIME_MISMATCH: return "TIME_MISMATCH";
        case RestartClockAdmissionStatus::ITERATION_OVERFLOW:
            return "ITERATION_OVERFLOW";
        case RestartClockAdmissionStatus::ITERATION_MISMATCH:
            return "ITERATION_MISMATCH";
        case RestartClockAdmissionStatus::PREVIOUS_DT_PRESENCE_MISMATCH:
            return "PREVIOUS_DT_PRESENCE_MISMATCH";
    }
    return "INVALID_STORE_RESULT";
}

struct RestartClockPublication
{
    double physicalTime = 0.0;
    int acceptedIteration = 0;
    std::optional<double> previousAcceptedDt;

    bool valid() const noexcept
    {
        return std::isfinite(physicalTime) && physicalTime >= 0.0
            && !std::signbit(physicalTime) && acceptedIteration >= 0
            && (!previousAcceptedDt.has_value()
                || (std::isfinite(*previousAcceptedDt)
                    && *previousAcceptedDt > 0.0));
    }

    bool publish(SolverBase& solver) const noexcept
    {
        if (!valid()) return false;
        try
        {
            solver.setRestartState(
                physicalTime, acceptedIteration, previousAcceptedDt);
            return true;
        }
        catch (...)
        {
            return false;
        }
    }
};

struct RestartClockAdmissionResult
{
    RestartClockAdmissionStatus status =
        RestartClockAdmissionStatus::INVALID_STORE_RESULT;
    std::optional<RestartClockPublication> publication;

    bool success() const noexcept
    {
        return status == RestartClockAdmissionStatus::SUCCESS
            && publication.has_value() && publication->valid();
    }
};

class RestartClockAdmission
{
public:
    static RestartClockAdmissionResult stage(
        const RestartGenerationStoreLoadResult& load,
        const RestartCheckpointClockObservation& checkpoint,
        bool requirePreviousAcceptedDt) noexcept
    {
        RestartClockAdmissionResult result;
        if (!load.success() || !load.candidate.has_value()) return result;

        const RestartGenerationState& generation =
            load.candidate->manifest.generationState;
        if (!std::isfinite(generation.physicalTime)
            || generation.physicalTime < 0.0
            || std::signbit(generation.physicalTime)
            || (generation.previousAcceptedDtPresent
                && (!std::isfinite(generation.previousAcceptedDt)
                    || generation.previousAcceptedDt <= 0.0))
            || (!generation.previousAcceptedDtPresent
                && (generation.previousAcceptedDt != 0.0
                    || std::signbit(generation.previousAcceptedDt))))
        {
            result.status =
                RestartClockAdmissionStatus::INVALID_GENERATION_CLOCK;
            return result;
        }
        if (generation.acceptedIteration
            > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        {
            result.status = RestartClockAdmissionStatus::ITERATION_OVERFLOW;
            return result;
        }
        if (!std::isfinite(checkpoint.physicalTime)
            || checkpoint.physicalTime < 0.0
            || std::signbit(checkpoint.physicalTime)
            || checkpoint.acceptedIteration < 0
            || checkpoint.threadCount <= 0)
        {
            result.status =
                RestartClockAdmissionStatus::INVALID_CHECKPOINT_CLOCK;
            return result;
        }
        if (checkpoint.physicalTime != generation.physicalTime)
        {
            result.status = RestartClockAdmissionStatus::TIME_MISMATCH;
            return result;
        }
        if (checkpoint.acceptedIteration
            != static_cast<int>(generation.acceptedIteration))
        {
            result.status = RestartClockAdmissionStatus::ITERATION_MISMATCH;
            return result;
        }
        if (generation.previousAcceptedDtPresent
            != requirePreviousAcceptedDt)
        {
            result.status =
                RestartClockAdmissionStatus::PREVIOUS_DT_PRESENCE_MISMATCH;
            return result;
        }

        RestartClockPublication publication;
        publication.physicalTime = generation.physicalTime;
        publication.acceptedIteration =
            static_cast<int>(generation.acceptedIteration);
        if (generation.previousAcceptedDtPresent)
            publication.previousAcceptedDt = generation.previousAcceptedDt;
        if (!publication.valid())
        {
            result.status =
                RestartClockAdmissionStatus::INVALID_GENERATION_CLOCK;
            return result;
        }

        result.status = RestartClockAdmissionStatus::SUCCESS;
        result.publication = publication;
        return result;
    }
};

} // End namespace Foam

#endif
