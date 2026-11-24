#include "ExperimentalCflAcousticScaler.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Foam
{

ExperimentalCflAcousticScaler::ExperimentalCflAcousticScaler()
:
    localMach_(1.0),
    acousticScaleSquared_(1.0),
    enabled_(false)
{
}

double ExperimentalCflAcousticScaler::computeAcousticScaleSquared(
    double mach) const
{
    if (!enabled_) return 1.0;

    const double floorSquared = MACH_FLOOR * MACH_FLOOR;

    if (mach < MACH_FLOOR)
    {
        return floorSquared;
    }
    if (mach < MACH_FULL_ACOUSTIC)
    {
        const double fraction = (mach - MACH_FLOOR)
            / (MACH_FULL_ACOUSTIC - MACH_FLOOR);
        return floorSquared * (1.0 - fraction) + fraction;
    }
    return 1.0;
}

void ExperimentalCflAcousticScaler::updateForMach(double mach)
{
    if (!std::isfinite(mach))
        throw std::invalid_argument(
            "ExperimentalCflAcousticScaler: non-finite Mach number");
    localMach_ = std::abs(mach);
    acousticScaleSquared_ = computeAcousticScaleSquared(localMach_);
}

double ExperimentalCflAcousticScaler::cflAcousticSpeed(
    double physicalSoundSpeed,
    const std::vector<double>& velocity
) const
{
    if (velocity.size() < 3)
        throw std::invalid_argument(
            "ExperimentalCflAcousticScaler: velocity requires 3 components");
    for (int component = 0; component < 3; ++component)
        if (!std::isfinite(velocity[component]))
            throw std::invalid_argument(
                "ExperimentalCflAcousticScaler: non-finite velocity");

    if (!std::isfinite(physicalSoundSpeed) || physicalSoundSpeed <= 0.0)
        throw std::invalid_argument(
            "ExperimentalCflAcousticScaler: invalid physical sound speed");
    if (!enabled_) return physicalSoundSpeed;

    const double velocityMagnitude = std::sqrt(
        velocity[0]*velocity[0] + velocity[1]*velocity[1]
        + velocity[2]*velocity[2]);
    const double acousticScale = std::sqrt(acousticScaleSquared_);
    const double scaled = acousticScale * physicalSoundSpeed
        + (1.0 - acousticScale) * velocityMagnitude;

    return std::max(scaled, velocityMagnitude + 1e-10);
}

double ExperimentalCflAcousticScaler::cflAcousticSpeed(
    const RealFluidState& state,
    const std::vector<double>& velocity
) const
{
    return cflAcousticSpeed(state.a(),velocity);
}

} // End namespace Foam
