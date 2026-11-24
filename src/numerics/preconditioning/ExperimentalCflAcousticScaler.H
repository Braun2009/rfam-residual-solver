#ifndef ExperimentalCflAcousticScaler_H
#define ExperimentalCflAcousticScaler_H

#include <vector>
#include "RealFluidState.H"

namespace Foam
{

// Experimental explicit-timestep heuristic only.
//
// This object does not transform the governing equations, flux residual, time
// derivative, or Jacobian. It must not be described as Turkel/Weiss all-Mach
// preconditioning and must not be used as evidence of low-Mach consistency.
class ExperimentalCflAcousticScaler
{
private:
    double localMach_;
    double acousticScaleSquared_;
    bool enabled_;

    static constexpr double MACH_FLOOR = 0.05;
    static constexpr double MACH_FULL_ACOUSTIC = 0.5;

    double computeAcousticScaleSquared(double mach) const;

public:
    ExperimentalCflAcousticScaler();

    void setEnabled(bool enable) { enabled_ = enable; }
    bool isEnabled() const { return enabled_; }

    void updateForMach(double mach);

    // Return a heuristic acoustic speed for an explicit CFL estimate. The
    // physical equilibrium sound speed is never modified.
    double cflAcousticSpeed(
        double physicalSoundSpeed,
        const std::vector<double>& velocity) const;

    // Transitional standalone overload for callers that do not own a Q cache.
    double cflAcousticSpeed(
        const RealFluidState& state,
        const std::vector<double>& velocity) const;

    double acousticScaleSquared() const { return acousticScaleSquared_; }
    double localMach() const { return localMach_; }

    static constexpr bool modifiesFluxResidual() { return false; }
    static constexpr bool modifiesTimeDerivative() { return false; }
    static constexpr bool modifiesJacobian() { return false; }
    static constexpr const char* capabilityLabel()
    {
        return "experimental-cfl-acoustic-scaler-only";
    }
};

} // End namespace Foam

#endif
