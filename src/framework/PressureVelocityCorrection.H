#ifndef PressureVelocityCorrection_H
#define PressureVelocityCorrection_H

#include "Domain.H"
#include "Operator.H"
#include <array>
#include <cmath>
#include <functional>
#include <vector>

namespace Foam
{

struct PressureVelocityCorrection
{
    std::vector<double> pressure;
    std::vector<std::array<double,3>> velocity;
    int linearSweeps=0;
    double linearResidualRatio=1.0;
    double maximumRelativePressure=0.0;
    double maximumVelocityMach=0.0;
    double minimumFacePseudoTime=0.0;
    double maximumFacePseudoTime=0.0;
    double maximumFaceReciprocityError=0.0;
    double maximumFaceNonOrthogonalityDegrees=0.0;
    int leastSquaresGradientApplications=0;
    int leastSquaresGradientFallbacks=0;
    int nonOrthogonalCorrections=0;
    double nonOrthogonalResidualRatio=1.0;
    int structuredKrylovIterations=0;
    double structuredResidualRatio=1.0;
    int structuredMomentumSolves=0;
    int structuredMomentumSweeps=0;
    double structuredMomentumMaximumResidualRatio=0.0;
    int structuredAcousticPreconditionerSolves=0;
    int structuredAcousticPreconditionerSweeps=0;
    double structuredAcousticPreconditionerMaximumResidualRatio=0.0;
    int structuredPressureDiffusionFaces=0;
    double structuredMinimumPressureDiffusionMobility=0.0;
    double structuredMaximumPressureDiffusionMobility=0.0;
    double structuredMaximumPressureDiffusionCoefficient=0.0;
    double structuredMaximumPressureDiffusionReciprocityError=0.0;
    int inletPressureDirichletFaces=0;
    int inletDirichletFaces=0;
    int inletNeumannFaces=0;
};

inline bool pressureSchurFlowResidualAdmitted(
    const std::vector<double>& base,
    const std::vector<double>& trial,
    double relativeNoise,
    double maximumMomentumComponentGrowth=0.01)
{
    if(base.size()<5||trial.size()<5
        ||!std::isfinite(relativeNoise)||relativeNoise<0.0
        ||!std::isfinite(maximumMomentumComponentGrowth)
        ||maximumMomentumComponentGrowth<0.0)
        return false;
    for(std::size_t component:{std::size_t(0),std::size_t(2),
        std::size_t(3),std::size_t(4)})
        if(!std::isfinite(base[component])||base[component]<0.0
            ||!std::isfinite(trial[component])||trial[component]<0.0)
            return false;
    // Density is the pressure-Schur continuity invariant.  A per-step noise
    // allowance creates a one-way ratchet: many individually admitted tiny
    // increases can accumulate into a material density-residual regression.
    // Keep the global residual noise margin separate and require continuity
    // itself to be non-increasing on every admitted correction.
    if(trial[0]>base[0]) return false;
    double baseSquared=base[0]*base[0];
    double trialSquared=trial[0]*trial[0];
    for(std::size_t component=2;component<=4;++component)
    {
        if(trial[component]>base[component]
            *(1.0+maximumMomentumComponentGrowth)) return false;
        baseSquared+=base[component]*base[component];
        trialSquared+=trial[component]*trial[component];
    }
    return trialSquared<baseSquared
        *(1.0-relativeNoise)*(1.0-relativeNoise);
}

using PressureVelocityCorrectionBuilder=std::function<bool(
    const Domain&,
    const std::vector<ConservedVars>&,
    const std::vector<unsigned char>&,
    const StructuredImplicitStencil&,
    double,
    PressureVelocityCorrection&)>;

} // End namespace Foam

#endif
