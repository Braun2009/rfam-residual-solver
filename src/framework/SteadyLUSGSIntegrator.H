#ifndef SteadyLUSGSIntegrator_H
#define SteadyLUSGSIntegrator_H

/*---------------------------------------------------------------------------*\
    SteadyLUSGSIntegrator - robust matrix-free steady approximate factorisation.

    The nonlinear residual is assembled once from the authoritative state. A
    positive scalar approximate Jacobian is then built from operator-advertised
    local spectral bounds. Operators may supply trustworthy six-neighbour
    couplings; black-box or non-smooth operators contribute diagonal damping
    only. No finite-difference J*v evaluation is used, so iterative thermo,
    limiters, and split chemistry cannot inject differencing noise.

    For A = D + L + U with negative off-diagonal M-matrix bounds, one LU-SGS
    correction solves the approximate factorisation

        (D + L) D^-1 (D + U) delta = R(U),

    using a forward and a backward structured sweep. The steady fixed point is
    still R(U)=0. Pseudo time never advances SolverBase physical time.
\*---------------------------------------------------------------------------*/

#include "TimeIntegrator.H"
#include "AcousticStateAdmission.H"
#include "Domain.H"
#include "Operator.H"
#include "PressureVelocityCorrection.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Foam
{

class SteadyLUSGSIntegrator : public TimeIntegrator
{
    double pseudoCfl_;
    double residualTolerance_;
    int maxBacktracks_;
    bool excludeStructuredHalo_;
    double maximumRelativeUpdate_;
    double minimumRelaxation_;
    double maximumRelaxation_;
    bool useAndersonAcceleration_;
    bool useSpeciesLimiter_=true;
    int andersonDepth_;
    double andersonRegularization_;
    double andersonMaximumUpdateRatio_;
    bool useResidualLineSearch_;
    bool useBidirectionalResidualLineSearch_;
    bool useHydroDiagonalBlock_;
    bool usePressureSchurFallback_;
    bool useSegregatedAcousticFallback_;
    bool useSegregatedFlowFallback_;
    bool useSegregatedEnergyFallback_;
    std::function<bool(Domain&)> nonlinearBoundaryHook_;
    PressureVelocityCorrectionBuilder pressureCorrectionBuilder_;
    double residualEnvelopeGrowth_;
    int residualEnvelopeMemory_;
    double residualEvaluationRelativeNoise_;
    bool allowResidualBootstrapStep_;
    double energyCorrectionScale_;
    mutable double initialResidual_=0.0;
    mutable double relaxationReferenceResidual_=
        std::numeric_limits<double>::infinity();
    mutable double previousResidual_=
        std::numeric_limits<double>::infinity();
    mutable double relaxation_=1.0;
    mutable double lastResidual_=0.0;
    mutable double lastResidualRatio_=1.0;
    mutable std::vector<double> lastComponentResiduals_;
    mutable int lastMaximumResidualCell_=-1;
    mutable int lastMaximumResidualComponent_=-1;
    mutable double lastMaximumResidual_=0.0;
    mutable std::vector<double> fixedReferenceScale_;
    mutable std::vector<double> currentAndersonState_;
    mutable std::vector<double> currentAndersonUpdate_;
    mutable std::vector<std::vector<double>> andersonStateHistory_;
    mutable std::vector<std::vector<double>> andersonUpdateHistory_;
    mutable int lastAndersonDepthUsed_=0;
    mutable long lastBlockFaceApplications_=0;
    mutable long lastBlockFaceFallbacks_=0;
    mutable long lastHydroBlockApplications_=0;
    mutable long lastHydroBlockFallbacks_=0;
    mutable long lastHydroBlockCommitFallbacks_=0;
    mutable long lastHydroBlockPivotFallbacks_=0;
    mutable long lastHydroBlockOutputFallbacks_=0;
    mutable long lastHydroBlockMatrixFallbacks_=0;
    mutable long lastHydroBlockRhsFallbacks_=0;
    mutable double lastHydroBlockMaximumNormalizedOffDiagonal_=0.0;
    mutable double lastRelativeUpdate_=0.0;
    mutable double lastDiagonalMinimum_=0.0;
    mutable double lastDiagonalMaximum_=0.0;
    mutable double lastMaximumFaceCoupling_=0.0;
    mutable double lastMaximumSpeciesProjection_=0.0;
    mutable int lastProjectedCells_=0;
    mutable int lastBacktracks_=0;
    mutable int lastResidualEvaluations_=0;
    mutable double lastTrialResidual_=0.0;
    mutable double lastFullTrialResidual_=0.0;
    mutable double lastBestTrialResidual_=0.0;
    mutable double bestResidual_=std::numeric_limits<double>::infinity();
    mutable std::vector<double> residualEnvelopeHistory_;
    mutable std::vector<double> lastFullTrialComponentResiduals_;
    mutable double lastAcceptedResidualScale_=0.0;
    mutable double lastAcceptedResidualDirection_=1.0;
    mutable int lastAcceptedResidualSubspace_=0;
    mutable int lastPressureSchurResidualEvaluations_=0;
    mutable double lastBestPressureSchurTrialResidual_=
        std::numeric_limits<double>::infinity();
    mutable std::vector<double> lastBestPressureSchurComponentResiduals_;
    mutable std::vector<double> lastAcceptedPressureSchurComponentResiduals_;
    mutable int lastPressureSchurLinearSweeps_=0;
    mutable double lastPressureSchurLinearResidualRatio_=1.0;
    mutable double lastPressureSchurMaximumRelativePressure_=0.0;
    mutable double lastPressureSchurMaximumVelocityMach_=0.0;
    mutable double lastPressureSchurMinimumFacePseudoTime_=0.0;
    mutable double lastPressureSchurMaximumFacePseudoTime_=0.0;
    mutable double lastPressureSchurMaximumFaceReciprocityError_=0.0;
    mutable double lastPressureSchurMaximumFaceNonOrthogonalityDegrees_=0.0;
    mutable int lastPressureSchurLeastSquaresGradientApplications_=0;
    mutable int lastPressureSchurLeastSquaresGradientFallbacks_=0;
    mutable int lastPressureSchurNonOrthogonalCorrections_=0;
    mutable double lastPressureSchurNonOrthogonalResidualRatio_=1.0;
    mutable int lastPressureSchurStructuredKrylovIterations_=0;
    mutable double lastPressureSchurStructuredResidualRatio_=1.0;
    mutable int lastPressureSchurStructuredMomentumSolves_=0;
    mutable int lastPressureSchurStructuredMomentumSweeps_=0;
    mutable double lastPressureSchurStructuredMomentumMaximumResidualRatio_=0.0;
    mutable int lastPressureSchurStructuredAcousticPreconditionerSolves_=0;
    mutable int lastPressureSchurStructuredAcousticPreconditionerSweeps_=0;
    mutable double lastPressureSchurStructuredAcousticPreconditionerMaximumResidualRatio_=0.0;
    mutable int lastPressureSchurStructuredPressureDiffusionFaces_=0;
    mutable double lastPressureSchurStructuredMinimumPressureDiffusionMobility_=0.0;
    mutable double lastPressureSchurStructuredMaximumPressureDiffusionMobility_=0.0;
    mutable double lastPressureSchurStructuredMaximumPressureDiffusionCoefficient_=0.0;
    mutable double lastPressureSchurStructuredMaximumPressureDiffusionReciprocityError_=0.0;
    mutable int lastPressureSchurInletPressureDirichletFaces_=0;
    mutable int lastPressureSchurInletDirichletFaces_=0;
    mutable int lastPressureSchurInletNeumannFaces_=0;
    mutable int lastAcousticResidualEvaluations_=0;
    mutable double lastBestAcousticTrialResidual_=
        std::numeric_limits<double>::infinity();
    mutable int lastFlowResidualEvaluations_=0;
    mutable double lastBestFlowTrialResidual_=
        std::numeric_limits<double>::infinity();
    mutable int lastEnergyResidualEvaluations_=0;
    mutable double lastBestEnergyTrialResidual_=
        std::numeric_limits<double>::infinity();
    mutable double lastResidualAcceptanceLimit_=0.0;
    mutable int acceptedResidualSafeguardSteps_=0;
    mutable bool converged_=false;

    static bool solveDenseSystem(
        std::vector<std::vector<double>>& matrix,
        std::vector<double>& rightHandSide)
    {
        const int count=static_cast<int>(rightHandSide.size());
        for(int column=0;column<count;++column)
        {
            int pivot=column;
            for(int row=column+1;row<count;++row)
                if(std::abs(matrix[row][column])
                    >std::abs(matrix[pivot][column])) pivot=row;
            if(std::abs(matrix[pivot][column])<1.0e-24) return false;
            if(pivot!=column)
            {
                std::swap(matrix[pivot],matrix[column]);
                std::swap(rightHandSide[pivot],rightHandSide[column]);
            }
            const double diagonal=matrix[column][column];
            for(int entry=column;entry<count;++entry)
                matrix[column][entry]/=diagonal;
            rightHandSide[column]/=diagonal;
            for(int row=0;row<count;++row)
            {
                if(row==column) continue;
                const double factor=matrix[row][column];
                for(int entry=column;entry<count;++entry)
                    matrix[row][entry]-=factor*matrix[column][entry];
                rightHandSide[row]-=factor*rightHandSide[column];
            }
        }
        return true;
    }

    static void assignScaledDifference(
        ConservedVars& result,
        const ConservedVars& candidate,
        const ConservedVars& current,
        double inverseStep)
    {
        result.rho=(candidate.rho-current.rho)*inverseStep;
        result.rhoE=(candidate.rhoE-current.rhoE)*inverseStep;
        for(std::size_t component=0;component<result.rhoU.size();++component)
            result.rhoU[component]=
                (candidate.rhoU[component]-current.rhoU[component])*inverseStep;
        for(std::size_t component=0;component<result.rhoY.size();++component)
            result.rhoY[component]=
                (candidate.rhoY[component]-current.rhoY[component])*inverseStep;
        result.rhoZ=(candidate.rhoZ-current.rhoZ)*inverseStep;
        result.rhoC=(candidate.rhoC-current.rhoC)*inverseStep;
    }

    static void divide(ConservedVars& value,double denominator)
    {
        value.rho/=denominator;
        value.rhoE/=denominator;
        for(double& component:value.rhoU) component/=denominator;
        for(double& component:value.rhoY) component/=denominator;
        value.rhoZ/=denominator;
        value.rhoC/=denominator;
    }

    static std::array<double,StructuredImplicitBlock::HYDRO_COMPONENTS>
    hydroScales(const ConservedVars& state,double acousticScale)
    {
        const double rhoScale=std::max(std::abs(state.rho),1.0e-12);
        const double momentumScale=rhoScale*std::max(acousticScale,1.0);
        const double energyScale=std::max(std::abs(state.rhoE),
            rhoScale*std::max(acousticScale,1.0)
                *std::max(acousticScale,1.0));
        return {rhoScale,momentumScale,momentumScale,momentumScale,
            energyScale};
    }

    // Solve only the hydro block after row/column normalization. This keeps
    // the matrix numerically well scaled across rho, momentum and energy, and
    // makes singular/indefinite blocks fail closed to the scalar LU-SGS path.
    static bool divideWithHydroBlock(
        ConservedVars& value,const StructuredImplicitBlock& block,
        double scalarDiagonal,
        const std::array<double,StructuredImplicitBlock::HYDRO_COMPONENTS>& scales,
        int* failureReason=nullptr)
    {
        auto fail=[&](int reason)
        {
            if(failureReason) *failureReason=reason;
            return false;
        };
        if(!block.finite()||!std::isfinite(scalarDiagonal)
            ||scalarDiagonal<=0.0) return fail(1);
        std::array<std::array<double,StructuredImplicitBlock::HYDRO_COMPONENTS>,
            StructuredImplicitBlock::HYDRO_COMPONENTS> matrix{};
        std::array<double,StructuredImplicitBlock::HYDRO_COMPONENTS> rhs{};
        double scalarNorm=0.0;
        for(int row=0;row<StructuredImplicitBlock::HYDRO_COMPONENTS;++row)
        {
            const double rowScale=scales[static_cast<std::size_t>(row)];
            if(!std::isfinite(rowScale)||rowScale<=0.0) return fail(1);
            rhs[static_cast<std::size_t>(row)]=
                (row==0?value.rho:(row<4
                    ?value.rhoU[static_cast<std::size_t>(row-1)]
                    :value.rhoE))/rowScale;
            if(!std::isfinite(rhs[static_cast<std::size_t>(row)])) return fail(4);
            scalarNorm=std::max(scalarNorm,
                std::abs(rhs[static_cast<std::size_t>(row)])
                    /scalarDiagonal);
            for(int column=0;
                column<StructuredImplicitBlock::HYDRO_COMPONENTS;++column)
            {
                const double columnScale=scales[
                    static_cast<std::size_t>(column)];
                const double scaleRatio=columnScale/rowScale;
                if(!std::isfinite(columnScale)||columnScale<=0.0
                    ||!std::isfinite(scaleRatio)) return fail(1);
                matrix[static_cast<std::size_t>(row)][
                    static_cast<std::size_t>(column)]=
                    block.value[row][column]*scaleRatio;
                if(!std::isfinite(matrix[static_cast<std::size_t>(row)][
                        static_cast<std::size_t>(column)])) return fail(1);
            }
        }
        for(int column=0;column<StructuredImplicitBlock::HYDRO_COMPONENTS;
            ++column)
        {
            int pivot=column;
            for(int row=column+1;
                row<StructuredImplicitBlock::HYDRO_COMPONENTS;++row)
                if(std::abs(matrix[static_cast<std::size_t>(row)][column])
                    >std::abs(matrix[static_cast<std::size_t>(pivot)][column]))
                    pivot=row;
            if(!std::isfinite(matrix[static_cast<std::size_t>(pivot)][column])
                ||std::abs(matrix[static_cast<std::size_t>(pivot)][column])
                    <1.0e-12) return fail(2);
            if(pivot!=column)
            {
                std::swap(matrix[static_cast<std::size_t>(pivot)],
                    matrix[static_cast<std::size_t>(column)]);
                std::swap(rhs[static_cast<std::size_t>(pivot)],
                    rhs[static_cast<std::size_t>(column)]);
            }
            const double diagonal=matrix[static_cast<std::size_t>(column)][column];
            for(int entry=column;entry<StructuredImplicitBlock::HYDRO_COMPONENTS;
                ++entry)
                matrix[static_cast<std::size_t>(column)][entry]/=diagonal;
            rhs[static_cast<std::size_t>(column)]/=diagonal;
            if(!std::isfinite(rhs[static_cast<std::size_t>(column)]))
                return fail(3);
            for(int row=0;row<StructuredImplicitBlock::HYDRO_COMPONENTS;++row)
            {
                if(row==column) continue;
                const double factor=matrix[static_cast<std::size_t>(row)][column];
                for(int entry=column;
                    entry<StructuredImplicitBlock::HYDRO_COMPONENTS;++entry)
                    matrix[static_cast<std::size_t>(row)][entry]-=factor
                        *matrix[static_cast<std::size_t>(column)][entry];
                rhs[static_cast<std::size_t>(row)]-=factor
                    *rhs[static_cast<std::size_t>(column)];
                if(!std::isfinite(rhs[static_cast<std::size_t>(row)]))
                    return fail(3);
            }
        }
        ConservedVars solved(value.rhoY.size());
        solved.rho=rhs[0]*scales[0];
        for(int direction=0;direction<3;++direction)
            solved.rhoU[static_cast<std::size_t>(direction)]=
                rhs[static_cast<std::size_t>(direction+1)]
                    *scales[static_cast<std::size_t>(direction+1)];
        solved.rhoE=rhs[4]*scales[4];
        for(std::size_t component=0;component<solved.rhoY.size();++component)
            solved.rhoY[component]=value.rhoY[component]/scalarDiagonal;
        solved.rhoZ=value.rhoZ/scalarDiagonal;
        solved.rhoC=value.rhoC/scalarDiagonal;
        double solvedNorm=0.0;
        solvedNorm=std::max(solvedNorm,
            std::abs(solved.rho)/scales[0]);
        for(int direction=0;direction<3;++direction)
            solvedNorm=std::max(solvedNorm,
                std::abs(solved.rhoU[static_cast<std::size_t>(direction)])
                    /scales[static_cast<std::size_t>(direction+1)]);
        solvedNorm=std::max(solvedNorm,std::abs(solved.rhoE)/scales[4]);
        if(!std::isfinite(solvedNorm)
            ||solvedNorm>8.0*std::max(scalarNorm,1.0e-12)) return fail(3);
        if(!std::isfinite(solved.rho)||!std::isfinite(solved.rhoE)
            ||!std::all_of(solved.rhoU.begin(),solved.rhoU.end(),
                [](double v){return std::isfinite(v);} )) return fail(3);
        value=std::move(solved);
        if(failureReason) *failureReason=0;
        return true;
    }

    static void addScaled(
        ConservedVars& target,const ConservedVars& source,double scale)
    {
        target.rho+=scale*source.rho;
        target.rhoE+=scale*source.rhoE;
        for(std::size_t component=0;component<target.rhoU.size();++component)
            target.rhoU[component]+=scale*source.rhoU[component];
        for(std::size_t component=0;component<target.rhoY.size();++component)
            target.rhoY[component]+=scale*source.rhoY[component];
        target.rhoZ+=scale*source.rhoZ;
        target.rhoC+=scale*source.rhoC;
    }

    static void assignCandidate(
        ConservedVars& candidate,const ConservedVars& current,
        const ConservedVars& correction,double scale)
    {
        candidate.rho=current.rho+scale*correction.rho;
        candidate.rhoE=current.rhoE+scale*correction.rhoE;
        for(std::size_t component=0;component<candidate.rhoU.size();++component)
            candidate.rhoU[component]=current.rhoU[component]
                +scale*correction.rhoU[component];
        for(std::size_t component=0;component<candidate.rhoY.size();++component)
            candidate.rhoY[component]=current.rhoY[component]
                +scale*correction.rhoY[component];
        candidate.rhoZ=current.rhoZ+scale*correction.rhoZ;
        candidate.rhoC=current.rhoC+scale*correction.rhoC;
    }

    // Positivity-preserving composition limiter. Split the species candidate
    // into a density-following part Y^n*rho_candidate and a zero-sum
    // composition correction. Scale only that composition correction by the
    // largest beta in [0,1] that keeps every rhoY nonnegative. Unlike an L2
    // simplex projection, this does not redistribute mass among unrelated
    // species; it preserves the LU-SGS composition direction whenever it is
    // admissible and collapses to frozen local composition only where an
    // absent species would otherwise become negative.
    static bool limitSpeciesComposition(
        const ConservedVars& current,
        ConservedVars& candidate,
        double& maximumRelativeChange)
    {
        maximumRelativeChange=0.0;
        if(!std::isfinite(candidate.rho)||candidate.rho<=0.0
            ||candidate.rhoY.empty()
            ||candidate.rhoY.size()!=current.rhoY.size()) return false;
        double currentSpeciesMass=0.0,candidateSpeciesMass=0.0;
        bool nonnegative=true;
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
        {
            if(!std::isfinite(current.rhoY[species])
                ||current.rhoY[species]<0.0
                ||!std::isfinite(candidate.rhoY[species])) return false;
            currentSpeciesMass+=current.rhoY[species];
            candidateSpeciesMass+=candidate.rhoY[species];
            if(candidate.rhoY[species]<0.0) nonnegative=false;
        }
        if(!std::isfinite(currentSpeciesMass)||currentSpeciesMass<=0.0
            ||!std::isfinite(candidateSpeciesMass)) return false;
        if(nonnegative
            &&std::abs(candidateSpeciesMass-candidate.rho)
                <=1.0e-12*candidate.rho)
        {
            const auto largest=std::max_element(
                candidate.rhoY.begin(),candidate.rhoY.end());
            if(largest==candidate.rhoY.end()) return false;
            *largest+=candidate.rho-candidateSpeciesMass;
            return std::isfinite(*largest)&&*largest>=0.0;
        }

        std::vector<double> base(candidate.rhoY.size(),0.0);
        std::vector<double> composition(candidate.rhoY.size(),0.0);
        double compositionSum=0.0;
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
        {
            const double massFraction=
                current.rhoY[species]/currentSpeciesMass;
            base[species]=massFraction*candidate.rho;
            composition[species]=candidate.rhoY[species]-base[species];
            compositionSum+=composition[species];
        }
        // Remove roundoff/nonconservative mass from the composition direction
        // along the current composition, leaving sum(composition)=0.
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
            composition[species]-=
                (current.rhoY[species]/currentSpeciesMass)*compositionSum;

        double beta=1.0;
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
            if(composition[species]<0.0)
                beta=std::min(beta,
                    base[species]/(-composition[species]));
        beta=std::clamp(beta,0.0,1.0);

        double limitedSum=0.0;
        std::size_t largest=0;
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
        {
            const double previous=candidate.rhoY[species];
            const double limited=std::max(
                base[species]+beta*composition[species],0.0);
            candidate.rhoY[species]=limited;
            limitedSum+=limited;
            if(limited>candidate.rhoY[largest]) largest=species;
            maximumRelativeChange=std::max(maximumRelativeChange,
                std::abs(limited-previous)/candidate.rho);
        }
        candidate.rhoY[largest]+=candidate.rho-limitedSum;
        return candidate.rhoY[largest]>=0.0
            &&std::isfinite(candidate.rhoY[largest]);
    }

public:
    // Diagnostics: disable the positivity composition limiter to probe
    // non-smoothness-induced residual plateaus.
    void setUseSpeciesLimiter(bool value){useSpeciesLimiter_=value;}

    explicit SteadyLUSGSIntegrator(
        double pseudoCfl=5.0,
        double residualTolerance=0.05,
        int maxBacktracks=24,
        bool excludeStructuredHalo=false,
        double maximumRelativeUpdate=0.05,
        double minimumRelaxation=1.0e-4,
        double maximumRelaxation=1.0,
        bool useAndersonAcceleration=false,
        int andersonDepth=4,
        double andersonRegularization=1.0e-10,
        double andersonMaximumUpdateRatio=2.0,
        bool useResidualLineSearch=false,
        std::function<bool(Domain&)> nonlinearBoundaryHook={},
        double energyCorrectionScale=1.0,
        double residualEnvelopeGrowth=1.0,
        bool allowResidualBootstrapStep=true,
        double residualEvaluationRelativeNoise=1.0e-5,
        int residualEnvelopeMemory=20,
        bool useBidirectionalResidualLineSearch=false,
        bool useHydroDiagonalBlock=false,
        bool useSegregatedEnergyFallback=false,
        bool useSegregatedFlowFallback=false,
        bool useSegregatedAcousticFallback=false,
        PressureVelocityCorrectionBuilder pressureCorrectionBuilder={})
        :pseudoCfl_(pseudoCfl),
         residualTolerance_(residualTolerance),
         maxBacktracks_(maxBacktracks),
         excludeStructuredHalo_(excludeStructuredHalo),
         maximumRelativeUpdate_(maximumRelativeUpdate),
         minimumRelaxation_(minimumRelaxation),
         maximumRelaxation_(maximumRelaxation),
         useAndersonAcceleration_(useAndersonAcceleration),
         andersonDepth_(std::max(1,andersonDepth)),
         andersonRegularization_(andersonRegularization),
         andersonMaximumUpdateRatio_(andersonMaximumUpdateRatio),
         useResidualLineSearch_(useResidualLineSearch),
         useBidirectionalResidualLineSearch_(
             useBidirectionalResidualLineSearch),
         useHydroDiagonalBlock_(useHydroDiagonalBlock),
         usePressureSchurFallback_(
             static_cast<bool>(pressureCorrectionBuilder)),
         useSegregatedAcousticFallback_(useSegregatedAcousticFallback),
         useSegregatedFlowFallback_(useSegregatedFlowFallback),
         useSegregatedEnergyFallback_(useSegregatedEnergyFallback),
         nonlinearBoundaryHook_(std::move(nonlinearBoundaryHook)),
         pressureCorrectionBuilder_(std::move(pressureCorrectionBuilder)),
         residualEnvelopeGrowth_(residualEnvelopeGrowth),
         residualEnvelopeMemory_(std::max(1,residualEnvelopeMemory)),
         residualEvaluationRelativeNoise_(
             residualEvaluationRelativeNoise),
         allowResidualBootstrapStep_(allowResidualBootstrapStep),
         energyCorrectionScale_(energyCorrectionScale)
    {
        if(!std::isfinite(pseudoCfl_)||pseudoCfl_<=0.0)
            throw std::invalid_argument("LU-SGS pseudo-CFL must be positive");
        if(!std::isfinite(residualTolerance_)
            ||residualTolerance_<=0.0||residualTolerance_>=1.0)
            throw std::invalid_argument(
                "LU-SGS residual tolerance must lie in (0,1)");
        if(maxBacktracks_<0)
            throw std::invalid_argument(
                "LU-SGS backtrack count must be nonnegative");
        if(!std::isfinite(maximumRelativeUpdate_)
            ||maximumRelativeUpdate_<=0.0)
            throw std::invalid_argument(
                "LU-SGS maximum relative update must be positive");
        if(!std::isfinite(minimumRelaxation_)
            ||!std::isfinite(maximumRelaxation_)
            ||minimumRelaxation_<=0.0
            ||maximumRelaxation_<minimumRelaxation_)
            throw std::invalid_argument(
                "LU-SGS relaxation bounds are invalid");
        if(!std::isfinite(andersonRegularization_)
            ||andersonRegularization_<0.0
            ||!std::isfinite(andersonMaximumUpdateRatio_)
            ||andersonMaximumUpdateRatio_<1.0)
            throw std::invalid_argument(
                "LU-SGS Anderson safeguards are invalid");
        if(useResidualLineSearch_&&!nonlinearBoundaryHook_)
            throw std::invalid_argument(
                "LU-SGS residual line search requires a boundary hook");
        if(useBidirectionalResidualLineSearch_&&!useResidualLineSearch_)
            throw std::invalid_argument(
                "bidirectional LU-SGS requires residual line search");
        if((usePressureSchurFallback_||useSegregatedAcousticFallback_
            ||useSegregatedFlowFallback_||useSegregatedEnergyFallback_)
            &&!useBidirectionalResidualLineSearch_)
            throw std::invalid_argument(
                "segregated fallback requires bidirectional line search");
        if(!std::isfinite(energyCorrectionScale_)
            ||energyCorrectionScale_<=0.0
            ||energyCorrectionScale_>1.0)
            throw std::invalid_argument(
                "LU-SGS energy correction scale must lie in (0,1]");
        if(!std::isfinite(residualEnvelopeGrowth_)
            ||residualEnvelopeGrowth_<1.0)
            throw std::invalid_argument(
                "LU-SGS residual envelope growth must be at least one");
        if(!std::isfinite(residualEvaluationRelativeNoise_)
            ||residualEvaluationRelativeNoise_<0.0
            ||residualEvaluationRelativeNoise_>1.0e-2)
            throw std::invalid_argument(
                "LU-SGS residual evaluation noise tolerance is invalid");
    }

    bool advancesPhysicalTime() const override { return false; }
    bool converged() const { return converged_; }
    double initialResidual() const { return initialResidual_; }
    double lastResidual() const { return lastResidual_; }
    double lastResidualRatio() const { return lastResidualRatio_; }
    const std::vector<double>& lastComponentResiduals() const
    { return lastComponentResiduals_; }
    int lastMaximumResidualCell() const
    { return lastMaximumResidualCell_; }
    int lastMaximumResidualComponent() const
    { return lastMaximumResidualComponent_; }
    double lastMaximumResidual() const
    { return lastMaximumResidual_; }
    int lastAndersonDepthUsed() const
    { return lastAndersonDepthUsed_; }
    long lastBlockFaceApplications() const
    { return lastBlockFaceApplications_; }
    long lastBlockFaceFallbacks() const
    { return lastBlockFaceFallbacks_; }
    long lastHydroBlockApplications() const
    { return lastHydroBlockApplications_; }
    long lastHydroBlockFallbacks() const
    { return lastHydroBlockFallbacks_; }
    long lastHydroBlockCommitFallbacks() const
    { return lastHydroBlockCommitFallbacks_; }
    long lastHydroBlockPivotFallbacks() const
    { return lastHydroBlockPivotFallbacks_; }
    long lastHydroBlockOutputFallbacks() const
    { return lastHydroBlockOutputFallbacks_; }
    long lastHydroBlockMatrixFallbacks() const
    { return lastHydroBlockMatrixFallbacks_; }
    long lastHydroBlockRhsFallbacks() const
    { return lastHydroBlockRhsFallbacks_; }
    double lastHydroBlockMaximumNormalizedOffDiagonal() const
    { return lastHydroBlockMaximumNormalizedOffDiagonal_; }
    double lastRelativeUpdate() const { return lastRelativeUpdate_; }
    double relaxation() const { return relaxation_; }
    double lastDiagonalMinimum() const { return lastDiagonalMinimum_; }
    double lastDiagonalMaximum() const { return lastDiagonalMaximum_; }
    double lastMaximumFaceCoupling() const
    { return lastMaximumFaceCoupling_; }
    double lastMaximumSpeciesProjection() const
    { return lastMaximumSpeciesProjection_; }
    int lastProjectedCells() const { return lastProjectedCells_; }
    int lastBacktracks() const { return lastBacktracks_; }
    int lastResidualEvaluations() const
    { return lastResidualEvaluations_; }
    double lastTrialResidual() const
    { return lastTrialResidual_; }
    double lastFullTrialResidual() const
    { return lastFullTrialResidual_; }
    double lastBestTrialResidual() const
    { return lastBestTrialResidual_; }
    const std::vector<double>& lastFullTrialComponentResiduals() const
    { return lastFullTrialComponentResiduals_; }
    double lastAcceptedResidualScale() const
    { return lastAcceptedResidualScale_; }
    double lastAcceptedResidualDirection() const
    { return lastAcceptedResidualDirection_; }
    int lastAcceptedResidualSubspace() const
    { return lastAcceptedResidualSubspace_; }
    int lastPressureSchurResidualEvaluations() const
    { return lastPressureSchurResidualEvaluations_; }
    double lastBestPressureSchurTrialResidual() const
    { return lastBestPressureSchurTrialResidual_; }
    const std::vector<double>& lastBestPressureSchurComponentResiduals() const
    { return lastBestPressureSchurComponentResiduals_; }
    const std::vector<double>& lastAcceptedPressureSchurComponentResiduals() const
    { return lastAcceptedPressureSchurComponentResiduals_; }
    int lastPressureSchurLinearSweeps() const
    { return lastPressureSchurLinearSweeps_; }
    double lastPressureSchurLinearResidualRatio() const
    { return lastPressureSchurLinearResidualRatio_; }
    double lastPressureSchurMaximumRelativePressure() const
    { return lastPressureSchurMaximumRelativePressure_; }
    double lastPressureSchurMaximumVelocityMach() const
    { return lastPressureSchurMaximumVelocityMach_; }
    double lastPressureSchurMinimumFacePseudoTime() const
    { return lastPressureSchurMinimumFacePseudoTime_; }
    double lastPressureSchurMaximumFacePseudoTime() const
    { return lastPressureSchurMaximumFacePseudoTime_; }
    double lastPressureSchurMaximumFaceReciprocityError() const
    { return lastPressureSchurMaximumFaceReciprocityError_; }
    double lastPressureSchurMaximumFaceNonOrthogonalityDegrees() const
    { return lastPressureSchurMaximumFaceNonOrthogonalityDegrees_; }
    int lastPressureSchurLeastSquaresGradientApplications() const
    { return lastPressureSchurLeastSquaresGradientApplications_; }
    int lastPressureSchurLeastSquaresGradientFallbacks() const
    { return lastPressureSchurLeastSquaresGradientFallbacks_; }
    int lastPressureSchurNonOrthogonalCorrections() const
    { return lastPressureSchurNonOrthogonalCorrections_; }
    double lastPressureSchurNonOrthogonalResidualRatio() const
    { return lastPressureSchurNonOrthogonalResidualRatio_; }
    int lastPressureSchurStructuredKrylovIterations() const
    { return lastPressureSchurStructuredKrylovIterations_; }
    double lastPressureSchurStructuredResidualRatio() const
    { return lastPressureSchurStructuredResidualRatio_; }
    int lastPressureSchurStructuredMomentumSolves() const
    { return lastPressureSchurStructuredMomentumSolves_; }
    int lastPressureSchurStructuredMomentumSweeps() const
    { return lastPressureSchurStructuredMomentumSweeps_; }
    double lastPressureSchurStructuredMomentumMaximumResidualRatio() const
    { return lastPressureSchurStructuredMomentumMaximumResidualRatio_; }
    int lastPressureSchurStructuredAcousticPreconditionerSolves() const
    { return lastPressureSchurStructuredAcousticPreconditionerSolves_; }
    int lastPressureSchurStructuredAcousticPreconditionerSweeps() const
    { return lastPressureSchurStructuredAcousticPreconditionerSweeps_; }
    double lastPressureSchurStructuredAcousticPreconditionerMaximumResidualRatio() const
    { return lastPressureSchurStructuredAcousticPreconditionerMaximumResidualRatio_; }
    int lastPressureSchurStructuredPressureDiffusionFaces() const
    { return lastPressureSchurStructuredPressureDiffusionFaces_; }
    double lastPressureSchurStructuredMinimumPressureDiffusionMobility() const
    { return lastPressureSchurStructuredMinimumPressureDiffusionMobility_; }
    double lastPressureSchurStructuredMaximumPressureDiffusionMobility() const
    { return lastPressureSchurStructuredMaximumPressureDiffusionMobility_; }
    double lastPressureSchurStructuredMaximumPressureDiffusionCoefficient() const
    { return lastPressureSchurStructuredMaximumPressureDiffusionCoefficient_; }
    double lastPressureSchurStructuredMaximumPressureDiffusionReciprocityError() const
    { return lastPressureSchurStructuredMaximumPressureDiffusionReciprocityError_; }
    int lastPressureSchurInletPressureDirichletFaces() const
    { return lastPressureSchurInletPressureDirichletFaces_; }
    int lastPressureSchurInletDirichletFaces() const
    { return lastPressureSchurInletDirichletFaces_; }
    int lastPressureSchurInletNeumannFaces() const
    { return lastPressureSchurInletNeumannFaces_; }
    int lastAcousticResidualEvaluations() const
    { return lastAcousticResidualEvaluations_; }
    double lastBestAcousticTrialResidual() const
    { return lastBestAcousticTrialResidual_; }
    int lastFlowResidualEvaluations() const
    { return lastFlowResidualEvaluations_; }
    double lastBestFlowTrialResidual() const
    { return lastBestFlowTrialResidual_; }
    int lastEnergyResidualEvaluations() const
    { return lastEnergyResidualEvaluations_; }
    double lastBestEnergyTrialResidual() const
    { return lastBestEnergyTrialResidual_; }
    double lastResidualAcceptanceLimit() const
    { return lastResidualAcceptanceLimit_; }

    bool step(
        Domain& d,
        const std::vector<std::unique_ptr<Operator>>& ops,
        double dtReference) const override
    {
        if(!std::isfinite(dtReference)||dtReference<=0.0)
            throw std::invalid_argument(
                "LU-SGS residual reference step must be positive");
        for(int cell=0;cell<d.N();++cell)
            if(!solverConsumerAdmitsAcousticState(d.Q(cell)))
            {
                d.recordFailure(Domain::FailureCode::STATE_RECOVERY_FAILED,
                    cell,name(),"solver acoustic state is not admitted");
                return false;
            }
        for(const auto& op:ops)
            if(op->ownsPersistentState())
            {
                d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    name(),"steady LU-SGS is not defined for operator-owned "
                        "state: "+op->name());
                return false;
            }

        const int total=d.N();
        const int nx=d.Nx(),ny=d.Ny(),nz=d.Nz();
        const int plane=nx*ny;
        std::vector<unsigned char> active(
            static_cast<std::size_t>(total),1U);
        if(excludeStructuredHalo_)
            for(int cell=0;cell<total;++cell)
            {
                const int i=cell%nx;
                const int j=(cell/nx)%ny;
                const int k=cell/plane;
                const bool isActive=i>0&&i<nx-1
                    &&(ny==1||(j>0&&j<ny-1))
                    &&(nz==1||(k>0&&k<nz-1));
                active[static_cast<std::size_t>(cell)]=isActive?1U:0U;
            }

        StructuredImplicitStencil stencil;
        stencil.ensureSize(static_cast<std::size_t>(total));
        for(const auto& op:ops) op->addStructuredImplicitStencil(d,stencil);

        std::vector<double> diagonal(static_cast<std::size_t>(total),0.0);
        lastDiagonalMinimum_=std::numeric_limits<double>::infinity();
        lastDiagonalMaximum_=0.0;
        lastMaximumFaceCoupling_=0.0;
        for(int cell=0;cell<total;++cell)
        {
            if(active[static_cast<std::size_t>(cell)]==0U) continue;
            const double spectral=stencil.diagonal[
                static_cast<std::size_t>(cell)];
            if(!std::isfinite(spectral)||spectral<0.0)
                throw std::runtime_error(
                    "LU-SGS operator stencil has an invalid diagonal");
            const double value=spectral>0.0
                ? spectral*(1.0+1.0/pseudoCfl_)
                : 1.0/dtReference;
            diagonal[static_cast<std::size_t>(cell)]=value;
            lastDiagonalMinimum_=std::min(lastDiagonalMinimum_,value);
            lastDiagonalMaximum_=std::max(lastDiagonalMaximum_,value);
            for(int face=0;face<StructuredImplicitStencil::FACE_COUNT;++face)
            {
                const double coupling=stencil.face[face][
                    static_cast<std::size_t>(cell)];
                if(!std::isfinite(coupling)||coupling<0.0)
                    throw std::runtime_error(
                        "LU-SGS operator stencil has an invalid face coupling");
                lastMaximumFaceCoupling_=std::max(
                    lastMaximumFaceCoupling_,coupling);
            }
        }
        if(!std::isfinite(lastDiagonalMinimum_))
        {
            lastDiagonalMinimum_=1.0/dtReference;
            lastDiagonalMaximum_=lastDiagonalMinimum_;
        }

        std::vector<StructuredImplicitBlock> hydroDiagonal(
            static_cast<std::size_t>(total));
        std::vector<unsigned char> hydroDiagonalReady(
            static_cast<std::size_t>(total),0U);
        lastHydroBlockApplications_=0;
        lastHydroBlockFallbacks_=0;
        lastHydroBlockCommitFallbacks_=0;
        lastHydroBlockPivotFallbacks_=0;
        lastHydroBlockOutputFallbacks_=0;
        lastHydroBlockMatrixFallbacks_=0;
        lastHydroBlockRhsFallbacks_=0;
        lastHydroBlockMaximumNormalizedOffDiagonal_=0.0;
        if(useHydroDiagonalBlock_)
        {
            for(int cell=0;cell<total;++cell)
            {
                if(active[static_cast<std::size_t>(cell)]==0U) continue;
                StructuredImplicitBlock& block=
                    hydroDiagonal[static_cast<std::size_t>(cell)];
                block.setScalarIdentity(
                    diagonal[static_cast<std::size_t>(cell)]);
                bool contributed=false;
                for(const auto& op:ops)
                {
                    if(!op->supportsStructuredDiagonalBlock()) continue;
                    StructuredImplicitBlock operatorBlock;
                    double scalarContribution=0.0;
                    if(!op->addStructuredImplicitDiagonalBlock(
                            d,cell,operatorBlock,scalarContribution)
                        ||!operatorBlock.finite()
                        ||!std::isfinite(scalarContribution)
                        ||scalarContribution<0.0)
                    {
                        ++lastHydroBlockFallbacks_;
                        continue;
                    }
                    for(int row=0;
                        row<StructuredImplicitBlock::HYDRO_COMPONENTS;++row)
                    for(int column=0;
                        column<StructuredImplicitBlock::HYDRO_COMPONENTS;
                        ++column)
                        block.add(row,column,
                            operatorBlock.value[row][column]
                            -(row==column?scalarContribution:0.0));
                    contributed=true;
                }
                if(contributed&&block.finite())
                {
                    const auto scales=hydroScales(
                        d.U(cell),d.Q(cell).acoustic.soundSpeed());
                    const double scalarDiagonal=diagonal[
                        static_cast<std::size_t>(cell)];
                    for(int row=0;
                        row<StructuredImplicitBlock::HYDRO_COMPONENTS;++row)
                    {
                        if(!std::isfinite(block.value[row][row])
                            ||block.value[row][row]<scalarDiagonal)
                            block.value[row][row]=scalarDiagonal;
                    }
                    // Frozen real-fluid derivatives can be poorly scaled near
                    // a dense LOX/fuel interface. Enforce normalized strict
                    // row diagonal dominance, not merely a bound on each
                    // individual entry; this prevents a formally finite block
                    // solve from producing an unbounded correction.
                    const double maximumAllowedRowSum=0.25;
                    double maximumNormalizedOffDiagonal=0.0;
                    for(int row=0;
                        row<StructuredImplicitBlock::HYDRO_COMPONENTS;++row)
                    {
                        const double diagonalMagnitude=std::max(
                            std::abs(block.value[row][row]),1.0e-300);
                        double normalizedRowSum=0.0;
                        for(int column=0;
                            column<StructuredImplicitBlock::HYDRO_COMPONENTS;
                            ++column)
                        {
                            if(row==column) continue;
                            normalizedRowSum+=std::abs(
                                block.value[row][column]
                                *scales[static_cast<std::size_t>(column)]
                                /scales[static_cast<std::size_t>(row)])
                                /diagonalMagnitude;
                        }
                        if(normalizedRowSum>maximumAllowedRowSum)
                        {
                            const double blend=maximumAllowedRowSum
                                /normalizedRowSum;
                            for(int column=0;
                                column<StructuredImplicitBlock::HYDRO_COMPONENTS;
                                ++column)
                                if(row!=column)
                                    block.value[row][column]*=blend;
                            normalizedRowSum=maximumAllowedRowSum;
                        }
                        maximumNormalizedOffDiagonal=std::max(
                            maximumNormalizedOffDiagonal,normalizedRowSum);
                    }
                    lastHydroBlockMaximumNormalizedOffDiagonal_=std::max(
                        lastHydroBlockMaximumNormalizedOffDiagonal_,
                        maximumNormalizedOffDiagonal);
                    hydroDiagonalReady[static_cast<std::size_t>(cell)]=1U;
                }
            }
        }

        detail::StepTransaction transaction(d,ops);
        if(!transaction.active()) return false;
        d.beginStep();
        if(!transaction.apply(dtReference)) return false;

        std::vector<ConservedVars> residual(
            static_cast<std::size_t>(total),ConservedVars(d.nSpecies()));
        std::vector<ConservedVars> forward(
            static_cast<std::size_t>(total),ConservedVars(d.nSpecies()));
        std::vector<ConservedVars> correction(
            static_cast<std::size_t>(total),ConservedVars(d.nSpecies()));
        const std::size_t componentsPerCell=static_cast<std::size_t>(
            7+d.nSpecies());
        if(fixedReferenceScale_.size()
            !=static_cast<std::size_t>(total)*componentsPerCell)
        {
            fixedReferenceScale_.assign(
                static_cast<std::size_t>(total)*componentsPerCell,1.0);
            for(int cell=0;cell<total;++cell)
            {
                const ConservedVars& current=d.U(cell);
                const double rhoScale=std::max(std::abs(current.rho),1e-12);
                const double acousticScale=std::max(
                    d.Q(cell).acoustic.soundSpeed(),1.0);
                const double momentumScale=rhoScale*acousticScale;
                const double energyScale=std::max(std::abs(current.rhoE),
                    rhoScale*acousticScale*acousticScale);
                std::size_t component=static_cast<std::size_t>(cell)
                    *componentsPerCell;
                fixedReferenceScale_[component++]=rhoScale;
                fixedReferenceScale_[component++]=energyScale;
                for(std::size_t direction=0;
                    direction<current.rhoU.size();++direction)
                    fixedReferenceScale_[component++]=momentumScale;
                for(std::size_t species=0;
                    species<current.rhoY.size();++species)
                    fixedReferenceScale_[component++]=rhoScale;
                fixedReferenceScale_[component++]=rhoScale;
                fixedReferenceScale_[component++]=rhoScale;
            }
        }
        std::vector<double> componentResidualSums(componentsPerCell,0.0);
        std::vector<std::size_t> componentResidualCounts(componentsPerCell,0U);
        // One nonlinear history must use one norm.  Rebuilding rho/a/rhoE
        // scales after every accepted pseudo-step can make an unchanged
        // physical residual appear to decrease in the trial transaction and
        // then rise when the next baseline is assembled.  The scale field is
        // frozen at the first admitted state of this integrator and shared by
        // the baseline and every safeguarded subspace trial.
        const std::vector<double>& residualReferenceScale=
            fixedReferenceScale_;
        lastMaximumResidualCell_=-1;
        lastMaximumResidualComponent_=-1;
        lastMaximumResidual_=0.0;
        double residualSum=0.0;
        std::size_t residualComponents=0;
        for(int cell=0;cell<total;++cell)
        {
            if(active[static_cast<std::size_t>(cell)]==0U) continue;
            const ConservedVars& current=d.U(cell);
            const ConservedVars& candidate=d.Unew(cell);
            ConservedVars& rate=residual[static_cast<std::size_t>(cell)];
            assignScaledDifference(rate,candidate,current,1.0/dtReference);

            std::size_t component=0;
            auto addResidual=[&](double value)
            {
                const double scale=residualReferenceScale[
                    static_cast<std::size_t>(cell)*componentsPerCell
                    +component];
                if(!std::isfinite(scale)||scale<=0.0)
                    throw std::runtime_error(
                        "LU-SGS residual metric has an invalid fixed scale");
                const double normalized=value/scale;
                const double squared=normalized*normalized;
                residualSum+=squared;
                ++residualComponents;
                componentResidualSums[component]+=squared;
                ++componentResidualCounts[component];
                if(std::abs(normalized)>lastMaximumResidual_)
                {
                    lastMaximumResidual_=std::abs(normalized);
                    lastMaximumResidualCell_=cell;
                    lastMaximumResidualComponent_=
                        static_cast<int>(component);
                }
                ++component;
            };
            addResidual(rate.rho);
            addResidual(rate.rhoE);
            for(double value:rate.rhoU) addResidual(value);
            for(double value:rate.rhoY) addResidual(value);
            addResidual(rate.rhoZ);
            addResidual(rate.rhoC);
        }
        lastResidual_=std::sqrt(residualSum/std::max<std::size_t>(
            residualComponents,1U));
        lastComponentResiduals_.assign(componentsPerCell,0.0);
        for(std::size_t component=0;component<componentsPerCell;++component)
            lastComponentResiduals_[component]=std::sqrt(
                componentResidualSums[component]
                /std::max<std::size_t>(
                    componentResidualCounts[component],1U));
        if(initialResidual_<=0.0)
            initialResidual_=std::max(lastResidual_,1e-300);
        lastResidualRatio_=lastResidual_/initialResidual_;
        if(useResidualLineSearch_)
            bestResidual_=std::min(bestResidual_,lastResidual_);
        const bool residualRising=std::isfinite(previousResidual_)
            &&lastResidual_>1.05*previousResidual_;
        if(useAndersonAcceleration_&&residualRising)
        {
            andersonStateHistory_.clear();
            andersonUpdateHistory_.clear();
            lastAndersonDepthUsed_=0;
        }
        if(!std::isfinite(relaxationReferenceResidual_))
            relaxationReferenceResidual_=lastResidual_;
        else if(lastResidual_>1.10*relaxationReferenceResidual_)
        {
            relaxation_=std::max(
                minimumRelaxation_,0.5*relaxation_);
            // Do not halve every subsequent step while recovering from the
            // same excursion. A new 10% rise from this state is required.
            relaxationReferenceResidual_=lastResidual_;
        }
        else
        {
            if(std::isfinite(previousResidual_)
                &&lastResidual_<=1.02*previousResidual_)
                relaxation_=std::min(
                    maximumRelaxation_,1.02*relaxation_);
            relaxationReferenceResidual_=std::min(
                relaxationReferenceResidual_,lastResidual_);
        }
        previousResidual_=lastResidual_;
        converged_=lastResidual_<residualTolerance_*initialResidual_;
        if(converged_)
        {
            d.beginStep();
            lastRelativeUpdate_=0.0;
            lastBacktracks_=0;
            transaction.accept();
            return true;
        }

        auto activeValue=[&](int cell)->bool
        {
            return cell>=0&&cell<total
                &&active[static_cast<std::size_t>(cell)]!=0U;
        };
        // An operator may write halo candidates even when the steady solve
        // deliberately excludes those cells from its residual. Such writes
        // are not an implicit unknown and must not leak into commit or line
        // search admission; boundary hooks own their authoritative state.
        for(int cell=0;cell<total;++cell)
            if(!activeValue(cell)) d.Unew(cell)=d.U(cell);
        lastBlockFaceApplications_=0;
        lastBlockFaceFallbacks_=0;
        std::vector<ConservedVars> scalarForward(
            static_cast<std::size_t>(total),ConservedVars(d.nSpecies()));
        std::vector<ConservedVars> scalarCorrection(
            static_cast<std::size_t>(total),ConservedVars(d.nSpecies()));
        for(int cell=0;cell<total;++cell)
        {
            if(!activeValue(cell)) continue;
            ConservedVars value=residual[static_cast<std::size_t>(cell)];
            auto lower=[&](int neighbor,int face)
            {
                if(!activeValue(neighbor)) return;
                addScaled(value,scalarForward[
                    static_cast<std::size_t>(neighbor)],
                    stencil.face[face][static_cast<std::size_t>(cell)]);
            };
            lower(cell-1,StructuredImplicitStencil::I_MINUS);
            if(ny>1) lower(cell-nx,StructuredImplicitStencil::J_MINUS);
            if(nz>1) lower(cell-plane,StructuredImplicitStencil::K_MINUS);
            divide(value,diagonal[static_cast<std::size_t>(cell)]);
            scalarForward[static_cast<std::size_t>(cell)]=std::move(value);
        }
        for(int cell=total-1;cell>=0;--cell)
        {
            if(!activeValue(cell)) continue;
            ConservedVars value=scalarForward[static_cast<std::size_t>(cell)];
            auto upper=[&](int neighbor,int face)
            {
                if(!activeValue(neighbor)) return;
                addScaled(value,scalarCorrection[
                    static_cast<std::size_t>(neighbor)],
                    stencil.face[face][static_cast<std::size_t>(cell)]
                        /diagonal[static_cast<std::size_t>(cell)]);
            };
            upper(cell+1,StructuredImplicitStencil::I_PLUS);
            if(ny>1) upper(cell+nx,StructuredImplicitStencil::J_PLUS);
            if(nz>1) upper(cell+plane,StructuredImplicitStencil::K_PLUS);
            scalarCorrection[static_cast<std::size_t>(cell)]=std::move(value);
        }
        auto addBlockNeighborDelta=[&](
            ConservedVars& value,int cell,int neighbor,int face,
            const ConservedVars& neighborCorrection)
        {
            for(const auto& op:ops)
            {
                if(!op->supportsStructuredBlockCoupling()) continue;
                ConservedVars contribution(d.nSpecies());
                if(op->addStructuredImplicitNeighborContribution(
                        d,cell,neighbor,face,neighborCorrection,contribution))
                {
                    addScaled(value,contribution,1.0);
                    ++lastBlockFaceApplications_;
                }
                else ++lastBlockFaceFallbacks_;
            }
        };
        for(int cell=0;cell<total;++cell)
        {
            if(!activeValue(cell)) continue;
            ConservedVars value=residual[static_cast<std::size_t>(cell)];
            auto lower=[&](int neighbor,int face)
            {
                if(!activeValue(neighbor)) return;
                const ConservedVars& neighborCorrection=
                    forward[static_cast<std::size_t>(neighbor)];
                addScaled(value,neighborCorrection,
                    stencil.face[face][static_cast<std::size_t>(cell)]);
                addBlockNeighborDelta(
                    value,cell,neighbor,face,neighborCorrection);
            };
            lower(cell-1,StructuredImplicitStencil::I_MINUS);
            if(ny>1) lower(cell-nx,StructuredImplicitStencil::J_MINUS);
            if(nz>1) lower(cell-plane,StructuredImplicitStencil::K_MINUS);
            bool blockApplied=false;
            if(hydroDiagonalReady[static_cast<std::size_t>(cell)]!=0U)
            {
                int blockFailure=0;
                blockApplied=divideWithHydroBlock(
                    value,hydroDiagonal[static_cast<std::size_t>(cell)],
                    diagonal[static_cast<std::size_t>(cell)],
                    hydroScales(d.U(cell),d.Q(cell).acoustic.soundSpeed()),
                    &blockFailure);
                if(blockApplied) ++lastHydroBlockApplications_;
                else
                {
                    ++lastHydroBlockFallbacks_;
                    if(blockFailure==1) ++lastHydroBlockMatrixFallbacks_;
                    else if(blockFailure==2) ++lastHydroBlockPivotFallbacks_;
                    else if(blockFailure==3) ++lastHydroBlockOutputFallbacks_;
                    else if(blockFailure==4) ++lastHydroBlockRhsFallbacks_;
                }
            }
            if(!blockApplied)
                divide(value,diagonal[static_cast<std::size_t>(cell)]);
            forward[static_cast<std::size_t>(cell)]=std::move(value);
        }
        for(int cell=total-1;cell>=0;--cell)
        {
            if(!activeValue(cell)) continue;
            ConservedVars value=forward[static_cast<std::size_t>(cell)];
            auto upper=[&](int neighbor,int face)
            {
                if(!activeValue(neighbor)) return;
                const ConservedVars& neighborCorrection=
                    correction[static_cast<std::size_t>(neighbor)];
                ConservedVars contribution(d.nSpecies());
                addScaled(contribution,neighborCorrection,
                    stencil.face[face][static_cast<std::size_t>(cell)]);
                addBlockNeighborDelta(
                    contribution,cell,neighbor,face,neighborCorrection);
                bool blockApplied=false;
                if(hydroDiagonalReady[static_cast<std::size_t>(cell)]!=0U)
                {
                    int blockFailure=0;
                    blockApplied=divideWithHydroBlock(
                        contribution,
                        hydroDiagonal[static_cast<std::size_t>(cell)],
                        diagonal[static_cast<std::size_t>(cell)],
                        hydroScales(
                            d.U(cell),d.Q(cell).acoustic.soundSpeed()),
                        &blockFailure);
                    if(blockApplied) ++lastHydroBlockApplications_;
                    else
                    {
                        ++lastHydroBlockFallbacks_;
                        if(blockFailure==1) ++lastHydroBlockMatrixFallbacks_;
                        else if(blockFailure==2) ++lastHydroBlockPivotFallbacks_;
                        else if(blockFailure==3)
                            ++lastHydroBlockOutputFallbacks_;
                        else if(blockFailure==4)
                            ++lastHydroBlockRhsFallbacks_;
                    }
                }
                if(!blockApplied)
                    divide(contribution,
                        diagonal[static_cast<std::size_t>(cell)]);
                addScaled(value,contribution,1.0);
            };
            upper(cell+1,StructuredImplicitStencil::I_PLUS);
            if(ny>1) upper(cell+nx,StructuredImplicitStencil::J_PLUS);
            if(nz>1) upper(cell+plane,StructuredImplicitStencil::K_PLUS);
            correction[static_cast<std::size_t>(cell)]=std::move(value);
        }

        auto correctionMaximum=[&](
            const std::vector<ConservedVars>& field)
        {
            double maximum=0.0;
            for(int cell=0;cell<total;++cell)
            {
                if(!activeValue(cell)) continue;
                const ConservedVars& current=d.U(cell);
                const ConservedVars& delta=field[
                    static_cast<std::size_t>(cell)];
                const double rhoScale=std::max(
                    std::abs(current.rho),1e-12);
                const double acousticScale=std::max(
                    d.Q(cell).acoustic.soundSpeed(),1.0);
                const double momentumScale=rhoScale*acousticScale;
                const double energyScale=std::max(std::abs(current.rhoE),
                    rhoScale*acousticScale*acousticScale);
                auto update=[&](double value,double scale)
                {
                    const double normalized=std::abs(value)
                        /std::max(scale,1e-300);
                    if(!std::isfinite(normalized))
                    {
                        maximum=std::numeric_limits<double>::infinity();
                        return;
                    }
                    maximum=std::max(maximum,normalized);
                };
                update(delta.rho,rhoScale);
                update(delta.rhoE,energyScale);
                for(double value:delta.rhoU) update(value,momentumScale);
                for(double value:delta.rhoY) update(value,rhoScale);
                update(delta.rhoZ,rhoScale);
                update(delta.rhoC,rhoScale);
            }
            return maximum;
        };
        double unscaledRelativeUpdateMaximum=correctionMaximum(correction);
        if(!std::isfinite(unscaledRelativeUpdateMaximum)
            &&useHydroDiagonalBlock_)
        {
            correction=scalarCorrection;
            unscaledRelativeUpdateMaximum=correctionMaximum(correction);
            lastHydroBlockFallbacks_+=static_cast<long>(2*std::count(
                active.begin(),active.end(),static_cast<unsigned char>(1U)));
            lastHydroBlockApplications_=0;
        }
        if(!std::isfinite(unscaledRelativeUpdateMaximum))
        {
            d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                name(),"LU-SGS correction norm is not finite");
            return false;
        }
        double correctionScale=relaxation_;
        if(correctionScale*unscaledRelativeUpdateMaximum
            >maximumRelativeUpdate_)
            correctionScale=maximumRelativeUpdate_
                /unscaledRelativeUpdateMaximum;

        double relativeUpdateMaximum=0.0;
        lastProjectedCells_=0;
        lastMaximumSpeciesProjection_=0.0;
        for(int cell=0;cell<total;++cell)
        {
            if(!activeValue(cell)) continue;
            const ConservedVars& current=d.U(cell);
            const ConservedVars& delta=correction[static_cast<std::size_t>(cell)];
            ConservedVars& candidate=d.Unew(cell);
            assignCandidate(candidate,current,delta,correctionScale);
            candidate.rhoE=current.rhoE
                +correctionScale*energyCorrectionScale_*delta.rhoE;
            double projection=0.0;
            if(limitSpeciesComposition(current,candidate,projection)
                &&projection>0.0)
            {
                ++lastProjectedCells_;
                lastMaximumSpeciesProjection_=std::max(
                    lastMaximumSpeciesProjection_,projection);
            }
            if(std::isfinite(candidate.rho)&&candidate.rho>0.0)
            {
                candidate.rhoZ=std::clamp(candidate.rhoZ,0.0,candidate.rho);
                candidate.rhoC=std::clamp(candidate.rhoC,0.0,candidate.rho);
            }
            const double rhoScale=std::max(std::abs(current.rho),1e-12);
            const double acousticScale=std::max(
                d.Q(cell).acoustic.soundSpeed(),1.0);
            const double momentumScale=rhoScale*acousticScale;
            const double energyScale=std::max(std::abs(current.rhoE),
                rhoScale*acousticScale*acousticScale);
            auto update=[&](double value,double scale)
            {
                relativeUpdateMaximum=std::max(relativeUpdateMaximum,
                    std::abs(value)/std::max(scale,1e-300));
            };
            update(candidate.rho-current.rho,rhoScale);
            update(candidate.rhoE-current.rhoE,energyScale);
            for(std::size_t component=0;
                component<candidate.rhoU.size();++component)
                update(candidate.rhoU[component]-current.rhoU[component],
                    momentumScale);
            for(std::size_t component=0;
                component<candidate.rhoY.size();++component)
                update(candidate.rhoY[component]-current.rhoY[component],
                    rhoScale);
            update(candidate.rhoZ-current.rhoZ,rhoScale);
            update(candidate.rhoC-current.rhoC,rhoScale);
        }
        lastAndersonDepthUsed_=0;
        if(useAndersonAcceleration_)
        {
            const std::size_t vectorSize=static_cast<std::size_t>(total)
                *componentsPerCell;
            currentAndersonState_.assign(vectorSize,0.0);
            currentAndersonUpdate_.assign(vectorSize,0.0);
            for(int cell=0;cell<total;++cell)
            {
                if(!activeValue(cell)) continue;
                const ConservedVars& current=d.U(cell);
                const ConservedVars& candidate=d.Unew(cell);
                std::size_t component=static_cast<std::size_t>(cell)
                    *componentsPerCell;
                auto store=[&](double state,double trial)
                {
                    const double scale=fixedReferenceScale_[component];
                    currentAndersonState_[component]=state/scale;
                    currentAndersonUpdate_[component]=(trial-state)/scale;
                    ++component;
                };
                store(current.rho,candidate.rho);
                store(current.rhoE,candidate.rhoE);
                for(std::size_t direction=0;
                    direction<current.rhoU.size();++direction)
                    store(current.rhoU[direction],candidate.rhoU[direction]);
                for(std::size_t species=0;
                    species<current.rhoY.size();++species)
                    store(current.rhoY[species],candidate.rhoY[species]);
                store(current.rhoZ,candidate.rhoZ);
                store(current.rhoC,candidate.rhoC);
            }
            andersonStateHistory_.push_back(currentAndersonState_);
            andersonUpdateHistory_.push_back(currentAndersonUpdate_);
            const std::size_t maximumHistory=static_cast<std::size_t>(
                andersonDepth_+1);
            if(andersonStateHistory_.size()>maximumHistory)
            {
                andersonStateHistory_.erase(andersonStateHistory_.begin());
                andersonUpdateHistory_.erase(andersonUpdateHistory_.begin());
            }
            const int depth=std::min(andersonDepth_,
                static_cast<int>(andersonUpdateHistory_.size())-1);
            if(depth>0)
            {
                const int start=static_cast<int>(
                    andersonUpdateHistory_.size())-depth-1;
                std::vector<std::vector<double>> normal(
                    depth,std::vector<double>(depth,0.0));
                std::vector<double> rightHandSide(depth,0.0);
                double diagonalSum=0.0;
                for(int column=0;column<depth;++column)
                for(int row=0;row<depth;++row)
                {
                    double dot=0.0;
                    const auto& column0=andersonUpdateHistory_[
                        static_cast<std::size_t>(start+column)];
                    const auto& column1=andersonUpdateHistory_[
                        static_cast<std::size_t>(start+column+1)];
                    const auto& row0=andersonUpdateHistory_[
                        static_cast<std::size_t>(start+row)];
                    const auto& row1=andersonUpdateHistory_[
                        static_cast<std::size_t>(start+row+1)];
                    for(std::size_t value=0;value<vectorSize;++value)
                        dot+=(column1[value]-column0[value])
                            *(row1[value]-row0[value]);
                    normal[column][row]=dot;
                    if(column==row) diagonalSum+=dot;
                }
                for(int column=0;column<depth;++column)
                {
                    double dot=0.0;
                    const auto& update0=andersonUpdateHistory_[
                        static_cast<std::size_t>(start+column)];
                    const auto& update1=andersonUpdateHistory_[
                        static_cast<std::size_t>(start+column+1)];
                    for(std::size_t value=0;value<vectorSize;++value)
                        dot+=(update1[value]-update0[value])
                            *currentAndersonUpdate_[value];
                    rightHandSide[column]=dot;
                }
                const double regularization=andersonRegularization_
                    *std::max(diagonalSum/std::max(depth,1),1.0e-30);
                for(int column=0;column<depth;++column)
                    normal[column][column]+=regularization;
                if(solveDenseSystem(normal,rightHandSide))
                {
                    std::vector<double> accelerated(vectorSize,0.0);
                    for(std::size_t value=0;value<vectorSize;++value)
                        accelerated[value]=currentAndersonState_[value]
                            +currentAndersonUpdate_[value];
                    for(int column=0;column<depth;++column)
                    {
                        const auto& state0=andersonStateHistory_[
                            static_cast<std::size_t>(start+column)];
                        const auto& state1=andersonStateHistory_[
                            static_cast<std::size_t>(start+column+1)];
                        const auto& update0=andersonUpdateHistory_[
                            static_cast<std::size_t>(start+column)];
                        const auto& update1=andersonUpdateHistory_[
                            static_cast<std::size_t>(start+column+1)];
                        for(std::size_t value=0;value<vectorSize;++value)
                            accelerated[value]-=rightHandSide[column]
                                *((state1[value]+update1[value])
                                  -(state0[value]+update0[value]));
                    }
                    // Anchor the accelerated state to the already safeguarded
                    // LU-SGS fixed-point candidate. Limit only the additional
                    // secant correction; a poor Anderson model therefore
                    // cannot replace the robust base direction outright.
                    double accelerationCorrectionMaximum=0.0;
                    for(std::size_t value=0;value<vectorSize;++value)
                    {
                        const double base=currentAndersonState_[value]
                            +currentAndersonUpdate_[value];
                        accelerationCorrectionMaximum=std::max(
                            accelerationCorrectionMaximum,
                            std::abs(accelerated[value]-base));
                    }
                    const double accelerationCorrectionLimit=
                        (andersonMaximumUpdateRatio_-1.0)
                        *std::max(relativeUpdateMaximum,1.0e-16);
                    const double correctionBlend=
                        accelerationCorrectionMaximum
                            >accelerationCorrectionLimit
                        ?accelerationCorrectionLimit
                            /accelerationCorrectionMaximum:1.0;
                    std::vector<double> finalState(vectorSize,0.0);
                    double finalUpdateMaximum=0.0;
                    for(std::size_t value=0;value<vectorSize;++value)
                    {
                        const double base=currentAndersonState_[value]
                            +currentAndersonUpdate_[value];
                        finalState[value]=base+correctionBlend
                            *(accelerated[value]-base);
                        finalUpdateMaximum=std::max(finalUpdateMaximum,
                            std::abs(finalState[value]
                                -currentAndersonState_[value]));
                    }
                    const double totalScale=finalUpdateMaximum
                        >maximumRelativeUpdate_
                        ?maximumRelativeUpdate_/finalUpdateMaximum:1.0;
                    relativeUpdateMaximum=0.0;
                    for(int cell=0;cell<total;++cell)
                    {
                        if(!activeValue(cell)) continue;
                        const ConservedVars& current=d.U(cell);
                        ConservedVars& candidate=d.Unew(cell);
                        std::size_t component=static_cast<std::size_t>(cell)
                            *componentsPerCell;
                        auto acceleratedValue=[&](double state)
                        {
                            const double normalized=currentAndersonState_[component]
                                +totalScale*(finalState[component]
                                    -currentAndersonState_[component]);
                            const double value=normalized
                                *fixedReferenceScale_[component];
                            ++component;
                            return value;
                        };
                        candidate.rho=acceleratedValue(current.rho);
                        candidate.rhoE=acceleratedValue(current.rhoE);
                        for(std::size_t direction=0;
                            direction<candidate.rhoU.size();++direction)
                            candidate.rhoU[direction]=acceleratedValue(
                                current.rhoU[direction]);
                        for(std::size_t species=0;
                            species<candidate.rhoY.size();++species)
                            candidate.rhoY[species]=acceleratedValue(
                                current.rhoY[species]);
                        candidate.rhoZ=acceleratedValue(current.rhoZ);
                        candidate.rhoC=acceleratedValue(current.rhoC);
                        double projection=0.0;
                        if(useSpeciesLimiter_&&limitSpeciesComposition(
                                current,candidate,projection)
                            &&projection>0.0)
                        {
                            ++lastProjectedCells_;
                            lastMaximumSpeciesProjection_=std::max(
                                lastMaximumSpeciesProjection_,projection);
                        }
                        if(std::isfinite(candidate.rho)&&candidate.rho>0.0)
                        {
                            candidate.rhoZ=std::clamp(
                                candidate.rhoZ,0.0,candidate.rho);
                            candidate.rhoC=std::clamp(
                                candidate.rhoC,0.0,candidate.rho);
                        }
                        const double rhoScale=std::max(
                            std::abs(current.rho),1.0e-12);
                        const double acousticScale=std::max(
                            d.Q(cell).acoustic.soundSpeed(),1.0);
                        const double momentumScale=rhoScale*acousticScale;
                        const double energyScale=std::max(
                            std::abs(current.rhoE),
                            rhoScale*acousticScale*acousticScale);
                        auto update=[&](double delta,double scale)
                        {
                            relativeUpdateMaximum=std::max(
                                relativeUpdateMaximum,
                                std::abs(delta)/std::max(scale,1.0e-300));
                        };
                        update(candidate.rho-current.rho,rhoScale);
                        update(candidate.rhoE-current.rhoE,energyScale);
                        for(std::size_t direction=0;
                            direction<candidate.rhoU.size();++direction)
                            update(candidate.rhoU[direction]
                                -current.rhoU[direction],momentumScale);
                        for(std::size_t species=0;
                            species<candidate.rhoY.size();++species)
                            update(candidate.rhoY[species]
                                -current.rhoY[species],rhoScale);
                        update(candidate.rhoZ-current.rhoZ,rhoScale);
                        update(candidate.rhoC-current.rhoC,rhoScale);
                    }
                    lastAndersonDepthUsed_=depth;
                }
            }
        }
        lastRelativeUpdate_=relativeUpdateMaximum;
        lastBacktracks_=0;
        bool committed=d.commit();
        if(!committed&&useHydroDiagonalBlock_)
        {
            // A block candidate that violates the conservative state contract
            // is never repaired by relaxing the contract. Rebuild the same
            // update from the scalar LU-SGS direction inside the transaction,
            // then let the normal residual line search audit that fallback.
            d.resetTrip();
            if(!d.restoreStepTransactionStartState()) return false;
            const double scalarMaximum=correctionMaximum(scalarCorrection);
            if(!std::isfinite(scalarMaximum))
            {
                d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    name(),"scalar LU-SGS fallback norm is not finite");
                return false;
            }
            double scalarScale=relaxation_;
            if(scalarScale*scalarMaximum>maximumRelativeUpdate_)
                scalarScale=maximumRelativeUpdate_/scalarMaximum;
            d.beginStep();
            for(int cell=0;cell<total;++cell)
            {
                if(!activeValue(cell)) continue;
                const ConservedVars& current=d.U(cell);
                const ConservedVars& delta=scalarCorrection[
                    static_cast<std::size_t>(cell)];
                ConservedVars& candidate=d.Unew(cell);
                assignCandidate(candidate,current,delta,scalarScale);
                candidate.rhoE=current.rhoE
                    +scalarScale*energyCorrectionScale_*delta.rhoE;
                double projection=0.0;
                if(useSpeciesLimiter_&&limitSpeciesComposition(
                    current,candidate,projection)
                    &&projection>0.0)
                {
                    ++lastProjectedCells_;
                    lastMaximumSpeciesProjection_=std::max(
                        lastMaximumSpeciesProjection_,projection);
                }
                if(std::isfinite(candidate.rho)&&candidate.rho>0.0)
                {
                    candidate.rhoZ=std::clamp(candidate.rhoZ,0.0,candidate.rho);
                    candidate.rhoC=std::clamp(candidate.rhoC,0.0,candidate.rho);
                }
            }
            committed=d.commit();
            if(committed)
            {
                const long activeCells=static_cast<long>(std::count(
                    active.begin(),active.end(),static_cast<unsigned char>(1U)));
                lastHydroBlockFallbacks_+=2*activeCells;
                ++lastHydroBlockCommitFallbacks_;
                lastHydroBlockApplications_=0;
            }
        }
        for(int backtrack=0;!committed&&backtrack<maxBacktracks_;++backtrack)
        {
            d.resetTrip();
            for(int cell=0;cell<total;++cell)
            {
                if(!activeValue(cell)) continue;
                const ConservedVars& current=d.U(cell);
                ConservedVars& candidate=d.Unew(cell);
                candidate.rho=current.rho+0.5*(candidate.rho-current.rho);
                candidate.rhoE=current.rhoE+0.5*(candidate.rhoE-current.rhoE);
                for(std::size_t component=0;
                    component<candidate.rhoU.size();++component)
                    candidate.rhoU[component]=current.rhoU[component]
                        +0.5*(candidate.rhoU[component]-current.rhoU[component]);
                for(std::size_t component=0;
                    component<candidate.rhoY.size();++component)
                    candidate.rhoY[component]=current.rhoY[component]
                        +0.5*(candidate.rhoY[component]-current.rhoY[component]);
                candidate.rhoZ=current.rhoZ
                    +0.5*(candidate.rhoZ-current.rhoZ);
                candidate.rhoC=current.rhoC
                    +0.5*(candidate.rhoC-current.rhoC);
            }
            ++lastBacktracks_;
            lastRelativeUpdate_*=0.5;
            committed=d.commit();
        }
        if(!committed)
        {
            andersonStateHistory_.clear();
            andersonUpdateHistory_.clear();
            lastAndersonDepthUsed_=0;
            return false;
        }

        lastResidualEvaluations_=0;
        lastTrialResidual_=lastResidual_;
        lastFullTrialResidual_=lastResidual_;
        lastBestTrialResidual_=lastResidual_;
        lastAcceptedResidualScale_=1.0;
        lastAcceptedResidualDirection_=1.0;
        lastAcceptedResidualSubspace_=0;
        lastPressureSchurResidualEvaluations_=0;
        lastBestPressureSchurTrialResidual_=
            std::numeric_limits<double>::infinity();
        lastBestPressureSchurComponentResiduals_.clear();
        lastAcceptedPressureSchurComponentResiduals_.clear();
        lastPressureSchurLinearSweeps_=0;
        lastPressureSchurLinearResidualRatio_=1.0;
        lastPressureSchurMaximumRelativePressure_=0.0;
        lastPressureSchurMaximumVelocityMach_=0.0;
        lastPressureSchurMinimumFacePseudoTime_=0.0;
        lastPressureSchurMaximumFacePseudoTime_=0.0;
        lastPressureSchurMaximumFaceReciprocityError_=0.0;
        lastPressureSchurMaximumFaceNonOrthogonalityDegrees_=0.0;
        lastPressureSchurLeastSquaresGradientApplications_=0;
        lastPressureSchurLeastSquaresGradientFallbacks_=0;
        lastPressureSchurNonOrthogonalCorrections_=0;
        lastPressureSchurNonOrthogonalResidualRatio_=1.0;
        lastPressureSchurStructuredKrylovIterations_=0;
        lastPressureSchurStructuredResidualRatio_=1.0;
        lastPressureSchurStructuredMomentumSolves_=0;
        lastPressureSchurStructuredMomentumSweeps_=0;
        lastPressureSchurStructuredMomentumMaximumResidualRatio_=0.0;
        lastPressureSchurStructuredAcousticPreconditionerSolves_=0;
        lastPressureSchurStructuredAcousticPreconditionerSweeps_=0;
        lastPressureSchurStructuredAcousticPreconditionerMaximumResidualRatio_=0.0;
        lastPressureSchurStructuredPressureDiffusionFaces_=0;
        lastPressureSchurStructuredMinimumPressureDiffusionMobility_=0.0;
        lastPressureSchurStructuredMaximumPressureDiffusionMobility_=0.0;
        lastPressureSchurStructuredMaximumPressureDiffusionCoefficient_=0.0;
        lastPressureSchurStructuredMaximumPressureDiffusionReciprocityError_=0.0;
        lastPressureSchurInletPressureDirichletFaces_=0;
        lastPressureSchurInletDirichletFaces_=0;
        lastPressureSchurInletNeumannFaces_=0;
        lastAcousticResidualEvaluations_=0;
        lastBestAcousticTrialResidual_=
            std::numeric_limits<double>::infinity();
        lastFlowResidualEvaluations_=0;
        lastBestFlowTrialResidual_=std::numeric_limits<double>::infinity();
        lastEnergyResidualEvaluations_=0;
        lastBestEnergyTrialResidual_=std::numeric_limits<double>::infinity();
        if(useResidualLineSearch_)
        {
            std::vector<ConservedVars> fullStepField;
            std::vector<PrimitiveVars> fullStepPrimitive;
            fullStepField.reserve(static_cast<std::size_t>(total));
            fullStepPrimitive.reserve(static_cast<std::size_t>(total));
            for(int cell=0;cell<total;++cell)
            {
                fullStepField.push_back(d.U(cell));
                fullStepPrimitive.push_back(d.W(cell));
            }
            bool residualAccepted=false;
            const bool bootstrap=!useBidirectionalResidualLineSearch_
                &&allowResidualBootstrapStep_
                &&residualEnvelopeHistory_.empty();
            double envelopeReference=lastResidual_;
            if(!residualEnvelopeHistory_.empty())
                envelopeReference=*std::max_element(
                    residualEnvelopeHistory_.begin(),
                    residualEnvelopeHistory_.end());
            lastResidualAcceptanceLimit_=useBidirectionalResidualLineSearch_
                ?lastResidual_*(1.0-residualEvaluationRelativeNoise_)
                :(bootstrap?std::numeric_limits<double>::infinity()
                    :residualEnvelopeGrowth_*envelopeReference);
            const double fullRelativeUpdate=lastRelativeUpdate_;
            const int directionCount=
                useBidirectionalResidualLineSearch_?2:1;
            for(int trial=0;
                trial<=maxBacktracks_&&!residualAccepted;++trial)
            {
                const double residualScale=std::ldexp(1.0,-trial);
                for(int directionIndex=0;
                    directionIndex<directionCount&&!residualAccepted;
                    ++directionIndex)
                {
                    const double residualDirection=
                        directionIndex==0?1.0:-1.0;
                    const bool existingPositiveFullStep=
                        trial==0&&directionIndex==0;
                    if(!existingPositiveFullStep)
                    {
                        d.resetTrip();
                        if(!d.restoreStepTransactionStartState()) break;
                        d.beginStep();
                        for(int cell=0;cell<total;++cell)
                        {
                            if(!activeValue(cell)) continue;
                            const ConservedVars& base=d.U(cell);
                            ConservedVars& candidate=d.Unew(cell);
                            const ConservedVars& full=
                                fullStepField[static_cast<std::size_t>(cell)];
                            const double scale=
                                residualDirection*residualScale;
                            candidate.rho=base.rho
                                +scale*(full.rho-base.rho);
                            candidate.rhoE=base.rhoE
                                +scale*(full.rhoE-base.rhoE);
                            for(std::size_t component=0;
                                component<candidate.rhoU.size();++component)
                                candidate.rhoU[component]=base.rhoU[component]
                                    +scale*(full.rhoU[component]
                                        -base.rhoU[component]);
                            for(std::size_t component=0;
                                component<candidate.rhoY.size();++component)
                                candidate.rhoY[component]=base.rhoY[component]
                                    +scale*(full.rhoY[component]
                                        -base.rhoY[component]);
                            candidate.rhoZ=base.rhoZ
                                +scale*(full.rhoZ-base.rhoZ);
                            candidate.rhoC=base.rhoC
                                +scale*(full.rhoC-base.rhoC);
                        }
                        ++lastBacktracks_;
                        lastRelativeUpdate_=
                            fullRelativeUpdate*residualScale;
                        if(!d.commit()) continue;
                    }
                    else lastRelativeUpdate_=fullRelativeUpdate;

                    if(!nonlinearBoundaryHook_(d)||d.tripped()) continue;
                    d.beginStep();
                    if(!transaction.apply(dtReference)||d.tripped()) continue;
                    std::vector<double> trialComponentSums(
                        componentsPerCell,0.0);
                    std::vector<std::size_t> trialComponentCounts(
                        componentsPerCell,0U);
                    double trialSum=0.0;
                    std::size_t trialComponents=0;
                    for(int cell=0;cell<total;++cell)
                    {
                        if(!activeValue(cell)) continue;
                        const ConservedVars& state=d.U(cell);
                        const ConservedVars& candidate=d.Unew(cell);
                        std::size_t component=0;
                        auto addTrial=[&](double delta)
                        {
                            const double scale=residualReferenceScale[
                                static_cast<std::size_t>(cell)*componentsPerCell
                                +component++];
                            const double normalized=delta
                                /(dtReference*std::max(scale,1.0e-300));
                            const double squared=normalized*normalized;
                            trialSum+=squared;
                            ++trialComponents;
                            trialComponentSums[component-1]+=squared;
                            ++trialComponentCounts[component-1];
                        };
                        addTrial(candidate.rho-state.rho);
                        addTrial(candidate.rhoE-state.rhoE);
                        for(std::size_t direction=0;
                            direction<candidate.rhoU.size();++direction)
                            addTrial(candidate.rhoU[direction]
                                -state.rhoU[direction]);
                        for(std::size_t species=0;
                            species<candidate.rhoY.size();++species)
                            addTrial(candidate.rhoY[species]
                                -state.rhoY[species]);
                        addTrial(candidate.rhoZ-state.rhoZ);
                        addTrial(candidate.rhoC-state.rhoC);
                    }
                    lastTrialResidual_=std::sqrt(trialSum
                        /std::max<std::size_t>(trialComponents,1U));
                    ++lastResidualEvaluations_;
                    if(trial==0&&directionIndex==0)
                    {
                        lastFullTrialResidual_=lastTrialResidual_;
                        lastFullTrialComponentResiduals_.assign(
                            componentsPerCell,0.0);
                        for(std::size_t component=0;
                            component<componentsPerCell;++component)
                            lastFullTrialComponentResiduals_[component]=
                                std::sqrt(trialComponentSums[component]
                                    /std::max<std::size_t>(
                                        trialComponentCounts[component],1U));
                    }
                    lastBestTrialResidual_=std::min(
                        lastBestTrialResidual_,lastTrialResidual_);
                    const bool accepted=std::isfinite(lastTrialResidual_)
                        &&(useBidirectionalResidualLineSearch_
                            ?lastTrialResidual_
                                <=lastResidualAcceptanceLimit_
                            :(bootstrap||lastTrialResidual_
                                <=lastResidualAcceptanceLimit_
                                    *(1.0
                                        +residualEvaluationRelativeNoise_)));
                    if(!accepted) continue;

                    residualAccepted=true;
                    lastAcceptedResidualScale_=residualScale;
                    lastAcceptedResidualDirection_=residualDirection;
                    bestResidual_=std::min(
                        bestResidual_,lastTrialResidual_);
                    residualEnvelopeHistory_.push_back(lastTrialResidual_);
                    if(residualEnvelopeHistory_.size()
                        >static_cast<std::size_t>(residualEnvelopeMemory_))
                        residualEnvelopeHistory_.erase(
                            residualEnvelopeHistory_.begin());
                    if(residualScale<1.0)
                        relaxation_=std::max(minimumRelaxation_,
                            relaxation_*std::sqrt(residualScale));
                    ++acceptedResidualSafeguardSteps_;
                    d.beginStep();
                }
            }
            if(!residualAccepted&&usePressureSchurFallback_)
            {
                d.resetTrip();
                if(!d.restoreStepTransactionStartState()) return false;
                PressureVelocityCorrection pressureCorrection;
                const bool correctionReady=pressureCorrectionBuilder_(
                    d,residual,active,stencil,pseudoCfl_,pressureCorrection);
                if(correctionReady
                    &&pressureCorrection.pressure.size()
                        ==static_cast<std::size_t>(total)
                    &&pressureCorrection.velocity.size()
                        ==static_cast<std::size_t>(total))
                {
                    lastPressureSchurLinearSweeps_=
                        pressureCorrection.linearSweeps;
                    lastPressureSchurLinearResidualRatio_=
                        pressureCorrection.linearResidualRatio;
                    lastPressureSchurMaximumRelativePressure_=
                        pressureCorrection.maximumRelativePressure;
                    lastPressureSchurMaximumVelocityMach_=
                        pressureCorrection.maximumVelocityMach;
                    lastPressureSchurMinimumFacePseudoTime_=
                        pressureCorrection.minimumFacePseudoTime;
                    lastPressureSchurMaximumFacePseudoTime_=
                        pressureCorrection.maximumFacePseudoTime;
                    lastPressureSchurMaximumFaceReciprocityError_=
                        pressureCorrection.maximumFaceReciprocityError;
                    lastPressureSchurMaximumFaceNonOrthogonalityDegrees_=
                        pressureCorrection.maximumFaceNonOrthogonalityDegrees;
                    lastPressureSchurLeastSquaresGradientApplications_=
                        pressureCorrection.leastSquaresGradientApplications;
                    lastPressureSchurLeastSquaresGradientFallbacks_=
                        pressureCorrection.leastSquaresGradientFallbacks;
                    lastPressureSchurNonOrthogonalCorrections_=
                        pressureCorrection.nonOrthogonalCorrections;
                    lastPressureSchurNonOrthogonalResidualRatio_=
                        pressureCorrection.nonOrthogonalResidualRatio;
                    lastPressureSchurStructuredKrylovIterations_=
                        pressureCorrection.structuredKrylovIterations;
                    lastPressureSchurStructuredResidualRatio_=
                        pressureCorrection.structuredResidualRatio;
                    lastPressureSchurStructuredMomentumSolves_=
                        pressureCorrection.structuredMomentumSolves;
                    lastPressureSchurStructuredMomentumSweeps_=
                        pressureCorrection.structuredMomentumSweeps;
                    lastPressureSchurStructuredMomentumMaximumResidualRatio_=
                        pressureCorrection.structuredMomentumMaximumResidualRatio;
                    lastPressureSchurStructuredAcousticPreconditionerSolves_=
                        pressureCorrection.structuredAcousticPreconditionerSolves;
                    lastPressureSchurStructuredAcousticPreconditionerSweeps_=
                        pressureCorrection.structuredAcousticPreconditionerSweeps;
                    lastPressureSchurStructuredAcousticPreconditionerMaximumResidualRatio_=
                        pressureCorrection.structuredAcousticPreconditionerMaximumResidualRatio;
                    lastPressureSchurStructuredPressureDiffusionFaces_=
                        pressureCorrection.structuredPressureDiffusionFaces;
                    lastPressureSchurStructuredMinimumPressureDiffusionMobility_=
                        pressureCorrection.structuredMinimumPressureDiffusionMobility;
                    lastPressureSchurStructuredMaximumPressureDiffusionMobility_=
                        pressureCorrection.structuredMaximumPressureDiffusionMobility;
                    lastPressureSchurStructuredMaximumPressureDiffusionCoefficient_=
                        pressureCorrection.structuredMaximumPressureDiffusionCoefficient;
                    lastPressureSchurStructuredMaximumPressureDiffusionReciprocityError_=
                        pressureCorrection.structuredMaximumPressureDiffusionReciprocityError;
                    lastPressureSchurInletPressureDirichletFaces_=
                        pressureCorrection.inletPressureDirichletFaces;
                    lastPressureSchurInletDirichletFaces_=
                        pressureCorrection.inletDirichletFaces;
                    lastPressureSchurInletNeumannFaces_=
                        pressureCorrection.inletNeumannFaces;
                    const double pressureVelocityRelativeUpdate=std::max(
                        pressureCorrection.maximumRelativePressure,
                        pressureCorrection.maximumVelocityMach);
                    for(int trial=0;
                        trial<=maxBacktracks_&&!residualAccepted;++trial)
                    {
                        const double residualScale=std::ldexp(1.0,-trial);
                        for(int directionIndex=0;
                            directionIndex<2&&!residualAccepted;
                            ++directionIndex)
                        {
                            const double residualDirection=
                                directionIndex==0?1.0:-1.0;
                            const double scale=
                                residualDirection*residualScale;
                            d.resetTrip();
                            if(!d.restoreStepTransactionStartState()) break;
                            std::vector<ConservedVars> candidateField;
                            candidateField.reserve(
                                static_cast<std::size_t>(total));
                            for(int cell=0;cell<total;++cell)
                                candidateField.push_back(d.U(cell));
                            bool trialValid=true;
                            for(int cell=0;cell<total;++cell)
                            {
                                if(!activeValue(cell)) continue;
                                const ConservedVars& base=d.U(cell);
                                const thermo::EquilibriumState& equilibrium=
                                    d.Q(cell);
                                const double sound=
                                    equilibrium.acoustic.soundSpeed();
                                const double densityDelta=scale
                                    *pressureCorrection.pressure[
                                        static_cast<std::size_t>(cell)]
                                    /(sound*sound);
                                ConservedVars& candidate=candidateField[
                                    static_cast<std::size_t>(cell)];
                                candidate.rho=base.rho+densityDelta;
                                if(!std::isfinite(candidate.rho)
                                    ||candidate.rho<=0.0)
                                { trialValid=false; break; }
                                double velocitySquared=0.0;
                                for(std::size_t direction=0;
                                    direction<candidate.rhoU.size();
                                    ++direction)
                                {
                                    candidate.rhoU[direction]=
                                        base.rhoU[direction]+scale*base.rho
                                            *pressureCorrection.velocity[
                                                static_cast<std::size_t>(cell)]
                                                [direction];
                                    const double velocity=
                                        candidate.rhoU[direction]/candidate.rho;
                                    if(!std::isfinite(velocity))
                                    { trialValid=false; break; }
                                    velocitySquared+=velocity*velocity;
                                }
                                if(!trialValid) break;
                                const double internalSpecific=equilibrium.e
                                    +equilibrium.p*densityDelta
                                        /(equilibrium.rho*equilibrium.rho);
                                candidate.rhoE=candidate.rho
                                    *(internalSpecific+0.5*velocitySquared);
                                const double densityRatio=
                                    candidate.rho/base.rho;
                                for(std::size_t species=0;
                                    species<candidate.rhoY.size();++species)
                                    candidate.rhoY[species]=
                                        base.rhoY[species]*densityRatio;
                                candidate.rhoZ=base.rhoZ*densityRatio;
                                candidate.rhoC=base.rhoC*densityRatio;
                                if(!std::isfinite(candidate.rhoE)
                                    ||!std::isfinite(internalSpecific))
                                { trialValid=false; break; }
                            }
                            if(!trialValid) continue;
                            if(!d.replaceConservedFieldInStepTransaction(
                                    candidateField)) continue;
                            ++lastBacktracks_;
                            lastRelativeUpdate_=
                                pressureVelocityRelativeUpdate*residualScale;
                            if(!nonlinearBoundaryHook_(d)||d.tripped())
                                continue;
                            d.beginStep();
                            if(!transaction.apply(dtReference)||d.tripped())
                                continue;
                            double trialSum=0.0;
                            std::size_t trialComponents=0;
                            std::vector<double> componentSums(
                                componentsPerCell,0.0);
                            std::vector<std::size_t> componentCounts(
                                componentsPerCell,0U);
                            for(int cell=0;cell<total;++cell)
                            {
                                if(!activeValue(cell)) continue;
                                const ConservedVars& state=d.U(cell);
                                const ConservedVars& candidate=d.Unew(cell);
                                std::size_t component=0;
                                auto addTrial=[&](double delta)
                                {
                                    const double reference=
                                        residualReferenceScale[
                                            static_cast<std::size_t>(cell)
                                                *componentsPerCell
                                            +component++];
                                    const double normalized=delta
                                        /(dtReference*std::max(
                                            reference,1.0e-300));
                                    const double squared=
                                        normalized*normalized;
                                    trialSum+=squared;
                                    ++trialComponents;
                                    componentSums[component-1]+=squared;
                                    ++componentCounts[component-1];
                                };
                                addTrial(candidate.rho-state.rho);
                                addTrial(candidate.rhoE-state.rhoE);
                                for(std::size_t direction=0;
                                    direction<candidate.rhoU.size();
                                    ++direction)
                                    addTrial(candidate.rhoU[direction]
                                        -state.rhoU[direction]);
                                for(std::size_t species=0;
                                    species<candidate.rhoY.size();++species)
                                    addTrial(candidate.rhoY[species]
                                        -state.rhoY[species]);
                                addTrial(candidate.rhoZ-state.rhoZ);
                                addTrial(candidate.rhoC-state.rhoC);
                            }
                            lastTrialResidual_=std::sqrt(trialSum
                                /std::max<std::size_t>(
                                    trialComponents,1U));
                            ++lastResidualEvaluations_;
                            ++lastPressureSchurResidualEvaluations_;
                            std::vector<double> trialComponentResiduals(
                                componentsPerCell,0.0);
                            for(std::size_t component=0;
                                component<componentsPerCell;++component)
                                trialComponentResiduals[component]=std::sqrt(
                                    componentSums[component]
                                    /std::max<std::size_t>(
                                        componentCounts[component],1U));
                            if(lastTrialResidual_
                                <lastBestPressureSchurTrialResidual_)
                            {
                                lastBestPressureSchurTrialResidual_=
                                    lastTrialResidual_;
                                lastBestPressureSchurComponentResiduals_=
                                    trialComponentResiduals;
                            }
                            const bool flowBlockAdmitted=
                                pressureSchurFlowResidualAdmitted(
                                    lastComponentResiduals_,
                                    trialComponentResiduals,
                                    residualEvaluationRelativeNoise_);
                            const bool accepted=
                                std::isfinite(lastTrialResidual_)
                                &&lastTrialResidual_
                                    <=lastResidualAcceptanceLimit_
                                &&flowBlockAdmitted;
                            if(!accepted) continue;
                            residualAccepted=true;
                            lastAcceptedResidualScale_=residualScale;
                            lastAcceptedResidualDirection_=
                                residualDirection;
                            lastAcceptedResidualSubspace_=4;
                            lastAcceptedPressureSchurComponentResiduals_=
                                trialComponentResiduals;
                            bestResidual_=std::min(
                                bestResidual_,lastTrialResidual_);
                            residualEnvelopeHistory_.push_back(
                                lastTrialResidual_);
                            if(residualEnvelopeHistory_.size()
                                >static_cast<std::size_t>(
                                    residualEnvelopeMemory_))
                                residualEnvelopeHistory_.erase(
                                    residualEnvelopeHistory_.begin());
                            d.beginStep();
                        }
                    }
                }
            }
            if(!residualAccepted&&useSegregatedAcousticFallback_)
            {
                d.resetTrip();
                if(!d.restoreStepTransactionStartState()) return false;
                double acousticRelativeUpdate=0.0;
                for(int cell=0;cell<total;++cell)
                {
                    if(!activeValue(cell)) continue;
                    const PrimitiveVars& base=d.W(cell);
                    const PrimitiveVars& full=fullStepPrimitive[
                        static_cast<std::size_t>(cell)];
                    const double pressureScale=std::max(
                        std::abs(base.p),d.Q(cell).rho
                            *d.Q(cell).acoustic.soundSpeed()
                            *d.Q(cell).acoustic.soundSpeed());
                    acousticRelativeUpdate=std::max(acousticRelativeUpdate,
                        std::abs(full.p-base.p)
                            /std::max(pressureScale,1.0e-300));
                    const double acousticScale=std::max(
                        d.Q(cell).acoustic.soundSpeed(),1.0);
                    for(std::size_t direction=0;
                        direction<base.U.size();++direction)
                        acousticRelativeUpdate=std::max(
                            acousticRelativeUpdate,
                            std::abs(full.U[direction]-base.U[direction])
                                /acousticScale);
                }
                for(int trial=0;
                    trial<=maxBacktracks_&&!residualAccepted;++trial)
                {
                    const double residualScale=std::ldexp(1.0,-trial);
                    for(int directionIndex=0;
                        directionIndex<2&&!residualAccepted;++directionIndex)
                    {
                        const double residualDirection=
                            directionIndex==0?1.0:-1.0;
                        d.resetTrip();
                        if(!d.restoreStepTransactionStartState()) break;
                        const double scale=
                            residualDirection*residualScale;
                        std::vector<Domain::CompleteCellReplacement>
                            replacements;
                        replacements.reserve(static_cast<std::size_t>(total));
                        bool trialValid=true;
                        for(int cell=0;cell<total;++cell)
                        {
                            if(!activeValue(cell)) continue;
                            PrimitiveVars primitive=d.W(cell);
                            const PrimitiveVars& full=fullStepPrimitive[
                                static_cast<std::size_t>(cell)];
                            primitive.p+=scale*(full.p-primitive.p);
                            for(std::size_t direction=0;
                                direction<primitive.U.size();++direction)
                                primitive.U[direction]+=scale
                                    *(full.U[direction]
                                        -primitive.U[direction]);
                            if(!std::isfinite(primitive.p)||primitive.p<=0.0)
                            { trialValid=false; break; }
                            thermo::EquilibriumResult initialized=
                                d.thermoModel().initializePTY(
                                    primitive.p,primitive.T,primitive.Y);
                            if(!initialized.success())
                            { trialValid=false; break; }
                            const int i=cell%nx;
                            const int j=(cell/nx)%ny;
                            const int k=cell/plane;
                            replacements.push_back(
                                Domain::CompleteCellReplacement{
                                    i,j,k,primitive,
                                    std::move(*initialized.state)});
                        }
                        if(!trialValid) continue;
                        if(!d.replaceCompleteCellsInStepTransaction(
                                replacements)) continue;
                        ++lastBacktracks_;
                        lastRelativeUpdate_=acousticRelativeUpdate
                            *residualScale;
                        if(!nonlinearBoundaryHook_(d)||d.tripped()) continue;
                        d.beginStep();
                        if(!transaction.apply(dtReference)||d.tripped())
                            continue;
                        double trialSum=0.0;
                        std::size_t trialComponents=0;
                        std::vector<double> acousticComponentSums(
                            componentsPerCell,0.0);
                        std::vector<std::size_t> acousticComponentCounts(
                            componentsPerCell,0U);
                        for(int cell=0;cell<total;++cell)
                        {
                            if(!activeValue(cell)) continue;
                            const ConservedVars& state=d.U(cell);
                            const ConservedVars& candidate=d.Unew(cell);
                            std::size_t component=0;
                            auto addTrial=[&](double delta)
                            {
                                const double reference=residualReferenceScale[
                                    static_cast<std::size_t>(cell)
                                        *componentsPerCell+component++];
                                const double normalized=delta
                                    /(dtReference*std::max(
                                        reference,1.0e-300));
                                const double squared=normalized*normalized;
                                trialSum+=squared;
                                ++trialComponents;
                                acousticComponentSums[component-1]+=squared;
                                ++acousticComponentCounts[component-1];
                            };
                            addTrial(candidate.rho-state.rho);
                            addTrial(candidate.rhoE-state.rhoE);
                            for(std::size_t direction=0;
                                direction<candidate.rhoU.size();++direction)
                                addTrial(candidate.rhoU[direction]
                                    -state.rhoU[direction]);
                            for(std::size_t species=0;
                                species<candidate.rhoY.size();++species)
                                addTrial(candidate.rhoY[species]
                                    -state.rhoY[species]);
                            addTrial(candidate.rhoZ-state.rhoZ);
                            addTrial(candidate.rhoC-state.rhoC);
                        }
                        lastTrialResidual_=std::sqrt(trialSum
                            /std::max<std::size_t>(trialComponents,1U));
                        ++lastResidualEvaluations_;
                        ++lastAcousticResidualEvaluations_;
                        lastBestAcousticTrialResidual_=std::min(
                            lastBestAcousticTrialResidual_,
                            lastTrialResidual_);
                        bool acousticComponentsAdmitted=true;
                        for(std::size_t component=0;
                            component<std::min<std::size_t>(
                                5,componentsPerCell);++component)
                        {
                            if(component==1) continue;
                            const double trialComponent=std::sqrt(
                                acousticComponentSums[component]
                                /std::max<std::size_t>(
                                    acousticComponentCounts[component],1U));
                            const double limit=lastComponentResiduals_[component]
                                *(1.0+residualEvaluationRelativeNoise_);
                            if(!std::isfinite(trialComponent)
                                ||trialComponent>limit)
                            {
                                acousticComponentsAdmitted=false;
                                break;
                            }
                        }
                        const bool accepted=std::isfinite(lastTrialResidual_)
                            &&lastTrialResidual_<=lastResidualAcceptanceLimit_
                            &&acousticComponentsAdmitted;
                        if(!accepted) continue;
                        residualAccepted=true;
                        lastAcceptedResidualScale_=residualScale;
                        lastAcceptedResidualDirection_=residualDirection;
                        lastAcceptedResidualSubspace_=3;
                        bestResidual_=std::min(bestResidual_,
                            lastTrialResidual_);
                        residualEnvelopeHistory_.push_back(
                            lastTrialResidual_);
                        if(residualEnvelopeHistory_.size()
                            >static_cast<std::size_t>(residualEnvelopeMemory_))
                            residualEnvelopeHistory_.erase(
                                residualEnvelopeHistory_.begin());
                        d.beginStep();
                    }
                }
            }
            if(!residualAccepted&&useSegregatedFlowFallback_)
            {
                d.resetTrip();
                if(!d.restoreStepTransactionStartState()) return false;
                double flowRelativeUpdate=0.0;
                for(int cell=0;cell<total;++cell)
                {
                    if(!activeValue(cell)) continue;
                    const ConservedVars& base=d.U(cell);
                    const ConservedVars& full=fullStepField[
                        static_cast<std::size_t>(cell)];
                    const double rhoScale=std::max(std::abs(base.rho),1.0e-12);
                    const double acousticScale=std::max(
                        d.Q(cell).acoustic.soundSpeed(),1.0);
                    const double momentumScale=rhoScale*acousticScale;
                    flowRelativeUpdate=std::max(flowRelativeUpdate,
                        std::abs(full.rho-base.rho)/rhoScale);
                    for(std::size_t direction=0;
                        direction<base.rhoU.size();++direction)
                        flowRelativeUpdate=std::max(flowRelativeUpdate,
                            std::abs(full.rhoU[direction]-base.rhoU[direction])
                                /momentumScale);
                }
                for(int trial=0;
                    trial<=maxBacktracks_&&!residualAccepted;++trial)
                {
                    const double residualScale=std::ldexp(1.0,-trial);
                    for(int directionIndex=0;
                        directionIndex<2&&!residualAccepted;++directionIndex)
                    {
                        const double residualDirection=
                            directionIndex==0?1.0:-1.0;
                        d.resetTrip();
                        if(!d.restoreStepTransactionStartState()) break;
                        d.beginStep();
                        bool trialValid=true;
                        for(int cell=0;cell<total;++cell)
                        {
                            if(!activeValue(cell)) continue;
                            const ConservedVars& base=d.U(cell);
                            const ConservedVars& full=fullStepField[
                                static_cast<std::size_t>(cell)];
                            ConservedVars& candidate=d.Unew(cell);
                            candidate=base;
                            const double scale=
                                residualDirection*residualScale;
                            candidate.rho=base.rho
                                +scale*(full.rho-base.rho);
                            if(!std::isfinite(candidate.rho)
                                ||candidate.rho<=0.0)
                            { trialValid=false; break; }
                            for(std::size_t direction=0;
                                direction<candidate.rhoU.size();++direction)
                                candidate.rhoU[direction]=base.rhoU[direction]
                                    +scale*(full.rhoU[direction]
                                        -base.rhoU[direction]);
                            const double densityRatio=candidate.rho/base.rho;
                            for(std::size_t species=0;
                                species<candidate.rhoY.size();++species)
                                candidate.rhoY[species]=
                                    base.rhoY[species]*densityRatio;
                            candidate.rhoZ=base.rhoZ*densityRatio;
                            candidate.rhoC=base.rhoC*densityRatio;
                            double baseMomentumSquared=0.0;
                            double candidateMomentumSquared=0.0;
                            for(std::size_t direction=0;
                                direction<candidate.rhoU.size();++direction)
                            {
                                baseMomentumSquared+=base.rhoU[direction]
                                    *base.rhoU[direction];
                                candidateMomentumSquared+=
                                    candidate.rhoU[direction]
                                    *candidate.rhoU[direction];
                            }
                            const double baseInternalSpecific=(base.rhoE
                                -0.5*baseMomentumSquared/base.rho)/base.rho;
                            candidate.rhoE=candidate.rho*baseInternalSpecific
                                +0.5*candidateMomentumSquared/candidate.rho;
                            if(!std::isfinite(candidate.rhoE))
                            { trialValid=false; break; }
                        }
                        if(!trialValid) continue;
                        ++lastBacktracks_;
                        lastRelativeUpdate_=flowRelativeUpdate*residualScale;
                        if(!d.commit()) continue;
                        if(!nonlinearBoundaryHook_(d)||d.tripped()) continue;
                        d.beginStep();
                        if(!transaction.apply(dtReference)||d.tripped()) continue;
                        double trialSum=0.0;
                        std::size_t trialComponents=0;
                        std::vector<double> flowComponentSums(
                            componentsPerCell,0.0);
                        std::vector<std::size_t> flowComponentCounts(
                            componentsPerCell,0U);
                        for(int cell=0;cell<total;++cell)
                        {
                            if(!activeValue(cell)) continue;
                            const ConservedVars& state=d.U(cell);
                            const ConservedVars& candidate=d.Unew(cell);
                            std::size_t component=0;
                            auto addTrial=[&](double delta)
                            {
                                const double reference=residualReferenceScale[
                                    static_cast<std::size_t>(cell)
                                        *componentsPerCell+component++];
                                const double normalized=delta
                                    /(dtReference*std::max(reference,1.0e-300));
                                const double squared=normalized*normalized;
                                trialSum+=squared;
                                ++trialComponents;
                                flowComponentSums[component-1]+=squared;
                                ++flowComponentCounts[component-1];
                            };
                            addTrial(candidate.rho-state.rho);
                            addTrial(candidate.rhoE-state.rhoE);
                            for(std::size_t direction=0;
                                direction<candidate.rhoU.size();++direction)
                                addTrial(candidate.rhoU[direction]
                                    -state.rhoU[direction]);
                            for(std::size_t species=0;
                                species<candidate.rhoY.size();++species)
                                addTrial(candidate.rhoY[species]
                                    -state.rhoY[species]);
                            addTrial(candidate.rhoZ-state.rhoZ);
                            addTrial(candidate.rhoC-state.rhoC);
                        }
                        lastTrialResidual_=std::sqrt(trialSum
                            /std::max<std::size_t>(trialComponents,1U));
                        ++lastResidualEvaluations_;
                        ++lastFlowResidualEvaluations_;
                        lastBestFlowTrialResidual_=std::min(
                            lastBestFlowTrialResidual_,lastTrialResidual_);
                        bool flowComponentsAdmitted=true;
                        for(std::size_t component=0;
                            component<std::min<std::size_t>(5,componentsPerCell);
                            ++component)
                        {
                            if(component==1) continue;
                            const double trialComponent=std::sqrt(
                                flowComponentSums[component]
                                /std::max<std::size_t>(
                                    flowComponentCounts[component],1U));
                            const double limit=lastComponentResiduals_[component]
                                *(1.0+residualEvaluationRelativeNoise_);
                            if(!std::isfinite(trialComponent)
                                ||trialComponent>limit)
                            {
                                flowComponentsAdmitted=false;
                                break;
                            }
                        }
                        const bool accepted=std::isfinite(lastTrialResidual_)
                            &&lastTrialResidual_<=lastResidualAcceptanceLimit_
                            &&flowComponentsAdmitted;
                        if(!accepted) continue;
                        residualAccepted=true;
                        lastAcceptedResidualScale_=residualScale;
                        lastAcceptedResidualDirection_=residualDirection;
                        lastAcceptedResidualSubspace_=2;
                        bestResidual_=std::min(bestResidual_,lastTrialResidual_);
                        residualEnvelopeHistory_.push_back(lastTrialResidual_);
                        if(residualEnvelopeHistory_.size()
                            >static_cast<std::size_t>(residualEnvelopeMemory_))
                            residualEnvelopeHistory_.erase(
                                residualEnvelopeHistory_.begin());
                        d.beginStep();
                    }
                }
            }
            if(!residualAccepted&&useSegregatedEnergyFallback_)
            {
                d.resetTrip();
                if(!d.restoreStepTransactionStartState()) return false;
                double energyRelativeUpdate=0.0;
                for(int cell=0;cell<total;++cell)
                {
                    if(!activeValue(cell)) continue;
                    const ConservedVars& base=d.U(cell);
                    const ConservedVars& full=fullStepField[
                        static_cast<std::size_t>(cell)];
                    const double rhoScale=std::max(std::abs(base.rho),1.0e-12);
                    const double acousticScale=std::max(
                        d.Q(cell).acoustic.soundSpeed(),1.0);
                    const double energyScale=std::max(std::abs(base.rhoE),
                        rhoScale*acousticScale*acousticScale);
                    energyRelativeUpdate=std::max(energyRelativeUpdate,
                        std::abs(full.rhoE-base.rhoE)
                            /std::max(energyScale,1.0e-300));
                }
                for(int trial=0;
                    trial<=maxBacktracks_&&!residualAccepted;++trial)
                {
                    const double residualScale=std::ldexp(1.0,-trial);
                    for(int directionIndex=0;
                        directionIndex<2&&!residualAccepted;++directionIndex)
                    {
                        const double residualDirection=
                            directionIndex==0?1.0:-1.0;
                        d.resetTrip();
                        if(!d.restoreStepTransactionStartState()) break;
                        d.beginStep();
                        for(int cell=0;cell<total;++cell)
                        {
                            if(!activeValue(cell)) continue;
                            const ConservedVars& base=d.U(cell);
                            const ConservedVars& full=fullStepField[
                                static_cast<std::size_t>(cell)];
                            ConservedVars& candidate=d.Unew(cell);
                            candidate=base;
                            candidate.rhoE=base.rhoE
                                +residualDirection*residualScale
                                    *(full.rhoE-base.rhoE);
                        }
                        ++lastBacktracks_;
                        lastRelativeUpdate_=energyRelativeUpdate
                            *residualScale;
                        if(!d.commit()) continue;
                        if(!nonlinearBoundaryHook_(d)||d.tripped()) continue;
                        d.beginStep();
                        if(!transaction.apply(dtReference)||d.tripped()) continue;
                        double trialSum=0.0;
                        std::size_t trialComponents=0;
                        for(int cell=0;cell<total;++cell)
                        {
                            if(!activeValue(cell)) continue;
                            const ConservedVars& state=d.U(cell);
                            const ConservedVars& candidate=d.Unew(cell);
                            std::size_t component=0;
                            auto addTrial=[&](double delta)
                            {
                                const double scale=residualReferenceScale[
                                    static_cast<std::size_t>(cell)
                                        *componentsPerCell+component++];
                                const double normalized=delta
                                    /(dtReference*std::max(scale,1.0e-300));
                                trialSum+=normalized*normalized;
                                ++trialComponents;
                            };
                            addTrial(candidate.rho-state.rho);
                            addTrial(candidate.rhoE-state.rhoE);
                            for(std::size_t direction=0;
                                direction<candidate.rhoU.size();++direction)
                                addTrial(candidate.rhoU[direction]
                                    -state.rhoU[direction]);
                            for(std::size_t species=0;
                                species<candidate.rhoY.size();++species)
                                addTrial(candidate.rhoY[species]
                                    -state.rhoY[species]);
                            addTrial(candidate.rhoZ-state.rhoZ);
                            addTrial(candidate.rhoC-state.rhoC);
                        }
                        lastTrialResidual_=std::sqrt(trialSum
                            /std::max<std::size_t>(trialComponents,1U));
                        ++lastResidualEvaluations_;
                        ++lastEnergyResidualEvaluations_;
                        lastBestEnergyTrialResidual_=std::min(
                            lastBestEnergyTrialResidual_,lastTrialResidual_);
                        const bool accepted=std::isfinite(lastTrialResidual_)
                            &&lastTrialResidual_<=lastResidualAcceptanceLimit_;
                        if(!accepted) continue;
                        residualAccepted=true;
                        lastAcceptedResidualScale_=residualScale;
                        lastAcceptedResidualDirection_=residualDirection;
                        lastAcceptedResidualSubspace_=1;
                        bestResidual_=std::min(bestResidual_,lastTrialResidual_);
                        residualEnvelopeHistory_.push_back(lastTrialResidual_);
                        if(residualEnvelopeHistory_.size()
                            >static_cast<std::size_t>(residualEnvelopeMemory_))
                            residualEnvelopeHistory_.erase(
                                residualEnvelopeHistory_.begin());
                        d.beginStep();
                    }
                }
            }
            if(!residualAccepted)
            {
                d.resetTrip();
                d.restoreStepTransactionStartState();
                andersonStateHistory_.clear();
                andersonUpdateHistory_.clear();
                lastAndersonDepthUsed_=0;
                d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    name(),"LU-SGS residual line search exhausted");
                return false;
            }
        }
        if(lastBacktracks_>0)
        {
            andersonStateHistory_.clear();
            andersonUpdateHistory_.clear();
            lastAndersonDepthUsed_=0;
        }
        transaction.accept();
        return true;
    }

    std::string name() const override
    { return "Steady LU-SGS (structured spectral bounds)"; }
};

} // End namespace Foam

#endif
