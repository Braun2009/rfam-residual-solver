#ifndef SectorCoarseCorrection_H
#define SectorCoarseCorrection_H

/*---------------------------------------------------------------------------*\
    SectorCoarseCorrection - aligned conservative coarse-grid defect transfer.

      U_f,new = U_f,base + alpha P(U_c,improved - U_c,base)

    This preserves the fine-grid base structure instead of replacing it with a
    prolonged coarse state. The coarse correction is prolonged piecewise
    constantly with exact volume scaling. Optional conservative pair smoothing
    acts on the correction, not on the fine base field.

    A local direction limiter protects density and bounded transported
    scalars. Species use a separate positivity-preserving projection of only
    the zero-sum composition direction, so an absent/tracer species cannot
    suppress the density/momentum/energy correction. A global line search then
    enforces complete rho-e-Y recoverability. Every limitation and inventory
    change is diagnosed because it modifies the requested coarse defect.
\*---------------------------------------------------------------------------*/

#include "Domain.H"
#include "SectorConservativeGridTransfer.H"
#include "SectorMeshMetrics.H"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace Foam
{

class SectorCoarseCorrection
{
public:
    struct Result
    {
        double acceptedGlobalScale=0.0;
        int localLimitedCells=0;
        int recoveryBacktracks=0;
        int compositionProjectedCells=0;
        int compositionFrozenCells=0;
        int compositionReconstructedCells=0;
        int compositionProjectedParents=0;
        int boundaryMaskedFineCells=0;
        double maximumLocalReduction=0.0;
        double maximumCompositionProjection=0.0;
        double maximumParentCompositionProjection=0.0;
        double maximumUnscaledRelativeCorrection=0.0;
        double trustRegionScale=1.0;
        double correctionDirection=1.0;
        double maximumInventoryDefectRelativeToBase=0.0;
        double maximumFlowInventoryDefectRelativeToBase=0.0;
        double maximumSpeciesInventoryDefectRelativeToBase=0.0;
        int firstRecoveryFailedCell=-1;
        int firstRecoveryEquilibriumStatus=0;
        std::string firstRecoveryBackendStatus;
        double firstFailureBaseDensity=0.0;
        double firstFailureCandidateDensity=0.0;
        double firstFailureBaseSpecificInternalEnergy=0.0;
        double firstFailureCandidateSpecificInternalEnergy=0.0;
        double firstFailureDensityCorrection=0.0;
        double firstFailureTotalEnergyCorrection=0.0;
        int firstCompositionFailureI=-1;
        int firstCompositionFailureJ=-1;
        int firstCompositionFailureK=-1;
        std::string firstCompositionFailureReason;
        double firstCompositionFailureValue=0.0;
        double firstCompositionFailureTolerance=0.0;
    };

private:
    static int index(const Domain& domain,int i,int j,int k)
    { return k*domain.Nx()*domain.Ny()+j*domain.Nx()+i; }

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

    static ConservedVars difference(
        const ConservedVars& improved,const ConservedVars& base)
    {
        ConservedVars result(base.rhoY.size());
        result.rho=improved.rho-base.rho;
        result.rhoE=improved.rhoE-base.rhoE;
        for(std::size_t component=0;component<result.rhoU.size();++component)
            result.rhoU[component]=
                improved.rhoU[component]-base.rhoU[component];
        for(std::size_t component=0;component<result.rhoY.size();++component)
            result.rhoY[component]=
                improved.rhoY[component]-base.rhoY[component];
        result.rhoZ=improved.rhoZ-base.rhoZ;
        result.rhoC=improved.rhoC-base.rhoC;
        return result;
    }

    static void scale(ConservedVars& value,double factor)
    {
        value.rho*=factor;
        value.rhoE*=factor;
        for(double& component:value.rhoU) component*=factor;
        for(double& component:value.rhoY) component*=factor;
        value.rhoZ*=factor;
        value.rhoC*=factor;
    }

    static void mixPair(
        ConservedVars& left,ConservedVars& right,
        double leftVolume,double rightVolume,double fraction)
    {
        const ConservedVars oldLeft=left;
        const ConservedVars oldRight=right;
        const double totalVolume=leftVolume+rightVolume;
        const double leftWeight=fraction*rightVolume/totalVolume;
        const double rightWeight=fraction*leftVolume/totalVolume;
        auto mix=[&](double a,double b,double& newA,double& newB)
        {
            newA=a+leftWeight*(b-a);
            newB=b+rightWeight*(a-b);
        };
        mix(oldLeft.rho,oldRight.rho,left.rho,right.rho);
        mix(oldLeft.rhoE,oldRight.rhoE,left.rhoE,right.rhoE);
        for(std::size_t component=0;component<left.rhoU.size();++component)
            mix(oldLeft.rhoU[component],oldRight.rhoU[component],
                left.rhoU[component],right.rhoU[component]);
        for(std::size_t component=0;component<left.rhoY.size();++component)
            mix(oldLeft.rhoY[component],oldRight.rhoY[component],
                left.rhoY[component],right.rhoY[component]);
        mix(oldLeft.rhoZ,oldRight.rhoZ,left.rhoZ,right.rhoZ);
        mix(oldLeft.rhoC,oldRight.rhoC,left.rhoC,right.rhoC);
    }

    static double specificInternalEnergy(const ConservedVars& state)
    {
        double velocitySquared=0.0;
        for(double momentum:state.rhoU)
        {
            const double velocity=momentum/state.rho;
            velocitySquared+=velocity*velocity;
        }
        return state.rhoE/state.rho-0.5*velocitySquared;
    }

    static double localAdmissibleScale(
        const ConservedVars& base,const ConservedVars& correction)
    {
        double result=1.0;
        auto lowerBound=[&](double state,double delta,double floor)
        {
            if(delta<0.0)
                result=std::min(result,
                    (state-floor)/(-delta));
        };
        const double rhoFloor=1.0e-12*std::max(base.rho,1.0);
        lowerBound(base.rho,correction.rho,rhoFloor);
        lowerBound(base.rhoZ,correction.rhoZ,0.0);
        lowerBound(base.rhoC,correction.rhoC,0.0);
        const double rhoCandidate=base.rho+correction.rho;
        const double zCandidate=base.rhoZ+correction.rhoZ;
        const double cCandidate=base.rhoC+correction.rhoC;
        if(zCandidate>rhoCandidate)
        {
            const double denominator=
                (correction.rhoZ-correction.rho);
            if(denominator>0.0)
                result=std::min(result,
                    (base.rho-base.rhoZ)/denominator);
        }
        if(cCandidate>rhoCandidate)
        {
            const double denominator=
                (correction.rhoC-correction.rho);
            if(denominator>0.0)
                result=std::min(result,
                    (base.rho-base.rhoC)/denominator);
        }
        return result<1.0
            ?std::clamp(0.999999*result,0.0,1.0):1.0;
    }

    // Split the species candidate into a density-following part based on the
    // fine-grid composition and a zero-sum composition correction. Limiting
    // only the latter retains the full density correction and does not let a
    // zero/tracer species collapse momentum or energy progress.
    static bool projectSpeciesComposition(
        const ConservedVars& base,
        ConservedVars& candidate,
        double& maximumRelativeProjection)
    {
        maximumRelativeProjection=0.0;
        if(!std::isfinite(candidate.rho)||candidate.rho<=0.0
            ||candidate.rhoY.empty()
            ||candidate.rhoY.size()!=base.rhoY.size()) return false;
        double baseSpeciesMass=0.0,candidateSpeciesMass=0.0;
        bool nonnegative=true;
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
        {
            if(!std::isfinite(base.rhoY[species])
                ||base.rhoY[species]<0.0
                ||!std::isfinite(candidate.rhoY[species])) return false;
            baseSpeciesMass+=base.rhoY[species];
            candidateSpeciesMass+=candidate.rhoY[species];
            if(candidate.rhoY[species]<0.0) nonnegative=false;
        }
        if(!std::isfinite(baseSpeciesMass)||baseSpeciesMass<=0.0
            ||!std::isfinite(candidateSpeciesMass)) return false;
        if(nonnegative
            &&std::abs(candidateSpeciesMass-candidate.rho)
                <=1.0e-12*candidate.rho)
            return true;

        std::vector<double> densityFollowing(candidate.rhoY.size(),0.0);
        std::vector<double> composition(candidate.rhoY.size(),0.0);
        double compositionSum=0.0;
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
        {
            const double massFraction=base.rhoY[species]/baseSpeciesMass;
            densityFollowing[species]=massFraction*candidate.rho;
            composition[species]=
                candidate.rhoY[species]-densityFollowing[species];
            compositionSum+=composition[species];
        }
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
            composition[species]-=
                (base.rhoY[species]/baseSpeciesMass)*compositionSum;

        double beta=1.0;
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
            if(composition[species]<0.0)
                beta=std::min(beta,
                    densityFollowing[species]/(-composition[species]));
        beta=std::clamp(beta,0.0,1.0);

        double limitedSum=0.0;
        std::size_t largest=0;
        for(std::size_t species=0;species<candidate.rhoY.size();++species)
        {
            const double previous=candidate.rhoY[species];
            const double limited=std::max(
                densityFollowing[species]+beta*composition[species],0.0);
            candidate.rhoY[species]=limited;
            limitedSum+=limited;
            if(limited>candidate.rhoY[largest]) largest=species;
            maximumRelativeProjection=std::max(maximumRelativeProjection,
                std::abs(limited-previous)/candidate.rho);
        }
        candidate.rhoY[largest]+=candidate.rho-limitedSum;
        return candidate.rhoY[largest]>=0.0
            &&std::isfinite(candidate.rhoY[largest]);
    }

    static bool projectParentSpeciesInventory(
        const std::vector<double>& baseInventory,
        double densityInventory,
        std::vector<double>& candidateInventory,
        double& maximumRelativeProjection)
    {
        maximumRelativeProjection=0.0;
        if(baseInventory.empty()
            ||candidateInventory.size()!=baseInventory.size()
            ||!std::isfinite(densityInventory)||densityInventory<=0.0)
            return false;
        double baseDensity=0.0,candidateDensity=0.0;
        for(std::size_t species=0;species<baseInventory.size();++species)
        {
            if(!std::isfinite(baseInventory[species])
                ||baseInventory[species]<0.0
                ||!std::isfinite(candidateInventory[species]))
                return false;
            baseDensity+=baseInventory[species];
            candidateDensity+=candidateInventory[species];
        }
        if(!std::isfinite(baseDensity)||baseDensity<=0.0
            ||!std::isfinite(candidateDensity)) return false;

        std::vector<double> densityFollowing(baseInventory.size(),0.0);
        std::vector<double> composition(baseInventory.size(),0.0);
        double compositionSum=0.0;
        for(std::size_t species=0;species<baseInventory.size();++species)
        {
            const double fraction=baseInventory[species]/baseDensity;
            densityFollowing[species]=fraction*densityInventory;
            composition[species]=
                candidateInventory[species]-densityFollowing[species];
            compositionSum+=composition[species];
        }
        for(std::size_t species=0;species<baseInventory.size();++species)
            composition[species]-=
                (baseInventory[species]/baseDensity)*compositionSum;

        double beta=1.0;
        for(std::size_t species=0;species<baseInventory.size();++species)
            if(composition[species]<0.0)
                beta=std::min(beta,
                    densityFollowing[species]/(-composition[species]));
        beta=std::clamp(beta,0.0,1.0);

        double projectedSum=0.0;
        std::size_t largest=0;
        for(std::size_t species=0;species<baseInventory.size();++species)
        {
            const double previous=candidateInventory[species];
            candidateInventory[species]=std::max(
                densityFollowing[species]+beta*composition[species],0.0);
            projectedSum+=candidateInventory[species];
            if(candidateInventory[species]>candidateInventory[largest])
                largest=species;
            maximumRelativeProjection=std::max(
                maximumRelativeProjection,
                std::abs(candidateInventory[species]-previous)
                    /densityInventory);
        }
        candidateInventory[largest]+=densityInventory-projectedSum;
        return candidateInventory[largest]>=0.0
            &&std::isfinite(candidateInventory[largest]);
    }

    static ConservedVars inventory(
        const Domain& domain,const SectorMeshMetrics& mesh,
        const std::vector<ConservedVars>* field=nullptr)
    {
        ConservedVars total(domain.nSpecies());
        for(int i=1;i<domain.Nx()-1;++i)
        for(int j=1;j<domain.Ny()-1;++j)
        for(int k=1;k<domain.Nz()-1;++k)
        {
            const int cell=index(domain,i,j,k);
            addScaled(total,
                field?(*field)[static_cast<std::size_t>(cell)]
                     :domain.U(i,j,k),
                mesh.volume(i,j,k));
        }
        return total;
    }

public:
    static Result apply(
        const Domain& fineBase,
        const SectorMeshMetrics& fineMesh,
        const Domain& coarseBase,
        const Domain& coarseImproved,
        const SectorMeshMetrics& coarseMesh,
        Domain& target,
        int correctionSmoothingSweeps=0,
        int maximumRecoveryBacktracks=24,
        double maximumRelativeCorrection=
            std::numeric_limits<double>::infinity(),
        bool excludeBoundaryParentLayer=false,
        double correctionDirection=1.0,
        bool freezeFineComposition=false,
        bool conservativeParentComposition=false)
    {
        if(!fineBase.stateReady()||fineBase.tripped()
            ||!coarseBase.stateReady()||coarseBase.tripped()
            ||!coarseImproved.stateReady()||coarseImproved.tripped()
            ||fineBase.stepTransactionActive()
            ||coarseBase.stepTransactionActive()
            ||coarseImproved.stepTransactionActive()
            ||target.stepTransactionActive())
            throw std::runtime_error(
                "coarse correction requires ready states outside transactions");
        if(fineBase.nSpecies()!=target.nSpecies()
            ||fineBase.speciesNames()!=target.speciesNames()
            ||fineBase.thermoModelId()!=target.thermoModelId()
            ||coarseBase.nSpecies()!=fineBase.nSpecies()
            ||coarseImproved.nSpecies()!=fineBase.nSpecies()
            ||coarseBase.speciesNames()!=fineBase.speciesNames()
            ||coarseImproved.speciesNames()!=fineBase.speciesNames()
            ||coarseBase.thermoModelId()!=fineBase.thermoModelId()
            ||coarseImproved.thermoModelId()!=fineBase.thermoModelId())
            throw std::runtime_error(
                "coarse correction state-contract mismatch");
        if(coarseBase.Nx()!=coarseImproved.Nx()
            ||coarseBase.Ny()!=coarseImproved.Ny()
            ||coarseBase.Nz()!=coarseImproved.Nz())
            throw std::runtime_error(
                "coarse correction grids do not match");
        if(!(maximumRelativeCorrection>0.0)
            ||std::isnan(maximumRelativeCorrection))
            throw std::runtime_error(
                "coarse correction trust region must be positive");
        if(!std::isfinite(correctionDirection)
            ||std::abs(std::abs(correctionDirection)-1.0)>1.0e-12)
            throw std::runtime_error(
                "coarse correction direction must be +1 or -1");
        if(freezeFineComposition&&conservativeParentComposition)
            throw std::runtime_error(
                "fine composition modes are mutually exclusive");
        const SectorConservativeGridTransfer::Plan transfer=
            SectorConservativeGridTransfer::plan(coarseBase,fineBase);
        if(transfer.direction
            !=SectorConservativeGridTransfer::Direction::PROLONG)
            throw std::runtime_error(
                "coarse correction requires coarse-to-fine aligned grids");

        Result result;
        result.correctionDirection=correctionDirection;
        auto boundaryParentMasked=[&](int i,int j)
        {
            if(!excludeBoundaryParentLayer) return false;
            const int ci=1+(i-1)/transfer.ratioI;
            const int cj=1+(j-1)/transfer.ratioJ;
            // i/j are physical inlet/outlet/axis/wall directions. k is
            // periodic and must not be treated as a correction boundary.
            return ci==1||ci==coarseBase.Nx()-2
                ||cj==1||cj==coarseBase.Ny()-2;
        };
        std::vector<ConservedVars> correction(
            static_cast<std::size_t>(fineBase.N()),
            ConservedVars(fineBase.nSpecies()));
        for(int i=1;i<fineBase.Nx()-1;++i)
        for(int j=1;j<fineBase.Ny()-1;++j)
        for(int k=1;k<fineBase.Nz()-1;++k)
        {
            const int ci=1+(i-1)/transfer.ratioI;
            const int cj=1+(j-1)/transfer.ratioJ;
            const int ck=1+(k-1)/transfer.ratioK;
            double childVolumeSum=0.0;
            const int beginI=1+(ci-1)*transfer.ratioI;
            const int beginJ=1+(cj-1)*transfer.ratioJ;
            const int beginK=1+(ck-1)*transfer.ratioK;
            for(int fi=beginI;fi<beginI+transfer.ratioI;++fi)
            for(int fj=beginJ;fj<beginJ+transfer.ratioJ;++fj)
            for(int fk=beginK;fk<beginK+transfer.ratioK;++fk)
                childVolumeSum+=fineMesh.volume(fi,fj,fk);
            ConservedVars delta=difference(
                coarseImproved.U(ci,cj,ck),coarseBase.U(ci,cj,ck));
            scale(delta,correctionDirection);
            if(boundaryParentMasked(i,j))
            {
                scale(delta,0.0);
                ++result.boundaryMaskedFineCells;
            }
            else
                scale(delta,coarseMesh.volume(ci,cj,ck)/childVolumeSum);
            correction[static_cast<std::size_t>(index(fineBase,i,j,k))]=
                std::move(delta);
        }

        if(correctionSmoothingSweeps>0)
        {
            const double fraction=0.25;
            for(int sweep=0;sweep<correctionSmoothingSweeps;++sweep)
            for(int direction=0;direction<3;++direction)
            for(int parity=0;parity<2;++parity)
            for(int i=1;i<fineBase.Nx()-1;++i)
            for(int j=1;j<fineBase.Ny()-1;++j)
            for(int k=1;k<fineBase.Nz()-1;++k)
            {
                const int coordinate=direction==0?i:(direction==1?j:k);
                if((coordinate-1)%2!=parity) continue;
                int ni=i,nj=j,nk=k;
                if(direction==0) ++ni;
                else if(direction==1) ++nj;
                else ++nk;
                if(ni>=fineBase.Nx()-1
                    ||nj>=fineBase.Ny()-1
                    ||nk>=fineBase.Nz()-1
                    ||boundaryParentMasked(i,j)
                    ||boundaryParentMasked(ni,nj)) continue;
                ConservedVars& left=correction[static_cast<std::size_t>(
                    index(fineBase,i,j,k))];
                ConservedVars& right=correction[static_cast<std::size_t>(
                    index(fineBase,ni,nj,nk))];
                mixPair(left,right,
                    fineMesh.volume(i,j,k),fineMesh.volume(ni,nj,nk),fraction);
            }
        }

        const ConservedVars baseInventory=inventory(fineBase,fineMesh);
        const ConservedVars requestedCorrectionInventory=
            inventory(fineBase,fineMesh,&correction);
        for(int i=1;i<fineBase.Nx()-1;++i)
        for(int j=1;j<fineBase.Ny()-1;++j)
        for(int k=1;k<fineBase.Nz()-1;++k)
        {
            const int cell=index(fineBase,i,j,k);
            const ConservedVars& base=fineBase.U(cell);
            const ConservedVars& delta=
                correction[static_cast<std::size_t>(cell)];
            const double rhoScale=std::max(std::abs(base.rho),1.0e-12);
            const double acousticScale=std::max(
                fineBase.Q(cell).acoustic.soundSpeed(),1.0);
            const double momentumScale=rhoScale*acousticScale;
            const double energyScale=std::max(
                std::abs(base.rhoE),
                rhoScale*acousticScale*acousticScale);
            auto update=[&](double value,double valueScale)
            {
                result.maximumUnscaledRelativeCorrection=std::max(
                    result.maximumUnscaledRelativeCorrection,
                    std::abs(value)/std::max(valueScale,1.0e-300));
            };
            update(delta.rho,rhoScale);
            update(delta.rhoE,energyScale);
            for(double value:delta.rhoU) update(value,momentumScale);
            for(double value:delta.rhoY) update(value,rhoScale);
            update(delta.rhoZ,rhoScale);
            update(delta.rhoC,rhoScale);
        }
        if(std::isfinite(maximumRelativeCorrection)
            &&result.maximumUnscaledRelativeCorrection
                >maximumRelativeCorrection)
            result.trustRegionScale=maximumRelativeCorrection
                /result.maximumUnscaledRelativeCorrection;
        std::vector<double> temperatureSeeds(
            static_cast<std::size_t>(fineBase.N()),300.0);
        for(int cell=0;cell<fineBase.N();++cell)
            temperatureSeeds[static_cast<std::size_t>(cell)]=fineBase.W(cell).T;

        double globalScale=result.trustRegionScale;
        std::vector<ConservedVars> candidate;
        std::vector<ConservedVars> appliedCorrection;
        for(int backtrack=0;backtrack<=maximumRecoveryBacktracks;++backtrack)
        {
            candidate.clear();
            appliedCorrection.assign(
                static_cast<std::size_t>(fineBase.N()),
                ConservedVars(fineBase.nSpecies()));
            candidate.reserve(static_cast<std::size_t>(fineBase.N()));
            int limitedCells=0;
            int projectedCells=0;
            int frozenCompositionCells=0;
            int reconstructedCompositionCells=0;
            int projectedCompositionParents=0;
            double maximumReduction=0.0;
            double maximumProjection=0.0;
            double maximumParentProjection=0.0;
            bool projectionFailed=false;
            for(int cell=0;cell<fineBase.N();++cell)
            {
                const ConservedVars& base=fineBase.U(cell);
                ConservedVars delta=correction[static_cast<std::size_t>(cell)];
                scale(delta,globalScale);
                double localScale=1.0;
                const int i=cell%fineBase.Nx();
                const int j=(cell/fineBase.Nx())%fineBase.Ny();
                const int k=cell/(fineBase.Nx()*fineBase.Ny());
                const bool active=i>0&&i<fineBase.Nx()-1
                    &&j>0&&j<fineBase.Ny()-1
                    &&k>0&&k<fineBase.Nz()-1;
                if(active)
                    localScale=localAdmissibleScale(base,delta);
                else
                {
                    // Halos are not coarse-correction unknowns. Their zero
                    // correction is retained and boundary operators will
                    // republish them before the next residual evaluation.
                    scale(delta,0.0);
                    localScale=1.0;
                }
                if(localScale<1.0)
                {
                    ++limitedCells;
                    maximumReduction=std::max(
                        maximumReduction,1.0-localScale);
                    scale(delta,localScale);
                }
                ConservedVars state=base;
                addScaled(state,delta,1.0);
                if(active)
                {
                    if(conservativeParentComposition)
                    {
                        // Keep the provisional signed species correction until
                        // every child in this parent aggregate is available.
                        // A conservative nonnegative reconstruction follows.
                    }
                    else if(freezeFineComposition)
                    {
                        double baseSpeciesMass=0.0;
                        for(double value:base.rhoY)
                            baseSpeciesMass+=value;
                        if(!std::isfinite(baseSpeciesMass)
                            ||baseSpeciesMass<=0.0
                            ||!std::isfinite(state.rho)||state.rho<=0.0)
                        {
                            projectionFailed=true;
                            break;
                        }
                        bool changed=false;
                        for(std::size_t species=0;
                            species<state.rhoY.size();++species)
                        {
                            const double frozen=state.rho
                                *base.rhoY[species]/baseSpeciesMass;
                            changed=changed||frozen!=state.rhoY[species];
                            state.rhoY[species]=frozen;
                        }
                        if(changed) ++frozenCompositionCells;
                    }
                    else
                    {
                        double projection=0.0;
                        if(!projectSpeciesComposition(base,state,projection))
                        {
                            projectionFailed=true;
                            break;
                        }
                        if(projection>0.0)
                        {
                            ++projectedCells;
                            maximumProjection=std::max(
                                maximumProjection,projection);
                        }
                    }
                    for(std::size_t species=0;
                        species<delta.rhoY.size();++species)
                        delta.rhoY[species]=
                            state.rhoY[species]-base.rhoY[species];
                }
                appliedCorrection[static_cast<std::size_t>(cell)]=delta;
                candidate.push_back(std::move(state));
            }
            if(!projectionFailed&&conservativeParentComposition)
            {
                for(int ci=1;ci<coarseBase.Nx()-1;++ci)
                for(int cj=1;cj<coarseBase.Ny()-1;++cj)
                for(int ck=1;ck<coarseBase.Nz()-1;++ck)
                {
                    auto recordCompositionFailure=[&](
                        const std::string& reason,double value,double tolerance)
                    {
                        if(result.firstCompositionFailureReason.empty())
                        {
                            result.firstCompositionFailureI=ci;
                            result.firstCompositionFailureJ=cj;
                            result.firstCompositionFailureK=ck;
                            result.firstCompositionFailureReason=reason;
                            result.firstCompositionFailureValue=value;
                            result.firstCompositionFailureTolerance=tolerance;
                        }
                    };
                    std::vector<double> speciesInventory(
                        static_cast<std::size_t>(fineBase.nSpecies()),0.0);
                    std::vector<double> baseSpeciesInventory(
                        static_cast<std::size_t>(fineBase.nSpecies()),0.0);
                    double densityInventory=0.0;
                    const int beginI=1+(ci-1)*transfer.ratioI;
                    const int beginJ=1+(cj-1)*transfer.ratioJ;
                    const int beginK=1+(ck-1)*transfer.ratioK;
                    for(int fi=beginI;fi<beginI+transfer.ratioI;++fi)
                    for(int fj=beginJ;fj<beginJ+transfer.ratioJ;++fj)
                    for(int fk=beginK;fk<beginK+transfer.ratioK;++fk)
                    {
                        const int cell=index(fineBase,fi,fj,fk);
                        const ConservedVars& state=
                            candidate[static_cast<std::size_t>(cell)];
                        const ConservedVars& base=fineBase.U(cell);
                        const double volume=fineMesh.volume(fi,fj,fk);
                        densityInventory+=volume*state.rho;
                        for(std::size_t species=0;
                            species<state.rhoY.size();++species)
                        {
                            speciesInventory[species]+=
                                volume*state.rhoY[species];
                            baseSpeciesInventory[species]+=
                                volume*base.rhoY[species];
                        }
                    }
                    if(!std::isfinite(densityInventory)
                        ||densityInventory<=0.0)
                    {
                        recordCompositionFailure(
                            "invalid parent density inventory",
                            densityInventory,0.0);
                        projectionFailed=true;
                        break;
                    }
                    double parentProjection=0.0;
                    if(!projectParentSpeciesInventory(
                            baseSpeciesInventory,densityInventory,
                            speciesInventory,parentProjection))
                    {
                        recordCompositionFailure(
                            "parent species inventory projection failed",
                            parentProjection,0.0);
                        projectionFailed=true;
                        break;
                    }
                    if(parentProjection>0.0)
                    {
                        ++projectedCompositionParents;
                        maximumParentProjection=std::max(
                            maximumParentProjection,parentProjection);
                    }
                    double speciesSum=0.0;
                    for(double value:speciesInventory)
                    {
                        if(!std::isfinite(value)||value<0.0)
                        {
                            recordCompositionFailure(
                                "projected parent species inventory is invalid",
                                value,0.0);
                            projectionFailed=true;
                            break;
                        }
                        speciesSum+=value;
                    }
                    const double inventoryMismatch=
                        std::abs(speciesSum-densityInventory);
                    const double inventoryTolerance=
                        1.0e-12*std::max(densityInventory,1.0e-300);
                    if(projectionFailed
                        ||inventoryMismatch>inventoryTolerance)
                    {
                        if(!projectionFailed)
                            recordCompositionFailure(
                                "projected parent species-density mismatch",
                                inventoryMismatch,inventoryTolerance);
                        projectionFailed=true;
                        break;
                    }
                    std::vector<int> childCells;
                    std::vector<double> childVolumes;
                    std::vector<double> rowTargets;
                    childCells.reserve(static_cast<std::size_t>(
                        transfer.ratioI*transfer.ratioJ*transfer.ratioK));
                    childVolumes.reserve(childCells.capacity());
                    rowTargets.reserve(childCells.capacity());
                    for(int fi=beginI;fi<beginI+transfer.ratioI;++fi)
                    for(int fj=beginJ;fj<beginJ+transfer.ratioJ;++fj)
                    for(int fk=beginK;fk<beginK+transfer.ratioK;++fk)
                    {
                        const int cell=index(fineBase,fi,fj,fk);
                        const double volume=fineMesh.volume(fi,fj,fk);
                        childCells.push_back(cell);
                        childVolumes.push_back(volume);
                        rowTargets.push_back(volume*candidate[
                            static_cast<std::size_t>(cell)].rho);
                    }
                    const std::size_t childCount=childCells.size();
                    const std::size_t speciesCount=speciesInventory.size();
                    std::vector<std::vector<double>> balanced(
                        childCount,std::vector<double>(speciesCount,0.0));
                    std::vector<double> column(speciesCount,0.0);
                    for(std::size_t child=0;child<childCount;++child)
                    {
                        const ConservedVars& base=fineBase.U(
                            childCells[child]);
                        double baseSpeciesMass=0.0;
                        for(double value:base.rhoY)
                            baseSpeciesMass+=value;
                        if(!std::isfinite(baseSpeciesMass)
                            ||baseSpeciesMass<=0.0)
                        {
                            recordCompositionFailure(
                                "invalid fine base species mass",
                                baseSpeciesMass,0.0);
                            projectionFailed=true;
                            break;
                        }
                        for(std::size_t species=0;
                            species<speciesCount;++species)
                        {
                            const double value=rowTargets[child]
                                *base.rhoY[species]/baseSpeciesMass;
                            balanced[child][species]=value;
                            column[species]+=value;
                        }
                    }
                    if(projectionFailed) break;

                    std::vector<double> margin(speciesCount,0.0);
                    for(std::size_t species=0;
                        species<speciesCount;++species)
                        margin[species]=std::max(
                            speciesInventory[species],0.0)-column[species];

                    // Move mass only within a child row. This preserves every
                    // child density exactly while matching the parent species
                    // columns. Donor availability is guaranteed by a negative
                    // column margin; the deterministic species/child order
                    // makes the reconstruction reproducible.
                    const double marginTolerance=
                        2.0e-13*std::max(densityInventory,1.0e-300);
                    for(std::size_t receiver=0;
                        receiver<speciesCount;++receiver)
                    {
                        if(margin[receiver]<=marginTolerance) continue;
                        for(std::size_t donor=0;
                            donor<speciesCount
                                &&margin[receiver]>marginTolerance;
                            ++donor)
                        {
                            if(margin[donor]>=-marginTolerance) continue;
                            for(std::size_t child=0;
                                child<childCount
                                    &&margin[receiver]>marginTolerance
                                    &&margin[donor]<-marginTolerance;
                                ++child)
                            {
                                const double transfer=std::min({
                                    margin[receiver],-margin[donor],
                                    balanced[child][donor]});
                                if(transfer<=0.0) continue;
                                balanced[child][donor]-=transfer;
                                balanced[child][receiver]+=transfer;
                                margin[receiver]-=transfer;
                                margin[donor]+=transfer;
                            }
                        }
                    }
                    double maximumMargin=0.0;
                    for(double value:margin)
                        maximumMargin=std::max(maximumMargin,std::abs(value));
                    if(!std::isfinite(maximumMargin)
                        ||maximumMargin>5.0*marginTolerance)
                    {
                        recordCompositionFailure(
                            "parent species margin transfer incomplete",
                            maximumMargin,5.0*marginTolerance);
                        projectionFailed=true;
                        break;
                    }
                    for(std::size_t child=0;child<childCount;++child)
                    {
                        const int cell=childCells[child];
                        ConservedVars& state=candidate[
                            static_cast<std::size_t>(cell)];
                        const ConservedVars& base=fineBase.U(cell);
                        ConservedVars& delta=appliedCorrection[
                            static_cast<std::size_t>(cell)];
                        double rowSum=0.0;
                        std::size_t largest=0;
                        for(std::size_t species=0;
                            species<speciesCount;++species)
                        {
                            state.rhoY[species]=
                                balanced[child][species]/childVolumes[child];
                            rowSum+=state.rhoY[species];
                            if(state.rhoY[species]>state.rhoY[largest])
                                largest=species;
                        }
                        state.rhoY[largest]+=state.rho-rowSum;
                        if(state.rhoY[largest]<0.0
                            ||!std::isfinite(state.rhoY[largest]))
                        {
                            recordCompositionFailure(
                                "fine row roundoff repair is negative",
                                state.rhoY[largest],0.0);
                            projectionFailed=true;
                            break;
                        }
                        for(std::size_t species=0;
                            species<speciesCount;++species)
                            delta.rhoY[species]=
                                state.rhoY[species]-base.rhoY[species];
                        ++reconstructedCompositionCells;
                    }
                    if(projectionFailed) break;
                }
            }
            if(projectionFailed)
            {
                globalScale*=0.5;
                continue;
            }
            if(target.tripped()) target.resetTrip();
            if(target.replaceConservedField(candidate,temperatureSeeds))
            {
                result.acceptedGlobalScale=globalScale;
                result.localLimitedCells=limitedCells;
                result.recoveryBacktracks=backtrack;
                result.compositionProjectedCells=projectedCells;
                result.compositionFrozenCells=frozenCompositionCells;
                result.compositionReconstructedCells=
                    reconstructedCompositionCells;
                result.compositionProjectedParents=
                    projectedCompositionParents;
                result.maximumLocalReduction=maximumReduction;
                result.maximumCompositionProjection=maximumProjection;
                result.maximumParentCompositionProjection=
                    maximumParentProjection;
                break;
            }
            if(result.firstRecoveryFailedCell<0)
            {
                const Domain::FailureInfo& failure=target.failure();
                result.firstRecoveryFailedCell=failure.cell;
                result.firstRecoveryEquilibriumStatus=
                    static_cast<int>(failure.equilibriumStatus);
                result.firstRecoveryBackendStatus=
                    failure.thermoDiagnostics.backendStatus;
                if(failure.cell>=0&&failure.cell<fineBase.N())
                {
                    const ConservedVars& base=fineBase.U(failure.cell);
                    const ConservedVars& failed=
                        candidate[static_cast<std::size_t>(failure.cell)];
                    result.firstFailureBaseDensity=base.rho;
                    result.firstFailureCandidateDensity=failed.rho;
                    result.firstFailureBaseSpecificInternalEnergy=
                        specificInternalEnergy(base);
                    result.firstFailureCandidateSpecificInternalEnergy=
                        specificInternalEnergy(failed);
                    result.firstFailureDensityCorrection=failed.rho-base.rho;
                    result.firstFailureTotalEnergyCorrection=
                        failed.rhoE-base.rhoE;
                }
            }
            globalScale*=0.5;
        }
        if(result.acceptedGlobalScale<=0.0)
            throw std::runtime_error(
                "coarse correction rho-e-Y recovery line search exhausted: "
                "cell="+std::to_string(result.firstRecoveryFailedCell)
                +" status="
                +std::to_string(result.firstRecoveryEquilibriumStatus)
                +" backend="+result.firstRecoveryBackendStatus
                +" rhoBase/candidate="
                +std::to_string(result.firstFailureBaseDensity)+"/"
                +std::to_string(result.firstFailureCandidateDensity)
                +" eBase/candidate="
                +std::to_string(
                    result.firstFailureBaseSpecificInternalEnergy)+"/"
                +std::to_string(
                    result.firstFailureCandidateSpecificInternalEnergy)
                +" drho/drhoE="
                +std::to_string(result.firstFailureDensityCorrection)+"/"
                +std::to_string(result.firstFailureTotalEnergyCorrection)
                +" compositionParent="
                +std::to_string(result.firstCompositionFailureI)+","
                +std::to_string(result.firstCompositionFailureJ)+","
                +std::to_string(result.firstCompositionFailureK)
                +" compositionReason="
                +result.firstCompositionFailureReason
                +" value/tolerance="
                +std::to_string(result.firstCompositionFailureValue)+"/"
                +std::to_string(result.firstCompositionFailureTolerance));

        const ConservedVars appliedCorrectionInventory=
            inventory(fineBase,fineMesh,&appliedCorrection);
        double maximumFlowDefect=0.0;
        double maximumSpeciesDefect=0.0;
        auto compare=[&](double requested,double applied,double base,
                         double& maximumDefect)
        {
            const double acceptedRequested=
                result.acceptedGlobalScale*requested;
            maximumDefect=std::max(maximumDefect,
                std::abs(acceptedRequested-applied)
                /std::max(std::abs(base),1.0e-300));
        };
        compare(requestedCorrectionInventory.rho,
            appliedCorrectionInventory.rho,baseInventory.rho,
            maximumFlowDefect);
        compare(requestedCorrectionInventory.rhoE,
            appliedCorrectionInventory.rhoE,baseInventory.rhoE,
            maximumFlowDefect);
        for(std::size_t component=0;
            component<baseInventory.rhoU.size();++component)
            compare(requestedCorrectionInventory.rhoU[component],
                appliedCorrectionInventory.rhoU[component],
                baseInventory.rhoU[component],maximumFlowDefect);
        for(std::size_t component=0;
            component<baseInventory.rhoY.size();++component)
            compare(requestedCorrectionInventory.rhoY[component],
                appliedCorrectionInventory.rhoY[component],
                baseInventory.rhoY[component],maximumSpeciesDefect);
        compare(requestedCorrectionInventory.rhoZ,
            appliedCorrectionInventory.rhoZ,baseInventory.rhoZ,
            maximumFlowDefect);
        compare(requestedCorrectionInventory.rhoC,
            appliedCorrectionInventory.rhoC,baseInventory.rhoC,
            maximumFlowDefect);
        result.maximumFlowInventoryDefectRelativeToBase=maximumFlowDefect;
        result.maximumSpeciesInventoryDefectRelativeToBase=
            maximumSpeciesDefect;
        result.maximumInventoryDefectRelativeToBase=
            std::max(maximumFlowDefect,maximumSpeciesDefect);
        return result;
    }
};

} // End namespace Foam

#endif
