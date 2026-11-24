#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../src/framework/SectorCoarseCorrection.H"
#include "../src/framework/SectorConservativeGridTransfer.H"
#include "../src/mesh/SectorMeshMetrics.H"
#include "../src/thermo/SpeciesRegistry.H"

using namespace Foam;

namespace
{

bool initialize(
    Domain& domain,const SpeciesRegistry& registry,
    bool heterogeneousOxygen=false)
{
    const int oxygen=registry.index("O2");
    const int nitrogen=registry.index("N2");
    const int carbonDioxide=registry.index("CO2");
    const int water=registry.index("H2O");
    for(int i=0;i<domain.Nx();++i)
    for(int j=0;j<domain.Ny();++j)
    for(int k=0;k<domain.Nz();++k)
    {
        PrimitiveVars& primitive=domain.stagedW(i,j,k);
        const double x=static_cast<double>(i)/std::max(domain.Nx()-1,1);
        const double r=static_cast<double>(j)/std::max(domain.Ny()-1,1);
        primitive.p=9.0e5+1.0e5*x;
        primitive.T=600.0+200.0*x+20.0*r;
        primitive.U={30.0+10.0*x,-1.0+2.0*r,0.5};
        primitive.Y.assign(static_cast<std::size_t>(registry.size()),0.0);
        const double oxygenFraction=heterogeneousOxygen
            ?(i%2==1?1.0e-14:0.4):0.2;
        primitive.Y[static_cast<std::size_t>(oxygen)]=oxygenFraction;
        primitive.Y[static_cast<std::size_t>(nitrogen)]=0.5-oxygenFraction;
        primitive.Y[static_cast<std::size_t>(carbonDioxide)]=0.25+0.05*x;
        primitive.Y[static_cast<std::size_t>(water)]=
            1.0-primitive.Y[static_cast<std::size_t>(oxygen)]
               -primitive.Y[static_cast<std::size_t>(nitrogen)]
               -primitive.Y[static_cast<std::size_t>(carbonDioxide)];
        primitive.Z=0.1+0.1*x;
        primitive.C=0.05+0.05*r;
    }
    return domain.primitiveToConserved();
}

std::vector<ConservedVars> capture(const Domain& domain)
{
    std::vector<ConservedVars> field;
    field.reserve(static_cast<std::size_t>(domain.N()));
    for(int cell=0;cell<domain.N();++cell) field.push_back(domain.U(cell));
    return field;
}

void scaleState(ConservedVars& state,double factor)
{
    state.rho*=factor;
    state.rhoE*=factor;
    for(double& component:state.rhoU) component*=factor;
    for(double& component:state.rhoY) component*=factor;
    state.rhoZ*=factor;
    state.rhoC*=factor;
}

double maximumFieldDifference(const Domain& left,const Domain& right)
{
    double maximum=0.0;
    auto compare=[&](double a,double b)
    {
        maximum=std::max(maximum,
            std::abs(a-b)/std::max({std::abs(a),std::abs(b),1.0e-300}));
    };
    for(int i=1;i<left.Nx()-1;++i)
    for(int j=1;j<left.Ny()-1;++j)
    for(int k=1;k<left.Nz()-1;++k)
    {
        const ConservedVars& a=left.U(i,j,k);
        const ConservedVars& b=right.U(i,j,k);
        compare(a.rho,b.rho);
        compare(a.rhoE,b.rhoE);
        for(std::size_t component=0;component<a.rhoU.size();++component)
            compare(a.rhoU[component],b.rhoU[component]);
        for(std::size_t component=0;component<a.rhoY.size();++component)
            compare(a.rhoY[component],b.rhoY[component]);
        compare(a.rhoZ,b.rhoZ);
        compare(a.rhoC,b.rhoC);
    }
    return maximum;
}

} // namespace

int main()
{
    SpeciesRegistry registry;
    const std::string profile=vortexProfilePath();
    SectorMeshMetrics fineMesh(profile,10,8,6,120.0);
    SectorMeshMetrics coarseMesh(profile,6,5,4,120.0);
    Domain fineBase(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    Domain fineCorrected(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    Domain coarseBase(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    Domain coarseImproved(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    Domain coarseCheck(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    for(Domain* domain:{&fineBase,&fineCorrected,&coarseBase,
                       &coarseImproved,&coarseCheck})
        domain->setSpecies(registry.names());
    if(!initialize(fineBase,registry))
    { std::cerr << "fine initialization failed\n"; return 1; }
    if(!fineCorrected.replaceConservedField(capture(fineBase)))
    { std::cerr << "fine clone failed\n"; return 1; }

    SectorConservativeGridTransfer::apply(
        fineBase,fineMesh,nullptr,coarseBase,coarseMesh,nullptr);
    std::vector<ConservedVars> improvedField=capture(coarseBase);
    std::vector<double> temperature(
        static_cast<std::size_t>(coarseBase.N()),700.0);
    for(int cell=0;cell<coarseBase.N();++cell)
    {
        // A spatially varying but composition/velocity/energy-consistent
        // density correction remains fully admissible.
        const int i=cell%coarseBase.Nx();
        const double factor=1.0+0.01
            *static_cast<double>(i)/std::max(coarseBase.Nx()-1,1);
        scaleState(improvedField[static_cast<std::size_t>(cell)],factor);
        temperature[static_cast<std::size_t>(cell)]=coarseBase.W(cell).T;
    }
    if(!coarseImproved.replaceConservedField(improvedField,temperature))
    { std::cerr << "coarse improved recovery failed\n"; return 1; }

    const SectorCoarseCorrection::Result correction=
        SectorCoarseCorrection::apply(
            fineBase,fineMesh,coarseBase,coarseImproved,coarseMesh,
            fineCorrected,0,24);
    if(correction.acceptedGlobalScale!=1.0
        ||correction.localLimitedCells!=0
        ||correction.recoveryBacktracks!=0
        ||correction.maximumInventoryDefectRelativeToBase>2.0e-12)
    {
        std::cerr << "correction gate failed scale/limited/backtrack/defect="
                  << correction.acceptedGlobalScale << "/"
                  << correction.localLimitedCells << "/"
                  << correction.recoveryBacktracks << "/"
                  << correction.maximumInventoryDefectRelativeToBase << "\n";
        return 1;
    }

    SectorConservativeGridTransfer::apply(
        fineCorrected,fineMesh,nullptr,coarseCheck,coarseMesh,nullptr);
    const double coarseRecoveryError=maximumFieldDifference(
        coarseImproved,coarseCheck);
    if(coarseRecoveryError>3.0e-12)
    {
        std::cerr << "coarse recovery error=" << coarseRecoveryError << "\n";
        return 1;
    }

    Domain fineTraceBase(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    Domain fineTraceCorrected(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    Domain coarseTraceBase(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    Domain coarseTraceImproved(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    for(Domain* domain:{&fineTraceBase,&fineTraceCorrected,
                       &coarseTraceBase,&coarseTraceImproved})
        domain->setSpecies(registry.names());
    if(!initialize(fineTraceBase,registry,true))
    { std::cerr << "trace initialization failed\n"; return 1; }
    if(!fineTraceCorrected.replaceConservedField(capture(fineTraceBase)))
    { std::cerr << "trace clone failed\n"; return 1; }
    SectorConservativeGridTransfer::apply(
        fineTraceBase,fineMesh,nullptr,
        coarseTraceBase,coarseMesh,nullptr);
    std::vector<ConservedVars> traceImproved=capture(coarseTraceBase);
    std::vector<double> traceTemperature(
        static_cast<std::size_t>(coarseTraceBase.N()),700.0);
    const int oxygen=registry.index("O2");
    const int nitrogen=registry.index("N2");
    for(int cell=0;cell<coarseTraceBase.N();++cell)
    {
        ConservedVars& state=traceImproved[static_cast<std::size_t>(cell)];
        const double shift=0.01*state.rho;
        state.rhoY[static_cast<std::size_t>(oxygen)]-=shift;
        state.rhoY[static_cast<std::size_t>(nitrogen)]+=shift;
        traceTemperature[static_cast<std::size_t>(cell)]=
            coarseTraceBase.W(cell).T;
    }
    if(!coarseTraceImproved.replaceConservedField(
            traceImproved,traceTemperature))
    { std::cerr << "trace coarse improved recovery failed\n"; return 1; }
    const SectorCoarseCorrection::Result traceCorrection=
        SectorCoarseCorrection::apply(
            fineTraceBase,fineMesh,coarseTraceBase,coarseTraceImproved,
            coarseMesh,fineTraceCorrected,0,24);
    if(traceCorrection.acceptedGlobalScale!=1.0
        ||traceCorrection.localLimitedCells!=0
        ||traceCorrection.recoveryBacktracks!=0
        ||traceCorrection.compositionProjectedCells<=0
        ||traceCorrection.maximumCompositionProjection<=0.0
        ||traceCorrection.maximumFlowInventoryDefectRelativeToBase>2.0e-12
        ||traceCorrection.maximumSpeciesInventoryDefectRelativeToBase<=0.0)
    {
        std::cerr << "trace projection gate failed scale/local/backtrack/projected/"
                     "projection/flow/species="
                  << traceCorrection.acceptedGlobalScale << "/"
                  << traceCorrection.localLimitedCells << "/"
                  << traceCorrection.recoveryBacktracks << "/"
                  << traceCorrection.compositionProjectedCells << "/"
                  << traceCorrection.maximumCompositionProjection << "/"
                  << traceCorrection.maximumFlowInventoryDefectRelativeToBase
                  << "/"
                  << traceCorrection.maximumSpeciesInventoryDefectRelativeToBase
                  << "\n";
        return 1;
    }
    for(int i=1;i<fineTraceCorrected.Nx()-1;++i)
    for(int j=1;j<fineTraceCorrected.Ny()-1;++j)
    for(int k=1;k<fineTraceCorrected.Nz()-1;++k)
    {
        const ConservedVars& state=fineTraceCorrected.U(i,j,k);
        double speciesMass=0.0;
        for(double value:state.rhoY)
        {
            if(value<0.0||!std::isfinite(value)) return 1;
            speciesMass+=value;
        }
        if(std::abs(speciesMass-state.rho)>2.0e-12*state.rho) return 1;
    }

    Domain coarseLargeImproved(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    Domain fineTrusted(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    coarseLargeImproved.setSpecies(registry.names());
    fineTrusted.setSpecies(registry.names());
    std::vector<ConservedVars> largeImprovement=capture(coarseBase);
    std::vector<double> largeTemperature(
        static_cast<std::size_t>(coarseBase.N()),700.0);
    for(int cell=0;cell<coarseBase.N();++cell)
    {
        scaleState(largeImprovement[static_cast<std::size_t>(cell)],2.0);
        largeTemperature[static_cast<std::size_t>(cell)]=coarseBase.W(cell).T;
    }
    if(!coarseLargeImproved.replaceConservedField(
            largeImprovement,largeTemperature)
        ||!fineTrusted.replaceConservedField(capture(fineBase)))
    { std::cerr << "trust-region setup failed\n"; return 1; }
    const SectorCoarseCorrection::Result trustedCorrection=
        SectorCoarseCorrection::apply(
            fineBase,fineMesh,coarseBase,coarseLargeImproved,coarseMesh,
            fineTrusted,4,24,0.05,true);
    if(!(trustedCorrection.maximumUnscaledRelativeCorrection>1.0)
        ||std::abs(trustedCorrection.trustRegionScale
                *trustedCorrection.maximumUnscaledRelativeCorrection-0.05)
            >2.0e-12
        ||std::abs(trustedCorrection.acceptedGlobalScale
                -trustedCorrection.trustRegionScale)>2.0e-12
        ||trustedCorrection.recoveryBacktracks!=0
        ||trustedCorrection.localLimitedCells!=0
        ||trustedCorrection.boundaryMaskedFineCells!=160
        ||trustedCorrection.maximumFlowInventoryDefectRelativeToBase>2.0e-12)
    {
        std::cerr << "trust gate failed unscaled/trust/accepted/backtrack/local/"
                     "masked/flow="
                  << trustedCorrection.maximumUnscaledRelativeCorrection << "/"
                  << trustedCorrection.trustRegionScale << "/"
                  << trustedCorrection.acceptedGlobalScale << "/"
                  << trustedCorrection.recoveryBacktracks << "/"
                  << trustedCorrection.localLimitedCells << "/"
                  << trustedCorrection.boundaryMaskedFineCells << "/"
                  << trustedCorrection.maximumFlowInventoryDefectRelativeToBase
                  << "\n";
        return 1;
    }
    bool interiorChanged=false;
    for(int i=1;i<fineTrusted.Nx()-1;++i)
    for(int j=1;j<fineTrusted.Ny()-1;++j)
    for(int k=1;k<fineTrusted.Nz()-1;++k)
    {
        const bool masked=i<=2||i>=fineTrusted.Nx()-3
            ||j<=2||j>=fineTrusted.Ny()-3;
        const double relativeDensityChange=std::abs(
            fineTrusted.U(i,j,k).rho-fineBase.U(i,j,k).rho)
            /fineBase.U(i,j,k).rho;
        if(masked&&relativeDensityChange>2.0e-12)
        {
            std::cerr << "masked fine cell changed\n";
            return 1;
        }
        if(!masked&&relativeDensityChange>1.0e-3)
            interiorChanged=true;
    }
    if(!interiorChanged)
    { std::cerr << "trusted interior correction was not applied\n"; return 1; }

    Domain coarseFrozenImproved(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    Domain fineCompositionFrozen(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    coarseFrozenImproved.setSpecies(registry.names());
    fineCompositionFrozen.setSpecies(registry.names());
    std::vector<ConservedVars> frozenImprovement=capture(coarseTraceBase);
    std::vector<double> frozenTemperature(
        static_cast<std::size_t>(coarseTraceBase.N()),700.0);
    for(int cell=0;cell<coarseTraceBase.N();++cell)
    {
        ConservedVars& state=frozenImprovement[static_cast<std::size_t>(cell)];
        scaleState(state,1.02);
        const double shift=0.01*state.rho;
        state.rhoY[static_cast<std::size_t>(oxygen)]-=shift;
        state.rhoY[static_cast<std::size_t>(nitrogen)]+=shift;
        frozenTemperature[static_cast<std::size_t>(cell)]=
            coarseTraceBase.W(cell).T;
    }
    if(!coarseFrozenImproved.replaceConservedField(
            frozenImprovement,frozenTemperature)
        ||!fineCompositionFrozen.replaceConservedField(capture(fineTraceBase)))
    { std::cerr << "frozen-composition setup failed\n"; return 1; }
    const SectorCoarseCorrection::Result frozenCorrection=
        SectorCoarseCorrection::apply(
            fineTraceBase,fineMesh,coarseTraceBase,coarseFrozenImproved,
            coarseMesh,fineCompositionFrozen,0,24,
            std::numeric_limits<double>::infinity(),false,1.0,true);
    bool frozenDensityChanged=false;
    for(int i=1;i<fineCompositionFrozen.Nx()-1;++i)
    for(int j=1;j<fineCompositionFrozen.Ny()-1;++j)
    for(int k=1;k<fineCompositionFrozen.Nz()-1;++k)
    {
        const ConservedVars& base=fineTraceBase.U(i,j,k);
        const ConservedVars& state=fineCompositionFrozen.U(i,j,k);
        if(std::abs(state.rho-base.rho)>1.0e-6*base.rho)
            frozenDensityChanged=true;
        for(std::size_t species=0;species<state.rhoY.size();++species)
        {
            const double baseFraction=base.rhoY[species]/base.rho;
            const double correctedFraction=state.rhoY[species]/state.rho;
            if(std::abs(baseFraction-correctedFraction)>2.0e-12)
            {
                std::cerr << "frozen fine composition changed\n";
                return 1;
            }
        }
    }
    if(frozenCorrection.acceptedGlobalScale!=1.0
        ||frozenCorrection.recoveryBacktracks!=0
        ||frozenCorrection.compositionFrozenCells<=0
        ||frozenCorrection.compositionProjectedCells!=0
        ||frozenCorrection.maximumFlowInventoryDefectRelativeToBase>2.0e-12
        ||frozenCorrection.maximumSpeciesInventoryDefectRelativeToBase<=0.0
        ||!frozenDensityChanged)
    {
        std::cerr << "frozen composition gate failed scale/backtrack/frozen/"
                     "projected/flow/species/density="
                  << frozenCorrection.acceptedGlobalScale << "/"
                  << frozenCorrection.recoveryBacktracks << "/"
                  << frozenCorrection.compositionFrozenCells << "/"
                  << frozenCorrection.compositionProjectedCells << "/"
                  << frozenCorrection.maximumFlowInventoryDefectRelativeToBase
                  << "/"
                  << frozenCorrection.maximumSpeciesInventoryDefectRelativeToBase
                  << "/" << frozenDensityChanged << "\n";
        return 1;
    }

    Domain fineParentComposition(10,8,6,fineMesh.length(),0.065,
        fineMesh.sectorAngle(),registry.size());
    Domain coarseParentCheck(6,5,4,coarseMesh.length(),0.065,
        coarseMesh.sectorAngle(),registry.size());
    fineParentComposition.setSpecies(registry.names());
    coarseParentCheck.setSpecies(registry.names());
    if(!fineParentComposition.replaceConservedField(capture(fineTraceBase)))
    { std::cerr << "parent composition clone failed\n"; return 1; }
    const SectorCoarseCorrection::Result parentCorrection=
        SectorCoarseCorrection::apply(
            fineTraceBase,fineMesh,coarseTraceBase,coarseFrozenImproved,
            coarseMesh,fineParentComposition,0,24,
            std::numeric_limits<double>::infinity(),false,1.0,false,true);
    SectorConservativeGridTransfer::apply(
        fineParentComposition,fineMesh,nullptr,
        coarseParentCheck,coarseMesh,nullptr);
    const double parentRecoveryError=maximumFieldDifference(
        coarseFrozenImproved,coarseParentCheck);
    for(int i=1;i<fineParentComposition.Nx()-1;++i)
    for(int j=1;j<fineParentComposition.Ny()-1;++j)
    for(int k=1;k<fineParentComposition.Nz()-1;++k)
    {
        const ConservedVars& state=fineParentComposition.U(i,j,k);
        double speciesMass=0.0;
        for(double value:state.rhoY)
        {
            if(!std::isfinite(value)||value<0.0)
            {
                std::cerr << "parent composition is not nonnegative\n";
                return 1;
            }
            speciesMass+=value;
        }
        if(std::abs(speciesMass-state.rho)>2.0e-12*state.rho)
        {
            std::cerr << "parent composition is not conservative\n";
            return 1;
        }
    }
    if(parentCorrection.acceptedGlobalScale!=1.0
        ||parentCorrection.recoveryBacktracks!=0
        ||parentCorrection.compositionReconstructedCells!=192
        ||parentCorrection.compositionProjectedCells!=0
        ||parentCorrection.compositionFrozenCells!=0
        ||parentCorrection.maximumInventoryDefectRelativeToBase>2.0e-12
        ||parentRecoveryError>3.0e-12)
    {
        std::cerr << "parent composition gate failed scale/backtrack/reconstructed/"
                     "projected/frozen/defect/recovery="
                  << parentCorrection.acceptedGlobalScale << "/"
                  << parentCorrection.recoveryBacktracks << "/"
                  << parentCorrection.compositionReconstructedCells << "/"
                  << parentCorrection.compositionProjectedCells << "/"
                  << parentCorrection.compositionFrozenCells << "/"
                  << parentCorrection.maximumInventoryDefectRelativeToBase
                  << "/" << parentRecoveryError << "\n";
        return 1;
    }

    std::cout << "coarse_recovery_error=" << coarseRecoveryError
              << " inventory_defect="
              << correction.maximumInventoryDefectRelativeToBase << "\n"
              << "trace_projected/max_projection/flow_defect/species_defect="
              << traceCorrection.compositionProjectedCells << "/"
              << traceCorrection.maximumCompositionProjection << "/"
              << traceCorrection.maximumFlowInventoryDefectRelativeToBase
              << "/"
              << traceCorrection.maximumSpeciesInventoryDefectRelativeToBase
              << "\ntrust_unscaled/scale/masked="
              << trustedCorrection.maximumUnscaledRelativeCorrection << "/"
              << trustedCorrection.trustRegionScale << "/"
              << trustedCorrection.boundaryMaskedFineCells
              << "\nfrozen_composition_cells/flow_defect/species_defect="
              << frozenCorrection.compositionFrozenCells << "/"
              << frozenCorrection.maximumFlowInventoryDefectRelativeToBase
              << "/"
              << frozenCorrection.maximumSpeciesInventoryDefectRelativeToBase
              << "\nparent_reconstructed/recovery_error/inventory_defect="
              << parentCorrection.compositionReconstructedCells << "/"
              << parentRecoveryError << "/"
              << parentCorrection.maximumInventoryDefectRelativeToBase
              << "\nSECTOR_COARSE_CORRECTION_PASS\n";
    return 0;
}
