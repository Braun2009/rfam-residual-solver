#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include "src/framework/Domain.H"
#include "src/framework/SectorCheckpoint.H"
#include "src/framework/SolverBase.H"
#include "src/mesh/SectorMeshMetrics.H"
#include "src/thermo/SpeciesRegistry.H"

using namespace Foam;

namespace
{

struct DomainSnapshot
{
    std::vector<ConservedVars> conserved;
    std::vector<ConservedVars> candidate;
    std::vector<PrimitiveVars> primitive;
    std::vector<thermo::EquilibriumState> equilibrium;
};

bool sameScalar(double a,double b)
{
    return a==b||(std::isnan(a)&&std::isnan(b));
}

bool sameVector(const std::vector<double>& a,const std::vector<double>& b)
{
    if(a.size()!=b.size()) return false;
    for(std::size_t i=0;i<a.size();++i)
        if(!sameScalar(a[i],b[i])) return false;
    return true;
}

bool samePhase(const thermo::EquilibriumPhaseProperties& a,
               const thermo::EquilibriumPhaseProperties& b)
{
    return sameVector(a.massFractions,b.massFractions)
        &&sameVector(a.moleFractions,b.moleFractions)
        &&sameScalar(a.molecularWeight,b.molecularWeight)
        &&sameScalar(a.rho,b.rho)&&sameScalar(a.h,b.h)
        &&sameScalar(a.e,b.e)&&sameScalar(a.s,b.s)
        &&sameScalar(a.compressibility,b.compressibility)
        &&a.entropyAvailable==b.entropyAvailable;
}

bool sameEquilibrium(const thermo::EquilibriumState& a,
                     const thermo::EquilibriumState& b)
{
    if(a.liquid.has_value()!=b.liquid.has_value()
       ||a.vapor.has_value()!=b.vapor.has_value()
       ||a.singlePhase.has_value()!=b.singlePhase.has_value()) return false;
    return a.region==b.region
        &&sameVector(a.overallMassFractions,b.overallMassFractions)
        &&sameScalar(a.p,b.p)&&sameScalar(a.T,b.T)
        &&sameScalar(a.rho,b.rho)&&sameScalar(a.h,b.h)
        &&sameScalar(a.e,b.e)&&sameScalar(a.s,b.s)
        &&sameScalar(a.vaporMoleFraction,b.vaporMoleFraction)
        &&sameScalar(a.vaporMassFraction,b.vaporMassFraction)
        &&sameScalar(a.vaporVolumeFraction,b.vaporVolumeFraction)
        &&a.entropyAvailable==b.entropyAvailable
        &&a.acoustic.closure==b.acoustic.closure
        &&sameScalar(a.acoustic.pressureDensityDerivative,
            b.acoustic.pressureDensityDerivative)
        &&(!a.liquid||samePhase(*a.liquid,*b.liquid))
        &&(!a.vapor||samePhase(*a.vapor,*b.vapor))
        &&(!a.singlePhase||samePhase(*a.singlePhase,*b.singlePhase));
}

bool sameConserved(const ConservedVars& a,const ConservedVars& b)
{
    return a.rho==b.rho&&a.rhoU==b.rhoU&&a.rhoE==b.rhoE
        &&a.rhoY==b.rhoY&&a.rhoZ==b.rhoZ&&a.rhoC==b.rhoC;
}

bool samePrimitive(const PrimitiveVars& a,const PrimitiveVars& b)
{
    return a.p==b.p&&a.T==b.T&&a.U==b.U&&a.Y==b.Y
        &&a.Z==b.Z&&a.C==b.C;
}

DomainSnapshot capture(const Domain& domain)
{
    DomainSnapshot snapshot;
    snapshot.conserved.reserve(domain.N());
    snapshot.candidate.reserve(domain.N());
    snapshot.primitive.reserve(domain.N());
    snapshot.equilibrium.reserve(domain.N());
    for(int i=0;i<domain.N();++i)
    {
        snapshot.conserved.push_back(domain.U(i));
        snapshot.candidate.push_back(domain.Unew(i));
        snapshot.primitive.push_back(domain.W(i));
        snapshot.equilibrium.push_back(domain.Q(i));
    }
    return snapshot;
}

bool unchanged(const Domain& domain,const DomainSnapshot& snapshot)
{
    if(domain.N()!=static_cast<int>(snapshot.conserved.size())) return false;
    for(int i=0;i<domain.N();++i)
        if(!sameConserved(domain.U(i),snapshot.conserved[i])
           ||!sameConserved(domain.Unew(i),snapshot.candidate[i])
           ||!samePrimitive(domain.W(i),snapshot.primitive[i])
           ||!sameEquilibrium(domain.Q(i),snapshot.equilibrium[i])) return false;
    return true;
}

void initialize(Domain& domain,const std::vector<double>& composition,
                double pressure,double temperature)
{
    for(int i=0;i<domain.Nx();++i) for(int j=0;j<domain.Ny();++j)
    for(int k=0;k<domain.Nz();++k)
    {
        PrimitiveVars& w=domain.stagedW(i,j,k);
        w.p=pressure+i;w.T=temperature+j;w.U={double(k),1.0,-1.0};
        w.Y=composition;
    }
    if(!domain.primitiveToConserved())
        throw std::runtime_error("checkpoint test initialization failed");
}

}

int main()
{
    const std::string path="test_sector_checkpoint.bin";
    const std::string noTurbulencePath="test_sector_checkpoint_no_turbulence.bin";
    const std::string missingSpeciesPath="test_sector_checkpoint_missing_species.bin";
    const std::string invalidWritePath="test_sector_checkpoint_invalid_write.bin";
    std::remove(path.c_str());
    std::remove(noTurbulencePath.c_str());
    std::remove(missingSpeciesPath.c_str());
    std::remove(invalidWritePath.c_str());
    SpeciesRegistry registry; const int nx=8,ny=6,nz=6;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,ny,nz,120.0);

    bool externalLayoutRejected=false;
    try
    {
        const auto emptyModel=thermo::makeLegacySinglePhaseThermoModel();
        Domain invalid(1,1.0,2,emptyModel);
    }
    catch(const std::invalid_argument&) { externalLayoutRejected=true; }
    if(!externalLayoutRejected) return 1;

    Domain missingSpecies(1,1.0,2);
    PrimitiveVars& unlabeled=missingSpecies.stagedW(0);
    unlabeled.p=4e6;unlabeled.T=700.0;unlabeled.Y={1.0,0.0};
    if(!missingSpecies.primitiveToConserved()) return 1;
    bool missingSpeciesRejected=false;
    try
    {
        SectorCheckpoint::write(
            missingSpeciesPath,missingSpecies,SectorCheckpointMetadata{});
    }
    catch(const std::runtime_error&) { missingSpeciesRejected=true; }
    if(!missingSpeciesRejected||std::filesystem::exists(missingSpeciesPath)) return 1;

    Domain source(nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    source.setSpecies(registry.names());
    std::vector<double> y(registry.size(),0.0); y[registry.index("O2")]=1.0;
    for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
    { auto& w=source.stagedW(i,j,k);w.p=4e6+i;w.T=700+j;w.U={double(k),0,0};w.Y=y; }
    if(!source.primitiveToConserved()) return 1;
    SolverBase restartAdmission(source);
    bool solverInvalidTimeRejected=false;
    try
    {
        restartAdmission.setRestartState(
            std::numeric_limits<double>::quiet_NaN(),42);
    }
    catch(const std::invalid_argument&) { solverInvalidTimeRejected=true; }
    bool solverInvalidIterationRejected=false;
    try { restartAdmission.setRestartState(1.0,-1); }
    catch(const std::invalid_argument&) { solverInvalidIterationRejected=true; }
    bool solverNegativeZeroRejected=false;
    try { restartAdmission.setRestartState(-0.0,42); }
    catch(const std::invalid_argument&) { solverNegativeZeroRejected=true; }
    if(!solverInvalidTimeRejected||!solverInvalidIterationRejected
       ||!solverNegativeZeroRejected||restartAdmission.time()!=0.0
       ||std::signbit(restartAdmission.time())
       ||restartAdmission.iteration()!=0) return 1;
    SectorTurbulenceOperator sourceTurbulence(mesh,registry);
    sourceTurbulence.apply(source,1e-10);
    SectorCheckpointMetadata metadata{1.25e-5,42,8,"rans-finite-rate","deadbeef",12345};
    SectorCheckpointMetadata invalidWriteMetadata=metadata;
    invalidWriteMetadata.time=-0.0;
    bool invalidWriteRejected=false;
    try
    {
        SectorCheckpoint::write(
            invalidWritePath,source,invalidWriteMetadata,&sourceTurbulence);
    }
    catch(const std::runtime_error&) { invalidWriteRejected=true; }
    bool invalidWriteCandidateExists=false;
    for(const auto& entry:std::filesystem::directory_iterator("."))
        if(entry.path().filename().string().find(invalidWritePath+".tmp.")==0)
            invalidWriteCandidateExists=true;
    if(!invalidWriteRejected||std::filesystem::exists(invalidWritePath)
       ||invalidWriteCandidateExists) return 1;
    SectorCheckpoint::write(path,source,metadata,&sourceTurbulence);
    SectorCheckpointMetadata replacement=metadata;
    replacement.iteration=43;
    SectorCheckpoint::write(path,source,replacement,&sourceTurbulence);
    const SectorCheckpointInfo replaced=SectorCheckpoint::inspect(path);
    if(replaced.metadata.iteration!=replacement.iteration) return 1;
    for(const auto& entry:std::filesystem::directory_iterator("."))
        if(entry.path().filename().string().find(path+".tmp.")==0) return 1;

    auto corruptClock = [&](const std::string& corruptedPath,
                            double time,int iteration)
    {
        std::filesystem::copy_file(
            path,corruptedPath,std::filesystem::copy_options::overwrite_existing);
        std::fstream stream(
            corruptedPath,std::ios::binary|std::ios::in|std::ios::out);
        // magic(8), version(4), nx/ny/nz/ns(16), then metadata time and iteration.
        stream.seekp(28);
        stream.write(reinterpret_cast<const char*>(&time),sizeof(time));
        stream.seekp(36);
        stream.write(reinterpret_cast<const char*>(&iteration),sizeof(iteration));
        if(!stream) throw std::runtime_error("checkpoint clock mutation failed");
    };
    const std::string invalidClockPath="test_sector_checkpoint_invalid_clock.bin";
    corruptClock(
        invalidClockPath,std::numeric_limits<double>::quiet_NaN(),
        replacement.iteration);
    bool invalidTimeRejected=false;
    try { SectorCheckpoint::inspect(invalidClockPath); }
    catch(const std::runtime_error&) { invalidTimeRejected=true; }
    std::remove(invalidClockPath.c_str());
    if(!invalidTimeRejected) return 1;

    corruptClock(invalidClockPath,-0.0,replacement.iteration);
    bool negativeZeroInspectRejected=false;
    try { SectorCheckpoint::inspect(invalidClockPath); }
    catch(const std::runtime_error&) { negativeZeroInspectRejected=true; }
    Domain negativeZeroTarget(
        nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    negativeZeroTarget.setSpecies(registry.names());
    initialize(negativeZeroTarget,y,3.5e6,900.0);
    const DomainSnapshot negativeZeroBefore=capture(negativeZeroTarget);
    bool negativeZeroReadRejected=false;
    try { SectorCheckpoint::read(invalidClockPath,negativeZeroTarget); }
    catch(const std::runtime_error&) { negativeZeroReadRejected=true; }
    std::remove(invalidClockPath.c_str());
    if(!negativeZeroInspectRejected||!negativeZeroReadRejected
       ||!unchanged(negativeZeroTarget,negativeZeroBefore)) return 1;

    corruptClock(invalidClockPath,replacement.time,-1);
    Domain invalidClockTarget(
        nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    invalidClockTarget.setSpecies(registry.names());
    initialize(invalidClockTarget,y,3.5e6,900.0);
    const DomainSnapshot invalidClockBefore=capture(invalidClockTarget);
    bool invalidIterationRejected=false;
    try { SectorCheckpoint::read(invalidClockPath,invalidClockTarget); }
    catch(const std::runtime_error&) { invalidIterationRejected=true; }
    std::remove(invalidClockPath.c_str());
    if(!invalidIterationRejected
       ||!unchanged(invalidClockTarget,invalidClockBefore)) return 1;

    Domain mismatchTarget(
        nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    mismatchTarget.setSpecies(registry.names());
    initialize(mismatchTarget,y,3.5e6,900.0);
    const DomainSnapshot mismatchBefore=capture(mismatchTarget);
    bool missingTurbulenceTargetRejected=false;
    try { SectorCheckpoint::read(path,mismatchTarget); }
    catch(const std::runtime_error&) { missingTurbulenceTargetRejected=true; }
    if(!missingTurbulenceTargetRejected
       ||!unchanged(mismatchTarget,mismatchBefore)) return 1;

    SectorCheckpoint::write(noTurbulencePath,source,replacement);
    Domain unexpectedTurbulenceTarget(
        nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    unexpectedTurbulenceTarget.setSpecies(registry.names());
    initialize(unexpectedTurbulenceTarget,y,3.25e6,850.0);
    const DomainSnapshot unexpectedBefore=capture(unexpectedTurbulenceTarget);
    SectorTurbulenceOperator unexpectedTurbulence(mesh,registry);
    const std::vector<double> turbulenceKBefore=unexpectedTurbulence.kField();
    const std::vector<double> turbulenceOmegaBefore=unexpectedTurbulence.omegaField();
    bool unexpectedTurbulenceTargetRejected=false;
    try
    {
        SectorCheckpoint::read(
            noTurbulencePath,unexpectedTurbulenceTarget,&unexpectedTurbulence);
    }
    catch(const std::runtime_error&) { unexpectedTurbulenceTargetRejected=true; }
    std::remove(noTurbulencePath.c_str());
    if(!unexpectedTurbulenceTargetRejected
       ||!unchanged(unexpectedTurbulenceTarget,unexpectedBefore)
       ||unexpectedTurbulence.kField()!=turbulenceKBefore
       ||unexpectedTurbulence.omegaField()!=turbulenceOmegaBefore) return 1;

    const std::string legacyPath="test_sector_checkpoint_unprovenanced.bin";
    std::filesystem::copy_file(path,legacyPath,
        std::filesystem::copy_options::overwrite_existing);
    {
        std::fstream legacy(legacyPath,
            std::ios::binary|std::ios::in|std::ios::out);
        const std::uint32_t legacyVersion=2;
        legacy.seekp(8);
        legacy.write(reinterpret_cast<const char*>(&legacyVersion),
            sizeof(legacyVersion));
        if(!legacy) return 1;
    }
    Domain legacyTarget(
        nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    legacyTarget.setSpecies(registry.names());
    bool legacyRejected=false;
    try { SectorCheckpoint::read(legacyPath,legacyTarget); }
    catch(const std::runtime_error&) { legacyRejected=true; }
    std::remove(legacyPath.c_str());
    if(!legacyRejected||legacyTarget.stateReady()) return 1;

    Domain restored(nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
    restored.setSpecies(registry.names()); SectorTurbulenceOperator restoredTurbulence(mesh,registry);
    const auto loaded=SectorCheckpoint::read(path,restored,&restoredTurbulence);
    std::remove(path.c_str());
    if(loaded.iteration!=replacement.iteration||loaded.time!=metadata.time
       ||loaded.mode!=metadata.mode||loaded.gitCommit!=metadata.gitCommit
       ||loaded.profileHash!=metadata.profileHash) return 1;
    for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
    {
        const auto& a=source.U(i,j,k);const auto& b=restored.U(i,j,k);
        if(a.rho!=b.rho||a.rhoE!=b.rhoE||a.rhoY!=b.rhoY||a.rhoU!=b.rhoU) return 1;
        const auto& wa=source.W(i,j,k);const auto& wb=restored.W(i,j,k);
        if(std::abs(wa.T-wb.T)>1e-8*std::max(1.0,std::abs(wa.T)))
        { std::cerr<<"checkpoint temperature branch mismatch\n"; return 1; }
    }
    if(sourceTurbulence.kField()!=restoredTurbulence.kField()
       ||sourceTurbulence.omegaField()!=restoredTurbulence.omegaField()) return 1;
    std::cout<<"SECTOR_CHECKPOINT_PASS\n";return 0;
}
