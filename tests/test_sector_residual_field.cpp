#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "../src/framework/SectorResidualField.H"

using namespace Foam;

namespace
{

bool initialize(Domain& domain)
{
    domain.setSpecies({"O2","C3H8"});
    for(int cell=0;cell<domain.N();++cell)
    {
        PrimitiveVars primitive(2);
        primitive.p=1.0e6;
        primitive.T=600.0;
        primitive.U={10.0,0.0,0.0};
        primitive.Y={0.7,0.3};
        domain.stagedW(cell)=primitive;
    }
    return domain.primitiveToConserved();
}

bool same(const ConservedVars& left,const ConservedVars& right)
{
    return left.rho==right.rho&&left.rhoU==right.rhoU
        &&left.rhoE==right.rhoE&&left.rhoY==right.rhoY
        &&left.rhoZ==right.rhoZ&&left.rhoC==right.rhoC;
}

} // namespace

int main()
{
    Domain domain(4,0.04,2);
    if(!initialize(domain)) return 1;
    SectorResidualField::Data source;
    source.metadata.kind=SectorResidualField::Kind::FAS_SOURCE;
    source.metadata.nx=domain.Nx();
    source.metadata.ny=domain.Ny();
    source.metadata.nz=domain.Nz();
    source.metadata.nSpecies=domain.nSpecies();
    source.metadata.referenceStep=1.0e-7;
    source.metadata.time=3.0e-4;
    source.metadata.iteration=17;
    source.metadata.profileHash=1234;
    source.metadata.stateCheckpointHash=5678;
    source.metadata.fineStateCheckpointHash=9012;
    source.metadata.mode="frozen-chemistry-les";
    source.metadata.gitCommit="test-commit";
    source.metadata.thermoModelId=domain.thermoModelId();
    source.metadata.species=domain.speciesNames();
    source.active={0U,1U,1U,0U};
    source.field.assign(
        static_cast<std::size_t>(domain.N()),
        ConservedVars(domain.nSpecies()));
    for(int cell=0;cell<domain.N();++cell)
    {
        ConservedVars& value=source.field[static_cast<std::size_t>(cell)];
        value.rho=-1.0-cell;
        value.rhoE=2.0+cell;
        value.rhoU={-3.0,4.0,-5.0};
        value.rhoY={-0.25,0.25};
        value.rhoZ=-6.0;
        value.rhoC=7.0;
    }

    const std::filesystem::path path=
        std::filesystem::temp_directory_path()
        /"rfam_sector_residual_field_test.bin";
    SectorResidualField::write(path.string(),source);
    const SectorResidualField::Data loaded=
        SectorResidualField::read(path.string());
    std::filesystem::remove(path);

    if(loaded.metadata.kind!=source.metadata.kind
        ||loaded.metadata.nx!=source.metadata.nx
        ||loaded.metadata.referenceStep!=source.metadata.referenceStep
        ||loaded.metadata.time!=source.metadata.time
        ||loaded.metadata.iteration!=source.metadata.iteration
        ||loaded.metadata.profileHash!=source.metadata.profileHash
        ||loaded.metadata.stateCheckpointHash
            !=source.metadata.stateCheckpointHash
        ||loaded.metadata.fineStateCheckpointHash
            !=source.metadata.fineStateCheckpointHash
        ||loaded.metadata.mode!=source.metadata.mode
        ||loaded.metadata.gitCommit!=source.metadata.gitCommit
        ||loaded.metadata.thermoModelId!=source.metadata.thermoModelId
        ||loaded.metadata.species!=source.metadata.species
        ||loaded.active!=source.active
        ||loaded.field.size()!=source.field.size())
        return 1;
    for(std::size_t cell=0;cell<source.field.size();++cell)
        if(!same(loaded.field[cell],source.field[cell])) return 1;

    SectorResidualField::requireCompatible(
        loaded,domain,1234,"frozen-chemistry-les",5678,
        SectorResidualField::Kind::FAS_SOURCE);
    bool mismatchRejected=false;
    try
    {
        SectorResidualField::requireCompatible(
            loaded,domain,1234,"frozen-chemistry-les",5679,
            SectorResidualField::Kind::FAS_SOURCE);
    }
    catch(const std::runtime_error&)
    {
        mismatchRejected=true;
    }
    if(!mismatchRejected) return 1;

    std::cout << "cells=" << loaded.field.size()
              << " negative_rho_rate=" << loaded.field[1].rho
              << "\nSECTOR_RESIDUAL_FIELD_PASS\n";
    return 0;
}
