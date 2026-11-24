#include "CondensedCarbonOutput.H"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace Foam
{
namespace
{

void validateOutputMesh(
    const Domain& domain,
    const SectorMeshMetrics& mesh,
    const SectorCondensedCarbonOperator& carbon,
    const char* outputKind)
{
    if(!carbon.meshMatches(domain)
       ||carbon.meshFingerprint()!=condensedCarbonSectorMeshFingerprint(mesh))
        throw std::runtime_error(
            std::string("condensed-carbon ")+outputKind+" output mesh mismatch");
}

} // namespace

void CondensedCarbonOutput::writeSectorFieldCsv(
    const std::string& path,
    const Domain& domain,
    const SectorMeshMetrics& mesh,
    const SectorCondensedCarbonOperator& carbon)
{
    validateOutputMesh(domain,mesh,carbon,"field CSV");
    std::ofstream out(path);
    if(!out)
        throw std::runtime_error(
            "cannot write condensed-carbon field: "+path);
    out<<"i,j,k,x_m,y_m,z_m,rho_soot_kg_m3,soot_mass_loading";
    if(carbon.numberDensityEnabled()) out<<",N_soot_particles_m3";
    out<<'\n'<<std::setprecision(12);
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const FaceAreaVector center=mesh.cellCenter(i,j,k);
        const double rhoSoot=carbon.rhoSoot(i,j,k);
        out<<i<<','<<j<<','<<k<<','<<center.x<<','<<center.y<<','
           <<center.z<<','<<rhoSoot<<','<<rhoSoot/domain.U(i,j,k).rho;
        if(carbon.numberDensityEnabled())
            out<<','<<carbon.numberDensity(i,j,k);
        out<<'\n';
    }
    if(!out)
        throw std::runtime_error(
            "condensed-carbon field write failed: "+path);
}

void CondensedCarbonOutput::writeWallInventoryCsv(
    const std::string& path,
    const Domain& domain,
    const SectorMeshMetrics& mesh,
    const SectorCondensedCarbonOperator& carbon)
{
    validateOutputMesh(domain,mesh,carbon,"wall CSV");
    std::ofstream out(path);
    if(!out)
        throw std::runtime_error(
            "cannot write wall carbon inventory: "+path);
    out<<"i,k,x_m,theta_rad,wall_carbon_kg_m2\n"
       <<std::setprecision(12);
    const int cellJ=domain.Ny()-2;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const FaceAreaVector center=mesh.cellCenter(i,cellJ,k);
        out<<i<<','<<k<<','<<center.x<<','<<mesh.thetaCenter(k)<<','
           <<carbon.wallCarbonMassPerArea(i,k)<<'\n';
    }
    if(!out)
        throw std::runtime_error(
            "wall carbon inventory write failed: "+path);
}

void CondensedCarbonOutput::writeSectorVtk(
    const std::string& path,
    const Domain& domain,
    const SectorMeshMetrics& mesh,
    const SectorCondensedCarbonOperator& carbon)
{
    validateOutputMesh(domain,mesh,carbon,"sector VTK");
    std::ofstream out(path);
    if(!out)
        throw std::runtime_error(
            "cannot write condensed-carbon sector VTK: "+path);

    const int cellNx=domain.Nx()-2;
    const int cellNr=domain.Ny()-2;
    const int cellNt=domain.Nz()-2;
    const int pointNx=cellNx+1;
    const int pointNr=cellNr+1;
    const int pointNt=cellNt+1;
    const std::size_t pointCount=static_cast<std::size_t>(pointNx)
        *pointNr*pointNt;
    const std::size_t cellCount=static_cast<std::size_t>(cellNx)
        *cellNr*cellNt;

    out<<"# vtk DataFile Version 3.0\n"
       <<"RFAM condensed carbon sector "<<carbon.meshFingerprint()<<"\n"
       <<"ASCII\n"
       <<"DATASET STRUCTURED_GRID\n"
       <<"DIMENSIONS "<<pointNx<<' '<<pointNr<<' '<<pointNt<<"\n"
       <<"POINTS "<<pointCount<<" double\n"
       <<std::setprecision(17);
    for(int k=1;k<domain.Nz();++k)
    for(int j=1;j<domain.Ny();++j)
    for(int i=1;i<domain.Nx();++i)
    {
        const double theta=mesh.thetaFace(k);
        const double radius=mesh.wallFace(i)*mesh.etaFace(j);
        out<<mesh.xFace(i)<<' '<<radius*std::cos(theta)<<' '
           <<radius*std::sin(theta)<<'\n';
    }

    out<<"CELL_DATA "<<cellCount<<"\n"
       <<"SCALARS rho_soot double 1\nLOOKUP_TABLE default\n";
    for(int k=1;k<domain.Nz()-1;++k)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int i=1;i<domain.Nx()-1;++i)
        out<<carbon.rhoSoot(i,j,k)<<'\n';

    if(carbon.numberDensityEnabled())
    {
        out<<"SCALARS N_soot double 1\nLOOKUP_TABLE default\n";
        for(int k=1;k<domain.Nz()-1;++k)
        for(int j=1;j<domain.Ny()-1;++j)
        for(int i=1;i<domain.Nx()-1;++i)
            out<<carbon.numberDensity(i,j,k)<<'\n';
    }
    if(!out)
        throw std::runtime_error(
            "condensed-carbon sector VTK write failed: "+path);
}

} // End namespace Foam
