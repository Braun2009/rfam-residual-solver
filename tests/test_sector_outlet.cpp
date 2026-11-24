#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "src/framework/Domain.H"
#include "src/framework/SectorNonReflectingOutlet.H"
#include "src/mesh/SectorMeshMetrics.H"

using namespace Foam;

namespace
{

class PrescribedAcousticModel final : public thermo::EquilibriumThermoModel
{
    std::vector<std::string> species_={"O2","C3H8"};
    thermo::LegacySinglePhaseThermoModel delegate_{species_};
    double derivative_;
    std::string id_="prescribed-sector-outlet-acoustic";

    thermo::EquilibriumResult prescribe(
        thermo::EquilibriumResult result) const
    {
        result.diagnostics.modelId=id_;
        if(!result.success()) return result;
        result.state->acoustic.pressureDensityDerivative=derivative_;
        result.state->thermodynamics.reset();
        if(!result.state->valid())
        {
            result.status=thermo::EquilibriumStatus::ACOUSTIC_CLOSURE_FAILED;
            result.state.reset();
        }
        return result;
    }

public:
    explicit PrescribedAcousticModel(double soundSpeed)
    : derivative_(soundSpeed*soundSpeed)
    {}

    thermo::EquilibriumResult recoverRhoE(
        double density,double internalEnergy,
        const std::vector<double>& massFractions,
        const thermo::EquilibriumSeed* seed) const override
    {
        return prescribe(delegate_.recoverRhoE(
            density,internalEnergy,massFractions,seed));
    }

    thermo::EquilibriumResult initializePTY(
        double pressure,double temperature,
        const std::vector<double>& massFractions) const override
    {
        return prescribe(delegate_.initializePTY(
            pressure,temperature,massFractions));
    }

    const std::string& modelId() const override { return id_; }
    const std::vector<std::string>& speciesNames() const override
    { return species_; }
    EnergyConvention energyConvention() const override
    {
        return EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
    }
    bool supportsPrimitivePTYInitialization() const override { return true; }
};

} // namespace

int main()
{
    const int nx=18,nr=10,nt=10;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
    Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2); w.p=4e6; w.T=700; w.U={0.0,0.0,0.0}; w.Y={0.8,0.2};
        d.stagedW(i,j,k)=w;
    }
    d.primitiveToConserved();
    SectorNonReflectingOutlet outlet(mesh,101325.0);
    outlet.apply(d,1e-8);
    for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        const auto& w=d.W(nx-1,j,k);
        if(std::abs(w.p-101325.0)>1e-8||w.U[0]<0.0||d.U(nx-1,j,k).rho<=0.0) return 1;
        PrimitiveVars flowing=w; flowing.p=2e5; flowing.U[0]=100.0;
        d.setCell(nx-1,j,k,flowing);
    }
    outlet.apply(d,1e-8);
    for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        const auto& w=d.W(nx-1,j,k);
        if(!std::isfinite(w.p)||!std::isfinite(w.T)||w.p<=0.0||d.U(nx-1,j,k).rho<=0.0) return 1;
    }

    constexpr double prescribedSoundSpeed=50.0;
    const auto acousticModel=
        std::make_shared<const PrescribedAcousticModel>(prescribedSoundSpeed);
    Domain qRouted(
        nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2,acousticModel);
    for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
    {
        PrimitiveVars w(2);
        w.p=2e5;w.T=700.0;w.U={10.0,0.0,0.0};w.Y={0.8,0.2};
        qRouted.stagedW(i,j,k)=w;
    }
    if(!qRouted.primitiveToConserved()) return 1;
    for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
    {
        PrimitiveVars boundary=qRouted.W(nx-1,j,k);
        boundary.U[0]=20.0;
        if(!qRouted.setCell(nx-1,j,k,boundary)) return 1;
    }
    const double dx=mesh.xFace(nx-1)-mesh.xFace(nx-2);
    const double dt=1e-8;
    const double rho=qRouted.Q(nx-1,nr/2,nt/2).rho;
    const double outgoingWave=(20.0+prescribedSoundSpeed)
        *rho*prescribedSoundSpeed*(10.0/dx);
    const double expectedPressure=2e5-0.5*dt*outgoingWave;
    SectorNonReflectingOutlet qOutlet(mesh,2e5,0.0);
    qOutlet.apply(qRouted,dt);
    for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
        if(std::abs(qRouted.W(nx-1,j,k).p-expectedPressure)
           >1e-10*std::max(1.0,std::abs(expectedPressure))) return 1;

    std::cout<<"SECTOR_OUTLET_PASS\n"; return 0;
}
