// Sample the exact legacy tetra thermodynamic configuration against an external
// CoolProp reference. CSV input: id,p_Pa,T_K,Y_O2. No replacement/fallback states.
#include "EquilibriumThermoModel.H"
#include "CoolPropThermoModel.H"
#include "GERGResidualTable.H"
#include "GERGSaturationTable.H"
#include "AcousticStateAdmission.H"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace Foam;
int main(int argc,char**argv)
{
    if(argc!=2&&argc!=3)return 2;
    const bool coolprop=argc==3&&std::string(argv[2])=="--coolprop";
    auto model=coolprop?thermo::makeCoolPropThermoModel():
        thermo::makeLegacySinglePhaseThermoModel({"O2","C3H8"});
    GERGResidualTable residual;GERGSaturationTable saturation;
    if(!residual.load(std::string(RFAM_SOURCE_DIR)+"/tests/reference/gerg_residual_table_v3.csv")
       ||!saturation.load(std::string(RFAM_SOURCE_DIR)+"/tests/reference/gerg_saturation_table_v3.csv"))return 2;
    model->setGERGTable(&residual);model->setSaturationTable(&saturation);
    std::ifstream input(argv[1]);std::string line;std::getline(input,line);
    std::cout<<std::setprecision(17)<<"id,status,rho,h,e,a,cp,cv,region,acoustic_admitted,roundtrip_status,roundtrip_dp_rel,roundtrip_de_rel\n";
    while(std::getline(input,line))
    {
        std::stringstream s(line);std::string id,p,t,y;
        std::getline(s,id,',');std::getline(s,p,',');std::getline(s,t,',');std::getline(s,y,',');
        const auto r=model->initializePTY(std::stod(p),std::stod(t),{std::stod(y),1-std::stod(y)});
        std::cout<<id<<","<<thermo::equilibriumStatusName(r.status);
        if(!r.success()){std::cout<<",nan,nan,nan,nan,nan,nan,UNKNOWN,0,NOT_RUN,nan,nan\n";continue;}
        const auto& q=*r.state;
        thermo::EquilibriumSeed seed;seed.temperature=q.T;seed.pressure=q.p;
        const auto rt=model->recoverRhoE(q.rho,q.e,q.overallMassFractions,&seed);
        std::cout<<","<<q.rho<<","<<q.h<<","<<q.e<<","<<q.acoustic.soundSpeed()
                 <<","<<(q.thermodynamics?q.thermodynamics->Cp:NAN)<<","<<(q.thermodynamics?q.thermodynamics->Cv:NAN)
                 <<","<<thermo::equilibriumRegionName(q.region)<<","<<solverConsumerAdmitsAcousticState(q)
                 <<","<<thermo::equilibriumStatusName(rt.status)<<","<<(rt.success()?(rt.state->p-q.p)/q.p:NAN)
                 <<","<<(rt.success()?(rt.state->e-q.e)/std::max(std::abs(q.e),1.):NAN)<<"\n";
    }
}
