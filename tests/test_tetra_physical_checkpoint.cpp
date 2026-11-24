#include "TetraPhysicalCheckpoint.H"
#include "TetraFiniteVolume.H"
#include "EquilibriumThermoModel.H"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace Foam;

int main()
{
    const auto root=std::filesystem::temp_directory_path()/
        ("rfam-tetra-checkpoint-"+std::to_string(
            std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(root);
    try
    {
        const auto path=root/"mesh.rfamtet";
        {std::ofstream out(path,std::ios::binary);auto put=[&](const auto& x){out.write(reinterpret_cast<const char*>(&x),sizeof(x));};
        out.write("RFAMTET1",8);std::uint64_t np=5,nc=2,nf=7;put(np);put(nc);put(nf);
        double points[15]={0,0,0, 1,0,0, 0,1,0, 0,0,1, 0,0,-1};put(points);
        std::int32_t cells[8]={0,1,2,3, 0,2,1,4};put(cells);
        double vols[2]={1./6,1./6};put(vols);
        std::int32_t faces[42]={0,2,1,0,1,0, 0,1,3,0,-1,0, 1,2,3,0,-1,0, 2,0,3,0,-1,0,
            0,2,4,1,-1,0, 2,1,4,1,-1,0, 1,0,4,1,-1,0};put(faces);}
        TetraMesh mesh;mesh.read(path.string());
        std::vector<std::string> names{"O2","C3H8"};
        auto model=thermo::makeLegacySinglePhaseThermoModel(names);
        TetraFiniteVolume original(mesh,names,std::make_unique<AUSMFluxAdapter>(names));
        PrimitiveVars w(2);w.p=5e5;w.T=450.;w.Y={.6,.4};w.U={30.,-10.,50.};
        for(std::size_t i=0;i<2;++i)
        {
            w.p+=1e5;w.T+=10.;w.U[0]+=4.;
            auto q=model->initializePTY(w.p,w.T,w.Y);
            if(!q.success())throw std::runtime_error("checkpoint test initialization failed");
            original.setPrimitive(i,w,*q.state);
        }
        TetraPhysicalCheckpoint::Metadata meta;
        meta.meshSha256=std::string(64,'a');meta.modelId=model->modelId();
        meta.boundaryIdentity="test-boundaries";meta.steps=3;meta.physicalTime=2e-9;meta.lastDt=1e-9;
        const auto checkpoint=root/"state.rfamtpc";
        TetraPhysicalCheckpoint::write(checkpoint.string(),original,meta);
        bool rejected=false;
        try{TetraPhysicalCheckpoint::write(checkpoint.string(),original,meta);}
        catch(const std::runtime_error&){rejected=true;}
        if(!rejected)throw std::runtime_error("accepted checkpoint overwritten");
        TetraFiniteVolume restored(mesh,names,std::make_unique<AUSMFluxAdapter>(names));
        const auto read=TetraPhysicalCheckpoint::read(checkpoint.string(),restored,meta);
        if(read.steps!=meta.steps||read.physicalTime!=meta.physicalTime)
            throw std::runtime_error("physical clock did not roundtrip");
        for(std::size_t i=0;i<2;++i)
        {
            if(original.U(i).rho!=restored.U(i).rho
               ||original.U(i).rhoE!=restored.U(i).rhoE
               ||original.U(i).rhoU!=restored.U(i).rhoU
               ||original.U(i).rhoY!=restored.U(i).rhoY
               ||original.W(i).p!=restored.W(i).p
               ||original.W(i).T!=restored.W(i).T
               ||original.W(i).U!=restored.W(i).U
               ||original.Q(i).h!=restored.Q(i).h
               ||original.Q(i).acoustic.pressureDensityDerivative
                    !=restored.Q(i).acoustic.pressureDensityDerivative)
                throw std::runtime_error("physical checkpoint state did not roundtrip exactly");
        }
        auto wrong=meta;wrong.boundaryIdentity="different-boundaries";
        rejected=false;
        try{TetraPhysicalCheckpoint::read(checkpoint.string(),restored,wrong);}
        catch(const std::runtime_error&){rejected=true;}
        if(!rejected||original.U(0).rho!=restored.U(0).rho)
            throw std::runtime_error("checkpoint identity mismatch changed live cells");
        const auto corrupt=root/"corrupt.rfamtpc";
        std::filesystem::copy_file(checkpoint,corrupt);
        {std::fstream file(corrupt,std::ios::binary|std::ios::in|std::ios::out);
         file.seekp(-20,std::ios::end);char byte=0;file.read(&byte,1);
         file.seekp(-20,std::ios::end);byte^=1;file.write(&byte,1);}
        rejected=false;
        try{TetraPhysicalCheckpoint::read(corrupt.string(),restored,meta);}
        catch(const std::runtime_error&){rejected=true;}
        if(!rejected||original.U(0).rho!=restored.U(0).rho)
            throw std::runtime_error("corrupt checkpoint changed live cells");
        std::filesystem::remove_all(root);
        std::cout<<"TETRA_PHYSICAL_CHECKPOINT_PASS\n";return 0;
    }
    catch(const std::exception& e)
    {
        std::filesystem::remove_all(root);
        std::cerr<<e.what()<<"\n";return 1;
    }
}
