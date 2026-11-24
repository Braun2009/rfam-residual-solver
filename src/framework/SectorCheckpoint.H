#ifndef SectorCheckpoint_H
#define SectorCheckpoint_H

#include "Domain.H"
#include "SectorTurbulenceOperator.H"
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace Foam
{

struct SectorCheckpointMetadata
{
    double time=0.0;
    int iteration=0,threadCount=1;
    std::string mode,gitCommit;
    std::uint64_t profileHash=0;
    std::string thermoModelId;
};

struct SectorCheckpointInfo
{
    std::uint32_t version=0;
    int nx=0,ny=0,nz=0,nSpecies=0;
    SectorCheckpointMetadata metadata;
    std::vector<std::string> species;
};

class SectorCheckpoint
{
    static std::filesystem::path temporaryWritePath(
        const std::filesystem::path& target)
    {
        static std::atomic<std::uint64_t> sequence{0};
        const std::uint64_t clock = static_cast<std::uint64_t>(
            std::chrono::high_resolution_clock::now()
                .time_since_epoch().count());
        std::random_device random;
        const std::uint64_t noise =
            (static_cast<std::uint64_t>(random()) << 32U)
            ^ static_cast<std::uint64_t>(random());
        std::filesystem::path temporary = target;
        temporary += ".tmp." + std::to_string(clock)
            + "." + std::to_string(noise)
            + "." + std::to_string(sequence.fetch_add(1));
        return temporary;
    }

    static const std::string& legacyThermoModelId()
    {
        static const std::string id="unprovenanced-legacy-checkpoint";
        return id;
    }

    static void replaceAtomically(
        const std::filesystem::path& source,
        const std::filesystem::path& target)
    {
#ifdef _WIN32
        if(!MoveFileExW(
                source.c_str(),target.c_str(),
                MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH))
        {
            const std::error_code error(
                static_cast<int>(GetLastError()),std::system_category());
            throw std::runtime_error(
                "checkpoint atomic replace failed: "+error.message());
        }
#else
        std::error_code error;
        std::filesystem::rename(source,target,error);
        if(error)
            throw std::runtime_error(
                "checkpoint atomic replace failed: "+error.message());
#endif
    }

    static bool readMagic(std::istream& in)
    {
        const char expected[8]={'R','F','A','M','C','P','1','\0'};
        char actual[8]={};
        in.read(actual,sizeof(actual));
        return static_cast<bool>(in)
            &&std::memcmp(actual,expected,sizeof(actual))==0;
    }

    template<class T> static void writeValue(std::ostream& out,const T& value)
    { out.write(reinterpret_cast<const char*>(&value),sizeof(T)); }
    template<class T> static void readValue(std::istream& in,T& value)
    { in.read(reinterpret_cast<char*>(&value),sizeof(T)); }
    static void writeString(std::ostream& out,const std::string& value)
    {
        const std::uint32_t length=static_cast<std::uint32_t>(value.size());
        writeValue(out,length); out.write(value.data(),length);
    }
    static std::string readString(std::istream& in)
    {
        std::uint32_t length=0; readValue(in,length);
        if(length>1024*1024) throw std::runtime_error("invalid checkpoint string length");
        std::string value(length,'\0'); in.read(value.data(),length); return value;
    }

    static void validateMetadata(const SectorCheckpointMetadata& metadata)
    {
        if(!std::isfinite(metadata.time)||metadata.time<0.0
           ||std::signbit(metadata.time)
           ||metadata.iteration<0
           ||metadata.threadCount<=0)
            throw std::runtime_error(
                "invalid checkpoint clock or execution metadata");
    }

public:
    static SectorCheckpointInfo inspect(const std::string& path)
    {
        std::ifstream in(path,std::ios::binary);
        if(!in) throw std::runtime_error("cannot inspect checkpoint: "+path);
        if(!readMagic(in)) throw std::runtime_error("invalid checkpoint magic");
        SectorCheckpointInfo info;
        readValue(in,info.version);readValue(in,info.nx);readValue(in,info.ny);
        readValue(in,info.nz);readValue(in,info.nSpecies);
        if(info.version!=1&&info.version!=2&&info.version!=3)
            throw std::runtime_error("unsupported checkpoint version");
        readValue(in,info.metadata.time);readValue(in,info.metadata.iteration);
        readValue(in,info.metadata.threadCount);readValue(in,info.metadata.profileHash);
        info.metadata.mode=readString(in);info.metadata.gitCommit=readString(in);
        if(info.version>=3)
        {
            info.metadata.thermoModelId=readString(in);
            if(info.metadata.thermoModelId.empty())
                throw std::runtime_error("checkpoint thermo model is missing");
        }
        else info.metadata.thermoModelId=legacyThermoModelId();
        validateMetadata(info.metadata);
        std::uint32_t speciesCount=0;readValue(in,speciesCount);
        if(speciesCount!=static_cast<std::uint32_t>(info.nSpecies))
            throw std::runtime_error("checkpoint species count mismatch");
        info.species.reserve(speciesCount);
        for(std::uint32_t s=0;s<speciesCount;++s) info.species.push_back(readString(in));
        if(!in) throw std::runtime_error("checkpoint header read failed: "+path);
        return info;
    }

    static std::uint64_t hashFile(const std::string& path)
    {
        std::ifstream in(path,std::ios::binary);
        if(!in) throw std::runtime_error("cannot hash profile: "+path);
        std::uint64_t hash=1469598103934665603ULL; char buffer[4096];
        while(in.read(buffer,sizeof(buffer))||in.gcount())
            for(std::streamsize i=0;i<in.gcount();++i)
            { hash^=static_cast<unsigned char>(buffer[i]); hash*=1099511628211ULL; }
        return hash;
    }

    static void write(const std::string& path,const Domain& d,
                      const SectorCheckpointMetadata& metadata,
                      const SectorTurbulenceOperator* turbulence=nullptr)
    {
        if(!d.stateReady())
            throw std::runtime_error("cannot checkpoint an uninitialized Domain");
        if(d.stepTransactionActive())
            throw std::runtime_error("cannot checkpoint during a step transaction");
        validateMetadata(metadata);
        const std::string& thermoModelId=d.thermoModelId();
        if(thermoModelId.empty())
            throw std::runtime_error("cannot checkpoint a Domain without thermo provenance");
        if(d.speciesNames().size()!=static_cast<std::size_t>(d.nSpecies()))
            throw std::runtime_error(
                "cannot checkpoint a Domain without an explicit species layout");
        const std::vector<double>* turbulenceK=nullptr;
        const std::vector<double>* turbulenceOmega=nullptr;
        if(turbulence)
        {
            if(turbulence->stepTransactionActive())
                throw std::runtime_error(
                    "cannot checkpoint turbulence during a step transaction");
            if(!turbulence->meshMatches(d))
                throw std::runtime_error("checkpoint turbulence mesh mismatch");
            turbulenceK=&turbulence->kField();
            turbulenceOmega=&turbulence->omegaField();
            if(!turbulence->fieldsAreRestorable(
                    *turbulenceK,*turbulenceOmega))
                throw std::runtime_error("invalid checkpoint turbulence fields");
        }
        const std::filesystem::path target(path);
        const std::filesystem::path temporary = temporaryWritePath(target);
        try
        {
            std::ofstream out(
                temporary,std::ios::binary|std::ios::trunc);
            if(!out)
                throw std::runtime_error(
                    "cannot write checkpoint candidate: "+temporary.string());
            const char magic[8]={'R','F','A','M','C','P','1','\0'};
            out.write(magic,sizeof(magic));
            const std::uint32_t version=3; writeValue(out,version);
            writeValue(out,d.Nx()); writeValue(out,d.Ny()); writeValue(out,d.Nz());
            writeValue(out,d.nSpecies()); writeValue(out,metadata.time);
            writeValue(out,metadata.iteration); writeValue(out,metadata.threadCount);
            writeValue(out,metadata.profileHash); writeString(out,metadata.mode);
            writeString(out,metadata.gitCommit);
            writeString(out,thermoModelId);
            const auto& species=d.speciesNames();
            const std::uint32_t speciesCount=
                static_cast<std::uint32_t>(species.size());
            writeValue(out,speciesCount);
            for(const auto& name:species) writeString(out,name);
            for(int i=0;i<d.Nx();++i) for(int j=0;j<d.Ny();++j)
            for(int k=0;k<d.Nz();++k)
            {
                const ConservedVars& u=d.U(i,j,k);
                writeValue(out,u.rho);
                for(double value:u.rhoU) writeValue(out,value);
                writeValue(out,u.rhoE);
                for(double value:u.rhoY) writeValue(out,value);
                writeValue(out,u.rhoZ); writeValue(out,u.rhoC);
                // Preserve the real-fluid C2P branch used by the uninterrupted
                // run. Conserved state alone is not a unique temperature seed
                // near strongly nonlinear or multi-root states.
                writeValue(out,d.W(i,j,k).T);
            }
            const std::uint8_t hasTurbulence=turbulence?1:0;
            writeValue(out,hasTurbulence);
            if(turbulence)
            {
                const std::uint64_t count=turbulenceK->size();
                writeValue(out,count);
                out.write(reinterpret_cast<const char*>(turbulenceK->data()),
                    count*sizeof(double));
                out.write(reinterpret_cast<const char*>(turbulenceOmega->data()),
                    count*sizeof(double));
            }
            out.flush();
            if(!out)
                throw std::runtime_error(
                    "checkpoint candidate write failed: "+temporary.string());
            out.close();
            if(!out)
                throw std::runtime_error(
                    "checkpoint candidate close failed: "+temporary.string());

            replaceAtomically(temporary,target);
        }
        catch(...)
        {
            std::error_code cleanupError;
            std::filesystem::remove(temporary,cleanupError);
            throw;
        }
    }

    static SectorCheckpointMetadata read(const std::string& path,Domain& d,
                                         SectorTurbulenceOperator* turbulence=nullptr)
    {
        if(d.tripped())
            throw std::runtime_error("cannot restore into a tripped Domain");
        if(d.stepTransactionActive())
            throw std::runtime_error("cannot restore during a step transaction");
        if(turbulence&&turbulence->stepTransactionActive())
            throw std::runtime_error(
                "cannot restore turbulence during a step transaction");
        std::ifstream in(path,std::ios::binary);
        if(!in) throw std::runtime_error("cannot read checkpoint: "+path);
        if(!readMagic(in)) throw std::runtime_error("invalid checkpoint magic");
        std::uint32_t version=0; int nx=0,ny=0,nz=0,ns=0; readValue(in,version);
        readValue(in,nx);readValue(in,ny);readValue(in,nz);readValue(in,ns);
        if(!in) throw std::runtime_error("checkpoint header read failed: "+path);
        if(version==1||version==2)
            throw std::runtime_error(
                "checkpoint predates thermodynamic provenance and cannot be restored");
        if(version!=3)
            throw std::runtime_error("unsupported checkpoint version");
        if(nx!=d.Nx()||ny!=d.Ny()||nz!=d.Nz()||ns!=d.nSpecies())
            throw std::runtime_error("checkpoint mesh/species dimensions mismatch");
        SectorCheckpointMetadata metadata; readValue(in,metadata.time);
        readValue(in,metadata.iteration);readValue(in,metadata.threadCount);
        readValue(in,metadata.profileHash);metadata.mode=readString(in);
        metadata.gitCommit=readString(in);
        metadata.thermoModelId=readString(in);
        validateMetadata(metadata);
        if(metadata.thermoModelId.empty())
            throw std::runtime_error("checkpoint thermo model is missing");
        if(metadata.thermoModelId!=d.thermoModelId())
            throw std::runtime_error("checkpoint thermo model mismatch");
        std::uint32_t speciesCount=0;readValue(in,speciesCount);
        if(speciesCount!=d.speciesNames().size())
            throw std::runtime_error("checkpoint species count mismatch");
        for(std::uint32_t s=0;s<speciesCount;++s)
            if(readString(in)!=d.speciesNames()[s])
                throw std::runtime_error("checkpoint species order mismatch");
        const int cellCount=nx*ny*nz;
        std::vector<ConservedVars> conserved(
            static_cast<std::size_t>(cellCount),ConservedVars(ns));
        std::vector<double> temperatureSeeds(
            static_cast<std::size_t>(cellCount),300.0);
        const auto cellIndex=[&](int i,int j,int k)
        { return static_cast<std::size_t>(k*nx*ny+j*nx+i); };
        for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
        {
            ConservedVars& u=conserved[cellIndex(i,j,k)]; readValue(in,u.rho);
            for(double& value:u.rhoU) readValue(in,value); readValue(in,u.rhoE);
            for(double& value:u.rhoY) readValue(in,value);
            readValue(in,u.rhoZ);readValue(in,u.rhoC);
            double& temperature=temperatureSeeds[cellIndex(i,j,k)];
            readValue(in,temperature);
            if(!std::isfinite(temperature)||temperature<=0.0)
                throw std::runtime_error("invalid checkpoint temperature seed");
        }
        std::uint8_t hasTurbulence=0;readValue(in,hasTurbulence);
        if(hasTurbulence>1)
            throw std::runtime_error("invalid checkpoint turbulence flag");
        if((hasTurbulence!=0)!=(turbulence!=nullptr))
            throw std::runtime_error(
                "checkpoint turbulence presence does not match the restore target");
        std::vector<double> turbulenceK,turbulenceOmega;
        if(hasTurbulence)
        {
            std::uint64_t count=0;readValue(in,count);
            if(count!=static_cast<std::uint64_t>(cellCount))
                throw std::runtime_error("sector turbulence checkpoint size mismatch");
            if(!turbulence->meshMatches(d))
                throw std::runtime_error("checkpoint turbulence mesh mismatch");
            turbulenceK.resize(static_cast<std::size_t>(count));
            turbulenceOmega.resize(static_cast<std::size_t>(count));
            in.read(reinterpret_cast<char*>(turbulenceK.data()),count*sizeof(double));
            in.read(reinterpret_cast<char*>(turbulenceOmega.data()),count*sizeof(double));
        }
        if(!in) throw std::runtime_error("checkpoint read failed: "+path);
        if(hasTurbulence
           &&!turbulence->fieldsAreRestorable(turbulenceK,turbulenceOmega))
            throw std::runtime_error("invalid checkpoint turbulence fields");
        if(!d.replaceConservedField(conserved,temperatureSeeds))
        {
            const Domain::FailureInfo& failure=d.failure();
            throw std::runtime_error(
                "checkpoint state recovery failed at cell "
                +std::to_string(failure.cell)+": "+failure.detail);
        }
        if(hasTurbulence)
            turbulence->restoreFields(
                std::move(turbulenceK),std::move(turbulenceOmega));
        return metadata;
    }
};

} // End namespace Foam
#endif
