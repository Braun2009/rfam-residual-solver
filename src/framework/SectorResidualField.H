#ifndef SectorResidualField_H
#define SectorResidualField_H

/*---------------------------------------------------------------------------*\
    SectorResidualField - provenance-bound steady residual/FAS source field.
\*---------------------------------------------------------------------------*/

#include "Domain.H"
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
#include <vector>

namespace Foam
{

class SectorResidualField
{
public:
    enum class Kind : std::uint32_t
    {
        STEADY_RESIDUAL=1,
        FAS_SOURCE=2
    };

    struct Metadata
    {
        Kind kind=Kind::STEADY_RESIDUAL;
        int nx=0,ny=0,nz=0,nSpecies=0;
        double referenceStep=0.0;
        double time=0.0;
        int iteration=0;
        std::uint64_t profileHash=0;
        std::uint64_t stateCheckpointHash=0;
        std::uint64_t fineStateCheckpointHash=0;
        std::string mode;
        std::string gitCommit;
        std::string thermoModelId;
        std::vector<std::string> species;
    };

    struct Data
    {
        Metadata metadata;
        std::vector<unsigned char> active;
        std::vector<ConservedVars> field;
    };

private:
    template<class T>
    static void writeValue(std::ostream& output,const T& value)
    {
        output.write(reinterpret_cast<const char*>(&value),sizeof(T));
    }

    template<class T>
    static void readValue(std::istream& input,T& value)
    {
        input.read(reinterpret_cast<char*>(&value),sizeof(T));
    }

    static void writeString(std::ostream& output,const std::string& value)
    {
        const std::uint32_t length=static_cast<std::uint32_t>(value.size());
        writeValue(output,length);
        output.write(value.data(),length);
    }

    static std::string readString(std::istream& input)
    {
        std::uint32_t length=0;
        readValue(input,length);
        if(length>1024U*1024U)
            throw std::runtime_error("invalid residual-field string length");
        std::string value(length,'\0');
        input.read(value.data(),length);
        return value;
    }

    static std::filesystem::path temporaryPath(
        const std::filesystem::path& target)
    {
        static std::atomic<std::uint64_t> sequence{0};
        const std::uint64_t clock=static_cast<std::uint64_t>(
            std::chrono::high_resolution_clock::now()
                .time_since_epoch().count());
        std::random_device random;
        std::filesystem::path temporary=target;
        temporary+=".tmp."+std::to_string(clock)+"."
            +std::to_string(random())+"."
            +std::to_string(sequence.fetch_add(1));
        return temporary;
    }

    static void validate(const Data& data)
    {
        const Metadata& metadata=data.metadata;
        if(metadata.kind!=Kind::STEADY_RESIDUAL
            &&metadata.kind!=Kind::FAS_SOURCE)
            throw std::runtime_error("invalid residual-field kind");
        if(metadata.nx<=0||metadata.ny<=0||metadata.nz<=0
            ||metadata.nSpecies<=0
            ||!std::isfinite(metadata.referenceStep)
            ||metadata.referenceStep<=0.0
            ||!std::isfinite(metadata.time)||metadata.time<0.0
            ||std::signbit(metadata.time)||metadata.iteration<0
            ||metadata.profileHash==0
            ||metadata.stateCheckpointHash==0
            ||metadata.mode.empty()||metadata.gitCommit.empty()
            ||metadata.thermoModelId.empty()
            ||metadata.species.size()
                !=static_cast<std::size_t>(metadata.nSpecies))
            throw std::runtime_error("invalid residual-field metadata");
        if(metadata.kind==Kind::FAS_SOURCE
            &&metadata.fineStateCheckpointHash==0)
            throw std::runtime_error(
                "FAS source is missing the fine-state checkpoint hash");
        const std::size_t count=static_cast<std::size_t>(metadata.nx)
            *static_cast<std::size_t>(metadata.ny)
            *static_cast<std::size_t>(metadata.nz);
        if(data.active.size()!=count||data.field.size()!=count)
            throw std::runtime_error("residual-field dimensions do not match");
        for(std::size_t cell=0;cell<count;++cell)
        {
            if(data.active[cell]>1U
                ||data.field[cell].rhoY.size()
                    !=static_cast<std::size_t>(metadata.nSpecies))
                throw std::runtime_error(
                    "invalid residual-field cell layout");
            const ConservedVars& value=data.field[cell];
            auto finite=[](double scalar){ return std::isfinite(scalar); };
            if(!finite(value.rho)||!finite(value.rhoE)
                ||!finite(value.rhoZ)||!finite(value.rhoC))
                throw std::runtime_error("nonfinite residual-field value");
            for(double scalar:value.rhoU)
                if(!finite(scalar))
                    throw std::runtime_error("nonfinite residual momentum");
            for(double scalar:value.rhoY)
                if(!finite(scalar))
                    throw std::runtime_error("nonfinite residual species");
        }
    }

public:
    static void write(const std::string& path,const Data& data)
    {
        validate(data);
        const std::filesystem::path target(path);
        const std::filesystem::path temporary=temporaryPath(target);
        try
        {
            std::ofstream output(
                temporary,std::ios::binary|std::ios::trunc);
            if(!output)
                throw std::runtime_error(
                    "cannot write residual field: "+temporary.string());
            const char magic[8]={'R','F','A','M','S','R','1','\0'};
            output.write(magic,sizeof(magic));
            const std::uint32_t version=1;
            writeValue(output,version);
            const Metadata& metadata=data.metadata;
            writeValue(output,metadata.kind);
            writeValue(output,metadata.nx);
            writeValue(output,metadata.ny);
            writeValue(output,metadata.nz);
            writeValue(output,metadata.nSpecies);
            writeValue(output,metadata.referenceStep);
            writeValue(output,metadata.time);
            writeValue(output,metadata.iteration);
            writeValue(output,metadata.profileHash);
            writeValue(output,metadata.stateCheckpointHash);
            writeValue(output,metadata.fineStateCheckpointHash);
            writeString(output,metadata.mode);
            writeString(output,metadata.gitCommit);
            writeString(output,metadata.thermoModelId);
            const std::uint32_t speciesCount=
                static_cast<std::uint32_t>(metadata.species.size());
            writeValue(output,speciesCount);
            for(const std::string& species:metadata.species)
                writeString(output,species);
            const std::uint64_t count=data.field.size();
            writeValue(output,count);
            output.write(
                reinterpret_cast<const char*>(data.active.data()),
                static_cast<std::streamsize>(data.active.size()));
            for(const ConservedVars& value:data.field)
            {
                writeValue(output,value.rho);
                for(double scalar:value.rhoU) writeValue(output,scalar);
                writeValue(output,value.rhoE);
                for(double scalar:value.rhoY) writeValue(output,scalar);
                writeValue(output,value.rhoZ);
                writeValue(output,value.rhoC);
            }
            output.flush();
            if(!output)
                throw std::runtime_error("residual-field write failed");
            output.close();
            std::error_code error;
            std::filesystem::rename(temporary,target,error);
            if(error)
            {
                std::filesystem::remove(target,error);
                error.clear();
                std::filesystem::rename(temporary,target,error);
            }
            if(error)
                throw std::runtime_error(
                    "residual-field atomic replace failed: "+error.message());
        }
        catch(...)
        {
            std::error_code error;
            std::filesystem::remove(temporary,error);
            throw;
        }
    }

    static Data read(const std::string& path)
    {
        std::ifstream input(path,std::ios::binary);
        if(!input)
            throw std::runtime_error("cannot read residual field: "+path);
        const char expected[8]={'R','F','A','M','S','R','1','\0'};
        char actual[8]={};
        input.read(actual,sizeof(actual));
        if(!input||std::memcmp(actual,expected,sizeof(expected))!=0)
            throw std::runtime_error("invalid residual-field magic");
        std::uint32_t version=0;
        readValue(input,version);
        if(version!=1)
            throw std::runtime_error("unsupported residual-field version");
        Data data;
        Metadata& metadata=data.metadata;
        readValue(input,metadata.kind);
        readValue(input,metadata.nx);
        readValue(input,metadata.ny);
        readValue(input,metadata.nz);
        readValue(input,metadata.nSpecies);
        readValue(input,metadata.referenceStep);
        readValue(input,metadata.time);
        readValue(input,metadata.iteration);
        readValue(input,metadata.profileHash);
        readValue(input,metadata.stateCheckpointHash);
        readValue(input,metadata.fineStateCheckpointHash);
        metadata.mode=readString(input);
        metadata.gitCommit=readString(input);
        metadata.thermoModelId=readString(input);
        std::uint32_t speciesCount=0;
        readValue(input,speciesCount);
        metadata.species.reserve(speciesCount);
        for(std::uint32_t species=0;species<speciesCount;++species)
            metadata.species.push_back(readString(input));
        std::uint64_t count=0;
        readValue(input,count);
        const std::uint64_t expectedCount=
            static_cast<std::uint64_t>(metadata.nx)
            *static_cast<std::uint64_t>(metadata.ny)
            *static_cast<std::uint64_t>(metadata.nz);
        if(count!=expectedCount||count>100000000ULL)
            throw std::runtime_error("invalid residual-field cell count");
        data.active.resize(static_cast<std::size_t>(count));
        input.read(reinterpret_cast<char*>(data.active.data()),
            static_cast<std::streamsize>(data.active.size()));
        data.field.assign(
            static_cast<std::size_t>(count),
            ConservedVars(metadata.nSpecies));
        for(ConservedVars& value:data.field)
        {
            readValue(input,value.rho);
            for(double& scalar:value.rhoU) readValue(input,scalar);
            readValue(input,value.rhoE);
            for(double& scalar:value.rhoY) readValue(input,scalar);
            readValue(input,value.rhoZ);
            readValue(input,value.rhoC);
        }
        if(!input)
            throw std::runtime_error("residual-field read failed: "+path);
        validate(data);
        return data;
    }

    static void requireCompatible(
        const Data& data,const Domain& domain,
        std::uint64_t profileHash,const std::string& mode,
        std::uint64_t stateCheckpointHash,Kind expectedKind)
    {
        const Metadata& metadata=data.metadata;
        if(metadata.kind!=expectedKind
            ||metadata.nx!=domain.Nx()||metadata.ny!=domain.Ny()
            ||metadata.nz!=domain.Nz()
            ||metadata.nSpecies!=domain.nSpecies()
            ||metadata.species!=domain.speciesNames()
            ||metadata.thermoModelId!=domain.thermoModelId()
            ||metadata.profileHash!=profileHash
            ||metadata.mode!=mode
            ||metadata.stateCheckpointHash!=stateCheckpointHash)
            throw std::runtime_error(
                "residual field is incompatible with the authoritative state");
    }
};

} // End namespace Foam

#endif
