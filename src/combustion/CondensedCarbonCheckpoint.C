#include "CondensedCarbonCheckpoint.H"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
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
namespace
{

template<class T>
void writeValue(std::ostream& out,const T& value)
{
    out.write(reinterpret_cast<const char*>(&value),sizeof(T));
}

template<class T>
void readValue(std::istream& in,T& value)
{
    in.read(reinterpret_cast<char*>(&value),sizeof(T));
}

void writeString(std::ostream& out,const std::string& value)
{
    if(value.size()>1024U*1024U)
        throw std::runtime_error(
            "condensed-carbon checkpoint string is too large");
    const std::uint32_t length=static_cast<std::uint32_t>(value.size());
    writeValue(out,length);
    out.write(value.data(),length);
}

std::string readString(std::istream& in)
{
    std::uint32_t length=0;
    readValue(in,length);
    if(length>1024U*1024U)
        throw std::runtime_error(
            "invalid condensed-carbon checkpoint string length");
    std::string value(length,'\0');
    in.read(value.data(),length);
    return value;
}

std::filesystem::path temporaryWritePath(
    const std::filesystem::path& target)
{
    static std::atomic<std::uint64_t> sequence{0};
    const std::uint64_t clock=static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::random_device random;
    const std::uint64_t noise=
        (static_cast<std::uint64_t>(random())<<32U)
        ^static_cast<std::uint64_t>(random());
    std::filesystem::path temporary=target;
    temporary+=".tmp."+std::to_string(clock)+"."+std::to_string(noise)
        +"."+std::to_string(sequence.fetch_add(1));
    return temporary;
}

void replaceAtomically(
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
            "condensed-carbon checkpoint atomic replace failed: "
            +error.message());
    }
#else
    std::error_code error;
    std::filesystem::rename(source,target,error);
    if(error)
        throw std::runtime_error(
            "condensed-carbon checkpoint atomic replace failed: "
            +error.message());
#endif
}

bool readMagic(std::istream& in)
{
    const char expected[8]={'R','F','A','M','C','C','1','\0'};
    char actual[8]={};
    in.read(actual,sizeof(actual));
    return static_cast<bool>(in)
        &&std::memcmp(actual,expected,sizeof(actual))==0;
}

} // namespace

void CondensedCarbonCheckpoint::write(
    const std::string& path,
    const Domain& domain,
    const SectorCondensedCarbonOperator& carbon,
    const CondensedCarbonCheckpointMetadata& metadata)
{
    if(!domain.stateReady()||domain.stepTransactionActive())
        throw std::runtime_error(
            "cannot checkpoint condensed carbon from an inactive fluid state");
    if(carbon.stepTransactionActive())
        throw std::runtime_error(
            "cannot checkpoint condensed carbon during a step transaction");
    if(!carbon.meshMatches(domain))
        throw std::runtime_error(
            "condensed-carbon checkpoint mesh mismatch");
    if(!std::isfinite(metadata.time)||metadata.time<0.0
       ||metadata.coupledFluidCheckpointId.empty())
        throw std::invalid_argument(
            "invalid condensed-carbon checkpoint metadata");
    if(domain.thermoModelId().empty()||carbon.modelFingerprint().empty())
        throw std::runtime_error(
            "condensed-carbon checkpoint provenance is missing");
    const std::vector<CondensedCarbonState>& state=carbon.stateField();
    const std::vector<WallCarbonInventory>& wall=carbon.wallField();
    if(!carbon.fieldsAreRestorable(
            state,wall,carbon.sootBoundaryMassNetIntoDomain(),
            carbon.sootBoundaryNumberNetIntoDomain()))
        throw std::runtime_error(
            "invalid condensed-carbon checkpoint fields");

    const std::filesystem::path target(path);
    const std::filesystem::path temporary=temporaryWritePath(target);
    try
    {
        std::ofstream out(temporary,std::ios::binary|std::ios::trunc);
        if(!out)
            throw std::runtime_error(
                "cannot write condensed-carbon checkpoint candidate");
        const char magic[8]={'R','F','A','M','C','C','1','\0'};
        out.write(magic,sizeof(magic));
        const std::uint32_t version=1;
        writeValue(out,version);
        writeValue(out,domain.Nx());
        writeValue(out,domain.Ny());
        writeValue(out,domain.Nz());
        const std::uint8_t mode=carbon.numberDensityEnabled()?2U:1U;
        writeValue(out,mode);
        writeString(out,carbon.meshFingerprint());
        writeValue(out,metadata.time);
        writeValue(out,metadata.iteration);
        writeString(out,metadata.coupledFluidCheckpointId);
        writeString(out,domain.thermoModelId());
        writeString(out,carbon.modelFingerprint());
        writeString(out,carbon.kineticsModelId());
        writeString(out,carbon.wallModelId());
        const std::uint64_t stateCount=state.size();
        writeValue(out,stateCount);
        for(const CondensedCarbonState& value:state)
        {
            writeValue(out,value.rhoSoot);
            writeValue(out,value.numberDensity);
        }
        const std::uint64_t wallCount=wall.size();
        writeValue(out,wallCount);
        for(const WallCarbonInventory& value:wall)
            writeValue(out,value.massPerArea);
        writeValue(out,carbon.sootBoundaryMassNetIntoDomain());
        writeValue(out,carbon.sootBoundaryNumberNetIntoDomain());
        out.flush();
        if(!out)
            throw std::runtime_error(
                "condensed-carbon checkpoint candidate write failed");
        out.close();
        if(!out)
            throw std::runtime_error(
                "condensed-carbon checkpoint candidate close failed");
        replaceAtomically(temporary,target);
    }
    catch(...)
    {
        std::error_code cleanupError;
        std::filesystem::remove(temporary,cleanupError);
        throw;
    }
}

CondensedCarbonCheckpointMetadata CondensedCarbonCheckpoint::read(
    const std::string& path,
    const Domain& domain,
    SectorCondensedCarbonOperator& carbon,
    const std::string& expectedFluidCheckpointId)
{
    if(!domain.stateReady()||domain.tripped()||domain.stepTransactionActive())
        throw std::runtime_error(
            "cannot restore condensed carbon into an inactive fluid state");
    if(carbon.stepTransactionActive())
        throw std::runtime_error(
            "cannot restore condensed carbon during a step transaction");
    if(expectedFluidCheckpointId.empty())
        throw std::invalid_argument(
            "expected fluid checkpoint id is empty");
    std::ifstream in(path,std::ios::binary);
    if(!in)
        throw std::runtime_error(
            "cannot read condensed-carbon checkpoint: "+path);
    if(!readMagic(in))
        throw std::runtime_error(
            "invalid condensed-carbon checkpoint magic");
    std::uint32_t version=0;
    int nx=0,ny=0,nz=0;
    std::uint8_t mode=0;
    readValue(in,version);
    readValue(in,nx);
    readValue(in,ny);
    readValue(in,nz);
    readValue(in,mode);
    if(!in||version!=1)
        throw std::runtime_error(
            "unsupported condensed-carbon checkpoint version");
    if(nx!=domain.Nx()||ny!=domain.Ny()||nz!=domain.Nz()
       ||!carbon.meshMatches(domain))
        throw std::runtime_error(
            "condensed-carbon checkpoint mesh mismatch");
    const std::uint8_t expectedMode=carbon.numberDensityEnabled()?2U:1U;
    if(mode!=expectedMode)
        throw std::runtime_error(
            "condensed-carbon checkpoint moment-model mismatch");
    const std::string meshFingerprint=readString(in);
    if(meshFingerprint!=carbon.meshFingerprint())
        throw std::runtime_error(
            "condensed-carbon checkpoint mesh-geometry mismatch");

    CondensedCarbonCheckpointMetadata metadata;
    readValue(in,metadata.time);
    readValue(in,metadata.iteration);
    metadata.coupledFluidCheckpointId=readString(in);
    const std::string thermoModelId=readString(in);
    const std::string modelFingerprint=readString(in);
    const std::string kineticsModelId=readString(in);
    const std::string wallModelId=readString(in);
    if(!std::isfinite(metadata.time)||metadata.time<0.0
       ||metadata.coupledFluidCheckpointId!=expectedFluidCheckpointId)
        throw std::runtime_error(
            "condensed-carbon checkpoint coupling metadata mismatch");
    if(thermoModelId!=domain.thermoModelId())
        throw std::runtime_error(
            "condensed-carbon checkpoint thermo model mismatch");
    if(modelFingerprint!=carbon.modelFingerprint()
       ||kineticsModelId!=carbon.kineticsModelId()
       ||wallModelId!=carbon.wallModelId())
        throw std::runtime_error(
            "condensed-carbon checkpoint model fingerprint mismatch");

    std::uint64_t stateCount=0;
    readValue(in,stateCount);
    const std::uint64_t expectedState=static_cast<std::uint64_t>(domain.N());
    if(stateCount!=expectedState)
        throw std::runtime_error(
            "condensed-carbon checkpoint state size mismatch");
    std::vector<CondensedCarbonState> state(
        static_cast<std::size_t>(stateCount));
    for(CondensedCarbonState& value:state)
    {
        readValue(in,value.rhoSoot);
        readValue(in,value.numberDensity);
    }
    std::uint64_t wallCount=0;
    readValue(in,wallCount);
    const std::uint64_t expectedWall=
        static_cast<std::uint64_t>(domain.Nx())*domain.Nz();
    if(wallCount!=expectedWall)
        throw std::runtime_error(
            "condensed-carbon checkpoint wall size mismatch");
    std::vector<WallCarbonInventory> wall(
        static_cast<std::size_t>(wallCount));
    for(WallCarbonInventory& value:wall)
        readValue(in,value.massPerArea);
    double boundaryMass=0.0,boundaryNumber=0.0;
    readValue(in,boundaryMass);
    readValue(in,boundaryNumber);
    if(!in||in.peek()!=std::char_traits<char>::eof())
        throw std::runtime_error(
            "condensed-carbon checkpoint payload is invalid");
    if(!carbon.fieldsAreRestorable(
            state,wall,boundaryMass,boundaryNumber))
        throw std::runtime_error(
            "invalid condensed-carbon checkpoint fields");
    carbon.restoreFields(
        std::move(state),std::move(wall),boundaryMass,boundaryNumber);
    return metadata;
}

} // End namespace Foam
