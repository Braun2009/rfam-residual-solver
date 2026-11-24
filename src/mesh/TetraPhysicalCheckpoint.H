#ifndef TetraPhysicalCheckpoint_H
#define TetraPhysicalCheckpoint_H

#include "TetraFiniteVolume.H"
#include "TetraAcousticAdmission.H"
#include "CompleteCellState.H"
#include "StateAcceptance.H"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace Foam
{
// Accepted-state only. Unlike RFAMTETF1, every conserved, primitive and EOS
// value is stored as a double; no flash or PT rounding occurs on restore.
class TetraPhysicalCheckpoint
{
public:
    struct Metadata
    {
        std::string meshSha256,modelId,boundaryIdentity;
        std::uint64_t steps=0;
        double physicalTime=0.,lastDt=0.;
    };
private:
    static constexpr std::uint64_t offset=14695981039346656037ULL;
    static constexpr std::uint64_t prime=1099511628211ULL;
    template<typename T> static void value(std::ostream& out,std::uint64_t& hash,const T& v)
    {
        static_assert(std::is_trivially_copyable<T>::value,"checkpoint scalar must be trivially copyable");
        const auto* bytes=reinterpret_cast<const unsigned char*>(&v);
        out.write(reinterpret_cast<const char*>(bytes),sizeof(T));
        if(!out)throw std::runtime_error("tetra checkpoint write failed");
        for(std::size_t i=0;i<sizeof(T);++i){hash^=bytes[i];hash*=prime;}
    }
    template<typename T> static T value(std::istream& in,std::uint64_t& hash)
    {
        static_assert(std::is_trivially_copyable<T>::value,"checkpoint scalar must be trivially copyable");
        T v{};auto* bytes=reinterpret_cast<unsigned char*>(&v);
        in.read(reinterpret_cast<char*>(bytes),sizeof(T));
        if(!in)throw std::runtime_error("truncated tetra checkpoint");
        for(std::size_t i=0;i<sizeof(T);++i){hash^=bytes[i];hash*=prime;}
        return v;
    }
    static void text(std::ostream& out,std::uint64_t& hash,const std::string& s)
    {
        if(s.empty()||s.size()>4096)throw std::runtime_error("invalid tetra checkpoint identity");
        value(out,hash,static_cast<std::uint32_t>(s.size()));
        for(unsigned char c:s)value(out,hash,c);
    }
    static std::string text(std::istream& in,std::uint64_t& hash)
    {
        const auto length=value<std::uint32_t>(in,hash);
        if(length==0||length>4096)throw std::runtime_error("invalid tetra checkpoint identity length");
        std::string s(length,'\0');
        for(char& c:s)c=static_cast<char>(value<unsigned char>(in,hash));
        return s;
    }
    static void vector(std::ostream& out,std::uint64_t& hash,const std::vector<double>& a,std::size_t limit)
    {
        if(a.size()>limit)throw std::runtime_error("invalid tetra checkpoint vector size");
        value(out,hash,static_cast<std::uint32_t>(a.size()));
        for(double v:a)value(out,hash,v);
    }
    static std::vector<double> vector(std::istream& in,std::uint64_t& hash,std::size_t limit)
    {
        const auto n=value<std::uint32_t>(in,hash);
        if(n>limit)throw std::runtime_error("invalid tetra checkpoint vector size");
        std::vector<double> a(n);
        for(double& v:a)v=value<double>(in,hash);
        return a;
    }
    static void phase(std::ostream& out,std::uint64_t& hash,
                      const std::optional<thermo::EquilibriumPhaseProperties>& p,std::size_t species)
    {
        value(out,hash,static_cast<std::uint8_t>(p.has_value()));
        if(!p)return;
        vector(out,hash,p->massFractions,species);
        vector(out,hash,p->moleFractions,species);
        for(double v:{p->molecularWeight,p->rho,p->h,p->e,p->s,p->compressibility})value(out,hash,v);
        value(out,hash,static_cast<std::uint8_t>(p->entropyAvailable));
    }
    static std::optional<thermo::EquilibriumPhaseProperties> phase(
        std::istream& in,std::uint64_t& hash,std::size_t species)
    {
        const auto present=value<std::uint8_t>(in,hash);
        if(present>1)throw std::runtime_error("invalid tetra checkpoint phase tag");
        if(!present)return std::nullopt;
        thermo::EquilibriumPhaseProperties p;
        p.massFractions=vector(in,hash,species);p.moleFractions=vector(in,hash,species);
        p.molecularWeight=value<double>(in,hash);p.rho=value<double>(in,hash);
        p.h=value<double>(in,hash);p.e=value<double>(in,hash);
        p.s=value<double>(in,hash);p.compressibility=value<double>(in,hash);
        const auto available=value<std::uint8_t>(in,hash);
        if(available>1)throw std::runtime_error("invalid tetra checkpoint entropy tag");
        p.entropyAvailable=available!=0;return p;
    }
    static void state(std::ostream& out,std::uint64_t& hash,const CompleteCellState& c,std::size_t species)
    {
        const auto& u=c.U;const auto& w=c.W;const auto& q=c.Q;
        value(out,hash,u.rho);vector(out,hash,u.rhoU,3);value(out,hash,u.rhoE);
        vector(out,hash,u.rhoY,species);value(out,hash,u.rhoZ);value(out,hash,u.rhoC);
        value(out,hash,w.p);value(out,hash,w.T);vector(out,hash,w.U,3);
        vector(out,hash,w.Y,species);value(out,hash,w.Z);value(out,hash,w.C);
        value(out,hash,static_cast<std::int32_t>(q.region));vector(out,hash,q.overallMassFractions,species);
        for(double v:{q.p,q.T,q.rho,q.h,q.e,q.s,q.vaporMoleFraction,q.vaporMassFraction,
                      q.vaporVolumeFraction,q.condensedMassFraction,q.densityProjectionDelta,
                      q.energyProjectionDelta})value(out,hash,v);
        value(out,hash,static_cast<std::uint8_t>(q.entropyAvailable));
        value(out,hash,static_cast<std::uint8_t>(q.conservedStateProjected));
        value(out,hash,static_cast<std::int32_t>(q.acoustic.closure));
        value(out,hash,static_cast<std::int32_t>(q.acoustic.evidence));
        value(out,hash,q.acoustic.pressureDensityDerivative);
        value(out,hash,static_cast<std::uint8_t>(q.thermodynamics.has_value()));
        if(q.thermodynamics)
        {
            const auto& p=*q.thermodynamics;
            value(out,hash,p.version);value(out,hash,static_cast<std::int32_t>(p.closure));
            for(double v:{p.Cp,p.Cv,p.dpdrho_T,p.dpdT_rho,p.dedrho_T,p.dedT_rho})value(out,hash,v);
        }
        phase(out,hash,q.liquid,species);phase(out,hash,q.vapor,species);
        phase(out,hash,q.singlePhase,species);phase(out,hash,q.condensed,species);
    }
    static CompleteCellState state(std::istream& in,std::uint64_t& hash,std::size_t species)
    {
        CompleteCellState c(static_cast<int>(species));auto& u=c.U;auto& w=c.W;auto& q=c.Q;
        u.rho=value<double>(in,hash);u.rhoU=vector(in,hash,3);u.rhoE=value<double>(in,hash);
        u.rhoY=vector(in,hash,species);u.rhoZ=value<double>(in,hash);u.rhoC=value<double>(in,hash);
        w.p=value<double>(in,hash);w.T=value<double>(in,hash);w.U=vector(in,hash,3);
        w.Y=vector(in,hash,species);w.Z=value<double>(in,hash);w.C=value<double>(in,hash);
        q.region=static_cast<thermo::EquilibriumRegion>(value<std::int32_t>(in,hash));
        q.overallMassFractions=vector(in,hash,species);
        q.p=value<double>(in,hash);q.T=value<double>(in,hash);q.rho=value<double>(in,hash);
        q.h=value<double>(in,hash);q.e=value<double>(in,hash);q.s=value<double>(in,hash);
        q.vaporMoleFraction=value<double>(in,hash);
        q.vaporMassFraction=value<double>(in,hash);
        q.vaporVolumeFraction=value<double>(in,hash);
        q.condensedMassFraction=value<double>(in,hash);
        q.densityProjectionDelta=value<double>(in,hash);
        q.energyProjectionDelta=value<double>(in,hash);
        const auto entropy=value<std::uint8_t>(in,hash);
        const auto projected=value<std::uint8_t>(in,hash);
        if(entropy>1||projected>1)throw std::runtime_error("invalid tetra checkpoint state flags");
        q.entropyAvailable=entropy!=0;q.conservedStateProjected=projected!=0;
        q.acoustic.closure=static_cast<thermo::AcousticClosure>(value<std::int32_t>(in,hash));
        q.acoustic.evidence=static_cast<thermo::AcousticEvidence>(value<std::int32_t>(in,hash));
        q.acoustic.pressureDensityDerivative=value<double>(in,hash);
        const auto thermo=value<std::uint8_t>(in,hash);
        if(thermo>1)throw std::runtime_error("invalid tetra checkpoint thermo tag");
        if(thermo)
        {
            thermo::EquilibriumThermodynamicPayload p;
            p.version=value<std::uint32_t>(in,hash);
            p.closure=static_cast<thermo::AcousticClosure>(value<std::int32_t>(in,hash));
            p.Cp=value<double>(in,hash);p.Cv=value<double>(in,hash);
            p.dpdrho_T=value<double>(in,hash);p.dpdT_rho=value<double>(in,hash);
            p.dedrho_T=value<double>(in,hash);p.dedT_rho=value<double>(in,hash);
            q.thermodynamics=p;
        }
        q.liquid=phase(in,hash,species);q.vapor=phase(in,hash,species);
        q.singlePhase=phase(in,hash,species);q.condensed=phase(in,hash,species);
        return c;
    }
    static void admitted(const CompleteCellState& c,std::size_t species)
    {
        if(!c.valid()||!tetraAdmitsAcoustics(c.Q)||c.U.rhoY.size()!=species
           ||c.W.Y.size()!=species||c.Q.overallMassFractions.size()!=species
           ||c.U.rhoU.size()!=3||c.W.U.size()!=3||c.Q.conservedStateProjected)
            throw std::runtime_error("tetra checkpoint cell contract failed");
        const double rho=c.U.rho;
        double speed2=0.;for(double m:c.U.rhoU)speed2+=(m/rho)*(m/rho);
        const double e=c.U.rhoE/rho-.5*speed2;
        if(std::abs(c.Q.rho-rho)>stateAcceptanceRel()*std::max(rho,1.)
           ||std::abs(c.Q.e-e)>stateAcceptanceRel()*std::max({std::abs(e),std::abs(c.Q.e),1.}))
            throw std::runtime_error("tetra checkpoint conserved/EOS closure failed");
    }
    static void metadataValid(const Metadata& m)
    {
        if(m.meshSha256.size()!=64||m.modelId.empty()||m.boundaryIdentity.empty()
           ||!std::isfinite(m.physicalTime)||m.physicalTime<0.
           ||!std::isfinite(m.lastDt)||m.lastDt<0.
           ||(m.steps==0&&m.physicalTime!=0.)||(m.steps>0&&(!(m.physicalTime>0.)||!(m.lastDt>0.))))
            throw std::runtime_error("invalid tetra checkpoint metadata");
    }
public:
    static void write(const std::string& path,const TetraFiniteVolume& fv,const Metadata& metadata)
    {
        metadataValid(metadata);
        if(std::filesystem::exists(path))
            throw std::runtime_error("refusing to overwrite accepted tetra checkpoint");
        std::filesystem::path temp=path+".tmp";
        if(std::filesystem::exists(temp))throw std::runtime_error("tetra checkpoint temp file exists");
        try
        {
            std::ofstream out(temp,std::ios::binary);
            if(!out)throw std::runtime_error("cannot open tetra checkpoint temp file");
            std::uint64_t hash=offset;
            for(char c:std::string("RFAMTPC1",8))value(out,hash,c);
            value(out,hash,std::uint32_t(1));
            value(out,hash,static_cast<std::uint64_t>(fv.cellCount()));
            value(out,hash,static_cast<std::uint32_t>(fv.species().size()));
            text(out,hash,metadata.meshSha256);text(out,hash,metadata.modelId);
            text(out,hash,metadata.boundaryIdentity);
            value(out,hash,metadata.steps);value(out,hash,metadata.physicalTime);value(out,hash,metadata.lastDt);
            for(std::size_t i=0;i<fv.cellCount();++i)
            {
                CompleteCellState c;c.U=fv.U(i);c.W=fv.W(i);c.Q=fv.Q(i);
                admitted(c,fv.species().size());state(out,hash,c,fv.species().size());
            }
            out.write(reinterpret_cast<const char*>(&hash),sizeof(hash));
            out.flush();if(!out)throw std::runtime_error("tetra checkpoint flush failed");
            out.close();if(!out)throw std::runtime_error("tetra checkpoint close failed");
            std::filesystem::rename(temp,path);
        }
        catch(...){std::filesystem::remove(temp);throw;}
    }
    static Metadata read(const std::string& path,TetraFiniteVolume& fv,const Metadata& expected)
    {
        metadataValid(expected);
        std::ifstream in(path,std::ios::binary|std::ios::ate);
        if(!in)throw std::runtime_error("cannot open tetra checkpoint");
        const auto size=in.tellg();in.seekg(0);
        if(size<64)throw std::runtime_error("truncated tetra checkpoint");
        std::uint64_t hash=offset;
        for(char c:std::string("RFAMTPC1",8))
            if(value<char>(in,hash)!=c)throw std::runtime_error("invalid tetra checkpoint magic");
        if(value<std::uint32_t>(in,hash)!=1)throw std::runtime_error("unsupported tetra checkpoint version");
        if(value<std::uint64_t>(in,hash)!=fv.cellCount()
           ||value<std::uint32_t>(in,hash)!=fv.species().size())
            throw std::runtime_error("tetra checkpoint mesh/species size mismatch");
        Metadata m;m.meshSha256=text(in,hash);m.modelId=text(in,hash);
        m.boundaryIdentity=text(in,hash);
        m.steps=value<std::uint64_t>(in,hash);
        m.physicalTime=value<double>(in,hash);m.lastDt=value<double>(in,hash);
        metadataValid(m);
        if(m.meshSha256!=expected.meshSha256||m.modelId!=expected.modelId
           ||m.boundaryIdentity!=expected.boundaryIdentity)
            throw std::runtime_error("tetra checkpoint identity mismatch");
        // Stage every cell and publish only after checksum, EOF and all state
        // contracts pass. A rejected restart never changes a live cell.
        std::vector<CompleteCellState> staged;
        staged.reserve(fv.cellCount());
        for(std::size_t i=0;i<fv.cellCount();++i)
        {
            auto c=state(in,hash,fv.species().size());
            admitted(c,fv.species().size());staged.push_back(std::move(c));
        }
        std::uint64_t saved=0;in.read(reinterpret_cast<char*>(&saved),sizeof(saved));
        if(!in||saved!=hash||in.peek()!=std::char_traits<char>::eof())
            throw std::runtime_error("tetra checkpoint checksum or length mismatch");
        for(std::size_t i=0;i<fv.cellCount();++i)
            fv.setState(i,staged[i].W,staged[i].U,staged[i].Q);
        return m;
    }
};
} // namespace Foam
#endif
