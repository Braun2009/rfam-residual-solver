#include "TetraFiniteVolume.H"
#include "CompleteCellState.H"
#include "TetraAcousticAdmission.H"
#include "TetraFaceInvariant.H"
#include "TetraMolecularSpectral.H"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>

namespace Foam
{
namespace
{
void addScaled(ConservedVars& a,const ConservedVars& b,double scale)
{
    a.rho+=scale*b.rho;a.rhoE+=scale*b.rhoE;
    a.rhoZ+=scale*b.rhoZ;a.rhoC+=scale*b.rhoC;
    for(int i=0;i<3;++i)a.rhoU[i]+=scale*b.rhoU[i];
    for(std::size_t s=0;s<a.rhoY.size();++s)a.rhoY[s]+=scale*b.rhoY[s];
}
void admit(const PrimitiveVars& w,const thermo::EquilibriumState& q)
{
    if(!tetraAdmitsAcoustics(q))
        throw std::runtime_error("tetra flux rejected invalid, projected or unadmitted acoustic state");
    if(!std::isfinite(w.p)||!std::isfinite(w.T)||w.U.size()!=3
       ||std::abs(w.p-q.p)>1e-8*std::max(w.p,1.)
       ||std::abs(w.T-q.T)>1e-8*std::max(w.T,1.)
       ||w.Y.size()!=q.overallMassFractions.size())
        throw std::runtime_error("tetra flux primitive/thermo mismatch");
    for(std::size_t s=0;s<w.Y.size();++s)
        if(!std::isfinite(w.Y[s])||std::abs(w.Y[s]-q.overallMassFractions[s])>1e-10)
            throw std::runtime_error("tetra flux species mismatch");
    for(double value:w.U)if(!std::isfinite(value))throw std::runtime_error("nonfinite tetra velocity");
    if(!std::isfinite(w.Z)||!std::isfinite(w.C))throw std::runtime_error("nonfinite tetra scalar");
}
}
TetraFiniteVolume::TetraFiniteVolume(const TetraMesh& mesh,std::vector<std::string> species,std::unique_ptr<FluxScheme> flux)
    :mesh_(mesh),flux_(std::move(flux)),species_(std::move(species))
{
    if(!flux_||species_.empty())throw std::invalid_argument("tetra kernel requires flux and species");
    W_.assign(mesh.cellCount(),PrimitiveVars(species_.size()));
    U_.assign(mesh.cellCount(),ConservedVars(species_.size()));Q_.resize(mesh.cellCount());
    areas_.resize(mesh.faceCount());normals_.resize(mesh.faceCount());
    for(std::size_t f=0;f<mesh.faceCount();++f)
    {
        const auto a=mesh.faceAreaVector(f);double mag=std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        if(!std::isfinite(mag)||mag<=0)throw std::runtime_error("invalid tetra face area");
        areas_[f]=mag;normals_[f]={a.x/mag,a.y/mag,a.z/mag};
    }
    faceFlux_.assign(mesh.faceCount(),ConservedVars(species_.size()));
}
void TetraFiniteVolume::setState(std::size_t c,const PrimitiveVars& w,const ConservedVars& u,const thermo::EquilibriumState& q)
{ W_.at(c)=w;U_.at(c)=u;Q_.at(c)=q; }
void TetraFiniteVolume::exchangeState(std::size_t c,PrimitiveVars& w,ConservedVars& u,thermo::EquilibriumState& q)
{ std::swap(W_.at(c),w);std::swap(U_.at(c),u);std::swap(Q_.at(c),q); }
void TetraFiniteVolume::setPrimitive(std::size_t c,const PrimitiveVars& w,const thermo::EquilibriumState& q)
{
    auto state=CompleteStateConversion::fromPrimitiveAndEquilibrium(w,q);
    if(!state.success())throw std::runtime_error("invalid tetra initialization W/Q");
    setState(c,state.state->W,state.state->U,state.state->Q);
}
void TetraFiniteVolume::setProgressCallback(std::function<void(std::size_t, std::size_t)> callback)
{
    progressCallback_ = std::move(callback);
}

void TetraFiniteVolume::setConsistentDensityHook(
    std::function<bool(double,double,double,double&)> hook)
{
    if(flux_)flux_->setConsistentDensityHook(std::move(hook));
}
std::size_t TetraFiniteVolume::consistentDensityFallbacks() const
{
    return flux_?flux_->consistentDensityFallbacks():0;
}
void TetraFiniteVolume::computeResidual(const TetraBoundaryProvider& boundary,std::vector<ConservedVars>& residual,TetraKernelFailure* failure) const
{
    if(failure)*failure=TetraKernelFailure{};
    // One flux per shared face, accumulated with opposite signs. A test with
    // unequal adjacent states is essential: uniform states hide wrong neighbours.
    for(std::size_t c=0;c<cellCount();++c)
    {
        try{admit(W_[c],Q_[c]);}
        catch(const std::exception& e)
        {
            residual.clear();
            if(failure){failure->failed=true;failure->cell=c;failure->detail=e.what();return;}
            throw;
        }
    }
    std::size_t bad=mesh_.faceCount();std::string detail;
    const std::size_t progressTotal=mesh_.faceCount();
    if(progressCallback_)progressCallback_(0,progressTotal);
    #pragma omp parallel for schedule(static)
    for(std::size_t f=0;f<mesh_.faceCount();++f)
    {
        if(progressCallback_&&(f&1023u)==0u)progressCallback_(f,progressTotal);
        try
        {
            const auto& face=mesh_.face(f);const auto& wl=W_[face.owner];const auto& ql=Q_[face.owner];
            PrimitiveVars wr(species_.size());thermo::EquilibriumState qr;
            if(face.neighbour>=0){wr=W_[face.neighbour];qr=Q_[face.neighbour];}
            else
            {
                if(!boundary)throw std::runtime_error("missing tetra boundary provider");
                boundary(face.owner,f,wr,qr);admit(wr,qr);
            }
            auto& result=faceFlux_[f];result=ConservedVars(species_.size());
            const bool prescribed=face.neighbour<0&&boundaryFlux_&&boundaryFlux_(face.owner,f,result);
            if(!prescribed)
                flux_->computeFluxCached(wl,wr,normals_[f],result,ql.rho,ql.acoustic.soundSpeed(),ql.h,qr.rho,qr.acoustic.soundSpeed(),qr.h);
            if(!prescribed)
            {
                const auto& upwind=result.rho>=0.?wl:wr;
                result.rhoZ=result.rho*upwind.Z;result.rhoC=result.rho*upwind.C;
            }
            if(!std::isfinite(result.rho)||!std::isfinite(result.rhoE)
               ||!std::isfinite(result.rhoZ)||!std::isfinite(result.rhoC))
                throw std::runtime_error("nonfinite tetra flux");
            for(double value:result.rhoU)if(!std::isfinite(value))throw std::runtime_error("nonfinite tetra momentum flux");
            for(double value:result.rhoY)if(!std::isfinite(value))throw std::runtime_error("nonfinite tetra species flux");
        }
        catch(const std::exception& e)
        {
            #pragma omp critical(rfam_tetra_flux_failure)
            { if(f<bad){bad=f;detail=e.what();} }
        }
    }
    if(bad<mesh_.faceCount())
    {
        residual.clear();
        if(failure){failure->failed=true;failure->cell=mesh_.face(bad).owner;failure->detail=detail;return;}
        throw std::runtime_error(detail);
    }
    residual.assign(cellCount(),ConservedVars(species_.size()));
    #pragma omp parallel for schedule(static)
    for(std::size_t c=0;c<cellCount();++c)
        for(auto f:mesh_.cellFaces(c))
            addScaled(residual[c],faceFlux_[f],mesh_.face(f).owner==int(c)?areas_[f]:-areas_[f]);
    if(molecularTransport_)
    {
        std::vector<TetraVec3> velocity(cellCount());
        std::vector<double> temperature(cellCount()),mu(cellCount()),lambda(cellCount());
        // Molecular-only transport, cached for the turbulence operator. It must
        // not see the eddy contribution it helped create: lastMu()/lastLambda()
        // feed nu = mu/rho into the SST wall omega and the nut/nu ceiling, and
        // returning mu + mu_t there is a positive feedback (nut ran to its
        // ceiling within two steps before this split).
        std::vector<double> muMolecular(cellCount()),lambdaMolecular(cellCount());
        std::vector<double> density(cellCount()),heatCapacity(cellCount());
        using TransportKey=std::tuple<double,double,std::vector<double>,int,double>;
        std::map<TransportKey,std::pair<double,double>> samples;
        for(std::size_t c=0;c<cellCount();++c)
        {
            const double quality=std::isfinite(Q_[c].vaporMassFraction)
                ?Q_[c].vaporMassFraction:-1.;
            const auto key=TransportKey{W_[c].p,W_[c].T,W_[c].Y,
                static_cast<int>(Q_[c].region),quality};
            auto it=samples.find(key);
            if(it==samples.end())
            {
                const auto properties=molecularTransport_(W_[c],Q_[c]);
                it=samples.emplace(key,properties).first;
            }
            const auto& properties=it->second;
            if(!(properties.first>0.)||!(properties.second>0.)
                ||!std::isfinite(properties.first)||!std::isfinite(properties.second))
                throw std::runtime_error("invalid tetra molecular transport");
            velocity[c]={W_[c].U[0],W_[c].U[1],W_[c].U[2]};
            temperature[c]=W_[c].T;
            mu[c]=properties.first;lambda[c]=properties.second;
            muMolecular[c]=properties.first;lambdaMolecular[c]=properties.second;
            if(turbulentTransport_)
            {
                const auto extra=turbulentTransport_(c,W_[c],Q_[c]);
                if(!std::isfinite(extra.first)||!std::isfinite(extra.second)
                   ||extra.first<0.)
                    throw std::runtime_error("invalid tetra turbulent transport");
                mu[c]+=extra.first;lambda[c]+=extra.second;
            }
            density[c]=Q_[c].rho;
            heatCapacity[c]=Q_[c].thermodynamics?Q_[c].thermodynamics->Cv:0.;
        }
        lastMu_=muMolecular; lastLambda_=lambdaMolecular;
        auto spectral=tetraMolecularSpectral(mesh_,density,heatCapacity,mu,lambda);
        addViscousResidual(mesh_,velocity,temperature,mu,lambda,residual,true);
        molecularSpectral_.swap(spectral);
    }
}
void TetraFiniteVolume::localTimeStep(double cfl,std::vector<double>& dt) const
{
    if(!std::isfinite(cfl)||cfl<=0)throw std::invalid_argument("CFL must be positive");
    dt.assign(cellCount(),0.);
    for(std::size_t c=0;c<cellCount();++c)
    {
        double sum=0;
        for(auto f:mesh_.cellFaces(c))
        {
            const auto& face=mesh_.face(f);int other=face.owner==int(c)?face.neighbour:face.owner;
            // Wave speed from the single owner (TetraFaceInvariant.H). The
            // ghost side of a boundary face has no state yet, so it
            // contributes zero here; the boundary-aware overload below is the
            // one that includes the trace. Both must use the same invariant or
            // the explicit and implicit step estimates disagree on stiff
            // boundary cells, which is the cell class the injector run stalls on.
            const auto speed=[&](std::size_t j)
            {return tetraFaceInvariant(W_[j],Q_[j],normals_[f]).spectralRadius;};
            sum+=areas_[f]*std::max(speed(c),other<0?0.:speed(other));
        }
        dt[c]=cfl*mesh_.volume(c)/(sum+molecularSpectral(c));
    }
}
void TetraFiniteVolume::localTimeStep(double cfl,std::vector<double>& dt,
                                     const TetraBoundaryProvider& boundary) const
{
    if(!std::isfinite(cfl)||cfl<=0.)throw std::invalid_argument("CFL must be positive");
    if(!boundary)throw std::invalid_argument("boundary-aware time step needs a boundary provider");
    dt.assign(cellCount(),0.);
    for(std::size_t c=0;c<cellCount();++c)
    {
        double sum=0.;
        for(auto f:mesh_.cellFaces(c))
        {
            const auto& face=mesh_.face(f);
            // Boundary-aware step: the boundary trace participates in the
            // spectral bound because it can be much faster than the adjacent
            // cell. Same invariant, same justification as above.
            const auto speed=[&](const PrimitiveVars& w,const thermo::EquilibriumState& q)
            {return tetraFaceInvariant(w,q,normals_[f]).spectralRadius;};
            double maximum=speed(W_[c],Q_[c]);
            if(face.neighbour>=0)
                maximum=std::max(maximum,speed(W_[face.neighbour],Q_[face.neighbour]));
            else
            {
                PrimitiveVars ghost(species_.size());thermo::EquilibriumState q;
                boundary(face.owner,f,ghost,q);admit(ghost,q);
                maximum=std::max(maximum,speed(ghost,q));
            }
            sum+=areas_[f]*maximum;
        }
        dt[c]=cfl*mesh_.volume(c)/(sum+molecularSpectral(c));
    }
}
double TetraFiniteVolume::faceMassFlow(std::size_t f,const PrimitiveVars& wl,const PrimitiveVars& wr,const thermo::EquilibriumState& ql,const thermo::EquilibriumState& qr) const
{
    admit(wl,ql);admit(wr,qr);ConservedVars result(species_.size());
    const auto& face=mesh_.face(f);
    if(!(face.neighbour<0&&boundaryFlux_&&boundaryFlux_(face.owner,f,result)))
        flux_->computeFluxCached(wl,wr,normals_[f],result,ql.rho,ql.acoustic.soundSpeed(),ql.h,qr.rho,qr.acoustic.soundSpeed(),qr.h);
    return result.rho*areas_[f];
}
} // namespace Foam
