#ifndef TETRA_TURBULENCE_OPERATOR_H
#define TETRA_TURBULENCE_OPERATOR_H

/*---------------------------------------------------------------------------*\
  TetraTurbulenceOperator - Menter k-omega SST for the tetrahedral path.

  Why this file exists
  --------------------
  src/turbulence/KOmegaSST.H is the closure only (constants, blending
  functions, source terms). It was wired into the structured sector path via
  SectorTurbulenceOperator; the unstructured tetra driver had no equivalent, so
  the tetra path stayed inviscid-slip with no boundary layer and could not
  answer a separation question. This operator fills the slot recorded in
  docs/PROJECT_TREE.md 1.7: the eddy viscosity enters through the same
  mu/lambda channel as the molecular transport in TetraFiniteVolume, so it also
  participates in the LU-SGS diagonal (molecularSpectral_).

  Documented limitations (do not widen silently)
  ----------------------------------------------
  1. One-way Boussinesq coupling: the closure reads (rho, nu, S) and returns
     nu_t. It is NOT variable-density or Favre-averaged.
  2. Pr_t = 0.85 fixed; lambda_t = mu_t*cp/Pr_t with cp = Cv + p/(rho*T).
  3. Operator split with the flow step.
  4. Wall distance via Dijkstra graph propagation (accuracy ~ one cell width).
\*---------------------------------------------------------------------------*/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

#include "TetraFiniteVolume.H"
#include "TetraMesh.H"
#include "KOmegaSST.H"

namespace Foam
{

class TetraTurbulenceOperator
{
public:
    struct Settings
    {
        double intensity=0.05;
        double lengthScale=1e-3;
        double prandtlT=0.85;
        double wallOmegaFactor=60.0;
        int wallZone=0;
    };

    TetraTurbulenceOperator(TetraFiniteVolume& fv)
        :TetraTurbulenceOperator(fv,Settings{}) {}

    TetraTurbulenceOperator(TetraFiniteVolume& fv,const Settings& settings)
        :fv_(fv),mesh_(fv.mesh()),settings_(settings)
    {
        if(!(settings_.intensity>0.)||!(settings_.intensity<1.))
            throw std::invalid_argument("turbulence intensity must be in (0,1)");
        if(!(settings_.lengthScale>0.))
            throw std::invalid_argument("turbulence length scale must be positive");
        const std::size_t n=mesh_.cellCount();
        k_.assign(n,KOmegaSST::k_min);
        omega_.assign(n,KOmegaSST::om_min);
        nut_.assign(n,0.);
        y_.assign(n,0.);
        strain_.assign(n,0.);
        buildWallDistance();
    }

    const std::vector<double>& kField() const {return k_;}
    const std::vector<double>& omegaField() const {return omega_;}
    const std::vector<double>& wallDistanceField() const {return y_;}
    // Roll the transported fields back to a snapshot. The physical-time path
    // rejects a trial step by restoring every cell state; the turbulence fields
    // are not part of that transaction, so they must be restored explicitly.
    void restoreFields(const std::vector<double>& kSnapshot,
                       const std::vector<double>& omegaSnapshot)
    {
        if(kSnapshot.size()!=k_.size()||omegaSnapshot.size()!=omega_.size())
            throw std::runtime_error("tetra turbulence: snapshot size mismatch");
        k_=kSnapshot;
        omega_=omegaSnapshot;
        applyWallBoundary();
        updateEddyViscosity();
    }
    double turbulentViscosity(std::size_t c) const {return nut_.at(c);}
    std::size_t wallAdjacentCells() const {return wallAdjacent_;}
    double maxNut() const
    {
        if(nut_.empty()) return 0.;
        // Restrict to cells that carry gas: in near-vacuum cells nu = mu/rho is
        // huge, the nut/nu ceiling stops binding, and the reported maximum
        // says nothing about the boundary layer.
        double m=0.;
        for(std::size_t c=0;c<nut_.size();++c)
            if(fv_.Q(c).rho>1.0) m=std::max(m,nut_[c]);
        return m;
    }

    // (mu_t, lambda_t) added to the molecular transport of every cell.
    std::pair<double,double> transportAddition(
        std::size_t c,const PrimitiveVars& w,const thermo::EquilibriumState& q) const
    {
        const double nut=nut_.at(c);
        if(!(nut>0.)||!std::isfinite(nut)) return {0.,0.};
        const double muT=q.rho*nut;
        double cp=0.;
        if(q.thermodynamics) cp=q.thermodynamics->Cv;
        if(std::isfinite(w.T)&&w.T>0.&&q.rho>0.) cp+=w.p/(q.rho*w.T);
        double lamT=0.;
        if(cp>0.&&std::isfinite(cp)) lamT=muT*cp/settings_.prandtlT;
        if(!std::isfinite(muT)||!std::isfinite(lamT)) return {0.,0.};
        return {muT,lamT};
    }

    void initialize()
    {
        const std::size_t n=mesh_.cellCount();
        for(std::size_t c=0;c<n;++c)
        {
            const auto& w=fv_.W(c);
            const double speed=std::sqrt(w.U[0]*w.U[0]+w.U[1]*w.U[1]+w.U[2]*w.U[2]);
            const double uRef=std::max(speed,1.0);
            const double kIn=1.5*std::pow(settings_.intensity*uRef,2);
            const double omIn=std::sqrt(std::max(kIn,KOmegaSST::k_min))
                /(std::pow(KOmegaSST::beta_star,0.25)*settings_.lengthScale);
            k_[c]=std::max(kIn,KOmegaSST::k_min);
            omega_[c]=std::max(omIn,KOmegaSST::om_min);
        }
        applyWallBoundary();
        updateEddyViscosity();
    }

    void advance(const std::vector<double>& dtCell)
    {
        const std::size_t n=mesh_.cellCount();
        if(dtCell.size()!=n) throw std::runtime_error("tetra turbulence: dt size mismatch");
        if(fv_.lastMu().size()!=n)
            throw std::runtime_error("tetra turbulence: molecular transport must be enabled");
        updateEddyViscosity();

        const auto& mu=fv_.lastMu();
        std::vector<double> nextK=k_,nextOmega=omega_;
        wallAdjacent_=0;
        for(std::size_t c=0;c<n;++c)
        {
            if(isWallAdjacent(c)){++wallAdjacent_;continue;}
            const double rho=std::max(fv_.Q(c).rho,1e-30);
            const double nu=mu[c]/rho;
            const double f2=KOmegaSST::F2(k_[c],omega_[c],y_[c],nu);
            const double nut=KOmegaSST::nu_t(k_[c],omega_[c],strain_[c],f2);
            const double pk=KOmegaSST::Pk(nut,strain_[c],k_[c],omega_[c]);
            const double cd=crossDiffusion(c);
            const double f1=KOmegaSST::F1(k_[c],omega_[c],y_[c],nu,std::max(cd,1e-20));
            const double beta=KOmegaSST::beta(f1);
            const double alpha=KOmegaSST::alpha(f1);

            double rateK=0.,rateOmega=0.;
            transportRate(c,k_,KOmegaSST::sigma_k(f1),rateK);
            transportRate(c,omega_,KOmegaSST::sigma_w(f1),rateOmega);

            // Point-implicit (backward Euler) update. The steady flow step is a
            // pseudo-time step, up to two orders larger than the physical step
            // the SST source terms tolerate, and the wall omega (60 nu/(beta dy^2),
            // order 1e10 1/s on the 0.2 um wall cells) makes an explicit update
            // diverge within one step. Backward Euler on both equations is
            // unconditionally stable and keeps each cell monotone:
            //   k:      (k*-k)/dt = Pk + rateK - beta* k* omega*
            //   omega:  (w*-w)/dt = alpha Pk/nut + CD + rateOmega - beta w*^2
            // The omega equation is a quadratic; the stable root form is used.
            const double dt=std::max(dtCell[c],0.);
            const double al=(nut>1e-30)?alpha*pk/nut:0.;
            double kk=k_[c], om=omega_[c];
            if(dt>0.)
            {
                const double sourceOmega=al+cd+rateOmega;
                const double R=om+dt*sourceOmega;
                const double disc=1.+4.*beta*dt*std::max(R,0.);
                const double wNew=2.*std::max(R,0.)/(1.+std::sqrt(disc));
                const double kNew=(kk+dt*(pk+rateK))/(1.+dt*KOmegaSST::beta_star*wNew);
                om=wNew;
                kk=kNew;
            }
            nextK[c]=std::min(std::max(kk,KOmegaSST::k_min),kMax_);
            nextOmega[c]=std::min(std::max(om,KOmegaSST::om_min),omMax_);
        }
        k_.swap(nextK);
        omega_.swap(nextOmega);
        applyWallBoundary();
        updateEddyViscosity();
    }

private:
    TetraFiniteVolume& fv_;
    const TetraMesh& mesh_;
    const Settings settings_;
    std::vector<double> k_,omega_,nut_,y_,strain_;
    std::size_t wallAdjacent_=0;
    // Guard rails: an explicit SST update on a coarse mesh can otherwise run
    // away through the nut = k/omega coupling. These bounds are far outside any
    // physical chamber value (k ~ 1e3 m2/s2, omega ~ 1e6 1/s, nut ~ 1 m2/s).
    static constexpr double kMax_=1e7;
    static constexpr double omMax_=1e12;
    static constexpr double nutMax_=1e3;
    static constexpr double maxNutRatio_=1e4;

    bool isWallAdjacent(std::size_t c) const
    {
        for(std::int32_t f:mesh_.cellFaces(c))
        {
            const auto& face=mesh_.face(f);
            if(face.neighbour<0&&face.boundaryZone==settings_.wallZone) return true;
        }
        return false;
    }

    // ---- Dijkstra wall-distance field -------------------------------------
    // Phase 1: seed wall-adjacent cells with Euclidean distance to nearest
    //   wall face centroid (exact for the adjacent layer).
    // Phase 2: propagate through the cell connectivity graph by relaxation.
    //   Accuracy: O(cell width), which is sufficient for SST F1/F2.
    //   Complexity: O(N log N), scales linearly with mesh size.
    void buildWallDistance()
    {
        std::cout<<"[turbulence] building wall-distance field (Dijkstra) for "
                 <<mesh_.cellCount()<<" cells...";
        std::cout.flush();

        const std::size_t n=mesh_.cellCount();
        const double INF=std::numeric_limits<double>::max();
        y_.assign(n,INF);

        // Phase 1: wall-adjacent cells
        bool anyWall=false;
        for(std::size_t f=0;f<mesh_.faceCount();++f)
        {
            const auto& face=mesh_.face(f);
            if(face.neighbour>=0) continue;
            if(face.boundaryZone!=settings_.wallZone) continue;
            const std::size_t c=static_cast<std::size_t>(face.owner);
            const auto fc=mesh_.faceCentroid(f);
            const auto cc=mesh_.cellCentroid(c);
            const double dx=fc.x-cc.x,dy=fc.y-cc.y,dz=fc.z-cc.z;
            const double d=std::sqrt(dx*dx+dy*dy+dz*dz);
            if(d<y_[c]) y_[c]=d;
            anyWall=true;
        }
        if(!anyWall)
            throw std::runtime_error(
                "tetra turbulence: mesh has no wall faces in zone "
                +std::to_string(settings_.wallZone));

        // Phase 2: Dijkstra propagation
        using P2=std::pair<double,std::size_t>;
        std::priority_queue<P2,std::vector<P2>,std::greater<P2>> pq;
        for(std::size_t c=0;c<n;++c)
            if(y_[c]<INF) pq.push({y_[c],c});

        while(!pq.empty())
        {
            auto [d,c]=pq.top(); pq.pop();
            if(d>y_[c]) continue;  // stale
            const auto cc=mesh_.cellCentroid(c);
            for(std::int32_t fi:mesh_.cellFaces(c))
            {
                const auto& face=mesh_.face(fi);
                if(face.neighbour<0) continue;
                const std::size_t nb=(static_cast<std::size_t>(face.owner)==c)
                    ?static_cast<std::size_t>(face.neighbour)
                    :static_cast<std::size_t>(face.owner);
                const auto nc=mesh_.cellCentroid(nb);
                const double ex=nc.x-cc.x,ey=nc.y-cc.y,ez=nc.z-cc.z;
                const double nd=d+std::sqrt(ex*ex+ey*ey+ez*ez);
                if(nd<y_[nb]){ y_[nb]=nd; pq.push({nd,nb}); }
            }
        }

        const double yMin=*std::min_element(y_.begin(),y_.end());
        const double yMax=*std::max_element(y_.begin(),y_.end());
        std::cout<<" done. y_min="<<yMin<<" y_max="<<yMax<<" m\n";
    }

    void applyWallBoundary()
    {
        const std::size_t n=mesh_.cellCount();
        const auto& mu=fv_.lastMu();
        for(std::size_t c=0;c<n;++c)
        {
            if(!isWallAdjacent(c)) continue;
            const double dy=std::max(y_[c],1e-9);
            const double rho=std::max(fv_.Q(c).rho,1e-30);
            const double nu=(mu.size()==n)?mu[c]/rho:0.;
            const double omW=(nu>0.)
                ?settings_.wallOmegaFactor*nu/(KOmegaSST::beta1*dy*dy)
                :KOmegaSST::om_min;
            k_[c]=KOmegaSST::k_min;
            omega_[c]=std::max(omW,KOmegaSST::om_min);
        }
    }

    void updateEddyViscosity()
    {
        const std::size_t n=mesh_.cellCount();
        const auto& mu=fv_.lastMu();
        for(std::size_t c=0;c<n;++c)
        {
            double g[9]={0,0,0,0,0,0,0,0,0};
            const double V=mesh_.volume(c);
            if(V>0.)
            {
                for(std::int32_t f:mesh_.cellFaces(c))
                {
                    const auto& face=mesh_.face(f);
                    if(face.neighbour<0) continue;
                    const std::size_t other=(static_cast<std::size_t>(face.owner)==c)
                        ?static_cast<std::size_t>(face.neighbour)
                        :static_cast<std::size_t>(face.owner);
                    const auto av=mesh_.faceAreaVector(f);
                    for(int m=0;m<3;++m)
                    {
                        const double phi=0.5*(fv_.W(c).U[m]+fv_.W(other).U[m]);
                        g[m*3+0]+=phi*av.x;
                        g[m*3+1]+=phi*av.y;
                        g[m*3+2]+=phi*av.z;
                    }
                }
                for(double& v:g) v/=V;
            }
            const double e12=0.5*(g[1]+g[3]),e13=0.5*(g[2]+g[6]),e23=0.5*(g[5]+g[7]);
            const double s2=2.0*(g[0]*g[0]+g[4]*g[4]+g[8]*g[8]
                                 +2.0*(e12*e12+e13*e13+e23*e23));
            strain_[c]=std::sqrt(std::max(s2,0.));

            const double rho=std::max(fv_.Q(c).rho,1e-30);
            if(rho<1e-6)
            {
                // Degenerate (sliver) cell: there is no meaningful kinematic
                // viscosity to scale against, and nut = k/omega would report a
                // spurious maximum. No eddy viscosity here.
                strain_[c]=0.;
                nut_[c]=0.;
                continue;
            }
            const double nu=(mu.size()==n)?mu[c]/rho:1e-30;
            const double f2=KOmegaSST::F2(k_[c],omega_[c],y_[c],nu);
            double nut=KOmegaSST::nu_t(k_[c],omega_[c],strain_[c],f2);
            if(!std::isfinite(nut)||nut<0.) nut=0.;
            // Practical eddy-viscosity ratio ceiling. SST has no intrinsic bound
            // when omega decays, and nut = k/omega then runs away; nut/nu above
            // 1e4 is outside any chamber measurement.
            nut=std::min(nut,std::min(nutMax_,maxNutRatio_*nu));
            nut_[c]=nut;
        }
    }

    void transportRate(std::size_t c,const std::vector<double>& field,
                       double sigma,double& rate) const
    {
        rate=0.;
        const double V=mesh_.volume(c);
        const double rho=std::max(fv_.Q(c).rho,1e-30);
        if(!(V>0.)) return;
        const auto& mu=fv_.lastMu();
        double acc=0.;
        for(std::int32_t f:mesh_.cellFaces(c))
        {
            const auto& face=mesh_.face(f);
            if(face.neighbour<0) continue;
            const std::size_t other=(static_cast<std::size_t>(face.owner)==c)
                ?static_cast<std::size_t>(face.neighbour)
                :static_cast<std::size_t>(face.owner);
            const auto av=mesh_.faceAreaVector(f);
            const double areaMag=std::sqrt(av.x*av.x+av.y*av.y+av.z*av.z);
            if(areaMag<=0.) continue;
            const auto a=mesh_.cellCentroid(c),b=mesh_.cellCentroid(other);
            const double dist=std::max(((b.x-a.x)*av.x+(b.y-a.y)*av.y
                                       +(b.z-a.z)*av.z)/areaMag,1e-12);
            const auto& wc=fv_.W(c);const auto& wn=fv_.W(other);
            const double rc=fv_.Q(c).rho,rn=fv_.Q(other).rho;
            const double massFlux=0.5*((rc*wc.U[0]+rn*wn.U[0])*av.x
                                      +(rc*wc.U[1]+rn*wn.U[1])*av.y
                                      +(rc*wc.U[2]+rn*wn.U[2])*av.z);
            const double faceValue=massFlux>=0.?field[c]:field[other];
            const double muFace=0.5*(mu[c]+mu[other]);
            const double muT=0.5*(rc*nut_[c]+rn*nut_[other]);
            const double diffusion=(muFace+sigma*muT)*(field[other]-field[c])
                *areaMag/dist;
            acc+=(-massFlux*(faceValue-field[c])+diffusion);
        }
        rate=acc/(rho*V);
    }

    double crossDiffusion(std::size_t c) const
    {
        const double V=mesh_.volume(c);
        if(!(V>0.)) return 0.;
        double gk[3]={0,0,0},go[3]={0,0,0};
        for(std::int32_t f:mesh_.cellFaces(c))
        {
            const auto& face=mesh_.face(f);
            if(face.neighbour<0) continue;
            const std::size_t other=(static_cast<std::size_t>(face.owner)==c)
                ?static_cast<std::size_t>(face.neighbour)
                :static_cast<std::size_t>(face.owner);
            const auto av=mesh_.faceAreaVector(f);
            const double phiK=0.5*(k_[c]+k_[other]);
            const double phiO=0.5*(omega_[c]+omega_[other]);
            gk[0]+=phiK*av.x;gk[1]+=phiK*av.y;gk[2]+=phiK*av.z;
            go[0]+=phiO*av.x;go[1]+=phiO*av.y;go[2]+=phiO*av.z;
        }
        for(int d=0;d<3;++d){gk[d]/=V;go[d]/=V;}
        return KOmegaSST::crossDiffusion(omega_[c],gk,go);
    }
};

} // namespace Foam
#endif
