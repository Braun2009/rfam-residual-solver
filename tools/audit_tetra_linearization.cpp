// Controlled reproduction of tetra LU-SGS composition/energy-reference errors.
// No CAD or boundary physics changes. The energy offset is a species-linear
// reference shift, so nonreacting physical results should be invariant.
#include "TetraFiniteVolume.H"
#include "TetraLUSGS.H"
#include "CoolPropThermoModel.H"
#include "CompleteCellState.H"
#include "FrozenEulerFluxJacobian.H"
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
using namespace Foam;
namespace {
class EnergyShiftModel final:public thermo::EquilibriumThermoModel
{
    std::shared_ptr<const thermo::EquilibriumThermoModel> base_;
    double shift_;
    std::string id_;
    thermo::EquilibriumResult translate(thermo::EquilibriumResult r)const
    {
        if(!r.success())return r;
        auto& q=*r.state;
        q.e+=shift_*q.overallMassFractions[1];q.h+=shift_*q.overallMassFractions[1];
        for(auto* phase:{&q.singlePhase,&q.liquid,&q.vapor})
            if(*phase){(**phase).e+=shift_*(**phase).massFractions[1];(**phase).h+=shift_*(**phase).massFractions[1];}
        r.diagnostics.modelId=id_;
        if(!q.valid())throw std::runtime_error("energy reference shift violated state contract");
        return r;
    }
public:
    EnergyShiftModel(std::shared_ptr<const thermo::EquilibriumThermoModel> base,double shift)
      :base_(base),shift_(shift),id_(base->modelId()+";diagnostic_energy_shift="+std::to_string(shift)){}
    thermo::EquilibriumResult initializePTY(double p,double t,const std::vector<double>& y)const override
    {return translate(base_->initializePTY(p,t,y));}
    thermo::EquilibriumResult recoverRhoE(double rho,double e,const std::vector<double>& y,
                                         const thermo::EquilibriumSeed* seed=nullptr)const override
    {return translate(base_->recoverRhoE(rho,e-shift_*y.at(1),y,seed));}
    bool stateValidForModel(const thermo::EquilibriumState& q)const override{return base_->stateValidForModel(q);}
    const std::string& modelId()const override{return id_;}
    const std::vector<std::string>& speciesNames()const override{return base_->speciesNames();}
    EnergyConvention energyConvention()const override{return base_->energyConvention();}
    bool supportsPrimitivePTYInitialization()const override{return true;}
};
CompleteCellState at(const thermo::EquilibriumThermoModel& model,double p,double t,double oxygen,double speed)
{
    PrimitiveVars w(2);w.p=p;w.T=t;w.Y={oxygen,1.-oxygen};w.U={0.,0.,speed};
    const auto r=CompleteStateConversion::fromPrimitive(w,model);
    if(!r.success())throw std::runtime_error(r.diagnostics.backendStatus);
    return *r.state;
}
ConservedVars difference(const ConservedVars& a,const ConservedVars& b,double denominator)
{
    ConservedVars d(2);d.rho=(a.rho-b.rho)/denominator;d.rhoE=(a.rhoE-b.rhoE)/denominator;
    for(int i=0;i<3;++i)d.rhoU[i]=(a.rhoU[i]-b.rhoU[i])/denominator;
    for(int i=0;i<2;++i)d.rhoY[i]=(a.rhoY[i]-b.rhoY[i])/denominator;
    d.rhoZ=(a.rhoZ-b.rhoZ)/denominator;d.rhoC=(a.rhoC-b.rhoC)/denominator;return d;
}
void writeMesh(const std::filesystem::path& path)
{
    std::ofstream out(path,std::ios::binary);
    auto put=[&](const auto& x){out.write(reinterpret_cast<const char*>(&x),sizeof(x));};
    out.write("RFAMTET1",8);std::uint64_t np=5,nc=2,nf=7;put(np);put(nc);put(nf);
    double points[15]={0,0,0, .01,0,0, 0,.01,0, 0,0,.01, 0,0,-.01};put(points);
    std::int32_t cells[8]={0,1,2,3,0,2,1,4};put(cells);
    double vols[2]={1e-6/6,1e-6/6};put(vols);
    std::int32_t faces[42]={0,2,1,0,1,0,0,1,3,0,-1,0,1,2,3,0,-1,0,2,0,3,0,-1,0,
                          0,2,4,1,-1,0,2,1,4,1,-1,0,1,0,4,1,-1,0};put(faces);
}
}
int main(int argc,char**argv)
{
    const bool check=argc==2&&std::string(argv[1])=="--check";
    const bool coupled=check||(argc==2&&std::string(argv[1])=="--coupled");
    const bool contact=argc==2&&std::string(argv[1])=="--contact";
    if(argc>2||(argc==2&&!coupled&&!contact)){std::cerr<<"usage: audit_tetra_linearization [--coupled|--check|--contact]\n";return 2;}
    if(!thermo::coolPropThermoAvailable()){std::cerr<<"CoolProp native backend required\n";return 77;}
    const auto path=std::filesystem::temp_directory_path()/
      ("rfam-linearization-"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count())+".rfamtet");
    try
    {
        auto base=thermo::makeCoolPropThermoModel();
        std::cout<<std::setprecision(17)<<"backend="<<base->modelId()<<"\ncomposition_coupling="<<coupled<<"\n";
        if(contact)
        {
            // Exact inviscid material contact: both sides have equal p and u.
            // Moving the unresolved contact across a FV cell gives this volume
            // average of conserved quantities. Exact continuum p stays fixed;
            // single homogeneous EOS recovery need not preserve it.
            const auto ambient=at(*base,4e6,400.,.784,20.);
            for(const auto& kind:{std::string("same_state"),std::string("fuel"),std::string("oxygen")})
            {
                const auto incoming=kind=="same_state"?ambient:kind=="fuel"?
                  at(*base,4e6,453.5,0.,20.):at(*base,4e6,133.5,1.,20.);
                for(double f:{0.,1e-4,1e-3,.01,.03,.1,.3,.5,.9,1.})
                {
                    auto u=ambient.U;
                    u.rho=(1-f)*ambient.U.rho+f*incoming.U.rho;
                    u.rhoE=(1-f)*ambient.U.rhoE+f*incoming.U.rhoE;
                    for(int i=0;i<3;++i)u.rhoU[i]=(1-f)*ambient.U.rhoU[i]+f*incoming.U.rhoU[i];
                    for(int i=0;i<2;++i)u.rhoY[i]=(1-f)*ambient.U.rhoY[i]+f*incoming.U.rhoY[i];
                    const auto r=recoverTetraState(u,ambient.W,ambient.Q,*base);
                    std::printf("MATERIAL_CONTACT_INPUT kind=%s volume_fraction=%.17g rho=%.17g e=%.17g yO2=%.17g\n",
                        kind.c_str(),f,u.rho,u.rhoE/u.rho-0.5*400.,u.rhoY[0]/u.rho);
                    std::cout<<"MATERIAL_CONTACT kind="<<kind<<" volume_fraction="<<f<<" success="<<r.success();
                    if(r.success())std::cout<<" p="<<r.state->Q.p<<" dp_Pa="<<r.state->Q.p-4e6
                      <<" T="<<r.state->Q.T<<" Y_O2="<<r.state->W.Y[0]
                      <<" region="<<thermo::equilibriumRegionName(r.state->Q.region);
                    else std::cout<<" status="<<r.diagnostics.backendStatus;
                    std::cout<<std::endl;
                }
            }
            return 0;
        }
        writeMesh(path);TetraMesh mesh;mesh.read(path.string());
        struct Physical{double p,t,y,u;};
        std::vector<Physical> explicitReference,implicitReference;
        double maxExplicitDifference=0.,maxImplicitDifference=0.;
        for(double shift:{0.,1e6,-1e6})
        {
            EnergyShiftModel model(base,shift);
            const auto state=at(model,4e6,400.,.6,0.);
            constexpr double dy=1e-5;
            const auto hi=at(model,4e6,400.,.6+dy,0.),lo=at(model,4e6,400.,.6-dy,0.);
            const auto d=difference(hi.U,lo.U,2.*dy);
            ConservedVars predicted(2);
            if(!FrozenEulerFluxJacobian::apply(state.W,state.Q,d,{0.,0.,1.},predicted))
                throw std::runtime_error("frozen Euler application failed");
            ConservedVars corrected(2);
            const auto coefficient=TetraBinaryEulerJacobian::prepare(state.Q,model);
            if(!TetraBinaryEulerJacobian::apply(state.W,state.Q,d,{0.,0.,1.},coefficient,corrected))
                throw std::runtime_error("composition-coupled Euler application failed");
            std::cout<<"COMPOSITION_TANGENT shift="<<shift
                     <<" dP_dY_exact=0 frozen_dP_dY="<<predicted.rhoU[2]
                     <<" coupled_dP_dY="<<corrected.rhoU[2]<<"\n";
            if(check&&std::abs(corrected.rhoU[2])>1e-6*state.Q.p)
                throw std::runtime_error("composition tangent pressure cancellation failed");
            // Complement the isobaric tangent with arbitrary p,T,Y changes
            // and check directional Euler momentum/energy flux derivatives.
            const auto moving=at(model,4e6,400.,.6,20.);
            const auto plus=at(model,4e6+10.,400.+.001,.6+dy,20.);
            const auto minus=at(model,4e6-10.,400.-.001,.6-dy,20.);
            const auto perturbation=difference(plus.U,minus.U,2.);
            const auto coefficients=TetraBinaryEulerJacobian::prepare(moving.Q,model);
            ConservedVars product(2);
            if(!TetraBinaryEulerJacobian::apply(moving.W,moving.Q,perturbation,{0.,0.,1.},coefficients,product))
                throw std::runtime_error("arbitrary tangent application failed");
            const double expectedMomentum=400.*perturbation.rho+10.;
            const double expectedEnergy=.5*(plus.U.rho*20.*(plus.Q.h+200.)-minus.U.rho*20.*(minus.Q.h+200.));
            if(check&&(std::abs(product.rhoU[2]-expectedMomentum)>1e-6*std::max(std::abs(expectedMomentum),1.)
              ||std::abs(product.rhoE-expectedEnergy)>1e-6*std::max(std::abs(expectedEnergy),1.)))
                throw std::runtime_error("coupled Euler finite difference mismatch");
            std::size_t comparison=0;
            for(double oxygenJump:{0.,.02,.4})
            {
                TetraFiniteVolume fv(mesh,base->speciesNames(),std::make_unique<AUSMFluxAdapter>(base->speciesNames()));
                for(int c=0;c<2;++c)
                {
                    auto initial=at(model,4e6,400.,.6+(c==0?-.5:.5)*oxygenJump,20.);
                    fv.setState(c,initial.W,initial.U,initial.Q);
                }
                const TetraBoundaryProvider ghost=[&](std::size_t c,std::size_t,PrimitiveVars& w,thermo::EquilibriumState& q)
                {w=fv.W(c);q=fv.Q(c);};
                std::vector<ConservedVars> residual;fv.computeResidual(ghost,residual);
                std::vector<double> dt;fv.localTimeStep(.5,dt);
                for(std::size_t c=0;c<2;++c)
                {
                    auto update=fv.U(c);double scale=.2*dt[c]/mesh.volume(c);
                    update.rho-=scale*residual[c].rho;update.rhoE-=scale*residual[c].rhoE;
                    for(int i=0;i<3;++i)update.rhoU[i]-=scale*residual[c].rhoU[i];
                    for(int i=0;i<2;++i)update.rhoY[i]-=scale*residual[c].rhoY[i];
                    auto explicitState=recoverTetraState(update,fv.W(c),fv.Q(c),model);
                    if(!explicitState.success())throw std::runtime_error("explicit reference failed: "+explicitState.diagnostics.backendStatus);
                    const auto& ew=explicitState.state->W;
                    if(shift==0.)explicitReference.push_back({ew.p,ew.T,ew.Y[0],ew.U[2]});
                    else
                    {
                        const auto ref=explicitReference.at(comparison+c);
                        maxExplicitDifference=std::max({maxExplicitDifference,std::abs(ew.p-ref.p)/4e6,
                          std::abs(ew.T-ref.t)/400.,std::abs(ew.Y[0]-ref.y),std::abs(ew.U[2]-ref.u)/20.});
                    }
                    std::cout<<"EXPLICIT shift="<<shift<<" jump="<<oxygenJump<<" cell="<<c
                             <<" p="<<explicitState.state->W.p<<" T="<<explicitState.state->W.T
                             <<" Y="<<explicitState.state->W.Y[0]<<" uz="<<explicitState.state->W.U[2]<<"\n";
                }
                TetraLUSGS lusgs(fv);lusgs.setCompositionCoupling(coupled);
                const auto report=lusgs.step(ghost,.5,1.,.2,model);
                if(check&&report.limited!=0)throw std::runtime_error("gauge test accidentally exercised limiting");
                for(std::size_t c=0;c<2;++c)
                {
                    const auto& w=fv.W(c);
                    if(shift==0.)implicitReference.push_back({w.p,w.T,w.Y[0],w.U[2]});
                    else
                    {
                        const auto ref=implicitReference.at(comparison+c);
                        maxImplicitDifference=std::max({maxImplicitDifference,std::abs(w.p-ref.p)/4e6,
                          std::abs(w.T-ref.t)/400.,std::abs(w.Y[0]-ref.y),std::abs(w.U[2]-ref.u)/20.});
                    }
                    std::cout<<"LUSGS shift="<<shift<<" jump="<<oxygenJump<<" cell="<<c
                             <<" p="<<w.p<<" T="<<w.T<<" Y="<<w.Y[0]
                             <<" uz="<<w.U[2]<<" limited="<<report.limited<<"\n";
                }
                comparison+=2;
            }
        }
        if(check&&(maxExplicitDifference>1e-9||maxImplicitDifference>1e-9))
            throw std::runtime_error("species energy reference changed physical update");
        std::cout<<"GAUGE_DIFFERENCE explicit="<<maxExplicitDifference<<" implicit="<<maxImplicitDifference<<"\n";
        if(check)std::cout<<"TETRA_BINARY_LINEARIZATION_PASS\n";
        std::filesystem::remove(path);return 0;
    }catch(const std::exception& e){std::filesystem::remove(path);std::cerr<<e.what()<<"\n";return 1;}
}
