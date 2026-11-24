// Diagnostic probe: reproduce the cell-7766 rho-e rejection recorded in
// work/coldflow_B_verified/native_last_accepted.rfamtetf.failure.txt and expose
// where the production recovery loses the state. Read-only; no solver state is
// modified. The verbose loops below are diagnostic copies of
// CoolPropThermoModel::recoverRhoE and cpdetail::recover with identical
// equations and tolerances; only logging differs.
#include "CoolPropThermoModel.H"
#include <algorithm>
#include "CoolPropEquilibrium.H"
#include "AbstractState.h"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>

namespace
{
using CPState=CoolProp::AbstractState;
void printPoint(const char* tag,CPState& s)
{
    std::printf("%s p=%.17g T=%.17g rho=%.17g u=%.17g h=%.17g s=%.17g phase=%d\n",
        tag,s.p(),s.T(),s.rhomass(),s.umass(),s.hmass(),s.smass(),int(s.phase()));
}
}

int main(int argc,char**argv)
{
    if(argc==5&&std::string(argv[1])=="--pt")
    {
        const double p=std::atof(argv[2]),t=std::atof(argv[3]),y0=std::atof(argv[4]);
        auto model=Foam::thermo::makeCoolPropThermoModel();
        const auto r=model->initializePTY(p,t,{y0,1.-y0});
        std::printf("PT_RESULT success=%d status=%s backend=%s\n",int(r.success()),
            Foam::thermo::equilibriumStatusName(r.status),r.diagnostics.backendStatus.c_str());
        if(r.state)std::printf("PT_STATE p=%.17g T=%.17g rho=%.17g e=%.17g s=%.17g region=%d\n",
            r.state->p,r.state->T,r.state->rho,r.state->e,r.state->s,int(r.state->region));
        return r.success()?0:1;
    }
    if(argc==5&&std::string(argv[1])=="--dt")
    {
        const double rho=std::atof(argv[2]),t=std::atof(argv[3]),y0=std::atof(argv[4]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        try
        {
            s->unspecify_phase();
            s->update(CoolProp::DmassT_INPUTS,rho,t);
            std::printf("DT_RESULT success=1 p=%.17g T=%.17g rho=%.17g u=%.17g s=%.17g phase=%d\n",
                s->p(),s->T(),s->rhomass(),s->umass(),s->smass(),int(s->phase()));
            return 0;
        }
        catch(const std::exception& e)
        {
            std::printf("DT_RESULT success=0 error=%s\n",e.what());
            return 1;
        }
    }
    if(argc==5&&std::string(argv[1])=="--du")
    {
        // Mixture-backend (rho,u) inversion: the same flash the pure path uses.
        const double rho=std::atof(argv[2]),u=std::atof(argv[3]),y0=std::atof(argv[4]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        try
        {
            s->unspecify_phase();
            s->update(CoolProp::DmassUmass_INPUTS,rho,u);
            std::printf("DU success=1 p=%.17g T=%.17g rho=%.17g u=%.17g s=%.17g Q=%.17g phase=%d\n",
                s->p(),s->T(),s->rhomass(),s->umass(),s->smass(),s->Q(),int(s->phase()));
            return 0;
        }
        catch(const std::exception& e){std::printf("DU success=0 error=%s\n",e.what());return 1;}
    }
    if(argc==5&&std::string(argv[1])=="--ptraw")
    {
        const double p=std::atof(argv[2]),t=std::atof(argv[3]),y0=std::atof(argv[4]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        try
        {
            s->unspecify_phase();s->update(CoolProp::PT_INPUTS,p,t);
            std::printf("PTRAW success=1 p=%.17g T=%.17g rho=%.17g u=%.17g s=%.17g Q=%.17g phase=%d\n",
                s->p(),s->T(),s->rhomass(),s->umass(),s->smass(),s->Q(),int(s->phase()));
            return 0;
        }
        catch(const std::exception& e){std::printf("PTRAW success=0 error=%s\n",e.what());return 1;}
    }
    if(argc==7&&std::string(argv[1])=="--rhos")
    {
        const double rho=std::atof(argv[2]),target=std::atof(argv[3]);
        const double p0=std::atof(argv[4]),t0=std::atof(argv[5]),y0=std::atof(argv[6]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        double p=p0,t=t0;bool converged=false;
        const double secondScale=std::max(std::abs(target),1.);
        for(int iter=0;iter<35;++iter)
        {
            const auto q=Foam::thermo::cpdetail::pt(*s,p,t);
            const double dr=rho-q.rho,dv=target-q.s;
            std::printf("  it=%d p=%.17g T=%.17g rho=%.17g s=%.17g dr=%.6g dv=%.6g\n",
                iter,p,t,q.rho,q.s,dr,dv);
            if(std::abs(dr)<=2e-11*std::max(rho,1.)&&std::abs(dv)<=2e-11*secondScale){converged=true;break;}
            const auto d=Foam::thermo::cpdetail::derivatives(*s,q);
            const double det=d.rP*d.sT-d.rT*d.sP;
            if(!std::isfinite(det)||std::abs(det)<1e-20){std::printf("  det=%.6g ABORT\n",det);break;}
            const double dp=(dr*d.sT-d.rT*dv)/det,dt=(d.rP*dv-dr*d.sP)/det;
            const double norm=std::abs(dr)/rho+std::abs(dv)/secondScale;
            bool advanced=false;
            for(double scale=1.;scale>=1./65536.;scale*=.5)
            {
                const double pn=p+scale*dp,tn=t+scale*dt;
                if(pn<=1.||tn<60.||tn>2000.)continue;
                try
                {
                    const auto trial=Foam::thermo::cpdetail::pt(*s,pn,tn);
                    const double err=std::abs(rho-trial.rho)/rho+std::abs(target-trial.s)/secondScale;
                    if(err<norm){p=pn;t=tn;advanced=true;break;}
                }catch(const std::exception&){}
            }
            if(!advanced){std::printf("  NO_ADVANCE\n");break;}
        }
        std::printf("RHOS_RESULT converged=%d p=%.17g T=%.17g\n",int(converged),p,t);
        return converged?0:1;
    }
    if(argc==5&&std::string(argv[1])=="--identity")
    {
        const double p=std::atof(argv[2]),t=std::atof(argv[3]),y0=std::atof(argv[4]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        try
        {
            s->unspecify_phase();s->update(CoolProp::PT_INPUTS,p,t);
            const auto base=Foam::thermo::cpdetail::current(*s);
            std::printf("IDENT base p=%.17g T=%.17g rho=%.17g s=%.17g twoPhase=%d\n",base.p,base.t,base.rho,base.s,int(base.twoPhase));
            for(double eps:{1e-3,5e-4,2.5e-4,1.25e-4,6.25e-5})
            {
                try
                {
                    const double a2=Foam::thermo::cpdetail::soundSpeedIdentitySquared(*s,base,eps);
                    std::printf("  eps=%.6g a2=%.9g a=%.6g\n",eps,a2,std::sqrt(a2));
                }
                catch(const std::exception& e){std::printf("  eps=%.6g threw %s\n",eps,e.what());}
            }
        }
        catch(const std::exception& e){std::printf("IDENT failed %s\n",e.what());}
        return 0;
    }
    if(argc==5&&std::string(argv[1])=="--slope")
    {
        const double p=std::atof(argv[2]),t=std::atof(argv[3]),y0=std::atof(argv[4]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        try
        {
            const double slope=Foam::thermo::cpdetail::entropySlopeT(*s,p,t);
            std::printf("SLOPE success=1 p=%.17g T=%.17g ds_dT=%.17g\n",p,t,slope);return 0;
        }
        catch(const std::exception& e){std::printf("SLOPE success=0 error=%s\n",e.what());return 1;}
    }
    if(argc==6&&std::string(argv[1])=="--ptseq")
    {
        const double p1=std::atof(argv[2]),p2=std::atof(argv[3]);
        const double t=std::atof(argv[4]),y0=std::atof(argv[5]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        for(double p:{p1,p2})
        {
            try
            {
                s->unspecify_phase();s->update(CoolProp::PT_INPUTS,p,t);
                std::printf("PTSEQ p=%.17g rho=%.17g s=%.17g Q=%.17g phase=%d\n",
                    p,s->rhomass(),s->smass(),s->Q(),int(s->phase()));
            }
            catch(const std::exception& e){std::printf("PTSEQ p=%.17g error=%s\n",p,e.what());}
        }
        return 0;
    }
    if(argc==5&&std::string(argv[1])=="--isent")
    {
        const double p=std::atof(argv[2]),t=std::atof(argv[3]),y0=std::atof(argv[4]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        try
        {
            s->unspecify_phase();s->update(CoolProp::PT_INPUTS,p,t);
            const auto base=Foam::thermo::cpdetail::current(*s);
            std::printf("ISENT base p=%.17g T=%.17g rho=%.17g s=%.17g twoPhase=%d\n",
                base.p,base.t,base.rho,base.s,int(base.twoPhase));
            for(double eps:{1e-3,5e-4,2.5e-4,1.25e-4,6.25e-5})
            {
                const double pLo=p*(1.-eps),pHi=p*(1.+eps);
                const auto q0=Foam::thermo::cpdetail::isentropePoint(*s,pLo,base.s,base.t);
                const auto q1=Foam::thermo::cpdetail::isentropePoint(*s,pHi,base.s,base.t);
                const double a2=(q0.twoPhase&&q1.twoPhase)?(pHi-pLo)/(q1.rho-q0.rho)
                    :q0.twoPhase?(p-pLo)/(q0.rho-base.rho)
                    :(pHi-p)/(q1.rho-base.rho);
                std::printf("  eps=%.6g q0(T=%.6g rho=%.6g s=%.6g tp=%d) q1(T=%.6g rho=%.6g s=%.6g tp=%d) a2=%.6g\n",
                    eps,q0.t,q0.rho,q0.s,int(q0.twoPhase),q1.t,q1.rho,q1.s,int(q1.twoPhase),a2);
            }
        }
        catch(const std::exception& e){std::printf("ISENT failed: %s\n",e.what());}
        return 0;
    }
    if(argc==7&&std::string(argv[1])=="--rhoe")
    {
        const double rho=std::atof(argv[2]),e=std::atof(argv[3]),y0=std::atof(argv[4]);
        const double p0=std::atof(argv[5]),t0=std::atof(argv[6]);
        auto model=Foam::thermo::makeCoolPropThermoModel();
        Foam::thermo::EquilibriumSeed seed;
        seed.pressure=p0;seed.temperature=t0;
        const auto r=model->recoverRhoE(rho,e,{y0,1.-y0},&seed);
        std::printf("RHOE rho=%.17g e=%.17g yO2=%.17g seed(p=%.17g,T=%.17g) success=%d status=%s\n",
            rho,e,y0,p0,t0,int(r.success()),r.diagnostics.backendStatus.c_str());
        if(r.success())
            std::printf("  -> p=%.10g T=%.10g rho=%.10g e=%.10g region=%d\n",
                r.state->p,r.state->T,r.state->rho,r.state->e,int(r.state->region));
        return r.success()?0:1;
    }
    if(argc==7&&std::string(argv[1])=="--rhoemin")
    {
        const double rho=std::atof(argv[2]),e=std::atof(argv[3]),y0=std::atof(argv[4]);
        const double p0=std::atof(argv[5]);
        const double t0=std::atof(argv[6]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        double href[2];
        for(int i=0;i<2;++i)
        {
            CPState* r=CPState::factory("HEOS",i==0?"Oxygen":"Propane");
            r->update(CoolProp::DmassT_INPUTS,1e-10,298.15);href[i]=r->hmass_idealgas();delete r;
        }
        const double offset=y0*href[0]+(1.-y0)*href[1];
        const double target=e+offset;
        double best=1e300,bp=0,bt=0;int count=0;
        for(int ip=0;ip<=40;++ip)
            for(int it=0;it<=40;++it)
            {
                const double pv=p0*(0.7+0.015*ip);
                const double tv=t0*(0.7+0.015*it);
                try
                {
                    const auto q=Foam::thermo::cpdetail::pt(*s,pv,tv);
                    const double norm=std::abs(q.rho-rho)/rho+std::abs(q.e-target)/2.2e5;
                    ++count;
                    if(norm<best){best=norm;bp=pv;bt=tv;}
                }
                catch(const std::exception&){}
            }
        std::printf("RHOEMIN samples=%d bestNorm=%.6g bestP=%.10g bestT=%.10g twoPhaseAtBest=%d\n",
            count,best,bp,bt,0);
        try
        {
            const auto q=Foam::thermo::cpdetail::pt(*s,bp,bt);
            std::printf("RHOEMIN_STATE rho=%.12g e_cp=%.12g target_cp=%.12g dr=%.6g dv=%.6g tp=%d\n",
                q.rho,q.e,target,q.rho-rho,q.e-target,int(q.twoPhase));
        }
        catch(const std::exception& ex){std::printf("RHOEMIN_STATE failed %s\n",ex.what());}
        return 0;
    }
    if(argc==11&&std::string(argv[1])=="--rhofine")
    {
        const double rho=std::atof(argv[2]),e=std::atof(argv[3]),y0=std::atof(argv[4]);
        const double pmin=std::atof(argv[5]),pmax=std::atof(argv[6]);
        const double tmin=std::atof(argv[7]),tmax=std::atof(argv[8]);
        const int nP=std::atoi(argv[9]),nT=std::atoi(argv[10]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        double href[2];
        for(int i=0;i<2;++i)
        {
            CPState* r=CPState::factory("HEOS",i==0?"Oxygen":"Propane");
            r->update(CoolProp::DmassT_INPUTS,1e-10,298.15);href[i]=r->hmass_idealgas();delete r;
        }
        const double offset=y0*href[0]+(1.-y0)*href[1];
        const double target=e+offset;
        struct Entry{double norm,p,t,rho,ee;int q;};
        std::vector<Entry> best;
        int fail=0,ok=0;
        for(int ip=0;ip<nP;++ip)
            for(int it=0;it<nT;++it)
            {
                const double pv=pmin+(pmax-pmin)*ip/std::max(1,nP-1);
                const double tv=tmin+(tmax-tmin)*it/std::max(1,nT-1);
                try
                {
                    const auto q=Foam::thermo::cpdetail::pt(*s,pv,tv);
                    ++ok;
                    const double norm=std::abs(q.rho-rho)/rho+std::abs(q.e-target)/std::max(std::abs(e),1.);
                    best.push_back({norm,pv,tv,q.rho,q.e,int(q.twoPhase)});
                }
                catch(const std::exception&){++fail;}
            }
        std::sort(best.begin(),best.end(),[](const Entry&a,const Entry&b){return a.norm<b.norm;});
        std::printf("RHO_FINE ok=%d fail=%d best=%zu\n",ok,fail,best.size());
        for(std::size_t i=0;i<best.size()&&i<5;++i)
            std::printf("  #%zu norm=%.6g p=%.8g T=%.8g rho=%.10g e=%.10g tp=%d\n",
                i+1,best[i].norm,best[i].p,best[i].t,best[i].rho,best[i].ee,best[i].q);
        return 0;
    }
    if(argc==7&&std::string(argv[1])=="--rhoev")
    {
        const double rho=std::atof(argv[2]),e=std::atof(argv[3]),y0=std::atof(argv[4]);
        const double p0=std::atof(argv[5]),t0=std::atof(argv[6]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        // energy offset exactly as the model computes it
        double href[2];
        for(int i=0;i<2;++i)
        {
            CPState* r=CPState::factory("HEOS",i==0?"Oxygen":"Propane");
            r->update(CoolProp::DmassT_INPUTS,1e-10,298.15);href[i]=r->hmass_idealgas();delete r;
        }
        const double offset=y0*href[0]+(1.-y0)*href[1];
        const double target=e+offset;
        double p=p0,t=t0;
        for(int iter=0;iter<50;++iter)
        {
            try
            {
                const auto q=Foam::thermo::cpdetail::pt(*s,p,t);
                const double dr=rho-q.rho,dv=target-q.e;
                std::printf("  it=%d p=%.17g T=%.17g rho=%.17g e=%.17g dr=%.6g dv=%.6g\n",
                    iter,p,t,q.rho,q.e,dr,dv);
                if(std::abs(dr)<=2e-11*std::max(rho,1.)&&std::abs(dv)<=2e-11*2.2e5){std::printf("CONVERGED\n");return 0;}
                const auto d=Foam::thermo::cpdetail::derivatives(*s,q,1e-5);
                const double det=d.rP*d.eT-d.rT*d.eP;
                std::printf("    rP=%.6g rT=%.6g eP=%.6g eT=%.6g det=%.6g\n",d.rP,d.rT,d.eP,d.eT,det);
                if(!std::isfinite(det)||std::abs(det)<1e-20){std::printf("  det fail\n");return 1;}
                const double dp=(dr*d.eT-d.rT*dv)/det,dt=(d.rP*dv-dr*d.eP)/det;
                std::printf("    dp=%.6g dt=%.6g\n",dp,dt);
                p+=dp;t+=dt;
                if(p<=1.||t<60.||t>2000.){std::printf("  out of range\n");return 1;}
            }
            catch(const std::exception& ex){std::printf("  EXCEPTION %s\n",ex.what());return 1;}
        }
        return 1;
    }
    if(argc==6&&std::string(argv[1])=="--pssolvev")
    {
        const double p=std::atof(argv[2]),target=std::atof(argv[3]);
        const double t0=std::atof(argv[4]),y0=std::atof(argv[5]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        double t=t0;
        for(int iter=0;iter<40;++iter)
        {
            try
            {
                const auto q=Foam::thermo::cpdetail::pt(*s,p,t);
                const double error=q.s-target;
                std::printf("  it=%d T=%.17g rho=%.17g s=%.17g twoPhase=%d error=%.6g\n",
                    iter,t,q.rho,q.s,int(q.twoPhase),error);
                if(std::abs(error)<=1e-10*std::max(std::abs(target),1.)){std::printf("CONVERGED\n");return 0;}
                const double slope=Foam::thermo::cpdetail::entropySlopeT(*s,p,t);
                const double step=std::clamp(-error/slope,-.1*t,.1*t);
                bool advanced=false;
                for(double scale=1.;scale>=1./1024.;scale*=.5)
                {
                    const double tn=t+scale*step;
                    if(tn<60.||tn>2000.)continue;
                    try
                    {
                        const auto trial=Foam::thermo::cpdetail::pt(*s,p,tn);
                        if(std::abs(trial.s-target)<std::abs(error)){t=tn;advanced=true;break;}
                    }
                    catch(const std::exception&){}
                }
                if(!advanced){std::printf("  NO_ADVANCE slope=%.6g step=%.6g\n",slope,step);return 1;}
            }
            catch(const std::exception& e){std::printf("  EXCEPTION %s\n",e.what());return 1;}
        }
        return 1;
    }
    if(argc==6&&std::string(argv[1])=="--pssolve")
    {
        const double p=std::atof(argv[2]),target=std::atof(argv[3]);
        const double t0=std::atof(argv[4]),y0=std::atof(argv[5]);
        auto s=CPState::factory("HEOS","Oxygen&Propane");
        s->set_mass_fractions({y0,1.-y0});
        try
        {
            const auto q=Foam::thermo::cpdetail::isentropePoint(*s,p,target,t0);
            std::printf("PSSOLVE success=1 p=%.17g T=%.17g rho=%.17g s=%.17g twoPhase=%d\n",
                q.p,q.t,q.rho,q.s,int(q.twoPhase));return 0;
        }
        catch(const std::exception& e){std::printf("PSSOLVE success=0 error=%s\n",e.what());return 1;}
    }
    if(argc==7&&std::string(argv[1])=="--recover2")
    {
        const double rho=std::atof(argv[2]),e=std::atof(argv[3]),y0=std::atof(argv[4]);
        const double sp=std::atof(argv[5]),st=std::atof(argv[6]);
        auto model=Foam::thermo::makeCoolPropThermoModel();
        Foam::thermo::EquilibriumSeed seed;
        seed.pressure=sp;seed.temperature=st;seed.previousTemperature=st;
        seed.previousRegion=Foam::thermo::EquilibriumRegion::TWO_PHASE;
        const auto r=model->recoverRhoE(rho,e,{y0,1.-y0},&seed);
        std::printf("RECOVER2 success=%d status=%s backend=%s\n",int(r.success()),
            Foam::thermo::equilibriumStatusName(r.status),r.diagnostics.backendStatus.c_str());
        if(r.state)std::printf("RECOVER2_STATE p=%.17g T=%.17g rho=%.17g e=%.17g region=%d\n",
            r.state->p,r.state->T,r.state->rho,r.state->e,int(r.state->region));
        return r.success()?0:1;
    }
    // Defaults are the recorded trial state and seed from the failure snapshot.
    double rho=44.31367393683184,e=29291.760254948,y0=0.71284999
          ,sp=4302281.193782107,st=411.25762958063888;
    if(argc==7){rho=std::atof(argv[1]);e=std::atof(argv[2]);y0=std::atof(argv[3]);
                sp=std::atof(argv[4]);st=std::atof(argv[5]);if(std::string(argv[6])!="run")return 2;}
    else if(argc!=1){std::fprintf(stderr,"usage: probe_rhoe_recovery [rho e yO2 seed_p seed_T run]\n");return 2;}
    std::vector<double> y{y0,1.-y0};

    std::printf("PROBE_INPUT rho=%.17g e_rfam=%.17g Y_O2=%.17g seed_p=%.17g seed_T=%.17g\n",rho,e,y0,sp,st);
    auto model=Foam::thermo::makeCoolPropThermoModel();
    std::printf("BACKEND %s\n",model->modelId().c_str());
    Foam::thermo::EquilibriumSeed seed;
    seed.pressure=sp;seed.temperature=st;seed.previousTemperature=st;
    seed.previousRegion=Foam::thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE;
    const auto result=model->recoverRhoE(rho,e,y,&seed);
    std::printf("MODEL_RESULT success=%d status=%s backend=%s\n",int(result.success()),
        Foam::thermo::equilibriumStatusName(result.status),result.diagnostics.backendStatus.c_str());
    if(result.state)
    {
        std::printf("MODEL_STATE p=%.17g T=%.17g rho=%.17g e=%.17g\n",
            result.state->p,result.state->T,result.state->rho,result.state->e);
        std::printf("MODEL_RESIDUALS density=%.6g energy=%.6g\n",
            result.diagnostics.densityResidual,result.diagnostics.energyResidual);
    }

    // Reference energy offset exactly as the model computes it.
    double href[2];CPState* ref[2];
    for(int i=0;i<2;++i)
    {
        ref[i]=CPState::factory("HEOS",i==0?"Oxygen":"Propane");
        ref[i]->update(CoolProp::DmassT_INPUTS,1e-10,298.15);
        href[i]=ref[i]->hmass_idealgas();delete ref[i];
    }
    const double offset=y[0]*href[0]+y[1]*href[1];
    const double target=e+offset;
    std::printf("ENERGY_OFFSET %.17g target_u=%.17g\n",offset,target);

    // Native CoolProp mixture rho-u flash (patched density bracket) as an
    // independent comparison to the production recovery path.
    {
        CPState* native=CPState::factory("HEOS","Oxygen&Propane");
        native->set_mass_fractions(y);
        try
        {
            native->update(CoolProp::DmassUmass_INPUTS,rho,target);
            std::printf("NATIVE_RHOE success=1 p=%.17g T=%.17g rho=%.17g u=%.17g phase=%d\n",
                native->p(),native->T(),native->rhomass(),native->umass(),int(native->phase()));
            std::printf("NATIVE_RESIDUALS density=%.6g energy=%.6g\n",
                native->rhomass()-rho,native->umass()-target);
        }
        catch(const std::exception& e)
        {
            std::printf("NATIVE_RHOE success=0 error=%s\n",e.what());
        }
        delete native;
    }

    CPState* s=CPState::factory("HEOS","Oxygen&Propane");
    s->set_mass_fractions(y);

    // Diagnostic copy of the production fast path (specified-gas DmassT loop).
    std::printf("FAST_PATH_BEGIN\n");
    try
    {
    s->specify_phase(CoolProp::iphase_gas);
    double trialT=st;
    for(int iter=0;iter<12;++iter)
    {
        s->update(CoolProp::DmassT_INPUTS,rho,trialT);
        const double residual=s->umass()-target;
        const double cv=s->cvmass();
        std::printf("  it=%d T=%.17g residual=%.6g cv=%.17g\n",iter,trialT,residual,cv);
        if(std::abs(residual)<=1e-10*std::max(std::abs(e),1.))break;
        if(!std::isfinite(cv)||cv<=0.)break;
        trialT+=std::clamp(-residual/cv,-.15*trialT,.15*trialT);
        if(trialT<60.||trialT>2000.)break;
    }
    s->update(CoolProp::DmassT_INPUTS,rho,trialT);
    const double candidateP=s->p();
    s->unspecify_phase();
    s->update(CoolProp::PT_INPUTS,candidateP,trialT);
    std::printf("FAST_PATH_VERIFY phase=%d d_rho=%.6g d_u=%.6g\n",int(s->phase()),
        s->rhomass()-rho,s->umass()-target);
    }
    catch(const std::exception& e)
    {
        std::printf("FAST_PATH_VERIFY threw=%s\n",e.what());
    }
    std::fflush(stdout);

    // Diagnostic copy of cpdetail::recover (same equations/tolerances).
    std::printf("NEWTON_BEGIN\n");
    double p=sp,t=st;
    const double secondScale=std::max(std::abs(target),1.);
    bool converged=false;
    for(int iter=0;iter<35;++iter)
    {
        const auto q=Foam::thermo::cpdetail::pt(*s,p,t);
        const double dr=rho-q.rho,dv=target-q.e;
        std::printf("  it=%d p=%.17g T=%.17g dr=%.6g dv=%.6g\n",iter,p,t,dr,dv);
        if(std::abs(dr)<=2e-11*std::max(rho,1.)&&std::abs(dv)<=2e-11*secondScale)
        {converged=true;break;}
        const auto d=Foam::thermo::cpdetail::derivatives(*s,q);
        const double det=d.rP*d.eT-d.rT*d.eP;
        if(!std::isfinite(det)||std::abs(det)<1e-20){std::printf("  det=%.6g ABORT\n",det);break;}
        const double dp=(dr*d.eT-d.rT*dv)/det,dt=(d.rP*dv-dr*d.eP)/det;
        const double norm=std::abs(dr)/rho+std::abs(dv)/secondScale;
        bool advanced=false;
        for(double scale=1.;scale>=1./65536.;scale*=.5)
        {
            const double pn=p+scale*dp,tn=t+scale*dt;
            if(pn<=1.||tn<60.||tn>2000.)continue;
            try
            {
                const auto trial=Foam::thermo::cpdetail::pt(*s,pn,tn);
                const double err=std::abs(rho-trial.rho)/rho+std::abs(target-trial.e)/secondScale;
                if(err<norm){p=pn;t=tn;advanced=true;std::printf("  accept scale=%.6g\n",scale);break;}
            }catch(const std::exception&){}
        }
        if(!advanced){std::printf("  NO_ADVANCE dp=%.6g dt=%.6g\n",dp,dt);break;}
    }
    std::printf("NEWTON_RESULT converged=%d final_p=%.17g final_T=%.17g\n",int(converged),p,t);
    printPoint("NEWTON_EOS",*s);

    // Diagnostic: bracketed equilibrium solve at fixed density (no imposed
    // phase) on the same state. This is a candidate replacement for the
    // gas-specified fast path when the target lies in the two-phase dome.
    try
    {
        std::fflush(stdout);
        CPState* b=CPState::factory("HEOS","Oxygen&Propane");
        b->set_mass_fractions(y);
        std::vector<std::pair<double,double>> samples;
        int usable=0;
        for(int i=0;i<=54;++i)
        {
            const double T=60.+10.*i;
            try
            {
                b->update(CoolProp::DmassT_INPUTS,rho,T);
                samples.push_back({T,b->umass()-target});++usable;
            }
            catch(const std::exception&){samples.push_back({T,std::nan("")});}
        }
        double lo=std::nan(""),hi=std::nan("");
        for(std::size_t i=1;i<samples.size();++i)
            if(std::isfinite(samples[i].second)&&std::isfinite(samples[i-1].second)
               &&samples[i-1].second*samples[i].second<=0.)
            {lo=samples[i-1].first;hi=samples[i].first;break;}
        if(std::isfinite(lo))
        {
            b->update(CoolProp::DmassT_INPUTS,rho,lo);
            const double flo=b->umass()-target;
            for(int i=0;i<80;++i)
            {
                const double mid=0.5*(lo+hi);
                b->update(CoolProp::DmassT_INPUTS,rho,mid);
                if((b->umass()-target)*flo>0.)lo=mid;else hi=mid;
            }
            const double T=0.5*(lo+hi);
            b->update(CoolProp::DmassT_INPUTS,rho,T);
            const double candidateP=b->p();const double candidateU=b->umass();
            b->unspecify_phase();b->update(CoolProp::PT_INPUTS,candidateP,T);
            std::printf("BRACKETED_DT success=1 T=%.17g p=%.17g phase=%d d_rho=%.6g d_u=%.6g\n",
                T,candidateP,int(b->phase()),b->rhomass()-rho,b->umass()-target);
        }
        else std::printf("BRACKETED_DT success=0 usable_samples=%d\n",usable);
        delete b;
    }
    catch(const std::exception& e)
    {
        std::printf("BRACKETED_DT success=0 error=%s\n",e.what());
    }
    std::fflush(stdout);
    delete s;
    return converged?0:1;
}
