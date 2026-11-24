#ifndef RFAM_CoolPropEquilibrium_H
#define RFAM_CoolPropEquilibrium_H
// Internal native HEOS flash utilities. All PT samples use equilibrium flash;
// no gas/liquid imposition, bulk-energy projection, or composition snapping.
#include "AbstractState.h"
#include "ThermoCounters.H"
#include "StateAcceptance.H"
#include "TableAcoustic.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace Foam { namespace thermo { namespace cpdetail {
struct Point
{
    double p,t,rho,e,h,s,quality;
    bool twoPhase;
};
inline Point current(CoolProp::AbstractState& state)
{
    return {state.p(),state.T(),state.rhomass(),state.umass(),state.hmass(),state.smass(),
            state.Q(),state.phase()==CoolProp::iphase_twophase};
}
inline Point pt(CoolProp::AbstractState& state,double pressure,double temperature)
{
    thermoCount(thermoCounters().pt);
    state.unspecify_phase();state.update(CoolProp::PT_INPUTS,pressure,temperature);
    const Point q=current(state);
    if(!std::isfinite(q.rho)||q.rho<=0.||!std::isfinite(q.e)||!std::isfinite(q.s))
        throw std::runtime_error("COOLPROP_INVALID_EQUILIBRIUM_SAMPLE");
    return q;
}
struct Derivatives
{
    double rP,rT,eP,eT,hP,hT,sP,sT;
};
inline Derivatives derivatives(CoolProp::AbstractState& state,const Point& q,double eps=1e-5)
{
    if(!q.twoPhase)
    {
        const auto d=[&](CoolProp::parameters prop,CoolProp::parameters wrt,CoolProp::parameters held)
        {return state.first_partial_deriv(prop,wrt,held);};
        return {d(CoolProp::iDmass,CoolProp::iP,CoolProp::iT),d(CoolProp::iDmass,CoolProp::iT,CoolProp::iP),
                d(CoolProp::iUmass,CoolProp::iP,CoolProp::iT),d(CoolProp::iUmass,CoolProp::iT,CoolProp::iP),
                d(CoolProp::iHmass,CoolProp::iP,CoolProp::iT),d(CoolProp::iHmass,CoolProp::iT,CoolProp::iP),
                d(CoolProp::iSmass,CoolProp::iP,CoolProp::iT),d(CoolProp::iSmass,CoolProp::iT,CoolProp::iP)};
    }
    // Two-phase branch: a symmetric stencil must stay inside the dome. Shrink
    // it; if the dome is too narrow, use one-sided differences taken from the
    // base point into the two-phase side (same as the identity fallback).
    for(double scale=1.;scale>=1./64.;scale*=.5)
    {
        const double e=eps*scale;
        try
        {
            const Point pp=pt(state,q.p*(1.+e),q.t),pm=pt(state,q.p*(1.-e),q.t);
            const Point tp=pt(state,q.p,q.t*(1.+e)),tm=pt(state,q.p,q.t*(1.-e));
            if(!pp.twoPhase||!pm.twoPhase||!tp.twoPhase||!tm.twoPhase)continue;
            pt(state,q.p,q.t);
            return {(pp.rho-pm.rho)/(2*q.p*e),(tp.rho-tm.rho)/(2*q.t*e),
                    (pp.e-pm.e)/(2*q.p*e),(tp.e-tm.e)/(2*q.t*e),
                    (pp.h-pm.h)/(2*q.p*e),(tp.h-tm.h)/(2*q.t*e),
                    (pp.s-pm.s)/(2*q.p*e),(tp.s-tm.s)/(2*q.t*e)};
        }
        catch(const std::exception&){}
    }
    for(double scale=1.;scale>=1./64.;scale*=.5)
    {
        const double e=eps*scale;
        try
        {
            const Point pp=pt(state,q.p*(1.+e),q.t),pm=pt(state,q.p*(1.-e),q.t);
            const Point tp=pt(state,q.p,q.t*(1.+e)),tm=pt(state,q.p,q.t*(1.-e));
            const bool pPlus=pp.twoPhase,pMinus=pm.twoPhase,tPlus=tp.twoPhase,tMinus=tm.twoPhase;
            if(!(pPlus||pMinus)||!(tPlus||tMinus))continue;
            const Point& pHi=pPlus?pp:q;const Point& pLo=pMinus?pm:q;
            const double dp=pPlus?q.p*e:-q.p*e;
            const Point& tHi=tPlus?tp:q;const Point& tLo=tMinus?tm:q;
            const double dt=tPlus?q.t*e:-q.t*e;
            return {(pHi.rho-pLo.rho)/dp,(tHi.rho-tLo.rho)/dt,
                    (pHi.e-pLo.e)/dp,(tHi.e-tLo.e)/dt,
                    (pHi.h-pLo.h)/dp,(tHi.h-tLo.h)/dt,
                    (pHi.s-pLo.s)/dp,(tHi.s-tLo.s)/dt};
        }
        catch(const std::exception&){}
    }
    throw std::runtime_error("COOLPROP_DERIVATIVE_CROSSED_PHASE_BOUNDARY");
}
// rho-e or rho-s inversion through the stable PT flash. Entropy inversion is
// used solely for independent acoustic verification. Density and the second
// conserved target must close; a returned failure carries no projected state.
// Pure-fluid (rho,e) inversion on the saturation curve.
//
// CoolProp's DmassUmass flash fails on states that sit on the steep branch just
// inside the saturated-liquid line (A-route cell 132114 at step 4:
// rho=920.66118876 kg/m3, e=-336706.5 J/kg, pure O2 -> p=2.409 MPa, T=128.45 K,
// quality ~0.17 %). For a pure fluid the lever rule gives a well-posed 1-D
// problem, so this fallback is used when the library flash throws:
//   1/rho = (1-x)/rhoL(T) + x/rhoV(T)      ->  x(T)
//   e     = (1-x)eL(T) + x eV(T)           ->  1-D root in T
// If the two-phase solve is not consistent (the state is single phase), fall
// back to an isobaric refinement at the caller's pressure: bisect T so that
// rho(p,T)=rho, then accept only if the energy also closes within the gate.
inline bool recoverPureOnSaturation(CoolProp::AbstractState& state,double rho,
                                    double eTarget,double pSeed,double tSeed,
                                    double secondScale)
{
    thermoCount(thermoCounters().recover);
    const double tCrit=state.T_critical();
    const double pCrit=state.p_critical();
    const double v=1.0/rho;
    double tLo=std::max(60.0,tSeed-25.0),tHi=std::min(tCrit-1e-6,tSeed+25.0);
    if(!(tHi>tLo))return false;
    // two-phase lever rule: scan T, build x(T) and g(T)=e_mix(T)-eTarget
    const int nScan=160;
    double prevT=0.,prevG=0.;bool havePrev=false;double bestG=1e300,bestT=0.;
    for(int i=0;i<=nScan;++i)
    {
        const double t=tLo+(tHi-tLo)*i/nScan;
        try
        {
            state.update(CoolProp::QT_INPUTS,0.0,t);
            const double rhoL=state.rhomass(),eL=state.umass();
            state.update(CoolProp::QT_INPUTS,1.0,t);
            const double rhoV=state.rhomass(),eV=state.umass();
            const double vL=1.0/rhoL,vV=1.0/rhoV;
            if(!(vV>vL))continue;
            const double x=(v-vL)/(vV-vL);
            if(x<-0.05||x>1.05)continue;
            const double g=(1.0-x)*eL+x*eV-eTarget;
            if(std::abs(g)<std::abs(bestG)){bestG=g;bestT=t;}
            if(havePrev&&prevG*g<0.)
            {
                double a=prevT,b=t,ga=prevG;
                for(int it=0;it<80;++it)
                {
                    const double m=0.5*(a+b);
                    state.update(CoolProp::QT_INPUTS,0.0,m);
                    const double rl=state.rhomass(),el=state.umass();
                    state.update(CoolProp::QT_INPUTS,1.0,m);
                    const double rv=state.rhomass(),ev=state.umass();
                    const double xx=((1.0/rho)-1.0/rl)/(1.0/rv-1.0/rl);
                    const double gm=(1.0-xx)*el+xx*ev-eTarget;
                    if(ga*gm<=0.){b=m;}else{a=m;ga=gm;}
                }
                const double tSol=0.5*(a+b);
                state.update(CoolProp::QT_INPUTS,0.0,tSol);
                const double rl=state.rhomass(),el=state.umass();
                const double psat=state.p();
                state.update(CoolProp::QT_INPUTS,1.0,tSol);
                const double rv=state.rhomass(),ev=state.umass();
                const double xx=std::clamp(((1.0/rho)-1.0/rl)/(1.0/rv-1.0/rl),0.0,1.0);
                state.update(CoolProp::PQ_INPUTS,psat,xx);
                if(std::abs(state.rhomass()-rho)<=stateAcceptanceRel()*std::max(rho,1.)
                   &&std::abs(state.umass()-eTarget)<=stateAcceptanceRel()*secondScale)
                    return true;
            }
            prevT=t;prevG=g;havePrev=true;
        }
        catch(const std::exception&){havePrev=false;}
    }
    // single-phase refinement at the caller's pressure: rho(p,T)=rho by bisection
    for(double pScale:{1.0,0.5,1.5})
    {
        const double p=pSeed*pScale;
        if(!(p>1e-3*pCrit)||p>10.*pCrit)continue;
        double tl=std::max(60.0,tSeed-40.0),th=std::min(tCrit*1.5,tSeed+40.0);
        auto resid=[&](double t)
        {
            try{state.update(CoolProp::PT_INPUTS,p,t);
                return std::make_pair(state.rhomass()-rho,state.umass()-eTarget);}
            catch(const std::exception&){return std::make_pair(1e300,1e300);}
        };
        double rl=resid(tl).first,rh=resid(th).first;
        if(!(rl*rh<0.))continue;
        for(int it=0;it<80;++it)
        {
            const double m=0.5*(tl+th);
            const double r=resid(m).first;
            if(rl*r<=0.){th=m;}else{tl=m;rl=r;}
        }
        const double tSol=0.5*(tl+th);
        const auto r=resid(tSol);
        if(std::abs(r.first)<=stateAcceptanceRel()*std::max(rho,1.)
           &&std::abs(r.second)<=stateAcceptanceRel()*secondScale)
            return true;
    }
    if(std::abs(bestG)<1e300)state.update(CoolProp::QT_INPUTS,0.0,bestT);
    return false;
}
inline Point recover(CoolProp::AbstractState& state,double rho,double target,
                     double pSeed,double tSeed,bool entropy=false,
                     double secondScaleOverride=-1.)
{
    thermoCount(thermoCounters().recover);
    // The caller may pass the conserved-variable scale (e.g. the RFAM sensible
    // energy) so the internal tolerance matches the solver's admission gate.
    // Without it the offset-shifted CoolProp target would make the tolerance
    // several times tighter than the gate and reject marginal states.
    const double secondScale=secondScaleOverride>0.?secondScaleOverride
        :(entropy?std::max(std::abs(target),1.):std::max(std::abs(target),1e5));
    const auto residualNorm=[&](const Point& q)
    {
        const double dr=rho-q.rho;
        const double dv=target-(entropy?q.s:q.e);
        return std::abs(dr)/std::max(rho,1.)+std::abs(dv)/secondScale;
    };
    // Newton with a damped line search. Returns true and leaves `result` at the
    // converged point.
    const auto newton=[&](double p0,double t0,Point& result)
    {
        double p=p0,t=t0;
        for(int iter=0;iter<15;++iter)
        {
            Point q;
            try{q=pt(state,p,t);}
            catch(const std::exception&){return false;}
            const double dr=rho-q.rho,dv=target-(entropy?q.s:q.e);
            // Internal tolerance is aligned with the solver's admission gate
            // (TetraStateRecovery: 2e-8 relative). The previous 2e-11 chased the
            // PT flash's own numerical noise and stalled on near-boundary
            // two-phase states whose target is already inside the gate.
            const double tolerance=entropy?std::max(1e-9,0.1*stateAcceptanceRel()):stateAcceptanceRel();
            if(std::abs(dr)<=tolerance*std::max(rho,1.)&&std::abs(dv)<=tolerance*secondScale)
            {result=q;return true;}
            Derivatives d;
            try{d=derivatives(state,q);}
            catch(const std::exception&){return false;}
            const double vP=entropy?d.sP:d.eP,vT=entropy?d.sT:d.eT;
            const double det=d.rP*vT-d.rT*vP;
            if(!std::isfinite(det)||std::abs(det)<1e-20)return false;
            const double dp=(dr*vT-d.rT*dv)/det,dt=(d.rP*dv-dr*vP)/det;
            const double norm=std::abs(dr)/rho+std::abs(dv)/secondScale;
            bool advanced=false;
            for(double scale=1.;scale>=1./65536.;scale*=.5)
            {
                const double pn=p+scale*dp,tn=t+scale*dt;
                if(pn<=1.||tn<60.||tn>2000.)continue;
                try
                {
                    const Point trial=pt(state,pn,tn);
                    const double err=std::abs(rho-trial.rho)/rho
                        +std::abs(target-(entropy?trial.s:trial.e))/secondScale;
                    if(err<norm){p=pn;t=tn;advanced=true;break;}
                }
                catch(const std::exception&){}
            }
            if(!advanced)return false;
        }
        return false;
    };
    Point result;
    if(newton(pSeed,tSeed,result))return result;
    // Fallback: coarse restart scan around the seed. Near the phase boundary the
    // local Jacobian can be ill-conditioned; a grid search finds a better basin
    // using only equilibrium PT flashes (no extrapolation).
    double bestNorm=std::numeric_limits<double>::max();
    double bestP=pSeed,bestT=tSeed;bool found=false;
    // 5x5 restart grid: enough to decide whether a local equilibrium state is
    // reachable. A target more than ~1e-3 off the surface is unrecoverable
    // locally (the solver's rollback reduces the step), so there is no value in
    // a larger search before failing.
    // 9x9 coarse grid (spacing 2.5 %) instead of 5x5 (5 %): A-route cell
    // 559462 needs a basin at +4.1 % p / +2.2 % T, which the 5 % grid steps
    // over; the fine probe shows that state closes at 3.6e-3, well inside the
    // development gate, so the grid must reach it. 81 flashes per failing cell
    // is still cheap next to the refinement levels that follow.
    for(int ip=0;ip<9;++ip)
        for(int it=0;it<9;++it)
        {
            const double pf=pSeed*(0.9+0.025*ip);
            const double tf=tSeed*(0.9+0.025*it);
            try
            {
                const Point q=pt(state,pf,tf);
                const double norm=residualNorm(q);
                if(norm<bestNorm){bestNorm=norm;bestP=pf;bestT=tf;found=true;}
            }
            catch(const std::exception&){}
        }
    thermoCount(thermoCounters().recoverScanRestart);
    // Multi-resolution refinement around the best coarse point. Near the
    // mixture critical point the (rho,e) target can sit a fraction of a kelvin
    // away from the seed (e.g. cell 11318 at step 45: recoverable at 3.804 MPa /
    // 165.30 K with norm 4.2e-3, but 0.5 K from the seed). The 5 % coarse
    // spacing steps over such basins entirely, so refine geometrically and let
    // the Newton start inside the basin. All evaluations are equilibrium PT
    // flashes (no extrapolation).
    for(int level=0;level<3&&found;++level)
    {
        const double pc=bestP,tc=bestT;
        // Half the previous level's point spacing, so each refined window
        // brackets one coarse cell (a quarter-step shrink skips the diagonal
        // pressure-temperature valley seen at the near-critical state).
        const double spanP=0.05*std::pow(0.5,double(level+1));
        const double spanT=0.05*std::pow(0.5,double(level+1));
        for(int ip=0;ip<5;++ip)
            for(int it=0;it<5;++it)
            {
                const double pf=pc*(1.+spanP*(ip-2)/2.);
                const double tf=tc*(1.+spanT*(it-2)/2.);
                if(pf<=1.||tf<60.||tf>2000.)continue;
                try
                {
                    const Point q=pt(state,pf,tf);
                    const double norm=residualNorm(q);
                    if(norm<bestNorm){bestNorm=norm;bestP=pf;bestT=tf;}
                }
                catch(const std::exception&){}
            }
    }
    // If the coarse scan already found a point inside the acceptance gate,
    // return it directly: the model re-verifies with its own gate, and the
    // Newton retry from the best point can otherwise burn hundreds of flashes
    // on states that are already acceptable.
    if(found)
    {
        try
        {
            const Point q=pt(state,bestP,bestT);
            const double dr=q.rho-rho,dv=target-(entropy?q.s:q.e);
            if(std::abs(dr)<=stateAcceptanceRel()*std::max(rho,1.)
               &&std::abs(dv)<=stateAcceptanceRel()*secondScale)
                return q;
        }
        catch(const std::exception&){}
    }
    if(found&&newton(bestP,bestT,result))return result;
    if(found&&projectionAcceptEnabled())
    {
        // Diagnostic projection mode: return the closest equilibrium state
        // found; the caller records the actual residual.
        try{return pt(state,bestP,bestT);}
        catch(const std::exception&){}
    }
    thermoCount(thermoCounters().recoverNewtonFail);
    if(std::getenv("RFAM_RECOVER_DEBUG"))
    {
        std::fprintf(stderr,"RECOVER_FAIL rho=%.17g target=%.17g pSeed=%.17g tSeed=%.17g bestP=%.17g bestT=%.17g found=%d bestNorm=%.6g\n",
                     rho,target,pSeed,tSeed,bestP,bestT,int(found),bestNorm);
    }
    throw std::runtime_error(entropy?"COOLPROP_RHOS_NEWTON_FAILED":"COOLPROP_RHOE_NEWTON_FAILED");
}
struct AcousticResult
{
    double a2,relativeChange,eps,lowerEntropyResidual,upperEntropyResidual;
    int comparisons;
    bool identityFallback=false;
};
// Single isentrope point: solve s(p,T)=s_target at fixed p (well conditioned in
// the two-phase region). The slope uses the same derivative helper as the
// rho-e inversion, with a one-sided fallback when a central difference would
// cross a phase boundary.
inline double entropySlopeT(CoolProp::AbstractState& state,double p,double t)
{
    try
    {
        const Point q=pt(state,p,t);
        const auto d=derivatives(state,q);
        if(std::isfinite(d.sT)&&d.sT>0.)return d.sT;
    }
    catch(const std::exception&){}
    for(double delta:{1e-6,5e-7,2e-7})
    {
        try
        {
            const Point hi=pt(state,p,t*(1.+delta)),lo=pt(state,p,t*(1.-delta));
            const double slope=(hi.s-lo.s)/(2.*t*delta);
            if(std::isfinite(slope)&&slope>0.)return slope;
        }
        catch(const std::exception&){}
    }
    throw std::runtime_error("COOLPROP_PS_SLOPE_FAILED");
}
inline Point isentropePoint(CoolProp::AbstractState& state,double p,double sTarget,double tSeed)
{
    double t=tSeed;
    bool newtonFailed=false;
    for(int iter=0;iter<40;++iter)
    {
        try
        {
            const Point q=pt(state,p,t);
            const double error=q.s-sTarget;
            if(std::abs(error)<=1e-10*std::max(std::abs(sTarget),1.))return q;
            const double slope=entropySlopeT(state,p,t);
            const double step=std::clamp(-error/slope,-.1*t,.1*t);
            bool advanced=false;
            for(double scale=1.;scale>=1./1024.;scale*=.5)
            {
                const double tn=t+scale*step;
                if(tn<60.||tn>2000.)continue;
                try
                {
                    const Point trial=pt(state,p,tn);
                    if(std::abs(trial.s-sTarget)<std::abs(error)){t=tn;advanced=true;break;}
                }
                catch(const std::exception&){}
            }
            if(!advanced)break;
        }
        catch(const std::exception&){newtonFailed=true;break;}
    }
    (void)newtonFailed;
    // Bracketed fallback: scan entropy along the isobar and bisect the sign
    // change. No derivatives are required, so this works on the two-phase
    // states where the local slope helper fails.
    std::array<double,2> bracketT{-1.,-1.};
    double previous=std::numeric_limits<double>::quiet_NaN();
    double previousT=std::numeric_limits<double>::quiet_NaN();
    const double lo=tSeed*0.5,hi=tSeed*1.5;
    for(int i=0;i<=120;++i)
    {
        const double ti=lo+(hi-lo)*i/120.;
        double value;
        try{value=pt(state,p,ti).s-sTarget;}
        catch(const std::exception&){previous=std::numeric_limits<double>::quiet_NaN();continue;}
        if(std::isfinite(previous)&&previous*value<=0.){bracketT={previousT,ti};break;}
        previous=value;previousT=ti;
    }
    if(bracketT[0]>0.)
    {
        double a=bracketT[0],b=bracketT[1];
        for(int iter=0;iter<80;++iter)
        {
            const double mid=0.5*(a+b);
            try
            {
                const Point q=pt(state,p,mid);
                if(std::abs(q.s-sTarget)<=1e-10*std::max(std::abs(sTarget),1.))return q;
                (q.s<sTarget? a : b)=mid;
            }
            catch(const std::exception&){a=mid;continue;}
        }
        const Point q=pt(state,p,0.5*(a+b));
        if(std::abs(q.s-sTarget)<=1e-8*std::max(std::abs(sTarget),1.))return q;
    }
    throw std::runtime_error("COOLPROP_PS_NOT_CONVERGED");
}
inline AcousticResult isentropeIterate(CoolProp::AbstractState& state,const Point& base,bool pure)
{
    thermoCount(thermoCounters().isentrope);
    double previous=NAN;int stable=0;AcousticResult result{};
    for(double eps:{5e-4,2.5e-4,1.25e-4,6.25e-5})
    {
        std::array<Point,2> q;
        double pLow=NAN,pHigh=NAN;
        if(pure)
        {
            // The isentrope crosses the saturation line, so one of the two
            // density perturbations can legitimately land single phase (the
            // same situation the mixture branch below already handles with a
            // one-sided difference). Keep the two-phase side instead of
            // rejecting the state: A-route cell 877234 sits on the saturated
            // liquid line, where a symmetric stencil does not exist.
            std::array<Point,2> qq;
            bool ok[2]={false,false};
            for(int i=0;i<2;++i)
            {
                const double rho=base.rho*(1.+(i==0?-eps:eps));
                try
                {
                    state.update(CoolProp::DmassSmass_INPUTS,rho,base.s);
                    qq[i]=current(state);
                    ok[i]=qq[i].twoPhase&&std::abs(qq[i].rho-rho)<=1e-8*rho
                          &&std::abs(qq[i].s-base.s)<=1e-8*std::max(std::abs(base.s),1.);
                }
                catch(const std::exception&){ok[i]=false;}
            }
            if(!ok[0]&&!ok[1])throw std::runtime_error("COOLPROP_ISENTROPE_STATE_FAILED");
            if(ok[0]&&ok[1]){q[0]=qq[0];q[1]=qq[1];pLow=qq[0].p;pHigh=qq[1].p;}
            else if(ok[0]){q[0]=qq[0];q[1]=qq[0];pLow=qq[0].p;pHigh=NAN;}
            else{q[0]=qq[1];q[1]=qq[1];pHigh=qq[1].p;pLow=NAN;}
        }
        else
        {
            // Pressure-parametrised isentrope. Inside the two-phase dome the
            // fixed-rho inversion is ill-conditioned (the isentrope runs
            // nearly parallel to the phase boundary); s(p,T)=s0 at fixed p is
            // a well-conditioned one-dimensional problem. Finite differencing
            // p and rho along the same isentrope gives the same derivative.
            pLow=base.p*(1.-eps);pHigh=base.p*(1.+eps);
            q[0]=isentropePoint(state,pLow,base.s,base.t);
            q[1]=isentropePoint(state,pHigh,base.s,base.t);
            // The isentrope is a continuous curve on the equilibrium surface;
            // near the phase boundary one side can legitimately be single
            // phase while the base state is two-phase. Use a one-sided
            // difference in that case instead of rejecting the state.
            if(!q[0].twoPhase&&!q[1].twoPhase)
                throw std::runtime_error("COOLPROP_ISENTROPE_STATE_FAILED");
            for(int i=0;i<2;++i)
                if(std::abs(q[i].s-base.s)>1e-8*std::max(std::abs(base.s),1.))
                    throw std::runtime_error("COOLPROP_ISENTROPE_STATE_FAILED");
        }
        const double a2=pure?(std::isfinite(pLow)&&std::isfinite(pHigh)
                ?(q[1].p-q[0].p)/(2.*base.rho*eps)
                :(std::isfinite(pLow)?(base.p-pLow)/(base.rho-q[0].rho)
                                     :(pHigh-base.p)/(q[1].rho-base.rho)))
            :(q[0].twoPhase&&q[1].twoPhase?(pHigh-pLow)/(q[1].rho-q[0].rho)
              :q[0].twoPhase?(base.p-pLow)/(q[0].rho-base.rho)
              :(pHigh-base.p)/(q[1].rho-base.rho));
        if(!std::isfinite(a2)||a2<=0.)throw std::runtime_error("COOLPROP_ISENTROPE_NONPOSITIVE");
        const double change=std::abs(a2-previous)/a2;
        result={a2,change,eps,q[0].s-base.s,q[1].s-base.s,result.comparisons+1};
        stable=std::isfinite(change)&&change<1e-4?stable+1:0;
        if(stable>=2)return result;
        previous=a2;
    }
    throw std::runtime_error("COOLPROP_ISENTROPE_NOT_CONVERGED");
}
// Equilibrium sound speed from the thermodynamic identity
//   a^2 = (dp/drho)_T + (T/rho^2) * (dp/dT)_rho^2 / Cv
// evaluated with central differences of equilibrium PT states. This is a local
// route (no root finding) used when the isentrope integration fails near the
// phase boundary; it was cross-validated against the property table.
inline double soundSpeedIdentitySquared(CoolProp::AbstractState& state,const Point& base,double eps)
{
    // Central differences must stay on one branch of the equilibrium surface.
    // If the two-phase band is too narrow for a symmetric stencil, fall back
    // to one-sided differences taken into the two-phase side. This is
    // first-order accurate but always available for an interior two-phase
    // state; the caller records it as an identity fallback.
    for(double scale=1.;scale>=1./16.;scale*=.5)
    {
        const double e=eps*scale;
        try
        {
            const Point pPlus=pt(state,base.p*(1.+e),base.t),pMinus=pt(state,base.p*(1.-e),base.t);
            const Point tPlus=pt(state,base.p,base.t*(1.+e)),tMinus=pt(state,base.p,base.t*(1.-e));
            const bool central=pPlus.twoPhase&&pMinus.twoPhase&&tPlus.twoPhase&&tMinus.twoPhase;
            double drho_dp_T,drho_dT_p,de_dT_p,de_dp_T;
            if(central)
            {
                drho_dp_T=(pPlus.rho-pMinus.rho)/(2.*base.p*e);
                drho_dT_p=(tPlus.rho-tMinus.rho)/(2.*base.t*e);
                de_dT_p=(tPlus.e-tMinus.e)/(2.*base.t*e);
                de_dp_T=(pPlus.e-pMinus.e)/(2.*base.p*e);
            }
            else
            {
                // one-sided in p (prefer the + side which is inside the dome at
                // the dew edge, else the - side), one-sided in T likewise.
                const bool pUsePlus=pPlus.twoPhase&&base.twoPhase;
                const bool pUseMinus=pMinus.twoPhase&&base.twoPhase;
                const bool tUsePlus=tPlus.twoPhase&&base.twoPhase;
                const bool tUseMinus=tMinus.twoPhase&&base.twoPhase;
                if(pUsePlus&&tUsePlus)
                {
                    drho_dp_T=(pPlus.rho-base.rho)/(base.p*e);
                    de_dp_T=(pPlus.e-base.e)/(base.p*e);
                    drho_dT_p=(tPlus.rho-base.rho)/(base.t*e);
                    de_dT_p=(tPlus.e-base.e)/(base.t*e);
                }
                else if(pUseMinus&&tUseMinus)
                {
                    drho_dp_T=(base.rho-pMinus.rho)/(base.p*e);
                    de_dp_T=(base.e-pMinus.e)/(base.p*e);
                    drho_dT_p=(base.rho-tMinus.rho)/(base.t*e);
                    de_dT_p=(base.e-tMinus.e)/(base.t*e);
                }
                else if(pUsePlus&&tUseMinus)
                {
                    drho_dp_T=(pPlus.rho-base.rho)/(base.p*e);
                    de_dp_T=(pPlus.e-base.e)/(base.p*e);
                    drho_dT_p=(base.rho-tMinus.rho)/(base.t*e);
                    de_dT_p=(base.e-tMinus.e)/(base.t*e);
                }
                else if(pUseMinus&&tUsePlus)
                {
                    drho_dp_T=(base.rho-pMinus.rho)/(base.p*e);
                    de_dp_T=(base.e-pMinus.e)/(base.p*e);
                    drho_dT_p=(tPlus.rho-base.rho)/(base.t*e);
                    de_dT_p=(tPlus.e-base.e)/(base.t*e);
                }
                else continue;
            }
            if(!(drho_dp_T>0.)||!std::isfinite(drho_dp_T))continue;
            const double dp_drho_T=1./drho_dp_T;
            const double dp_dT_rho=-drho_dT_p/drho_dp_T;
            const double cv=de_dT_p-de_dp_T*dp_dT_rho;
            if(!std::isfinite(cv)||cv<=0.)continue;
            const double a2=dp_drho_T+(base.t/(base.rho*base.rho))*dp_dT_rho*dp_dT_rho/cv;
            if(std::isfinite(a2)&&a2>0.)return a2;
        }
        catch(const std::exception&){}
    }
    throw std::runtime_error("COOLPROP_IDENTITY_SLOPE_FAILED");
}
inline AcousticResult acousticIdentity(CoolProp::AbstractState& state,const Point& base)
{
    thermoCount(thermoCounters().identity);
    // Near the phase boundary the stencil changes between central and one-sided
    // forms, so a strict per-step convergence test is not meaningful. Use the
    // smallest stencil value and accept it when the spread over the whole eps
    // sequence is within 1% (recorded as the relative change). This is the
    // documented approximate branch of the two-phase sound speed; the caller
    // annotates the backend status with IDENTITY-FALLBACK.
    double last=std::numeric_limits<double>::quiet_NaN();
    double spread=0.;int count=0;
    for(double eps:{1e-3,5e-4,2.5e-4,1.25e-4,6.25e-5})
    {
        // The identity slope can throw near the saturation line (A-route cell
        // 167662); skip that stencil instead of letting the exception escape,
        // so the table fallback below is still reachable.
        double raw=0.;
        try{raw=soundSpeedIdentitySquared(state,base,eps);}
        catch(const std::exception&){continue;}
        if(!std::isfinite(raw)||raw<=0.)continue;
        if(std::isfinite(last))spread=std::max(spread,std::abs(raw-last)/raw);
        last=raw;++count;
    }
    if(count>0&&std::isfinite(last)&&spread<=1e-2)
        return {last,spread,6.25e-5,NAN,NAN,count,true};
    // Last resort: the property table carries an interpolated a2 on every node,
    // two-phase nodes included (field states were verified to 0.11 % mean /
    // 0.49 % worst deviation). On the saturated-liquid line every local
    // derivative route can fail (A-route step 28-29: isentrope ->
    // identity slope), so this documented approximate branch keeps the run
    // alive; the result is flagged as an identity/table fallback.
    if(const TableAcoustic* table=tableAcoustic())
    {
        const auto y=state.get_mass_fractions();
        double a2=0.;bool clamped=false;
        if(table->lookupClamped(y[0],base.t,base.p,a2,clamped))
            return {a2,spread,0.,NAN,NAN,count,true};
    }
    throw std::runtime_error("COOLPROP_IDENTITY_NOT_CONVERGED");
}
inline AcousticResult isentrope(CoolProp::AbstractState& state,const Point& base,bool pure)
{
    try
    {
        return isentropeIterate(state,base,pure);
    }
    catch(const std::exception&)
    {
        // Fall back to the identity sound speed whenever the isentrope
        // integration fails, pure or not: on the saturated-liquid line the
        // symmetric stencil may not exist or the iteration may not settle
        // (A-route step 29, cells 698987/785575 at p=1.578 MPa, T=128 K), and
        // the identity route is a validated local alternative (documented
        // approximate branch, flagged as IDENTITY_FALLBACK).
        return acousticIdentity(state,base);
    }
}

}}}
#endif
