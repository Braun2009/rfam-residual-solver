#include "CoolPropThermoModel.H"
#include "StateAcceptance.H"
#include "TableAcoustic.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#ifdef RFAM_HAVE_COOLPROP
#include "AbstractState.h"
#include "CoolProp.h"
#include "CoolPropEquilibrium.H"
#include <mutex>
#endif

namespace Foam { namespace thermo {
#ifdef RFAM_HAVE_COOLPROP
namespace {
using CPState=CoolProp::AbstractState;
std::mutex factoryMutex;
std::unique_ptr<CPState> factory(const std::string& name)
{
    std::lock_guard<std::mutex> lock(factoryMutex);
    return std::unique_ptr<CPState>(CPState::factory("HEOS",name));
}
struct ThreadStates
{
    std::unique_ptr<CPState> oxygen=factory("Oxygen"),fuel=factory("Propane"),mixture=factory("Oxygen&Propane");
};
bool validY(const std::vector<double>& y)
{return y.size()==2&&std::isfinite(y[0])&&std::isfinite(y[1])&&y[0]>=0.&&y[1]>=0.&&y[0]<=1.&&y[1]<=1.&&std::abs(y[0]+y[1]-1.)<=1e-12;}
CPState& stateFor(const std::vector<double>& y)
{
    thread_local ThreadStates states;
    CPState& s=y[0]==1.?*states.oxygen:y[1]==1.?*states.fuel:*states.mixture;
    s.unspecify_phase();if(y[0]!=1.&&y[1]!=1.)s.set_mass_fractions(y);
    return s;
}
class Model final:public EquilibriumThermoModel
{
    std::vector<std::string> names_{"O2","C3H8"};
    std::vector<double> mw_,href_;
    std::string id_;
    EquilibriumResult failure(EquilibriumStatus status,const std::string& text)const
    {EquilibriumResult r;r.status=status;r.diagnostics.modelId=id_;r.diagnostics.backendStatus=text;return r;}
    double offset(const std::vector<double>& y)const{return y[0]*href_[0]+y[1]*href_[1];}
    EquilibriumResult publishTwoPhase(CPState& state,const std::vector<double>& y)const
    {
        const bool pure=y[0]==1.||y[1]==1.;
        thermoCount(thermoCounters().twoPhasePublish);
        const auto base=cpdetail::current(state);
        EquilibriumState q;q.region=EquilibriumRegion::TWO_PHASE;
        q.p=base.p;q.T=base.t;q.rho=base.rho;q.e=base.e-offset(y);q.h=base.h-offset(y);
        q.s=base.s;q.entropyAvailable=true;q.overallMassFractions=y;
        const auto makePhase=[&](bool vapor)
        {
            EquilibriumPhaseProperties ph;
            ph.moleFractions=pure?y:(vapor?state.mole_fractions_vapor_double():state.mole_fractions_liquid_double());
            ph.molecularWeight=ph.moleFractions[0]*mw_[0]+ph.moleFractions[1]*mw_[1];
            ph.massFractions={ph.moleFractions[0]*mw_[0]/ph.molecularWeight,ph.moleFractions[1]*mw_[1]/ph.molecularWeight};
            const auto value=[&](CoolProp::parameters key){return vapor?state.saturated_vapor_keyed_output(key):state.saturated_liquid_keyed_output(key);};
            ph.rho=value(CoolProp::iDmass);ph.e=value(CoolProp::iUmass)-offset(ph.massFractions);
            ph.h=value(CoolProp::iHmass)-offset(ph.massFractions);ph.s=value(CoolProp::iSmass);ph.entropyAvailable=true;
            ph.compressibility=q.p*ph.molecularWeight/(ph.rho*8.314462618*q.T);return ph;
        };
        q.liquid=makePhase(false);q.vapor=makePhase(true);
        q.vaporMoleFraction=base.quality;
        const double meanMW=(1.-base.quality)*q.liquid->molecularWeight+base.quality*q.vapor->molecularWeight;
        q.vaporMassFraction=base.quality*q.vapor->molecularWeight/meanMW;
        q.vaporVolumeFraction=q.vaporMassFraction*q.rho/q.vapor->rho;
        double material=0.;
        for(int i=0;i<2;++i)material=std::max(material,std::abs(y[i]-(1.-q.vaporMassFraction)*q.liquid->massFractions[i]-q.vaporMassFraction*q.vapor->massFractions[i]));
        const double leverRho=1./((1.-q.vaporMassFraction)/q.liquid->rho+q.vaporMassFraction/q.vapor->rho);
        const double leverE=(1.-q.vaporMassFraction)*q.liquid->e+q.vaporMassFraction*q.vapor->e;
        if(material>1e-9||std::abs(q.rho-leverRho)>1e-8*q.rho||std::abs(q.e-leverE)>1e-8*std::max(std::abs(q.e),1.))
            return failure(EquilibriumStatus::RESIDUAL_CONTRACT_FAILED,"COOLPROP_PHASE_LEVER_CLOSURE_FAILED");
        // Restore the working state to the base point (mixtures only): the
        // saturation queries above change the state cache/history that the
        // mixture isentrope inversion relies on. Pure fluids use the
        // DmassSmass input directly and must keep the state the rho-e flash
        // produced.
        if(!pure)
        {
            state.unspecify_phase();
            state.update(CoolProp::PT_INPUTS, base.p, base.t);
        }
        // The converged rho-s perturbations are independent equilibrium flashes.
        // Leave the single-phase derivative payload absent: mixture derivatives
        // across phase transfer are not frozen single-phase Cv/Cp quantities.
        // Optional table-backed acoustic path (RFAM_TETRA_TABLE): the a^2 is
        // interpolated from a verified RFAMTAB2 table; otherwise the converged
        // EOS isentrope is computed as before.
        double tableA2=-1.;bool fromTable=false;
        if(const TableAcoustic* table=tableAcoustic())
            fromTable=table->lookup(y[0],base.t,base.p,tableA2);
        if(fromTable)thermoCount(thermoCounters().tableAcousticHit);
        const auto acoustic=fromTable?cpdetail::AcousticResult{}:cpdetail::isentrope(state,base,pure);
        q.acoustic.closure=AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
        q.acoustic.evidence=fromTable?AcousticEvidence::TABLE_INTERPOLATED
                                     :AcousticEvidence::CONVERGED_RHOE_ISENTROPE;
        q.acoustic.pressureDensityDerivative=fromTable?tableA2:acoustic.a2;
        if(!stateValidForModel(q))return failure(EquilibriumStatus::RESIDUAL_CONTRACT_FAILED,"COOLPROP_TWO_PHASE_CONTRACT_FAILED");
        EquilibriumResult r;r.status=EquilibriumStatus::SUCCESS;r.state=std::move(q);
        auto& d=r.diagnostics;d.modelId=id_;
        d.backendStatus=fromTable?"COOLPROP_HEOS_EQUILIBRIUM_TABLE_INTERPOLATED"
            :acoustic.identityFallback
                ?"COOLPROP_HEOS_EQUILIBRIUM_ISENTROPE+IDENTITY_FALLBACK"
                :"COOLPROP_HEOS_EQUILIBRIUM_ISENTROPE";
        d.materialResidual=material;d.energyIdentityResidual=r.state->e-(r.state->h-r.state->p/r.state->rho);
        d.acousticLowerEntropyResidual=acoustic.lowerEntropyResidual;d.acousticUpperEntropyResidual=acoustic.upperEntropyResidual;
        d.acousticUsedRelativeDensityPerturbation=acoustic.eps;d.acousticDerivativeRelativeChange=acoustic.relativeChange;
        d.acousticDerivativeConvergenceComparisons=acoustic.comparisons;return r;
    }
    EquilibriumResult publish(CPState& s,const std::vector<double>& y)const
    {
        if(s.phase()==CoolProp::iphase_twophase)return publishTwoPhase(s,y);
        EquilibriumState q;
        // CoolProp may label a supercritical binary PT root as liquid. Keep
        // mixture phase identity unclassified instead of inventing a label.
        q.region=EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE;
        if(y[0]==1.||y[1]==1.)
        {
            switch(s.phase())
            {
                case CoolProp::iphase_liquid:case CoolProp::iphase_supercritical_liquid:q.region=EquilibriumRegion::SINGLE_LIQUID;break;
                case CoolProp::iphase_gas:case CoolProp::iphase_supercritical_gas:q.region=EquilibriumRegion::SINGLE_VAPOR;break;
                case CoolProp::iphase_supercritical:q.region=EquilibriumRegion::SUPERCRITICAL;break;
                default:break;
            }
        }
        q.overallMassFractions=y;q.p=s.p();q.T=s.T();q.rho=s.rhomass();
        q.h=s.hmass()-offset(y);q.e=s.umass()-offset(y);q.s=s.smass();q.entropyAvailable=true;
        q.acoustic.closure=AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
        q.acoustic.evidence=AcousticEvidence::ANALYTIC_FROZEN_SINGLE_PHASE;
        q.acoustic.pressureDensityDerivative=std::pow(s.speed_sound(),2);
        // (d a / d p)_s for the pressure outlet second-order acoustic
        // increment. CoolProp supplies it directly; a backend failure leaves
        // NaN and the outlet keeps the exact isentrope path.
        // (d a / d p)_s by the chain rule on CoolProp keyed derivatives:
        //   (da/dp)_s = (da/dp)_T + (da/dT)_p * (dT/dp)_s
        // CoolProp has no direct (A, P, S) keyed derivative (iA is the
        // Helmholtz energy, not the sound speed; the sound speed is
        // ispeed_sound). A backend that cannot supply the pieces leaves NaN
        // and the outlet keeps the first-order shortcut with its tight band.
        try
        {
            const double da_dp_T=
                s.first_partial_deriv(CoolProp::ispeed_sound,CoolProp::iP,CoolProp::iT);
            const double da_dT_p=
                s.first_partial_deriv(CoolProp::ispeed_sound,CoolProp::iT,CoolProp::iP);
            const double dT_dp_s=
                s.first_partial_deriv(CoolProp::iT,CoolProp::iP,CoolProp::iSmass);
            q.acoustic.soundSpeedPressureDerivative=da_dp_T+da_dT_p*dT_dp_s;
            if(!std::isfinite(q.acoustic.soundSpeedPressureDerivative))
                q.acoustic.soundSpeedPressureDerivative=
                    std::numeric_limits<double>::quiet_NaN();
        }
        catch(const std::exception&){q.acoustic.soundSpeedPressureDerivative=
            std::numeric_limits<double>::quiet_NaN();}
        EquilibriumThermodynamicPayload td;
        td.closure=q.acoustic.closure;td.Cp=s.cpmass();td.Cv=s.cvmass();
        td.dpdrho_T=s.first_partial_deriv(CoolProp::iP,CoolProp::iDmass,CoolProp::iT);
        td.dpdT_rho=s.first_partial_deriv(CoolProp::iP,CoolProp::iT,CoolProp::iDmass);
        td.dedrho_T=s.first_partial_deriv(CoolProp::iUmass,CoolProp::iDmass,CoolProp::iT);
        td.dedT_rho=s.first_partial_deriv(CoolProp::iUmass,CoolProp::iT,CoolProp::iDmass);q.thermodynamics=td;
        EquilibriumPhaseProperties ph;ph.massFractions=y;
        ph.molecularWeight=1./(y[0]/mw_[0]+y[1]/mw_[1]);
        ph.moleFractions={y[0]*ph.molecularWeight/mw_[0],y[1]*ph.molecularWeight/mw_[1]};
        ph.rho=q.rho;ph.h=q.h;ph.e=q.e;ph.s=q.s;ph.entropyAvailable=true;
        ph.compressibility=q.p*ph.molecularWeight/(q.rho*8.314462618*q.T);q.singlePhase=ph;
        if(!stateValidForModel(q))return failure(EquilibriumStatus::RESIDUAL_CONTRACT_FAILED,"COOLPROP_STATE_CONTRACT_FAILED");
        EquilibriumResult r;r.status=EquilibriumStatus::SUCCESS;r.state=std::move(q);
        r.diagnostics.modelId=id_;r.diagnostics.backendStatus="COOLPROP_HEOS_SINGLE_PHASE";
        r.diagnostics.energyIdentityResidual=r.state->e-(r.state->h-r.state->p/r.state->rho);return r;
    }
public:
    Model()
    {
        if(CoolProp::get_global_param_string("version")!="7.2.0")throw std::runtime_error("native CoolProp must be 7.2.0");
        for(const auto& name:{"Oxygen","Propane"})
        {
            auto s=factory(name);s->update(CoolProp::DmassT_INPUTS,1e-10,298.15);
            href_.push_back(s->hmass_idealgas());mw_.push_back(s->molar_mass());
        }
        id_="coolprop-7.2.0-HEOS-Oxygen-Propane-h0-298.15-isentrope-v1:git="
            +CoolProp::get_global_param_string("gitrevision");
#ifdef RFAM_COOLPROP_PROVENANCE
        id_+=";source=" RFAM_COOLPROP_PROVENANCE;
#endif
    }
    const std::string& modelId()const override{return id_;}
    const std::vector<std::string>& speciesNames()const override{return names_;}
    EnergyConvention energyConvention()const override{return EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;}
    bool supportsPrimitivePTYInitialization()const override{return true;}
    bool stateValidForModel(const EquilibriumState& q)const override
    {return q.valid()&&validateEquilibriumSpeciesBasis(q,mw_);}
    EquilibriumResult initializePTY(double p,double t,const std::vector<double>& y)const override
    {
        if(!validY(y)||!std::isfinite(p)||p<=0.||!std::isfinite(t)||t<=0.)return failure(EquilibriumStatus::INVALID_INPUT,"COOLPROP_INVALID_PTY");
        try{auto& s=stateFor(y);s.update(CoolProp::PT_INPUTS,p,t);return publish(s,y);}
        catch(const std::exception& e){return failure(EquilibriumStatus::NO_EQUILIBRIUM_STATE,e.what());}
    }
    EquilibriumResult initializePSY(double p,double entropy,const std::vector<double>& y,const EquilibriumSeed* seed)const override
    {
        if(!validY(y)||!std::isfinite(p)||p<=0.||!std::isfinite(entropy))
            return failure(EquilibriumStatus::INVALID_INPUT,"COOLPROP_INVALID_PSY");
        try
        {
            auto& s=stateFor(y);
            if(y[0]==1.||y[1]==1.)
            {
                try{s.update(CoolProp::PSmass_INPUTS,p,entropy);}
                catch(const std::exception&)
                {
                    // CoolProp's pure-fluid PS flash can fail close to the O2
                    // critical pressure even for an ordinary 130 K liquid.
                    // Only a supplied PT seed may initiate this single-phase
                    // bracket; no off-entropy state is ever returned.
                    if(!seed||!std::isfinite(seed->temperature))throw;
                    const double t0=seed->temperature;
                    double belowT=0.,aboveT=0.;bool bracketed=false;
                    for(int level=0;level<17&&!bracketed;++level)
                    {
                        const double step=.005*t0*std::pow(1.5,level);
                        for(double candidate:{t0-step,t0+step})
                        {
                            if(candidate<60.||candidate>2000.)continue;
                            try
                            {
                                const auto point=cpdetail::pt(s,p,candidate);
                                if(point.quality>=0.)continue;
                                if(point.s<entropy)belowT=std::max(belowT,candidate);
                                if(point.s>entropy&&(aboveT==0.||candidate<aboveT))aboveT=candidate;
                            }
                            catch(const std::exception&){}
                        }
                        bracketed=belowT>0.&&aboveT>belowT;
                    }
                    if(!bracketed)throw;
                    bool closed=false;
                    for(int iter=0;iter<65;++iter)
                    {
                        const double t=.5*(belowT+aboveT);
                        const auto point=cpdetail::pt(s,p,t);
                        if(point.quality>=0.)throw std::runtime_error("pure PS PT bracket crossed saturation");
                        const double error=point.s-entropy;
                        if(std::abs(error)<=1e-10*std::max(std::abs(entropy),1.))
                        {closed=true;break;}
                        if(error<0.)belowT=t;else aboveT=t;
                    }
                    if(!closed)throw std::runtime_error("pure PS PT bracket did not close entropy");
                }
            }
            else
            {
                double t=seed&&std::isfinite(seed->temperature)?seed->temperature:400.;
                bool converged=false;
                for(int iter=0;iter<30;++iter)
                {
                    const auto q=cpdetail::pt(s,p,t);const double error=q.s-entropy;
                    if(std::abs(error)<1e-10*std::max(std::abs(entropy),1.)){converged=true;break;}
                    const auto d=cpdetail::derivatives(s,q);
                    if(!std::isfinite(d.sT)||d.sT<=0.)break;
                    const double correction=std::clamp(-error/d.sT,-.1*t,.1*t);
                    bool advanced=false;
                    for(double scale=1.;scale>=1./1024.;scale*=.5)
                    {
                        try
                        {
                            const double next=t+scale*correction;
                            if(next<60.||next>2000.)continue;
                            if(std::abs(cpdetail::pt(s,p,next).s-entropy)<std::abs(error))
                            {t=next;advanced=true;break;}
                        }catch(const std::exception&){}
                    }
                    if(!advanced)break;
                }
                if(!converged)return failure(EquilibriumStatus::NOT_CONVERGED,"COOLPROP_PSY_NEWTON_FAILED");
            }
            auto result=publish(s,y);
            if(result.success()&&std::abs(result.state->s-entropy)>1e-8*std::max(std::abs(entropy),1.))
                return failure(EquilibriumStatus::RESIDUAL_CONTRACT_FAILED,"COOLPROP_PSY_ENTROPY_MISMATCH");
            return result;
        }catch(const std::exception& e){return failure(EquilibriumStatus::NO_EQUILIBRIUM_STATE,e.what());}
    }
    EquilibriumResult recoverRhoE(double rho,double e,const std::vector<double>& y,const EquilibriumSeed* seed)const override
    {
        if(!validY(y)||!std::isfinite(rho)||rho<=0.||!std::isfinite(e))return failure(EquilibriumStatus::INVALID_INPUT,"COOLPROP_INVALID_RHOE");
        try
        {
            // Trace-species snap: the FV update leaves ~1e-14 mass fractions of
            // the minor species on a nominally pure stream (A-route cell 132114:
            // rhoY_C3H8 = 2.7e-11 against rho=920.66). Routing that through the
            // mixture flash asks CoolProp to solve a pure-fluid state on the
            // unstable side of the mixture phase split and fails, although the
            // pure-fluid (rho,u) inversion exists. Compositions whose minor
            // species are below 1e-8 are snapped to pure for the EOS; the
            // conserved mass is untouched. The trace grows as the mixing front
            // advances (3e-14 at step 4 -> 1.1e-12 at step 10), so the snap is
            // an order below the production acceptance gate (2e-8).
            // Snap tolerance derives from the acceptance gate: a composition
            // snap shifts the energy reference by (1-max y)*dHref (dHref ~ 362
            // kJ/kg between the pure components), so allowing half the gate
            // gives the largest defensible threshold. The cap keeps a physically
            // meaningful mixture (>= 100 ppm minor species) on the mixture path.
            const double dHref=std::abs(offset({{1.,0.}})-offset({{0.,1.}}));
            const double snapY=std::clamp(0.5*stateAcceptanceRel()*std::max(std::abs(e),1.)
                                          /std::max(dHref,1.),1e-12,1e-4);
            std::vector<double> yEff=y;
            {
                double mx=0.;for(double v:y)mx=std::max(mx,v);
                if(mx>=1.0-snapY)
                    for(auto& v:yEff)v=(v==mx?1.0:0.0);
            }
            auto& s=stateFor(yEff);const double target=e+offset(yEff);
            // State used for the final publish; the last-resort pure retry
            // rebinds it to the pure-fluid backend (stateFor returns separate
            // Oxygen/Propane/mixture backends, and set_mass_fractions on the
            // mixture backend does NOT leave the mixture critical-point code).
            CPState* publishState=&s;
            if(yEff[0]==1.||yEff[1]==1.)
            {
                // Pure fluid: the library's (rho,u) flash is the fast path, but
                // it fails on the steep branch just inside the saturated-liquid
                // line (A-route cell 132114). Fall back to the lever-rule
                // inversion on the saturation curve before giving up.
                bool pureOk=true;
                try{s.update(CoolProp::DmassUmass_INPUTS,rho,target);}
                catch(const std::exception&){pureOk=false;}
                if(!pureOk)
                {
                    const double p0=seed&&std::isfinite(seed->pressure)?seed->pressure:1e6;
                    const double t0=seed&&std::isfinite(seed->temperature)?seed->temperature:400.;
                    if(!cpdetail::recoverPureOnSaturation(s,rho,target,p0,t0,
                                                          std::max(std::abs(e),1.)))
                        throw std::runtime_error("COOLPROP_RHOE_NEWTON_FAILED");
                }
            }
            else
            {
                const double p=seed&&std::isfinite(seed->pressure)?seed->pressure:1e6;
                const double t=seed&&std::isfinite(seed->temperature)?seed->temperature:400.;
                // Fast single-phase candidate at fixed density. The imposed
                // phase here enables only local EOS evaluation; it is removed
                // before the authoritative PT stability/phase flash. The
                // candidate cannot be published until that flash agrees.
                bool verified=false;
                // Single-phase fast branch. The phase prior matters: forcing
                // iphase_gas on a compressed-liquid LOX cell made CoolProp chase a
                // metastable vapour root, so the rho-e Newton below never saw the
                // liquid branch and the whole recovery fell through to the flash
                // scan, which then fails inside the near-saturation envelope
                // (measured 2026-10-16: p=2.85 MPa, rho=922, Y_O2=0.784 ->
                // COOLPROP_RHOE_NEWTON_FAILED at step 8-11 of the B-grid steady run).
                // Try the composition-aware phase first, then the other branch; each
                // attempt is still fully verified before it can be published.
                const auto trySpecifiedPhase=[&](CoolProp::phases phase,double tSeed)->bool
                {
                    bool ok=false;
                    try
                    {
                        s.specify_phase(phase);
                        double trialT=tSeed;
                        for(int iter=0;iter<12;++iter)
                        {
                            s.update(CoolProp::DmassT_INPUTS,rho,trialT);
                            const double residual=s.umass()-target;
                            if(std::abs(residual)<=stateAcceptanceRel()*std::max(std::abs(e),1.))
                            {
                                const double candidateP=s.p();
                                // Table-assisted verification: when the 3x3x3
                                // table block around the state is entirely
                                // single phase, the full mixture PT stability
                                // flash (the dominant per-cell cost) can be
                                // skipped. Near the phase boundary the strict
                                // flash still runs.
                                bool tableConfident=false;
                                if(const TableAcoustic* table=tableAcoustic())
                                    tableConfident=table->singlePhaseBlock(y[0],trialT,candidateP);
                                if(tableConfident)
                                {
                                    // Safety check: a metastable single-phase
                                    // candidate (true equilibrium two-phase)
                                    // has a non-physical frozen sound speed.
                                    // Reject it and use the strict flash.
                                    try
                                    {
                                        const double a=s.speed_sound();
                                        if(!std::isfinite(a)||!(a>0.))tableConfident=false;
                                    }
                                    catch(const std::exception&){tableConfident=false;}
                                }
                                if(tableConfident)thermoCount(thermoCounters().tableAcousticHit);
                                if(tableConfident){ok=true;verified=true;break;}
                                s.unspecify_phase();s.update(CoolProp::PT_INPUTS,candidateP,trialT);
                                ok=s.phase()!=CoolProp::iphase_twophase
                                    &&std::abs(s.rhomass()-rho)<=1e-10*std::max(rho,1.)
                                    &&std::abs(s.umass()-target)<=stateAcceptanceRel()*std::max(std::abs(e),1.);
                                verified=ok;
                                break;
                            }
                            const double cv=s.cvmass();
                            if(!std::isfinite(cv)||cv<=0.)break;
                            trialT+=std::clamp(-residual/cv,-.15*trialT,.15*trialT);
                            if(trialT<60.||trialT>2000.)break;
                        }
                        s.unspecify_phase();
                    }
                    catch(const std::exception&){s.unspecify_phase();}
                    return ok;
                };
                if(seed&&seed->previousRegion!=EquilibriumRegion::TWO_PHASE)
                {
                    // Phase order: the region the cell is coming from, then the
                    // other single-phase branch. A dense state (rho above the
                    // mixture critical density estimate) starts on the liquid
                    // branch because the label is UNCLASSIFIED for mixtures.
                    const bool dense=rho>0.75*(1.0/(y[0]/436.14+(1.0-y[0])/220.48));
                    CoolProp::phases first=CoolProp::iphase_gas;
                    if(seed->previousRegion==EquilibriumRegion::SINGLE_LIQUID)first=CoolProp::iphase_liquid;
                    else if(seed->previousRegion==EquilibriumRegion::SUPERCRITICAL)first=CoolProp::iphase_supercritical;
                    else if(seed->previousRegion==EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE)
                        first=dense?CoolProp::iphase_liquid:CoolProp::iphase_gas;
                    const CoolProp::phases second=(first==CoolProp::iphase_liquid)
                        ?CoolProp::iphase_gas:CoolProp::iphase_liquid;
                    verified=trySpecifiedPhase(first,t);
                    if(!verified&&second!=first)verified=trySpecifiedPhase(second,t);
                }
                if(!verified)
                {
                    // Deterministic starting point. The per-thread CoolProp
                    // state still holds the previous cell's root; without
                    // re-specifying the seed the backward flash starts from a
                    // history-dependent guess (A-route cell 559462 recovers in
                    // isolation but not inside the run, where recovFail grows).
                    try{s.unspecify_phase();s.update(CoolProp::PT_INPUTS,p,t);}
                    catch(const std::exception&){}
                    try{cpdetail::recover(s,rho,target,p,t,false,std::max(std::abs(e),1.));}
                    catch(const std::exception&)
                    {
                        // Last resort: when the mixture flash fails on a stream
                        // that is >=99.9 % one component (the LOX channels carry
                        // a growing trace of the other propellant as the mixing
                        // front advances), retry as a pure fluid. The snap and
                        // its size are recorded per step; the conserved mass is
                        // never modified.
                        double mx=0.;for(double v:y)mx=std::max(mx,v);
                        if(mx<1.0-1e-3)throw;
                        thermoCount(thermoCounters().pureFallbackEvents);
                        thermoMax(thermoCounters().maxPureSnapFraction,1.0-mx);
                        for(auto& v:yEff)v=(v==mx?1.0:0.0);
                        publishState=&stateFor(yEff);   // real pure-fluid backend
                        if(!cpdetail::recoverPureOnSaturation(*publishState,rho,e+offset(yEff),p,t,
                                                              std::max(std::abs(e),1.)))
                            throw;
                    }
                }
            }
            EquilibriumResult r;
            try{r=publish(s,yEff);}
            catch(const std::exception& ex)
            {
                // The mixture stability/critical-point routines can throw at
                // trace compositions; retry as a pure fluid when the stream is
                // essentially one component (recorded, mass untouched).
                double mx=0.;for(double v:yEff)mx=std::max(mx,v);
                if(mx<1.0-1e-3)return failure(EquilibriumStatus::NO_EQUILIBRIUM_STATE,ex.what());
                thermoCount(thermoCounters().pureFallbackEvents);
                thermoMax(thermoCounters().maxPureSnapFraction,1.0-mx);
                for(auto& v:yEff)v=(v==mx?1.0:0.0);
                publishState=&stateFor(yEff);   // real pure-fluid backend
                const double p0=seed&&std::isfinite(seed->pressure)?seed->pressure:1e6;
                const double t0=seed&&std::isfinite(seed->temperature)?seed->temperature:400.;
                bool ok=true;
                try{publishState->update(CoolProp::DmassUmass_INPUTS,rho,e+offset(yEff));}
                catch(const std::exception&){ok=false;}
                if(!ok&&!cpdetail::recoverPureOnSaturation(*publishState,rho,e+offset(yEff),p0,t0,
                                                           std::max(std::abs(e),1.)))
                    return failure(EquilibriumStatus::NO_EQUILIBRIUM_STATE,ex.what());
                r=publish(*publishState,yEff);
            }
            if(!r.success())return r;
            r.diagnostics.densityResidual=r.state->rho-rho;r.diagnostics.energyResidual=r.state->e-e;
            thermoMax(thermoCounters().maxRelRho,
                      std::abs(r.diagnostics.densityResidual)/std::max(rho,1.));
            thermoMax(thermoCounters().maxRelEnergy,
                      std::abs(r.diagnostics.energyResidual)/std::max(std::abs(e),1.));
            if(!projectionAcceptEnabled()
               &&(std::abs(r.diagnostics.densityResidual)>stateAcceptanceRel()*std::max(rho,1.)
                  ||std::abs(r.diagnostics.energyResidual)>stateAcceptanceRel()*std::max(std::abs(e),1.)))
                return failure(EquilibriumStatus::RESIDUAL_CONTRACT_FAILED,"COOLPROP_RHOE_CLOSURE_FAILED");
            return r;
        }catch(const std::exception& ex){return failure(EquilibriumStatus::NO_EQUILIBRIUM_STATE,ex.what());}
    }
};
}
std::shared_ptr<const EquilibriumThermoModel> makeCoolPropThermoModel(){return std::make_shared<Model>();}
bool coolPropThermoAvailable(){return true;}
#else
std::shared_ptr<const EquilibriumThermoModel> makeCoolPropThermoModel()
{throw std::runtime_error("CoolProp backend unavailable: configure RFAM_HAVE_COOLPROP with the native 7.2.0 library");}
bool coolPropThermoAvailable(){return false;}
#endif
}}
