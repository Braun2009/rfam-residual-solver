#ifndef RFAM_TETRA_CHEMISTRY_CVODE_H
#define RFAM_TETRA_CHEMISTRY_CVODE_H
#include "ChemistryModel.H"
#include "SandiegoNasa7Data.H"
#include "TetraStateRecovery.H"
#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

namespace Foam
{
// A reusable BDF workspace, owned by one worker; no ideal-gas EOS fallback.
class TetraChemistryCVODE
{
public:
    struct Result
    {
        std::optional<CompleteCellState> state;
        long steps=0,rhsCalls=0;
        double reachedTime=0.,compositionCorrection=0.,energyDefect=0.,elementDefect=0.;
        std::string reason;
        bool success() const {return state.has_value();}
    };
private:
    const ChemistryModel& chemistry_;
    const thermo::EquilibriumThermoModel& model_;
    std::vector<std::string> names_;
    std::vector<std::array<double,6>> atoms_;
    std::vector<double> formation_,Y_,omega_;
    SUNContext context_=nullptr;
    N_Vector values_=nullptr,atol_=nullptr,constraints_=nullptr;
    SUNMatrix matrix_=nullptr;
    SUNLinearSolver solver_=nullptr;
    void* memory_=nullptr;
    bool initialized_=false;
    double rho_=0.;
    std::string lastRhsFailure_;
    long failedRhs_=0;
    long eosFailures_=0;
    double lastT_=0.,lastMinimumY_=0.;
    thermo::EquilibriumSeed seed_;
    static constexpr double rtol_=1e-7,speciesAtol_=1e-16,energyAtol_=1e-4;
    static bool gasState(const thermo::EquilibriumState& q)
    {
        const bool region=q.region==thermo::EquilibriumRegion::SINGLE_VAPOR
            ||q.region==thermo::EquilibriumRegion::SUPERCRITICAL
            ||(q.region==thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE&&q.singlePhase);
        return region&&q.T>=500.&&q.thermodynamics&&!q.condensed
            &&!q.conservedStateProjected&&tetraAdmitsAcoustics(q);
    }
    void freeResources() noexcept
    {
        if(memory_)CVodeFree(&memory_);
        if(solver_)SUNLinSolFree(solver_);
        if(matrix_)SUNMatDestroy(matrix_);
        if(constraints_)N_VDestroy(constraints_);
        if(atol_)N_VDestroy(atol_);
        if(values_)N_VDestroy(values_);
        if(context_)SUNContext_Free(&context_);
    }
    bool prepare(const sunrealtype* raw)
    {
        double sum=0.;
        for(std::size_t s=0;s<names_.size();++s)
        {
            if(!std::isfinite(raw[s])||raw[s]<-10.*speciesAtol_)return false;
            Y_[s]=std::max(0.,double(raw[s]));sum+=Y_[s];
        }
        if(!std::isfinite(sum)||sum<=0.)return false;
        for(double& y:Y_)y/=sum;
        return true;
    }
    static int rhs(sunrealtype,N_Vector y,N_Vector derivative,void* data) noexcept
    {
        auto& self=*static_cast<TetraChemistryCVODE*>(data);
        try
        {
            const auto* raw=N_VGetArrayPointer(y);
            if(!self.prepare(raw)||!std::isfinite(raw[self.names_.size()]))
            {
                ++self.failedRhs_;self.lastRhsFailure_="composition/nonfinite energy";
                self.lastMinimumY_=*std::min_element(raw,raw+self.names_.size());return 1;
            }
            const auto q=self.model_.recoverRhoE(self.rho_,raw[self.names_.size()],self.Y_,&self.seed_);
            if(!q.success()||!gasState(*q.state))
            {
                ++self.failedRhs_;self.lastRhsFailure_=q.diagnostics.backendStatus;
                ++self.eosFailures_;
                if(self.eosFailures_<=1&&std::getenv("RFAM_CVODE_RHS_DIAG"))
                {
                    std::cerr<<std::setprecision(17)<<"CVODE_EOS_REQUEST rho="<<self.rho_
                        <<" e="<<raw[self.names_.size()]<<" Tseed="<<self.seed_.temperature
                        <<" status="<<q.diagnostics.backendStatus<<" Y=";
                    for(double y:self.Y_)std::cerr<<y<<",";
                    std::cerr<<std::endl;
                }
                if(q.success())self.lastRhsFailure_+=" inadmitted_region";return 1;
            }
            self.lastT_=q.state->T;
            double heat=0.;
            self.chemistry_.computeSourceTerms(q.state->T,q.state->p,self.rho_,self.Y_,self.omega_,heat);
            if(self.omega_.size()!=self.names_.size()||!std::isfinite(heat))return -1;
            auto* out=N_VGetArrayPointer(derivative);
            for(std::size_t s=0;s<self.names_.size();++s)
            {
                if(!std::isfinite(self.omega_[s]))return -1;
                out[s]=self.omega_[s]/self.rho_;
            }
            out[self.names_.size()]=heat/self.rho_;return 0;
        }
        catch(const std::exception& e)
        {++self.failedRhs_;self.lastRhsFailure_=e.what();return 1;}
        catch(...){++self.failedRhs_;self.lastRhsFailure_="unknown RHS exception";return 1;}
    }
public:
    TetraChemistryCVODE(const ChemistryModel& chemistry,const thermo::EquilibriumThermoModel& model);
    ~TetraChemistryCVODE(){freeResources();}
    TetraChemistryCVODE(const TetraChemistryCVODE&)=delete;
    TetraChemistryCVODE& operator=(const TetraChemistryCVODE&)=delete;
    Result advance(const PrimitiveVars& w,const ConservedVars& u,
                   const thermo::EquilibriumState& q,double dt,long maxSteps);
};
inline TetraChemistryCVODE::TetraChemistryCVODE(
    const ChemistryModel& chemistry,const thermo::EquilibriumThermoModel& model)
    :chemistry_(chemistry),model_(model),names_(model.speciesNames()),
     atoms_(names_.size()),formation_(names_.size()),Y_(names_.size()),omega_(names_.size())
{
    if(names_!=chemistry.speciesNames()||names_.size()!=SandiegoNasa7::nSpecies
       ||model.energyConvention()!=EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT)
        throw std::invalid_argument("CVODE chemistry/thermo identity mismatch");
    const auto& mw=chemistry.molecularWeights();
    if(mw.size()!=names_.size())throw std::invalid_argument("CVODE molecular weight dimension mismatch");
    for(std::size_t s=0;s<names_.size();++s)
    {
        const auto& entry=SandiegoNasa7::species[s];
        if(names_[s]!=entry.name||!std::isfinite(mw[s])||mw[s]<=0.
           ||std::abs(mw[s]*1e-3-entry.mw)>1e-12*entry.mw)
            throw std::invalid_argument("CVODE San Diego species identity mismatch");
        const double m=mw[s]*1e-3;
        atoms_[s]={entry.nC/m,entry.nH/m,entry.nO/m,entry.nN/m,
            (names_[s]=="AR"?1.:0.)/m,(names_[s]=="HE"?1.:0.)/m};
        formation_[s]=SandiegoNasa7::hMole(s,298.15)/m;
    }
    try
    {
        if(SUNContext_Create(SUN_COMM_NULL,&context_))throw std::runtime_error("SUNContext_Create failed");
        const auto n=static_cast<sunindextype>(names_.size()+1);
        values_=N_VNew_Serial(n,context_);atol_=N_VNew_Serial(n,context_);
        constraints_=N_VNew_Serial(n,context_);
        if(!values_||!atol_||!constraints_)throw std::runtime_error("CVODE vector allocation failed");
        for(std::size_t s=0;s<names_.size();++s)
        {NV_Ith_S(atol_,s)=speciesAtol_;NV_Ith_S(constraints_,s)=1.;}
        NV_Ith_S(atol_,names_.size())=energyAtol_;NV_Ith_S(constraints_,names_.size())=0.;
        memory_=CVodeCreate(CV_BDF,context_);matrix_=SUNDenseMatrix(n,n,context_);
        solver_=SUNLinSol_Dense(values_,matrix_,context_);
        if(!memory_||!matrix_||!solver_)throw std::runtime_error("CVODE solver allocation failed");
    }
    catch(...){freeResources();throw;}
}
inline TetraChemistryCVODE::Result TetraChemistryCVODE::advance(
    const PrimitiveVars& w,const ConservedVars& u,const thermo::EquilibriumState& q,
    double dt,long maxSteps)
{
    Result result;
    if(!std::isfinite(dt)||dt<=0.||maxSteps<=0||w.Y.size()!=names_.size()
       ||u.rhoY.size()!=names_.size()||!gasState(q)
       ||!std::isfinite(u.rho)||u.rho<=0.)
    {result.reason="CVODE invalid initial state/time";return result;}
    rho_=u.rho;seed_.pressure=q.p;seed_.temperature=q.T;
    failedRhs_=0;eosFailures_=0;lastRhsFailure_.clear();lastT_=q.T;lastMinimumY_=0.;
    seed_.previousTemperature=q.T;seed_.previousRegion=q.region;
    double kinetic=0.;for(double velocity:w.U)kinetic+=.5*velocity*velocity;
    const double initialEnergy=u.rhoE/u.rho-kinetic;
    std::vector<double> initialY(names_.size());
    for(std::size_t s=0;s<names_.size();++s)
    {initialY[s]=u.rhoY[s]/u.rho;NV_Ith_S(values_,s)=initialY[s];}
    NV_Ith_S(values_,names_.size())=initialEnergy;
    if(!prepare(N_VGetArrayPointer(values_)))
    {result.reason="CVODE invalid initial composition";return result;}
    int flag=initialized_?CVodeReInit(memory_,0.,values_):CVodeInit(memory_,rhs,0.,values_);
    if(flag!=CV_SUCCESS){result.reason="CVODE initialization failed";return result;}
    initialized_=true;
    if(CVodeSVtolerances(memory_,rtol_,atol_)!=CV_SUCCESS
       ||CVodeSetUserData(memory_,this)!=CV_SUCCESS
       ||CVodeSetLinearSolver(memory_,solver_,matrix_)!=CV_SUCCESS
       ||CVodeSetConstraints(memory_,constraints_)!=CV_SUCCESS
       ||CVodeSetMaxNumSteps(memory_,maxSteps)!=CV_SUCCESS
       ||CVodeSetMaxStep(memory_,dt)!=CV_SUCCESS||CVodeSetStopTime(memory_,dt)!=CV_SUCCESS)
    {result.reason="CVODE configuration failed";return result;}
    sunrealtype reached=0.;flag=CVode(memory_,dt,values_,&reached,CV_NORMAL);
    result.reachedTime=reached;CVodeGetNumSteps(memory_,&result.steps);
    CVodeGetNumRhsEvals(memory_,&result.rhsCalls);
    if((flag!=CV_SUCCESS&&flag!=CV_TSTOP_RETURN)
       ||std::abs(double(reached)-dt)>32.*std::numeric_limits<double>::epsilon()*dt)
    {
        result.reason="CVODE did not reach requested time flag="+std::to_string(flag)
            +" steps="+std::to_string(result.steps)+" rhs="+std::to_string(result.rhsCalls)
            +" rhs_failures="+std::to_string(failedRhs_)+" last_T="+std::to_string(lastT_)
            +" min_Y="+std::to_string(lastMinimumY_)+" last_rhs_reason="+lastRhsFailure_;
        return result;
    }
    const auto* raw=N_VGetArrayPointer(values_);double sum=0.;
    for(std::size_t s=0;s<names_.size();++s)
    {
        if(!std::isfinite(raw[s])||raw[s]<-10.*speciesAtol_)
        {result.reason="CVODE invalid final composition";return result;}
        sum+=std::max(0.,double(raw[s]));
    }
    if(!std::isfinite(sum)||std::abs(sum-1.)>1e-8||!prepare(raw))
    {result.reason="CVODE final mass mismatch";return result;}
    std::array<double,6> before{},after{};double oldFormation=0.,newFormation=0.;
    for(std::size_t s=0;s<names_.size();++s)
    {
        result.compositionCorrection=std::max(result.compositionCorrection,std::abs(Y_[s]-raw[s]));
        oldFormation+=formation_[s]*initialY[s];newFormation+=formation_[s]*Y_[s];
        for(std::size_t e=0;e<before.size();++e)
        {before[e]+=initialY[s]*atoms_[s][e];after[e]+=Y_[s]*atoms_[s][e];}
    }
    for(std::size_t e=0;e<before.size();++e)
        result.elementDefect=std::max(result.elementDefect,
            std::abs(before[e]-after[e])/std::max(std::abs(before[e]),1.));
    if(result.elementDefect>1e-7)
    {result.reason="CVODE element balance failed";return result;}
    const double finalEnergy=raw[names_.size()];
    result.energyDefect=std::abs(finalEnergy-initialEnergy+newFormation-oldFormation);
    const double energyScale=std::max({std::abs(initialEnergy),std::abs(finalEnergy),
        std::abs(oldFormation),std::abs(newFormation),1.});
    if(!std::isfinite(finalEnergy)||result.energyDefect>100.*(energyAtol_+rtol_*energyScale))
    {result.reason="CVODE sensible/formation energy balance failed";return result;}
    ConservedVars trial=u;trial.rhoE=u.rho*(finalEnergy+kinetic);
    for(std::size_t s=0;s<names_.size();++s)trial.rhoY[s]=u.rho*Y_[s];
    const auto recovered=recoverTetraState(trial,w,q,model_);
    if(!recovered.success()||!gasState(recovered.state->Q))
    {result.reason="CVODE final EOS commit failed: "+recovered.diagnostics.backendStatus;return result;}
    result.state=*recovered.state;return result;
}
}
#endif
