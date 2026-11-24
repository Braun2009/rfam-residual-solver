#include "CoolPropThermoModel.H"
#include "CompleteCellState.H"
#include "AcousticStateAdmission.H"
#include "TetraAcousticAdmission.H"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace Foam;
int main()
{
    if(!thermo::coolPropThermoAvailable()){std::cerr<<"native CoolProp unavailable\n";return 77;}
    try
    {
        auto model=thermo::makeCoolPropThermoModel();
        struct Case{double y,p,t,rho,a;};
        // Frozen Python CoolProp 7.2.0 references. Two-phase a from independent
        // scipy root solves for rho +/- 1.25e-4 at fixed equilibrium entropy.
        const std::vector<Case> cases{
            {1.,6e6,133.5,909.2854934453625,568.8966254585362},
            {1.,4e6,133.5,892.9560579394473,537.5793502804447},
            {0.,4.5e6,453.5,65.51777992025248,258.65954707402983},
            {.73,4e6,400.,42.17299386884054,346.201034436359},
            {.9,4e6,170.,128.54147156832138,199.16837500145996},
            {.99,4e6,150.,222.62395513215998,131.6227952890004},
            {.73,4e6,250.,81.88388252297813,238.44833670861885},
            {.1,4e6,250.,335.59659173069804,143.2219375353223}};
        for(const auto& c:cases)
        {
            const auto result=model->initializePTY(c.p,c.t,{c.y,1-c.y});
            if(!result.success())throw std::runtime_error("initialize: "+result.diagnostics.backendStatus);
            const auto& q=*result.state;
            if(std::abs(q.rho/c.rho-1.)>1e-8||std::abs(q.acoustic.soundSpeed()/c.a-1.)>2e-6)
                throw std::runtime_error("CoolProp frozen reference mismatch");
            if(q.region==thermo::EquilibriumRegion::TWO_PHASE)
            {
                if(q.acoustic.evidence!=thermo::AcousticEvidence::CONVERGED_RHOE_ISENTROPE
                   ||q.conservedStateProjected||q.thermodynamics
                   ||result.diagnostics.acousticDerivativeConvergenceComparisons<3)
                    throw std::runtime_error("invalid equilibrium acoustic evidence");
                // This backend proof does not globally enable unfinished sector
                // phase boundary/transport/chemistry capabilities.
                if(solverConsumerAdmitsAcousticState(q))throw std::runtime_error("global phase readiness was weakened");
                if(!tetraAdmitsAcoustics(q))throw std::runtime_error("verified HEM state rejected by tetra prototype");
                auto invalid=q;invalid.acoustic.evidence=thermo::AcousticEvidence::UNVERIFIED;
                if(tetraAdmitsAcoustics(invalid))throw std::runtime_error("placeholder HEM sound accepted");
                invalid=q;invalid.conservedStateProjected=true;
                if(tetraAdmitsAcoustics(invalid))throw std::runtime_error("projected HEM state accepted");
            }
            PrimitiveVars w(2);w.p=q.p;w.T=q.T;w.Y=q.overallMassFractions;w.U={2.,3.,40.};
            auto initial=CompleteStateConversion::fromPrimitiveAndEquilibrium(w,q,*model);
            if(!initial.success())throw std::runtime_error("CoolProp complete initialization failed");
            thermo::EquilibriumSeed seed;seed.pressure=c.p*1.001;seed.temperature=c.t*1.001;
            auto roundtrip=model->recoverRhoE(q.rho,q.e,q.overallMassFractions,&seed);
            if(!roundtrip.success())throw std::runtime_error("roundtrip: "+roundtrip.diagnostics.backendStatus);
            if(std::abs(roundtrip.state->p/c.p-1.)>1e-6||std::abs(roundtrip.state->T/c.t-1.)>1e-6)
                throw std::runtime_error("CoolProp PT-rhoe roundtrip mismatch");
            auto converted=CompleteStateConversion::fromConserved(initial.state->U,w,*model,&seed);
            if(!converted.success())throw std::runtime_error("complete recovery: "+converted.diagnostics.backendStatus);
            std::cout<<"case Y="<<c.y<<" T="<<c.t<<" p="<<c.p<<" a="<<q.acoustic.soundSpeed()<<" PASS\n";
        }
        // CFD regression: a ~6e-8 K change previously crossed a numerical
        // stationary-point classification boundary in CoolProp 7.2.0. The
        // density-root bracket patch must preserve the target and full flash.
        thermo::EquilibriumSeed warmSeed;warmSeed.pressure=4302281.193782107;
        warmSeed.temperature=411.25762958063888;
        const auto warm=model->recoverRhoE(44.31367393683184,29291.760254948385,
            {.712849989822318,1.-.712849989822318},&warmSeed);
        if(!warm.success())throw std::runtime_error("cell7766 density bracket regression: "+warm.diagnostics.backendStatus);
        if(std::abs(warm.state->T-411.25762964204574)>1e-6
           ||std::abs(warm.state->p-4302281.195420209)>0.1)
            throw std::runtime_error("density bracket changed the CFD target");
        // A-route step 4 (cell 132114): near-saturated pure liquid oxygen.
        // CoolProp's DmassUmass flash fails on this steep branch; the
        // saturation lever-rule fallback in recoverRhoE must close it.
        thermo::EquilibriumSeed satSeed;satSeed.pressure=2409150.4709759345;
        satSeed.temperature=128.45339349058102;
        satSeed.previousRegion=thermo::EquilibriumRegion::SINGLE_LIQUID;
        // exactly the A-route composition: trace propane (2.7e-11 kg/m3 -> 3e-14 mass
        // fraction) must not force the mixture flash on a pure stream
        const double traceC3=2.7249785543015885e-11/920.66118876027952;
        const auto nearsat=model->recoverRhoE(920.66118876027952,-336706.5,{1.-traceC3,traceC3},&satSeed);
        if(!nearsat.success())
            throw std::runtime_error("near-saturated LOX recovery: "+nearsat.diagnostics.backendStatus);
        if(nearsat.state->T<126.||nearsat.state->T>131.)
            throw std::runtime_error("near-saturated LOX recovery landed at the wrong temperature");
        std::cout<<"near-saturated LOX: T="<<nearsat.state->T<<" p="<<nearsat.state->p
                 <<" rho="<<nearsat.state->rho<<" e="<<nearsat.state->e<<" PASS\n";
        // A-route step 10 (cell 867814): the mixing front has raised the trace
        // propane to ~1.1e-12 mass fraction; the pure-fluid two-phase lever rule
        // must still close the state.
        const double trace2=9.7507170918775226e-10/917.53016241659714;
        const auto traceState=model->recoverRhoE(917.53016241659714,-336725.2,{1.-trace2,trace2},&satSeed);
        if(!traceState.success())
            throw std::runtime_error("trace-composition two-phase recovery: "+traceState.diagnostics.backendStatus);
        std::cout<<"trace-composition recovery: T="<<traceState.state->T<<" p="<<traceState.state->p
                 <<" rho="<<traceState.state->rho<<" PASS\n";
        // A-route step 11 (cell 312054): the trace has grown to 12 ppm propane
        // in LOX; the mixture flash fails there and the last-resort pure retry
        // must close the state (mirrors the run's species-repair behaviour).
        const double trace3=0.01106524398624753/915.95909510918091;
        const auto traceState3=model->recoverRhoE(915.95909510918091,-336691.2,{1.-trace3,trace3},&satSeed);
        if(!traceState3.success())
            throw std::runtime_error("12 ppm trace recovery: "+traceState3.diagnostics.backendStatus);
        std::cout<<"12 ppm trace recovery: T="<<traceState3.state->T<<" p="<<traceState3.state->p
                 <<" rho="<<traceState3.state->rho<<" PASS\n";
        // Pure O2 saturation: rho-u uniquely fixes quality, although PT alone
        // does not. Reference at T=120 K, vapor mass fraction 0.35.
        const auto saturated=model->recoverRhoE(104.47714814463698,-300168.9399434588,{1.,0.});
        if(!saturated.success()||std::abs(saturated.state->T-120.)>1e-6
           ||std::abs(saturated.state->vaporMassFraction-.35)>1e-8
           ||std::abs(saturated.state->acoustic.soundSpeed()-94.08377530088498)>1e-4)
            throw std::runtime_error("pure O2 two-phase rho-e flash mismatch");
        // Query order must not leak the preceding cell's phase/composition.
        int failed=0;
        #pragma omp parallel for reduction(+:failed)
        for(int i=0;i<16;++i)
        {
            const auto& c=cases[i%4];auto r=model->initializePTY(c.p,c.t,{c.y,1-c.y});
            if(!r.success()||std::abs(r.state->rho/c.rho-1.)>1e-8)++failed;
        }
        if(failed)throw std::runtime_error("thread-local CoolProp state contamination");
        if(model->initializePTY(4e6,400.,{-1e-9,1+1e-9}).success()
           ||model->recoverRhoE(-1.,1.,{.5,.5}).success())throw std::runtime_error("invalid input accepted");
        std::cout<<"COOLPROP_BACKEND_PASS\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
