#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/framework/SteadyResidualEvaluator.H"

using namespace Foam;

namespace
{

struct Snapshot
{
    std::vector<ConservedVars> current;
    std::vector<ConservedVars> candidate;
    std::vector<PrimitiveVars> primitive;
    std::vector<thermo::EquilibriumState> equilibrium;
};

bool sameScalar(double left,double right)
{ return left==right||(std::isnan(left)&&std::isnan(right)); }

bool sameVector(const std::vector<double>& left,const std::vector<double>& right)
{
    if(left.size()!=right.size()) return false;
    for(std::size_t entry=0;entry<left.size();++entry)
        if(!sameScalar(left[entry],right[entry])) return false;
    return true;
}

bool same(const ConservedVars& left,const ConservedVars& right)
{
    return sameScalar(left.rho,right.rho)
        &&sameVector(left.rhoU,right.rhoU)
        &&sameScalar(left.rhoE,right.rhoE)
        &&sameVector(left.rhoY,right.rhoY)
        &&sameScalar(left.rhoZ,right.rhoZ)
        &&sameScalar(left.rhoC,right.rhoC);
}

bool same(const PrimitiveVars& left,const PrimitiveVars& right)
{
    return sameScalar(left.p,right.p)
        &&sameScalar(left.T,right.T)
        &&sameVector(left.U,right.U)
        &&sameVector(left.Y,right.Y)
        &&sameScalar(left.Z,right.Z)
        &&sameScalar(left.C,right.C);
}

bool same(const thermo::EquilibriumPhaseProperties& left,
          const thermo::EquilibriumPhaseProperties& right)
{
    return sameVector(left.massFractions,right.massFractions)
        &&sameVector(left.moleFractions,right.moleFractions)
        &&sameScalar(left.molecularWeight,right.molecularWeight)
        &&sameScalar(left.rho,right.rho)
        &&sameScalar(left.h,right.h)
        &&sameScalar(left.e,right.e)
        &&sameScalar(left.s,right.s)
        &&sameScalar(left.compressibility,right.compressibility)
        &&left.entropyAvailable==right.entropyAvailable;
}

bool same(const thermo::EquilibriumThermodynamicPayload& left,
          const thermo::EquilibriumThermodynamicPayload& right)
{
    return left.version==right.version
        &&left.closure==right.closure
        &&sameScalar(left.Cp,right.Cp)
        &&sameScalar(left.Cv,right.Cv)
        &&sameScalar(left.dpdrho_T,right.dpdrho_T)
        &&sameScalar(left.dpdT_rho,right.dpdT_rho)
        &&sameScalar(left.dedrho_T,right.dedrho_T)
        &&sameScalar(left.dedT_rho,right.dedT_rho);
}

bool same(const thermo::EquilibriumState& left,
          const thermo::EquilibriumState& right)
{
    if(left.liquid.has_value()!=right.liquid.has_value()
        ||left.vapor.has_value()!=right.vapor.has_value()
        ||left.singlePhase.has_value()!=right.singlePhase.has_value()
        ||left.thermodynamics.has_value()!=right.thermodynamics.has_value())
        return false;
    return left.region==right.region
        &&sameVector(left.overallMassFractions,right.overallMassFractions)
        &&sameScalar(left.p,right.p)
        &&sameScalar(left.T,right.T)
        &&sameScalar(left.rho,right.rho)
        &&sameScalar(left.h,right.h)
        &&sameScalar(left.e,right.e)
        &&sameScalar(left.s,right.s)
        &&sameScalar(left.vaporMoleFraction,right.vaporMoleFraction)
        &&sameScalar(left.vaporMassFraction,right.vaporMassFraction)
        &&sameScalar(left.vaporVolumeFraction,right.vaporVolumeFraction)
        &&left.acoustic.closure==right.acoustic.closure
        &&sameScalar(left.acoustic.pressureDensityDerivative,
            right.acoustic.pressureDensityDerivative)
        &&left.entropyAvailable==right.entropyAvailable
        &&(!left.thermodynamics
            ||same(*left.thermodynamics,*right.thermodynamics))
        &&(!left.liquid||same(*left.liquid,*right.liquid))
        &&(!left.vapor||same(*left.vapor,*right.vapor))
        &&(!left.singlePhase||same(*left.singlePhase,*right.singlePhase));
}

Snapshot capture(const Domain& domain)
{
    Snapshot result;
    result.current.reserve(static_cast<std::size_t>(domain.N()));
    result.candidate.reserve(static_cast<std::size_t>(domain.N()));
    result.primitive.reserve(static_cast<std::size_t>(domain.N()));
    result.equilibrium.reserve(static_cast<std::size_t>(domain.N()));
    for(int cell=0;cell<domain.N();++cell)
    {
        result.current.push_back(domain.U(cell));
        result.candidate.push_back(domain.Unew(cell));
        result.primitive.push_back(domain.W(cell));
        result.equilibrium.push_back(domain.Q(cell));
    }
    return result;
}

bool sameState(const Domain& domain,const Snapshot& snapshot)
{
    if(domain.N()!=static_cast<int>(snapshot.current.size())) return false;
    for(int cell=0;cell<domain.N();++cell)
        if(!same(domain.U(cell),snapshot.current[static_cast<std::size_t>(cell)])
            ||!same(domain.Unew(cell),
                snapshot.candidate[static_cast<std::size_t>(cell)])
            ||!same(domain.W(cell),
                snapshot.primitive[static_cast<std::size_t>(cell)])
            ||!same(domain.Q(cell),
                snapshot.equilibrium[static_cast<std::size_t>(cell)]))
            return false;
    return true;
}

bool initialize(Domain& domain)
{
    for(int cell=0;cell<domain.N();++cell)
    {
        PrimitiveVars primitive(2);
        primitive.p=1.0e6;
        primitive.T=600.0+static_cast<double>(cell);
        primitive.U={50.0,0.0,0.0};
        primitive.Y={0.8,0.2};
        primitive.Z=0.2;
        primitive.C=0.1;
        domain.stagedW(cell)=primitive;
    }
    return domain.primitiveToConserved();
}

bool check(bool condition,const std::string& message)
{
    if(!condition) std::cerr << message << "\n";
    return condition;
}

class StatefulBoundarySensitiveOperator : public Operator
{
    mutable double state_=7.0;
    mutable double saved_=7.0;

public:
    bool ownsPersistentState() const override { return true; }
    void beginStepTransaction(const Domain&) const override
    { saved_=state_; }
    void rollbackStepTransaction() const noexcept override
    { state_=saved_; }
    void acceptStepTransaction() const noexcept override {}

    void apply(Domain& domain,double dt) const override
    {
        state_+=3.0;
        const double rate=1.0e-3*domain.W(0).T;
        for(int cell=0;cell<domain.N();++cell)
        {
            domain.Unew(cell).rhoZ+=dt*rate;
            domain.Unew(cell).rhoC+=dt*(rate+state_);
        }
    }

    double state() const { return state_; }
    std::string name() const override
    { return "stateful boundary-sensitive residual"; }
};

class StatefulThrowOperator : public Operator
{
    mutable int state_=4;
    mutable int saved_=4;

public:
    bool ownsPersistentState() const override { return true; }
    void beginStepTransaction(const Domain&) const override
    { saved_=state_; }
    void rollbackStepTransaction() const noexcept override
    { state_=saved_; }
    void apply(Domain& domain,double) const override
    {
        state_=99;
        domain.Unew(2).rhoE+=123.0;
        throw std::runtime_error("injected residual failure");
    }
    int state() const { return state_; }
    std::string name() const override
    { return "stateful throwing residual"; }
};

bool sameResidual(
    const SteadyResidualEvaluator::Result& left,
    const SteadyResidualEvaluator::Result& right)
{
    if(left.referenceStep!=right.referenceStep
        ||left.active!=right.active
        ||left.field.size()!=right.field.size()) return false;
    for(std::size_t cell=0;cell<left.field.size();++cell)
        if(!same(left.field[cell],right.field[cell])) return false;
    return true;
}

} // namespace

int main()
{
    Domain domain(6,0.06,2);
    if(!initialize(domain)) return 1;
    const Snapshot before=capture(domain);
    const std::uint64_t acceptedBefore=domain.acceptedSteps();
    const std::uint64_t rejectedBefore=domain.rejectedSteps();

    std::vector<std::unique_ptr<Operator>> operators;
    auto stateful=std::make_unique<StatefulBoundarySensitiveOperator>();
    StatefulBoundarySensitiveOperator* statefulPtr=stateful.get();
    operators.push_back(std::move(stateful));
    const double operatorStateBefore=statefulPtr->state();
    int boundaryCalls=0;
    auto boundaryHook=[&](Domain& trial)
    {
        ++boundaryCalls;
        PrimitiveVars primitive=trial.W(0);
        primitive.T=750.0;
        const thermo::EquilibriumResult initialized=
            trial.thermoModel().initializePTY(
                primitive.p,primitive.T,primitive.Y);
        if(!initialized.success()) return false;
        std::vector<Domain::CompleteCellReplacement> replacements;
        replacements.push_back(Domain::CompleteCellReplacement{
            0,0,0,primitive,*initialized.state});
        return trial.replaceCompleteCellsInStepTransaction(replacements);
    };

    SteadyResidualEvaluator::Result first,second;
    const bool firstOk=SteadyResidualEvaluator::evaluate(
        domain,operators,1.0e-6,first,boundaryHook,true);
    const bool firstRestored=sameState(domain,before)
        &&statefulPtr->state()==operatorStateBefore;
    const bool secondOk=SteadyResidualEvaluator::evaluate(
        domain,operators,1.0e-6,second,boundaryHook,true);
    const bool repeated=sameResidual(first,second);
    const bool boundaryUsed=first.field.size()==6
        &&first.active==std::vector<unsigned char>({0U,1U,1U,1U,1U,0U})
        &&std::abs(first.field[1].rhoZ-0.75)<1.0e-9
        &&first.field[0].rhoZ==0.0
        &&first.field[5].rhoZ==0.0;
    if(!boundaryUsed)
    {
        std::cerr << "boundary residual diagnostics size/active/rhoZ="
                  << first.field.size() << "/";
        for(unsigned char value:first.active)
            std::cerr << static_cast<int>(value);
        std::cerr << "/";
        for(const ConservedVars& value:first.field)
            std::cerr << value.rhoZ << ",";
        std::cerr << "\n";
    }
    if(!check(firstOk&&secondOk,"residual evaluation failed")
        ||!check(firstRestored&&sameState(domain,before),
            "successful residual probe did not restore U/Unew/W/Q/operator state")
        ||!check(repeated,"repeated residual evaluations are not identical")
        ||!check(boundaryCalls==2&&boundaryUsed,
            "transactional boundary state was not used by residual assembly")
        ||!check(domain.acceptedSteps()==acceptedBefore
                &&domain.rejectedSteps()==rejectedBefore,
            "successful residual probes polluted the physical step ledger"))
        return 1;

    std::vector<std::unique_ptr<Operator>> throwingOperators;
    auto throwing=std::make_unique<StatefulThrowOperator>();
    StatefulThrowOperator* throwingPtr=throwing.get();
    throwingOperators.push_back(std::move(throwing));
    SteadyResidualEvaluator::Result failed;
    const bool rejected=!SteadyResidualEvaluator::evaluate(
        domain,throwingOperators,1.0e-6,failed,{},true);
    const bool failureMetadata=domain.failure().code
        ==Domain::FailureCode::OPERATOR_EXCEPTION
        &&domain.failure().component=="stateful throwing residual";
    const bool failureRestored=sameState(domain,before)
        &&throwingPtr->state()==4;
    const bool failureLedgerClean=
        domain.acceptedSteps()==acceptedBefore
        &&domain.rejectedSteps()==rejectedBefore;
    domain.resetTrip();

    if(!check(rejected,"failing residual evaluation was accepted")
        ||!check(failureMetadata,
            "failing residual evaluation lost failure metadata")
        ||!check(failureRestored,
            "failing residual evaluation did not restore state atomically")
        ||!check(failureLedgerClean,
            "failing residual evaluation polluted the physical step ledger")
        ||!check(sameState(domain,before),
            "failure diagnostic reset changed the authoritative state"))
        return 1;

    std::cout << "boundary_calls=" << boundaryCalls
              << " active_residual_rhoZ=" << first.field[1].rhoZ
              << " operator_state=" << statefulPtr->state() << "\n"
              << "STEADY_RESIDUAL_EVALUATOR_PASS\n";
    return 0;
}
