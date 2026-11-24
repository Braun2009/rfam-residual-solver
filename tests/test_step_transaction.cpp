#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/cht/CHTOperator2D.H"
#include "../src/framework/PointImplicit.H"
#include "../src/framework/ProfiledOperator.H"
#include "../src/framework/SolverBase.H"
#include "../src/turbulence/TurbulenceOperator.H"

using namespace Foam;

namespace
{

struct DomainState
{
    std::vector<ConservedVars> current;
    std::vector<ConservedVars> candidate;
    std::vector<PrimitiveVars> primitive;
    std::vector<thermo::EquilibriumState> equilibrium;
};

bool sameScalar(double a, double b)
{
    return a == b || (std::isnan(a) && std::isnan(b));
}

bool sameVector(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!sameScalar(a[i], b[i])) return false;
    return true;
}

bool same(
    const thermo::EquilibriumPhaseProperties& a,
    const thermo::EquilibriumPhaseProperties& b)
{
    return sameVector(a.massFractions, b.massFractions)
        && sameVector(a.moleFractions, b.moleFractions)
        && sameScalar(a.molecularWeight, b.molecularWeight)
        && sameScalar(a.rho, b.rho)
        && sameScalar(a.h, b.h)
        && sameScalar(a.e, b.e)
        && sameScalar(a.s, b.s)
        && sameScalar(a.compressibility, b.compressibility)
        && a.entropyAvailable == b.entropyAvailable;
}

bool same(
    const thermo::EquilibriumThermodynamicPayload& a,
    const thermo::EquilibriumThermodynamicPayload& b)
{
    return a.version == b.version
        && a.closure == b.closure
        && sameScalar(a.Cp, b.Cp)
        && sameScalar(a.Cv, b.Cv)
        && sameScalar(a.dpdrho_T, b.dpdrho_T)
        && sameScalar(a.dpdT_rho, b.dpdT_rho)
        && sameScalar(a.dedrho_T, b.dedrho_T)
        && sameScalar(a.dedT_rho, b.dedT_rho);
}

bool same(const thermo::EquilibriumState& a, const thermo::EquilibriumState& b)
{
    const bool phases = a.liquid.has_value() == b.liquid.has_value()
        && a.vapor.has_value() == b.vapor.has_value()
        && a.singlePhase.has_value() == b.singlePhase.has_value()
        && a.thermodynamics.has_value() == b.thermodynamics.has_value();
    if (!phases) return false;
    return a.region == b.region
        && sameVector(a.overallMassFractions, b.overallMassFractions)
        && sameScalar(a.p, b.p)
        && sameScalar(a.T, b.T)
        && sameScalar(a.rho, b.rho)
        && sameScalar(a.h, b.h)
        && sameScalar(a.e, b.e)
        && sameScalar(a.s, b.s)
        && sameScalar(a.vaporMoleFraction, b.vaporMoleFraction)
        && sameScalar(a.vaporMassFraction, b.vaporMassFraction)
        && sameScalar(a.vaporVolumeFraction, b.vaporVolumeFraction)
        && a.acoustic.closure == b.acoustic.closure
        && sameScalar(a.acoustic.pressureDensityDerivative,
            b.acoustic.pressureDensityDerivative)
        && a.entropyAvailable == b.entropyAvailable
        && (!a.thermodynamics
            || same(*a.thermodynamics, *b.thermodynamics))
        && (!a.liquid || same(*a.liquid, *b.liquid))
        && (!a.vapor || same(*a.vapor, *b.vapor))
        && (!a.singlePhase || same(*a.singlePhase, *b.singlePhase));
}

bool same(const ConservedVars& a, const ConservedVars& b)
{
    return a.rho == b.rho && a.rhoU == b.rhoU && a.rhoE == b.rhoE
        && a.rhoY == b.rhoY && a.rhoZ == b.rhoZ && a.rhoC == b.rhoC;
}

bool same(const PrimitiveVars& a, const PrimitiveVars& b)
{
    return a.p == b.p && a.T == b.T && a.U == b.U && a.Y == b.Y
        && a.Z == b.Z && a.C == b.C;
}

DomainState capture(const Domain& domain)
{
    DomainState state;
    state.current.reserve(domain.N());
    state.candidate.reserve(domain.N());
    state.primitive.reserve(domain.N());
    state.equilibrium.reserve(domain.N());
    for (int i = 0; i < domain.N(); ++i)
    {
        state.current.push_back(domain.U(i));
        state.candidate.push_back(domain.Unew(i));
        state.primitive.push_back(domain.W(i));
        state.equilibrium.push_back(domain.Q(i));
    }
    return state;
}

bool sameState(const Domain& domain, const DomainState& expected)
{
    if (domain.N() != static_cast<int>(expected.current.size())) return false;
    for (int i = 0; i < domain.N(); ++i)
        if (!same(domain.U(i), expected.current[i])
            || !same(domain.Unew(i), expected.candidate[i])
            || !same(domain.W(i), expected.primitive[i])
            || !same(domain.Q(i), expected.equilibrium[i]))
            return false;
    return true;
}

bool initialize(Domain& domain, double temperature = 600.0)
{
    for (int i = 0; i < domain.N(); ++i)
    {
        PrimitiveVars state(2);
        state.p = 1.0e6;
        state.T = temperature;
        state.U = {100.0, 0.0, 0.0};
        state.Y = {0.8, 0.2};
        state.Z = 0.2;
        state.C = 0.1;
        domain.stagedW(i) = state;
    }
    return domain.primitiveToConserved();
}

bool check(bool condition, const std::string& message)
{
    if (!condition) std::cerr << message << '\n';
    return condition;
}

class InvalidCandidateOperator : public Operator
{
public:
    void apply(Domain& domain, double) const override
    {
        domain.Unew(0).rho = -1.0;
    }

    std::string name() const override { return "invalid candidate"; }
};

class MutateThenThrowOperator : public Operator
{
public:
    void apply(Domain& domain, double) const override
    {
        domain.Unew(0).rhoE += 123.0;
        throw std::runtime_error("injected operator failure");
    }

    std::string name() const override { return "mutate then throw"; }
};

class StageTwoRejectOperator : public Operator
{
    mutable int calls_ = 0;

public:
    void apply(Domain& domain, double dt) const override
    {
        ++calls_;
        if (calls_ == 1)
        {
            domain.Unew(0).rhoZ += dt;
            domain.Unew(0).rhoC += 2.0*dt;
        }
        else
        {
            domain.Unew(0).rho = -1.0;
        }
    }

    int calls() const { return calls_; }
    std::string name() const override { return "RK2 stage-two reject"; }
};

class ScalarRateOperator : public Operator
{
    double zRate_;
    double cRate_;

public:
    ScalarRateOperator(double zRate, double cRate)
        : zRate_(zRate), cRate_(cRate)
    {}

    void apply(Domain& domain, double dt) const override
    {
        domain.Unew(0).rhoZ += dt*zRate_;
        domain.Unew(0).rhoC += dt*cRate_;
    }

    std::string name() const override { return "scalar rate"; }
};

class SwitchableThrowOperator : public Operator
{
    bool enabled_ = false;

public:
    void setEnabled(bool enabled) { enabled_ = enabled; }

    void apply(Domain&, double) const override
    {
        if (enabled_) throw std::runtime_error("injected downstream failure");
    }

    std::string name() const override { return "switchable downstream failure"; }
};

class PersistentNoopOperator : public Operator
{
public:
    bool ownsPersistentState() const override { return true; }
    void apply(Domain&, double) const override {}
    std::string name() const override { return "persistent no-op"; }
};

class AverageThrowOperator : public Operator
{
public:
    bool ownsPersistentState() const override { return true; }
    void apply(Domain&, double) const override {}
    void averageWithStepStartState() const override
    {
        throw std::runtime_error("injected RK2 averaging failure");
    }
    std::string name() const override { return "average throw"; }
};

bool testEulerCommitRejection()
{
    Domain domain(4, 0.04, 2);
    if (!initialize(domain)) return check(false, "Euler setup failed");
    const DomainState before = capture(domain);

    SolverBase solver(domain);
    solver.addOperator(std::unique_ptr<Operator>(new InvalidCandidateOperator()));

    const bool accepted = solver.step(1.0e-6);
    return check(!accepted, "Euler accepted an invalid candidate")
        && check(solver.time() == 0.0 && solver.iteration() == 0,
            "Euler rejection advanced solver time or iteration")
        && check(domain.acceptedSteps() == 0 && domain.rejectedSteps() == 1,
            "Euler rejection ledger is incorrect")
        && check(domain.failure().code
                == Domain::FailureCode::INVALID_CONSERVED_STATE,
            "Euler rejection did not preserve the failure code")
        && check(sameState(domain, before),
            "Euler rejection did not restore Domain state");
}

bool testOperatorExceptionRollback()
{
    Domain domain(4, 0.04, 2);
    if (!initialize(domain)) return check(false, "exception setup failed");
    const DomainState before = capture(domain);

    SolverBase solver(domain);
    solver.addOperator(std::unique_ptr<Operator>(new MutateThenThrowOperator()));

    bool threw = false;
    try
    {
        (void)solver.step(1.0e-6);
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }

    return check(threw, "operator exception was not propagated")
        && check(solver.time() == 0.0 && solver.iteration() == 0,
            "operator exception advanced solver time or iteration")
        && check(domain.failure().code == Domain::FailureCode::OPERATOR_EXCEPTION,
            "operator exception failure metadata is incorrect")
        && check(domain.rejectedSteps() == 1,
            "operator exception was not recorded as a rejected step")
        && check(sameState(domain, before),
            "operator exception did not restore Domain state");
}

bool testRK2StageTwoRejection()
{
    Domain domain(4, 0.04, 2);
    if (!initialize(domain)) return check(false, "RK2 rejection setup failed");
    const DomainState before = capture(domain);

    std::unique_ptr<StageTwoRejectOperator> rejector(
        new StageTwoRejectOperator());
    StageTwoRejectOperator* rejectorPtr = rejector.get();

    SolverBase solver(domain);
    solver.setIntegrator(std::unique_ptr<TimeIntegrator>(new RK2Integrator()));
    solver.addOperator(std::move(rejector));

    const bool accepted = solver.step(1.0e-4);
    return check(!accepted, "RK2 accepted an invalid stage-two candidate")
        && check(rejectorPtr->calls() == 2,
            "RK2 rejection was not injected in stage two")
        && check(solver.time() == 0.0 && solver.iteration() == 0,
            "RK2 stage-two rejection advanced the clock")
        && check(domain.rejectedSteps() == 1 && domain.acceptedSteps() == 0,
            "RK2 rejection ledger is incorrect")
        && check(sameState(domain, before),
            "RK2 stage-two rejection did not restore the pre-step state");
}

bool testRK2ScalarAverage()
{
    Domain domain(4, 0.04, 2);
    if (!initialize(domain)) return check(false, "RK2 scalar setup failed");

    const double z0 = domain.U(0).rhoZ;
    const double c0 = domain.U(0).rhoC;
    const double dt = 1.0e-4;
    const double zRate = 2.0;
    const double cRate = 3.0;

    SolverBase solver(domain);
    solver.setIntegrator(std::unique_ptr<TimeIntegrator>(new RK2Integrator()));
    solver.addOperator(std::unique_ptr<Operator>(
        new ScalarRateOperator(zRate, cRate)));

    const bool accepted = solver.step(dt);
    const double tolerance = 1.0e-12;
    return check(accepted, "RK2 scalar step was rejected")
        && check(std::abs(domain.U(0).rhoZ - (z0 + dt*zRate)) < tolerance,
            "RK2 rhoZ did not use the transaction-start average")
        && check(std::abs(domain.U(0).rhoC - (c0 + dt*cRate)) < tolerance,
            "RK2 rhoC did not use the transaction-start average")
        && check(domain.acceptedSteps() == 1 && domain.rejectedSteps() == 0,
            "RK2 accepted-step ledger is incorrect")
        && check(std::abs(solver.time() - dt) < tolerance
                 && solver.iteration() == 1,
            "RK2 accepted step did not advance the clock once")
        && check(solver.hasPreviousAcceptedDt()
                 && solver.previousAcceptedDt().has_value()
                 && *solver.previousAcceptedDt() == dt,
            "RK2 accepted step did not publish its previous accepted dt");
}

bool testRestartClockAtomicAdmission()
{
    Domain domain(4, 0.04, 2);
    if (!initialize(domain)) return check(false, "restart clock setup failed");

    SolverBase solver(domain);
    solver.setRestartState(1.25e-5, 42, 1.0e-8);
    if (!check(solver.time() == 1.25e-5 && solver.iteration() == 42
                 && solver.hasPreviousAcceptedDt()
                 && solver.previousAcceptedDt().has_value()
                 && *solver.previousAcceptedDt() == 1.0e-8,
            "valid restart clock publishes all fields"))
        return false;

    const double timeBefore = solver.time();
    const int iterationBefore = solver.iteration();
    const std::optional<double> dtBefore = solver.previousAcceptedDt();
    auto unchanged = [&]()
    {
        return solver.time() == timeBefore
            && solver.iteration() == iterationBefore
            && solver.previousAcceptedDt() == dtBefore;
    };

    bool threw = false;
    try { solver.setRestartState(std::numeric_limits<double>::quiet_NaN(), 9, 1.0e-8); }
    catch (const std::invalid_argument&) { threw = true; }
    if (!check(threw && unchanged(), "NaN restart time is atomic fail-close"))
        return false;

    threw = false;
    try { solver.setRestartState(2.0, -1, 1.0e-8); }
    catch (const std::invalid_argument&) { threw = true; }
    if (!check(threw && unchanged(), "negative restart iteration is atomic fail-close"))
        return false;

    threw = false;
    try { solver.setRestartState(2.0, 43, 0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    if (!check(threw && unchanged(), "nonpositive previous dt is atomic fail-close"))
        return false;

    for (double invalidDt : {
            0.0,
            -1.0,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::infinity()})
    {
        threw = false;
        try { (void)solver.step(invalidDt); }
        catch (const std::invalid_argument&) { threw = true; }
        if (!check(threw && unchanged(),
                "invalid proposed dt is rejected before clock publication"))
            return false;
    }

    solver.addOperator(std::unique_ptr<Operator>(
        new InvalidCandidateOperator()));
    if (!check(!solver.step(1.0e-6) && unchanged(),
            "rejected physical step preserves previous accepted dt"))
        return false;

    solver.setRestartState(0.0, 0);
    return check(!solver.hasPreviousAcceptedDt()
                 && !solver.previousAcceptedDt().has_value(),
        "legacy restart admission explicitly clears absent previous dt");
}

bool testPointImplicitFailClose()
{
    Domain domain(4, 0.04, 2);
    if (!initialize(domain))
        return check(false, "point-implicit setup failed");
    const DomainState before = capture(domain);

    SolverBase solver(domain);
    solver.setIntegrator(std::unique_ptr<TimeIntegrator>(
        new PointImplicitIntegrator()));
    solver.addOperator(std::unique_ptr<Operator>(new PersistentNoopOperator()));

    const bool accepted = solver.step(1.0e-4);
    return check(!accepted,
            "point-implicit accepted operator-owned persistent state")
        && check(solver.time() == 0.0 && solver.iteration() == 0,
            "point-implicit fail-close advanced the solver clock")
        && check(domain.failure().code
                == Domain::FailureCode::TRANSACTION_PROTOCOL,
            "point-implicit fail-close metadata is incorrect")
        && check(domain.acceptedSteps() == 0 && domain.rejectedSteps() == 1,
            "point-implicit fail-close ledger is incorrect")
        && check(sameState(domain, before),
            "point-implicit fail-close changed Domain state");
}

bool testRK2AverageExceptionMetadata()
{
    Domain domain(4, 0.04, 2);
    if (!initialize(domain))
        return check(false, "RK2 average-exception setup failed");
    const DomainState before = capture(domain);

    SolverBase solver(domain);
    solver.setIntegrator(std::unique_ptr<TimeIntegrator>(new RK2Integrator()));
    solver.addOperator(std::unique_ptr<Operator>(new AverageThrowOperator()));

    bool threw = false;
    try
    {
        (void)solver.step(1.0e-4);
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }

    return check(threw, "RK2 state-average exception was not propagated")
        && check(domain.failure().code
                == Domain::FailureCode::TRANSACTION_PROTOCOL
                 && domain.failure().component == "average throw",
            "RK2 state-average exception metadata is incorrect")
        && check(domain.rejectedSteps() == 1 && domain.acceptedSteps() == 0,
            "RK2 state-average exception ledger is incorrect")
        && check(solver.time() == 0.0 && solver.iteration() == 0,
            "RK2 state-average exception advanced the clock")
        && check(sameState(domain, before),
            "RK2 state-average exception did not restore Domain state");
}

bool testTurbulenceRollbackThroughProfile()
{
    Domain domain(5, 5, 0.05, 0.05, 2);
    if (!initialize(domain)) return check(false, "SST setup failed");

    std::unique_ptr<TurbulenceOperator> turbulence(
        new TurbulenceOperator(1.0e-5));
    TurbulenceOperator* turbulencePtr = turbulence.get();
    std::unique_ptr<Operator> profiled(
        new ProfiledOperator(std::move(turbulence)));
    std::unique_ptr<SwitchableThrowOperator> thrower(
        new SwitchableThrowOperator());
    SwitchableThrowOperator* throwerPtr = thrower.get();

    SolverBase solver(domain);
    solver.addOperator(std::move(profiled));
    solver.addOperator(std::move(thrower));

    const double dt = 1.0e-6;
    if (!solver.step(dt)) return check(false, "SST baseline step was rejected");

    const double kBefore = turbulencePtr->k(2, 2, domain.Nx());
    const double omegaBefore = turbulencePtr->omega(2, 2, domain.Nx());
    if (!check(kBefore != 37.5, "SST baseline did not evolve k")) return false;
    const DomainState domainBefore = capture(domain);
    const double timeBefore = solver.time();
    const int iterationBefore = solver.iteration();

    throwerPtr->setEnabled(true);
    bool threw = false;
    try
    {
        (void)solver.step(dt);
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }

    return check(threw, "downstream SST failure was not propagated")
        && check(turbulencePtr->k(2, 2, domain.Nx()) == kBefore
                 && turbulencePtr->omega(2, 2, domain.Nx()) == omegaBefore,
            "SST k/omega advanced across a rejected step")
        && check(solver.time() == timeBefore
                 && solver.iteration() == iterationBefore,
            "SST rejected step advanced the solver clock")
        && check(domain.acceptedSteps() == 1 && domain.rejectedSteps() == 1,
            "SST transaction ledger is incorrect")
        && check(sameState(domain, domainBefore),
            "SST downstream failure did not restore Domain state");
}

bool testTransactionalBoundaryRefreshRollback()
{
    Domain domain(4,0.04,2);
    if(!initialize(domain))
        return check(false,"transactional boundary setup failed");
    const DomainState before=capture(domain);
    if(!domain.beginStepTransaction())
        return check(false,"transactional boundary begin failed");

    PrimitiveVars replacement=domain.W(0);
    replacement.T+=25.0;
    const thermo::EquilibriumResult initialized=
        domain.thermoModel().initializePTY(
            replacement.p,replacement.T,replacement.Y);
    if(!initialized.success())
        return check(false,"transactional boundary PTY setup failed");
    std::vector<Domain::CompleteCellReplacement> replacements;
    replacements.push_back(Domain::CompleteCellReplacement{
        0,0,0,replacement,*initialized.state});
    const bool published=
        domain.replaceCompleteCellsInStepTransaction(replacements);
    const bool changed=published&&domain.W(0).T==replacement.T;
    domain.rollbackStepTransaction();

    return check(published&&changed,
            "transactional boundary refresh was not published")
        &&check(domain.acceptedSteps()==0&&domain.rejectedSteps()==1,
            "transactional boundary rollback ledger is incorrect")
        &&check(sameState(domain,before),
            "transactional boundary rollback did not restore U/W/Q");
}

bool testTransactionalConservedFieldReplacement()
{
    Domain domain(4,0.04,2);
    if(!initialize(domain))
        return check(false,"transactional field setup failed");
    const DomainState before=capture(domain);

    std::vector<ConservedVars> trial;
    trial.reserve(static_cast<std::size_t>(domain.N()));
    for(int cell=0;cell<domain.N();++cell) trial.push_back(domain.U(cell));
    trial[1].rhoU[0]+=0.01*trial[1].rho;

    if(!domain.beginStepTransaction())
        return check(false,"transactional field begin failed");
    const bool published=
        domain.replaceConservedFieldInStepTransaction(trial);
    const bool changed=published
        &&domain.U(1).rhoU[0]==trial[1].rhoU[0]
        &&domain.W(1).U[0]!=before.primitive[1].U[0];
    const bool restored=domain.restoreStepTransactionStartState();
    const bool transactionRemainsActive=domain.stepTransactionActive();
    const bool exactRestore=sameState(domain,before);
    domain.rollbackStepTransaction(false);
    if(!check(published&&changed,
            "transactional conserved trial was not published")
        ||!check(restored&&transactionRemainsActive&&exactRestore,
            "in-transaction trial restore did not preserve the open snapshot")
        ||!check(domain.acceptedSteps()==0&&domain.rejectedSteps()==0,
            "matrix-free trial rollback polluted the physical step ledger")
        ||!check(sameState(domain,before),
            "transactional conserved rollback did not restore U/Unew/W/Q"))
        return false;

    trial[2].rho=-1.0;
    if(!domain.beginStepTransaction())
        return check(false,"invalid transactional field begin failed");
    const bool rejected=
        !domain.replaceConservedFieldInStepTransaction(trial);
    const bool atomic=sameState(domain,before);
    const bool metadata=domain.failure().code
        ==Domain::FailureCode::STATE_RECOVERY_FAILED
        &&domain.failure().cell==2;
    domain.rollbackStepTransaction(false);
    domain.resetTrip();

    return check(rejected&&atomic,
            "invalid transactional conserved field was partially published")
        &&check(metadata,
            "invalid transactional conserved field metadata is incorrect")
        &&check(domain.acceptedSteps()==0&&domain.rejectedSteps()==0,
            "transactional conserved trial polluted the physical ledger")
        &&check(sameState(domain,before),
            "invalid transactional conserved rollback did not restore state");
}

bool testCHTRollback()
{
    Domain domain(4, 3, 0.04, 0.03, 2);
    if (!initialize(domain, 800.0)) return check(false, "CHT setup failed");

    SpeciesRegistry registry;
    std::unique_ptr<CHTOperator2D> cht(new CHTOperator2D(
        6, 0.003, 8000.0, 500.0, 15.0, 5000.0,
        SolidBC::fixedT(300.0), registry, 300.0));
    CHTOperator2D* chtPtr = cht.get();
    std::unique_ptr<SwitchableThrowOperator> thrower(
        new SwitchableThrowOperator());
    SwitchableThrowOperator* throwerPtr = thrower.get();

    SolverBase solver(domain);
    solver.addOperator(std::move(cht));
    solver.addOperator(std::move(thrower));

    const double dt = 1.0e-7;
    if (!solver.step(dt)) return check(false, "CHT baseline step was rejected");

    const std::vector<double> wallBefore = chtPtr->wallTProfile();
    const double heatInBefore = chtPtr->heatToSolid();
    const double heatOutBefore = chtPtr->heatOutOuter();
    const DomainState domainBefore = capture(domain);
    const double timeBefore = solver.time();
    const int iterationBefore = solver.iteration();

    if (!check(!wallBefore.empty() && heatInBefore != 0.0,
        "CHT baseline did not initialize its persistent state"))
        return false;

    throwerPtr->setEnabled(true);
    bool threw = false;
    try
    {
        (void)solver.step(dt);
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }

    return check(threw, "downstream CHT failure was not propagated")
        && check(chtPtr->wallTProfile() == wallBefore,
            "CHT solid temperature advanced across a rejected step")
        && check(chtPtr->heatToSolid() == heatInBefore
                 && chtPtr->heatOutOuter() == heatOutBefore,
            "CHT heat ledger advanced across a rejected step")
        && check(solver.time() == timeBefore
                 && solver.iteration() == iterationBefore,
            "CHT rejected step advanced the solver clock")
        && check(domain.acceptedSteps() == 1 && domain.rejectedSteps() == 1,
            "CHT transaction ledger is incorrect")
        && check(sameState(domain, domainBefore),
            "CHT downstream failure did not restore Domain state");
}

} // namespace

int main()
{
    if (!testEulerCommitRejection()) return 1;
    if (!testOperatorExceptionRollback()) return 1;
    if (!testRK2StageTwoRejection()) return 1;
    if (!testRK2ScalarAverage()) return 1;
    if (!testRestartClockAtomicAdmission()) return 1;
    if (!testPointImplicitFailClose()) return 1;
    if (!testRK2AverageExceptionMetadata()) return 1;
    if (!testTurbulenceRollbackThroughProfile()) return 1;
    if (!testTransactionalBoundaryRefreshRollback()) return 1;
    if (!testTransactionalConservedFieldReplacement()) return 1;
    if (!testCHTRollback()) return 1;

    std::cout << "STEP_TRANSACTION_PASS\n";
    return 0;
}
