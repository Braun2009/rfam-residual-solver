#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "FiniteRateReactor.H"
#include "FiniteRate/RealFluidCVODEWrapper.H"
#include "RealFluidState.H"
#include "SpeciesRegistry.H"
#include "TwoStepCO.H"

using namespace Foam;

namespace
{

int failures = 0;

void check(const std::string& name, bool ok)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name << '\n';
    if (!ok) ++failures;
}

struct ReactorState
{
    std::vector<double> Y;
    double rho;
    double energy;
    double temperature;
    double pressure;
};

bool advance(
    FiniteRateReactor& reactor,
    ReactorState& state,
    double totalTime,
    int steps)
{
    const double dt = totalTime/steps;
    for (int step = 0; step < steps; ++step)
    {
        double deltaRhoE = 0.0;
        if (!reactor.reactConservative(
            state.rho,state.energy,state.Y,dt,deltaRhoE,
            state.temperature,state.pressure))
        {
            std::cout << "  reactor failure code=" << reactor.lastFailureCode()
                      << " rejectedChem=" << reactor.lastRejectedChemIndex()
                      << " before=" << reactor.lastRejectedBefore()
                      << " after=" << reactor.lastRejectedAfter() << '\n';
            return false;
        }
        state.energy += deltaRhoE/state.rho;
    }
    return true;
}

double stateDifference(const ReactorState& lhs, const ReactorState& rhs)
{
    double difference = std::abs(lhs.energy-rhs.energy)
        /std::max({std::abs(lhs.energy),std::abs(rhs.energy),1.0});
    for (size_t i = 0; i < lhs.Y.size(); ++i)
        difference = std::max(
            difference,std::abs(lhs.Y[i]-rhs.Y[i]));
    return difference;
}

bool advanceWithTolerances(
    ChemistryModel* chemistry,
    const SpeciesRegistry& registry,
    const SpeciesMapper& mapper,
    const ReactorState& initial,
    double totalTime,
    double relativeTolerance,
    double speciesAbsoluteTolerance,
    double energyAbsoluteTolerance,
    ReactorState& result)
{
    RealFluidCVODEWrapper solver(chemistry,registry,mapper);
    solver.setTolerances(
        relativeTolerance,
        speciesAbsoluteTolerance,
        energyAbsoluteTolerance);
    result=initial;
    solver.init(
        0.0,result.rho,result.Y,result.energy,result.temperature);
    return solver.advanceTo(
        totalTime,result.Y,result.energy,result.temperature,result.pressure);
}

} // namespace

int main()
{
    std::cout << "=== Conservative real-fluid chemistry energy coupling ===\n";
    SpeciesRegistry registry;
    TwoStepCO chemistry;
    FiniteRateReactor reactor(&chemistry,registry);
    bool invalidEnergyToleranceRejected=false;
    try
    {
        RealFluidCVODEWrapper invalidToleranceSolver(
            &chemistry,registry,reactor.mapper());
        invalidToleranceSolver.setTolerances(1.0e-6,1.0e-11,0.0);
    }
    catch(const std::invalid_argument&)
    {
        invalidEnergyToleranceRejected=true;
    }
    check("invalid CVODE energy tolerance is rejected",
        invalidEnergyToleranceRejected);

    std::vector<double> mixture(registry.size(),0.0);
    mixture[registry.index("O2")]=0.784;
    mixture[registry.index("C3H8")]=0.216;

    RealFluidState thermo;
    thermo.setSpecies(registry.names());
    const auto initialStatus=thermo.updatePTY(4.0e6,1200.0,mixture);
    check("initial real-fluid state succeeds",
        initialStatus==RealFluidState::Status::SUCCESS);
    if(initialStatus!=RealFluidState::Status::SUCCESS) return 1;

    ReactorState initial{
        mixture,thermo.rho(),thermo.e(),thermo.T(),thermo.p()};
    const double sourceDt=reactor.estimateConservativeMaxDt(
        initial.rho,initial.pressure,initial.temperature,initial.Y,0.1);
    check("active chemistry has finite source timestep",
        std::isfinite(sourceDt)&&sourceDt>0.0&&sourceDt<1.0e29);

    const double totalTime=std::min(1.0e-7,0.5*sourceDt);
    std::cout << "  sourceDt=" << sourceDt
              << " totalTime=" << totalTime << '\n';
    ReactorState one=initial,two=initial,four=initial,eight=initial;
    const bool oneOk=advance(reactor,one,totalTime,1);
    const int oneSubcycles=reactor.lastSubcycles();
    const bool twoOk=advance(reactor,two,totalTime,2);
    const int twoSubcycles=reactor.lastSubcycles();
    const bool fourOk=advance(reactor,four,totalTime,4);
    const int fourSubcycles=reactor.lastSubcycles();
    const bool eightOk=advance(reactor,eight,totalTime,8);
    const int eightSubcycles=reactor.lastSubcycles();
    const int fuelIndex=registry.index("C3H8");
    std::cout << "  one: Yf=" << one.Y[fuelIndex]
              << " e=" << one.energy << " T=" << one.temperature << '\n'
              << "  two: Yf=" << two.Y[fuelIndex]
              << " e=" << two.energy << " T=" << two.temperature << '\n'
              << "  four: Yf=" << four.Y[fuelIndex]
              << " e=" << four.energy << " T=" << four.temperature << '\n'
              << "  eight: Yf=" << eight.Y[fuelIndex]
              << " e=" << eight.energy << " T=" << eight.temperature << '\n'
              << "  accepted N: " << oneSubcycles << ' '
              << twoSubcycles << ' ' << fourSubcycles << ' '
              << eightSubcycles << '\n';
    check("one-step conservative reactor succeeds",oneOk);
    check("two-half-step conservative reactor succeeds",twoOk);
    check("four-quarter-step conservative reactor succeeds",fourOk);
    check("eight-eighth-step conservative reactor succeeds",eightOk);

    if(oneOk&&twoOk&&fourOk&&eightOk)
    {
        ReactorState tightReference=initial;
        const bool tightReferenceOk=advanceWithTolerances(
            &chemistry,registry,reactor.mapper(),initial,totalTime,
            1.0e-11,1.0e-15,1.0e-6,tightReference);
        check("tight CVODE reference integration succeeds",tightReferenceOk);

        ReactorState tolerance1e6=initial;
        ReactorState tolerance1e7=initial;
        ReactorState tolerance1e8=initial;
        const bool tolerance1e6Ok=advanceWithTolerances(
            &chemistry,registry,reactor.mapper(),initial,totalTime,
            1.0e-6,1.0e-11,1.0e-2,tolerance1e6);
        const bool tolerance1e7Ok=advanceWithTolerances(
            &chemistry,registry,reactor.mapper(),initial,totalTime,
            1.0e-7,1.0e-12,1.0e-3,tolerance1e7);
        const bool tolerance1e8Ok=advanceWithTolerances(
            &chemistry,registry,reactor.mapper(),initial,totalTime,
            1.0e-8,1.0e-13,1.0e-4,tolerance1e8);
        check("production-tolerance reference integration succeeds",
            tolerance1e6Ok);
        check("medium-tolerance reference integration succeeds",
            tolerance1e7Ok);
        check("fine-tolerance reference integration succeeds",
            tolerance1e8Ok);

        if(tightReferenceOk&&tolerance1e6Ok&&tolerance1e7Ok&&tolerance1e8Ok)
        {
            const double productionError=
                stateDifference(tolerance1e6,tightReference);
            const double mediumError=
                stateDifference(tolerance1e7,tightReference);
            const double fineError=
                stateDifference(tolerance1e8,tightReference);
            const double fineConfiguredFloor = std::max(
                1.0e-13,
                1.0e-4/std::max(std::abs(tightReference.energy),1.0));
            std::cout << "  tolerance errors production=" << productionError
                      << " medium=" << mediumError
                      << " fine=" << fineError
                      << " fineFloor=" << fineConfiguredFloor << '\n';
            check("production CVODE tolerance meets the accuracy budget",
                productionError<=3.0e-6);
            check("CVODE error decreases at 1e-7 relative tolerance",
                mediumError<=0.5*productionError);
            check("CVODE error decreases at 1e-8 or reaches configured floor",
                fineError < mediumError
                    && (fineError<=0.5*mediumError
                        || fineError<=fineConfiguredFloor));

            const double oneError=stateDifference(one,tightReference);
            const double twoError=stateDifference(two,tightReference);
            const double fourError=stateDifference(four,tightReference);
            const double eightError=stateDifference(eight,tightReference);
            std::cout << "  split errors one=" << oneError
                      << " two=" << twoError
                      << " four=" << fourError
                      << " eight=" << eightError << '\n';
            constexpr double splitConsistencyTolerance=3.0e-6;
            check("production chemistry is consistent across 1/2/4/8 splits",
                oneError<=splitConsistencyTolerance &&
                twoError<=splitConsistencyTolerance &&
                fourError<=splitConsistencyTolerance &&
                eightError<=splitConsistencyTolerance);
        }
        check("chemistry changes composition",
            four.Y[registry.index("C3H8")]
                <initial.Y[registry.index("C3H8")]);

        RealFluidState finalState;
        finalState.setSpecies(registry.names());
        check("final rho-e state is invertible",
            finalState.updateRhoEY(
                four.rho,four.energy,four.Y,four.temperature)
                ==RealFluidState::Status::SUCCESS);
        check("reported final temperature follows rho-e state",
            std::abs(finalState.T()-four.temperature)
                <=1.0e-7*std::max(four.temperature,1.0));

        const double fixedCvEnergy=
            four.rho*1500.0*(four.temperature-initial.temperature);
        const double conservativeEnergy=
            four.rho*(four.energy-initial.energy);
        check("production energy is not fixed-Cv reconstruction",
            std::abs(conservativeEnergy-fixedCvEnergy)
                >1.0e-4*std::max(std::abs(conservativeEnergy),1.0));
    }

    std::vector<double> inert(registry.size(),0.0);
    inert[registry.index("O2")]=1.0;
    RealFluidState inertState;
    inertState.setSpecies(registry.names());
    check("inert state succeeds",
        inertState.updatePTY(1.0e6,500.0,inert)
            ==RealFluidState::Status::SUCCESS);
    const double inertDt=reactor.estimateConservativeMaxDt(
        inertState.rho(),inertState.p(),inertState.T(),inert,0.1);
    check("inert chemistry permits unbounded source timestep",inertDt>=1.0e29);

    std::cout << (failures==0
        ? "REAL_FLUID_CHEMISTRY_ENERGY_PASS\n"
        : "REAL_FLUID_CHEMISTRY_ENERGY_FAIL\n");
    return failures==0?0:1;
}
