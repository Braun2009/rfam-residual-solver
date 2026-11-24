// Copyright 2026 Braun
// SPDX-License-Identifier: Apache-2.0
#include "framework/FixedConservedSourceOperator.H"
#include "framework/MachResidualIntegrator.H"
#include "cht/SolidEnergyEqn.H"
#include <precice/precice.hpp>

#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace
{
using namespace Foam;
constexpr int cells = 8;
constexpr int species = 2;
constexpr int wallCell = cells/2;
constexpr double area = 1.0; // m^2; exchange integrated power, not flux density
constexpr double initialSolidTemperature = 300.0;

struct Trial
{
    double output;
    double heat;
    double balanceError;
};

class Fluid
{
    Domain domain_{cells, 0.02, species};
    MachResidualIntegrator integrator_;
    std::vector<std::unique_ptr<Operator>> ops_;

    double energy() const
    {
        double total = 0.0;
        for (int c = 1; c < cells-1; ++c) total += domain_.U(c).rhoE;
        return total*domain_.dx()*area;
    }

public:
    using State = Domain;
    explicit Fluid(double machFloor) : integrator_(machFloor)
    {
        for (int c = 0; c < cells; ++c)
        {
            auto& w = domain_.stagedW(c);
            w.p = 1e5;
            w.T = 500.0;
            w.U[0] = 2.0;
            w.Y[0] = 1.0;
        }
        if (!domain_.primitiveToConserved())
            throw std::runtime_error("fluid initial state: " + domain_.failure().detail);
        ops_.push_back(std::make_unique<ConvectionOperator>(std::make_unique<AUSMFluxAdapter>()));
    }
    const State& state() const { return domain_; }
    int innerIterations() const { return integrator_.iterations(); }
    double innerDefect() const { return integrator_.defect(); }
    std::uint64_t acceptedSteps() const { return domain_.acceptedSteps(); }
    void restore(const State& state) { domain_ = state; }
    Trial advance(double dt, double wallTemperature)
    {
        const double before = energy();
        // Hold the film power fixed during this local trial. preCICE repeats
        // the full window with a new wall temperature until it converges.
        const double power = 150.0*area*(domain_.W(wallCell).T - wallTemperature);
        std::vector<ConservedVars> source(cells, ConservedVars(species));
        source[wallCell].rhoE = -power/(domain_.dx()*area);
        ops_.push_back(std::make_unique<FixedConservedSourceOperator>(std::move(source)));
        const bool ok = integrator_.step(domain_, ops_, dt);
        ops_.pop_back();
        if (!ok) throw std::runtime_error("fluid step: " + domain_.failure().detail);

        // Fixed ghost cells supply boundary fluxes. Account for that energy
        // separately from the power exchanged with the solid.
        ResidualField flux;
        zeroResidual(flux, cells, species);
        ops_.front()->residual(domain_, flux);
        double outward = 0.0;
        for (int c = 1; c < cells-1; ++c) outward += flux[c].rhoE;
        const double heat = power*dt;
        return {power, heat, energy() - before + outward*domain_.dx()*area*dt + heat};
    }
};

class Solid
{
    static constexpr double rhoCp = 2600.0*900.0;
    SolidEnergyEqn solid_{5, 0.005, 2600.0, 900.0, 100.0, initialSolidTemperature};

public:
    using State = std::vector<double>;
    const State& state() const { return solid_.field(); }
    int innerIterations() const { return 0; }
    double innerDefect() const { return 0.0; }
    void restore(const State& state) { solid_.setField(state); }
    Trial advance(double dt, double power)
    {
        const State before = state();
        const double flux = power/area;
        solid_.step(dt, SolidBC::fixedFlux(flux), SolidBC::adiabatic());
        const double surfaceTemperature = solid_.T(0) + flux*0.5*solid_.dx()/solid_.k();
        double heatStored = 0.0;
        for (int c = 0; c < solid_.N(); ++c) heatStored += solid_.T(c)-before[c];
        return {surfaceTemperature, power*dt, heatStored*rhoCp*solid_.dx()*area - power*dt};
    }
};

template<class Model>
void couple(Model& model, const char* participant, const std::string& config,
            const std::string& mesh, const std::string& input,
            const std::string& output, double initialOutput)
{
    precice::Participant coupling(participant, config, 0, 1);
    const std::array<double, 2> position{0.0, 0.0};
    const std::array<precice::VertexID, 1> vertex{coupling.setMeshVertex(mesh, position)};
    if (coupling.requiresInitialData())
        coupling.writeData(mesh, output, vertex, std::array<double, 1>{initialOutput});
    coupling.initialize();

    std::optional<typename Model::State> checkpoint;
    int windows = 0, restores = 0, maxInner = 0;
    double heat = 0.0, maxError = 0.0, maxDefect = 0.0, physicalTime = 0.0;
    while (coupling.isCouplingOngoing())
    {
        if (coupling.requiresWritingCheckpoint()) checkpoint.emplace(model.state());
        const double dt = coupling.getMaxTimeStepSize();
        std::array<double, 1> incoming{};
        coupling.readData(mesh, input, vertex, dt, incoming);
        if (!std::isfinite(incoming[0])) throw std::runtime_error("non-finite coupling input");
        const auto trial = model.advance(dt, incoming[0]);
        if (!std::isfinite(trial.output) || !std::isfinite(trial.balanceError))
            throw std::runtime_error("non-finite coupling trial");
        maxInner = std::max(maxInner, model.innerIterations());
        maxDefect = std::max(maxDefect, model.innerDefect());
        coupling.writeData(mesh, output, vertex, std::array<double, 1>{trial.output});
        coupling.advance(dt);
        if (coupling.requiresReadingCheckpoint())
        {
            if (!checkpoint) throw std::logic_error("preCICE requested a missing checkpoint");
            model.restore(*checkpoint);
            ++restores;
        }
        else if (coupling.isTimeWindowComplete())
        {
            ++windows;
            physicalTime += dt;
            heat += trial.heat;
            maxError = std::max(maxError, std::abs(trial.balanceError));
            checkpoint.reset();
        }
        else throw std::logic_error("this example requires one local step per window");
    }
    coupling.finalize();
    std::cout << std::setprecision(16) << "RFAM_PRECICE " << participant
              << " windows=" << windows << " restores=" << restores
              << " heat_J=" << heat << " max_balance_error_J=" << maxError
              << " physical_time_s=" << physicalTime
              << " max_inner_iterations=" << maxInner << " max_inner_defect=" << maxDefect;
    if constexpr (std::is_same_v<Model, Fluid>)
        std::cout << " solver_accepted_steps=" << model.acceptedSteps();
    std::cout << '\n';
}
} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc < 3 || argc > 4)
            throw std::invalid_argument("usage: rfam_precice_thermal Fluid|Solid config.xml [fluid-mach-floor]");
        const std::string role = argv[1];
        if (role == "Fluid")
        {
            Fluid fluid(argc == 4 ? std::stod(argv[3]) : 0.2);
            couple(fluid, "Fluid", argv[2], "Fluid-Wall", "Wall-Temperature", "Heat-Power", 0.0);
        }
        else if (role == "Solid" && argc == 3)
        {
            Solid solid;
            couple(solid, "Solid", argv[2], "Solid-Wall", "Heat-Power", "Wall-Temperature", initialSolidTemperature);
        }
        else throw std::invalid_argument("unexpected participant or arguments");
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
