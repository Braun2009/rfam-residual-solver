// Copyright 2026 Braun
// SPDX-License-Identifier: Apache-2.0
#include "framework/FixedConservedSourceOperator.H"
#include "framework/MachResidualIntegrator.H"

#include <array>
#include <iomanip>
#include <iostream>

using namespace Foam;

int main(int argc, char** argv)
{
    try
    {
        const bool stress = argc == 2 && std::string(argv[1]) == "--stress";
        if (argc != 1 && !stress)
            throw std::invalid_argument("usage: test_mach_residual_integrator [--stress]");
        const auto require = [](bool test, const char* msg) {
            if (!test) throw std::runtime_error(msg);
        };
        Domain start(8, 0.02, 2);
        for (int c = 0; c < start.N(); ++c)
        {
            auto& w = start.stagedW(c);
            w.p = 1e5;
            w.T = 500.0;
            w.U[0] = 2.0;
            w.Y = {1.0, 0.0};
        }
        if (!start.primitiveToConserved()) throw std::runtime_error(start.failure().detail);
        const auto& q = start.Q(3);
        const auto& w = start.W(3);
        constexpr double floor = 0.2;
        MachPressurePreconditioner matrix(floor);
        auto mode = start.U(3);
        mode.rho = 1.0;
        mode.rhoE = q.h + 0.5*w.U[0]*w.U[0];
        mode.rhoU = w.U;
        mode.rhoY = w.Y;
        for (double sigma : {0.0, 1e-4, 0.15, 1.0, std::numeric_limits<double>::infinity()})
        {
            auto scaled = mode;
            matrix.apply(w,q,scaled,sigma);
            const double expected = std::isinf(sigma) ? floor*floor
                : (1.0 + sigma)/(1.0 + sigma/(floor*floor));
            require(std::abs(scaled.rho - expected) < 1e-14,
                    "matrix pressure mode does not match the physical/pseudo mass inverse");
        }
        auto velocityMode = start.U(3);
        velocityMode.rho = 0.0;
        velocityMode.rhoU = {1.0, 0.0, 0.0};
        velocityMode.rhoE = w.U[0]*velocityMode.rhoU[0];
        velocityMode.rhoY = {0.0, 0.0};
        const auto originalMode = velocityMode;
        matrix.apply(w,q,velocityMode,0.15);
        require(velocityMode.rho == originalMode.rho
                && velocityMode.rhoU == originalMode.rhoU
                && velocityMode.rhoE == originalMode.rhoE,
                "preconditioner modified a pressure-neutral velocity correction");
        auto supersonic = w;
        supersonic.U[0] = 2.0*q.acoustic.soundSpeed();
        auto unchanged = start.U(3);
        matrix.apply(supersonic,q,unchanged,0.15);
        require(unchanged.rho == start.U(3).rho && unchanged.rhoE == start.U(3).rhoE
                && unchanged.rhoU == start.U(3).rhoU && unchanged.rhoY == start.U(3).rhoY,
                "matrix is not identity at supersonic Mach");

        std::vector<std::unique_ptr<Operator>> ops;
        ops.push_back(std::make_unique<ConvectionOperator>(std::make_unique<AUSMFluxAdapter>()));
        ResidualField source(start.N(), ConservedVars(2));
        source[3].rhoE = -1.2e7;
        ops.push_back(std::make_unique<FixedConservedSourceOperator>(std::move(source)));
        constexpr double dt = 1e-6;
        auto plain = start, conditioned = start, rejected = start;
        const MachResidualIntegrator off(1.0), on(floor), tooFew(floor, 1);
        require(off.step(plain,ops,dt), "unconditioned implicit step failed");
        require(on.step(conditioned,ops,dt), "preconditioned implicit step failed");
        for (int c = 0; c < start.N(); ++c)
        {
            const auto& a = plain.U(c);
            const auto& b = conditioned.U(c);
            const double rhoDiff = std::abs(a.rho-b.rho)/a.rho;
            const double energyDiff = std::abs(a.rhoE-b.rhoE)/std::max(1.0,std::abs(a.rhoE));
            // Use an acoustic momentum scale, not the small advective
            // momentum. The EOS inversion is itself finite-tolerance.
            const double momentumScale = a.rho*plain.Q(c).acoustic.soundSpeed();
            const double momentumDiff = std::abs(a.rhoU[0]-b.rhoU[0])/momentumScale;
            if (!(rhoDiff < 1e-9 && energyDiff < 1e-9 && momentumDiff < 1e-9))
            {
                std::cerr << "cell=" << c << " rho=" << rhoDiff << " energy_rel="
                          << energyDiff << " momentum=" << momentumDiff << '\n';
                throw std::runtime_error("preconditioning changed the accepted physical solution");
            }
        }
        require(!tooFew.step(rejected,ops,dt), "unconverged correction was accepted");
        require(rejected.acceptedSteps() == start.acceptedSteps()
                && !rejected.stepTransactionActive(), "failed step was not rolled back");
        for (int c = 0; c < start.N(); ++c)
            require(start.U(c).rhoE == rejected.U(c).rhoE && start.W(c).p == rejected.W(c).p,
                    "rejected inner correction changed the window start");
        std::cout << "MACH_RESIDUAL_PASS off=" << off.iterations()
                  << " on=" << on.iterations() << " defect=" << on.defect() << '\n';

        // Hold the EOS, source, geometry and nonlinear tolerance fixed.
        // Sweep the initial Mach and a dimensionless acoustic step size.
        const double a0 = q.acoustic.soundSpeed();
        const std::vector<double> stepSizes = stress
            ? std::vector<double>{0.6, 0.8, 1.0, 1.3, 1.6}
            : std::vector<double>{0.025, 0.18, 0.4};
        for (const double mach : std::array<double,5>{0.005, 0.05, 0.3, 0.9, 1.2})
        {
            auto initial = start;
            for (int c = 0; c < initial.N(); ++c) initial.stagedW(c).U[0] = mach*a0;
            require(initial.primitiveToConserved(), "could not initialize Mach sweep state");
            for (const double cfl : stepSizes)
            {
                auto baseline = initial, tested = initial;
                const double step = cfl*initial.dx()/(a0*(1.0 + mach));
                const bool offOK = off.step(baseline,ops,step);
                const bool onOK = on.step(tested,ops,step);
                std::cout << std::setprecision(8) << "MACH_SWEEP mach=" << mach
                          << " acoustic_cfl=" << cfl << " off=" << off.iterations()
                          << " on=" << on.iterations() << " defect=" << on.defect()
                          << " off_accepted=" << offOK << " on_accepted=" << onOK << '\n';
                if (stress && !(offOK && onOK))
                {
                    const auto verifyRejected = [&](const Domain& rejectedStep)
                    {
                        require(!rejectedStep.stepTransactionActive()
                                && rejectedStep.acceptedSteps() == initial.acceptedSteps(),
                                "stress: rejected step changed accepted-step state");
                        for (int c = 0; c < initial.N(); ++c)
                        {
                            const auto& a = initial.U(c);
                            const auto& b = rejectedStep.U(c);
                            require(a.rho == b.rho && a.rhoE == b.rhoE && a.rhoU == b.rhoU
                                    && a.rhoY == b.rhoY && a.rhoZ == b.rhoZ && a.rhoC == b.rhoC
                                    && initial.W(c).p == rejectedStep.W(c).p,
                                    "stress: rejected step changed accepted fields");
                        }
                    };
                    if (!offOK) verifyRejected(baseline);
                    if (!onOK) verifyRejected(tested);
                    std::cout << "MACH_SWEEP_REJECT off_detail=" << baseline.failure().detail
                              << " | on_detail=" << tested.failure().detail << '\n';
                    continue;
                }
                if (!offOK) throw std::runtime_error("unscaled sweep case: " + baseline.failure().detail);
                if (!onOK) throw std::runtime_error("scaled sweep case: " + tested.failure().detail);
                require(off.defect() <= 1e-12 && on.defect() <= 1e-12,
                        "sweep case did not meet the same unscaled tolerance");
                for (int c = 0; c < start.N(); ++c)
                {
                    const auto& b = baseline.U(c);
                    const auto& u = tested.U(c);
                    require(std::abs(u.rho-b.rho)/b.rho < 2e-9
                            && std::abs(u.rhoE-b.rhoE)/std::max(1.0,std::abs(b.rhoE)) < 2e-9
                            && std::abs(u.rhoU[0]-b.rhoU[0])/(b.rho*a0) < 2e-9,
                            "sweep: preconditioning and baseline converged to different states");
                }
                if (mach > 1.0)
                    require(off.iterations() == on.iterations() && off.defect() == on.defect(),
                            "supersonic step should be identical with scaling on or off");
            }
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
