// Copyright 2026 Braun
// SPDX-License-Identifier: Apache-2.0
#include "framework/FixedConservedSourceOperator.H"
#include "framework/TimeIntegrator.H"

#include <iostream>

using namespace Foam;

namespace
{
Domain initialState()
{
    Domain d(32, 0.032, 2);
    for (int c = 0; c < d.N(); ++c)
    {
        auto& w = d.stagedW(c);
        const double wave = std::sin(2.0*std::acos(-1.0)*d.x(c)/d.L());
        w.p = 1e5*(1.0 + 0.1*wave);
        w.T = 300.0*(1.0 + 0.05*wave);
        w.U[0] = 50.0;
        w.Y = {0.23, 0.77};
    }
    if (!d.primitiveToConserved()) throw std::runtime_error(d.failure().detail);
    return d;
}

std::vector<std::unique_ptr<Operator>> operators(const Domain& d, bool sourceFirst)
{
    ResidualField source(d.N(), ConservedVars(d.nSpecies()));
    for (int c = d.N()/3; c < 2*d.N()/3; ++c) source[c].rhoE = 4e8;
    source[d.N()/2].rhoY = {1e-3, -1e-3};
    std::vector<std::unique_ptr<Operator>> ops;
    ops.push_back(std::make_unique<ConvectionOperator>(std::make_unique<AUSMFluxAdapter>()));
    ops.push_back(std::make_unique<FixedConservedSourceOperator>(std::move(source)));
    if (sourceFirst) std::swap(ops[0], ops[1]);
    return ops;
}

double difference(const ConservedVars& a, const ConservedVars& b)
{
    const auto relative = [](double x, double y)
    {
        return std::abs(x-y)/std::max({1.0, std::abs(x), std::abs(y)});
    };
    double result = std::max({relative(a.rho,b.rho), relative(a.rhoE,b.rhoE),
        relative(a.rhoZ,b.rhoZ), relative(a.rhoC,b.rhoC)});
    for (std::size_t i = 0; i < a.rhoU.size(); ++i)
        result = std::max(result,relative(a.rhoU[i],b.rhoU[i]));
    for (std::size_t i = 0; i < a.rhoY.size(); ++i)
        result = std::max(result,relative(a.rhoY[i],b.rhoY[i]));
    return result;
}

void require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}

class ApplyOnly : public Operator
{
public:
    void apply(Domain&, double) const override {}
    std::string name() const override { return "apply-only"; }
};
}

int main()
{
    try
    {
        constexpr double dt = 2e-7;
        const auto start = initialState();
        auto additive = start, residual = start, reordered = start;
        const auto ops = operators(start, false);
        const auto reverse = operators(start, true);
        const OperatorSplittingIntegrator legacy;
        const ResidualIntegrator euler;
        for (int n = 0; n < 5; ++n)
        {
            require(legacy.step(additive, ops, dt), "legacy step rejected");
            require(euler.step(residual, ops, dt), "residual step rejected");
            require(euler.step(reordered, reverse, dt), "reordered step rejected");
        }
        for (int c = 0; c < start.N(); ++c)
        {
            require(difference(additive.U(c), residual.U(c)) < 1e-12,
                    "additive and residual trajectories differ");
            require(difference(residual.U(c), reordered.U(c)) < 1e-14,
                    "operator order changed the residual trajectory");
        }

        ResidualField r;
        zeroResidual(r, start.N(), start.nSpecies());
        for (const auto& op : ops) op->residual(start, r);
        for (double step : {dt, 2*dt})
        {
            auto d = start;
            require(euler.step(d, ops, step), "dt test step rejected");
            for (int c = 0; c < start.N(); ++c)
            {
                auto expected = start.U(c);
                subtractScaledRate(expected, r[c], step);
                require(difference(d.U(c), expected) < 1e-14,
                        "Euler update does not use dt exactly once");
            }
        }

        auto rejected = start;
        std::vector<std::unique_ptr<Operator>> unsupported;
        unsupported.push_back(std::make_unique<ApplyOnly>());
        require(!euler.step(rejected, unsupported, dt), "apply-only operator was silently accepted");
        require(rejected.failure().component == "apply-only", "missing operator diagnostic");
        require(!rejected.stepTransactionActive(), "rollback did not close transaction");
        for (int c = 0; c < start.N(); ++c)
            require(difference(rejected.U(c),start.U(c)) == 0.0, "failed step changed accepted state");
        std::cout << "RESIDUAL_FORM_PASS\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
