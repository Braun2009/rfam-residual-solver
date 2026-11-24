#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"

using namespace Foam;

static std::vector<ConservedVars> runCase(int threads)
{
#ifdef _OPENMP
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
    Domain d(12, 8, 6, 0.12, 0.08, 0.06, 2);
    for (int i = 0; i < d.Nx(); ++i)
    for (int j = 0; j < d.Ny(); ++j)
    for (int k = 0; k < d.Nz(); ++k)
    {
        PrimitiveVars w(2);
        w.p = 1.0e6 + 1000.0 * i;
        w.T = 450.0 + 2.0 * j;
        w.U = {20.0 + i, 0.5 * j, 0.25 * k};
        w.Y[0] = 0.7;
        w.Y[1] = 0.3;
        d.stagedW(i,j,k) = w;
    }
    d.primitiveToConserved();
    d.beginStep();
    ConvectionOperator op(std::unique_ptr<FluxScheme>(new AUSMFluxAdapter()));
    op.apply(d, 1e-8);
    d.commit();

    std::vector<ConservedVars> out;
    out.reserve(d.N());
    for (int n = 0; n < d.N(); ++n) out.push_back(d.U(n));
    return out;
}

int main()
{
    const auto serial = runCase(1);
    const auto parallel = runCase(4);
    if (serial.size() != parallel.size()) return 1;

    for (size_t i = 0; i < serial.size(); ++i)
    {
        if (serial[i].rho != parallel[i].rho ||
            serial[i].rhoE != parallel[i].rhoE)
        {
            std::cerr << "OpenMP mismatch at cell " << i << "\n";
            return 1;
        }
        for (int d = 0; d < 3; ++d)
            if (serial[i].rhoU[d] != parallel[i].rhoU[d]) return 1;
        for (size_t s = 0; s < serial[i].rhoY.size(); ++s)
            if (serial[i].rhoY[s] != parallel[i].rhoY[s]) return 1;
    }

    std::cout << "OPENMP_CONSISTENCY_PASS\n";
    return 0;
}
