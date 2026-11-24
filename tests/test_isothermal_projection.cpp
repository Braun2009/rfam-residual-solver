#include <cmath>
#include <iostream>

#include "src/framework/Domain.H"
#include "src/framework/IsothermalProjection.H"

using namespace Foam;

int main()
{
    Domain d(20, 0.2, 2);
    for (int i = 0; i < d.N(); ++i)
    {
        PrimitiveVars w(2);
        w.p = 1.0e6 + 2.0e4 * i;
        w.T = 350.0 + 10.0 * i;
        w.U = {30.0, 0.0, 0.0};
        w.Y[0] = 0.8;
        w.Y[1] = 0.2;
        d.stagedW(i) = w;
    }
    d.primitiveToConserved();
    d.beginStep();

    IsothermalProjectionOperator projection(700.0);
    projection.apply(d, 0.0);
    d.commit();
    if (d.tripped()) return 1;

    for (int i = 0; i < d.N(); ++i)
    {
        if (!std::isfinite(d.W(i).T) || std::abs(d.W(i).T - 700.0) > 1e-5)
        {
            std::cerr << "isothermal projection mismatch at " << i
                      << ": " << d.W(i).T << " K\n";
            return 1;
        }
    }
    std::cout << "ISOTHERMAL_PROJECTION_PASS\n";
    return 0;
}
