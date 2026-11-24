#include <cmath>
#include <iostream>

#include "src/framework/Domain.H"

using namespace Foam;

static void scaleConserved(ConservedVars& state, double scale)
{
    state.rho *= scale;
    state.rhoE *= scale;
    state.rhoZ *= scale;
    state.rhoC *= scale;
    for (double& value : state.rhoU) value *= scale;
    for (double& value : state.rhoY) value *= scale;
}

int main()
{
    Domain d(4, 3, 3, 0.4, 0.3, 0.3, 2);

    for (int i = 0; i < d.Nx(); ++i)
    for (int j = 0; j < d.Ny(); ++j)
    for (int k = 0; k < d.Nz(); ++k)
    {
        PrimitiveVars w(2);
        w.p = 1.0e6;
        w.T = 500.0;
        w.U = {10.0, 0.0, 0.0};
        w.Y[0] = 1.0;
        w.Y[1] = 0.0;
        d.stagedW(i, j, k) = w;
    }
    d.primitiveToConserved();

    const int i = d.Nx() - 2;
    const int j = d.Ny() - 2;
    const int k = d.Nz() - 1;
    const double rho0 = d.U(i, j, k).rho;

    d.beginStep();
    scaleConserved(d.Unew(i, j, k), 1.1);
    d.commit();
    if (d.tripped() || std::abs(d.U(i, j, k).rho - rho0 * 1.1) > 1e-12)
    {
        std::cerr << "3D commit skipped the final azimuthal layer\n";
        return 1;
    }

    const double rho1 = d.U(i, j, k).rho;
    d.saveStateForRK2();
    d.beginStep();
    scaleConserved(d.Unew(i, j, k), 1.2);
    d.averageWithSavedState();
    d.commit();

    const double expected = rho1 * 1.1;
    if (d.tripped() || std::abs(d.U(i, j, k).rho - expected) > 1e-12)
    {
        std::cerr << "3D RK2 averaging skipped the final azimuthal layer\n";
        return 1;
    }

    std::cout << "DOMAIN3D_LIFECYCLE_PASS\n";
    return 0;
}
