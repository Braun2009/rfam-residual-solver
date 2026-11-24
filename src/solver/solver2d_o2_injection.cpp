#include <iostream>
#include <vector>
#include <cmath>
#include "../mesh/Mesh2D.H"
#include "../variables/ConservedVars.H"
#include "../numerics/flux/AUSMPlusUpFlux.H"
#include "../thermo/RealFluidState.H"

using namespace Foam;
using namespace std;

class Solver2D
{
private:
    Mesh2D mesh_;
    vector<vector<ConservedVars>> U_;
    vector<vector<PrimitiveVars>> W_;
    AUSMPlusUpFlux fluxScheme_;
    double t_;
    int iter_;

public:
    Solver2D(int Nx, int Ny, double Lx, double Ly)
        : mesh_(Nx, Ny, Lx, Ly), t_(0.0), iter_(0)
    {
        U_.resize(Nx);
        W_.resize(Nx);
        for (int i = 0; i < Nx; i++)
        {
            U_[i].resize(Ny);
            W_[i].resize(Ny);
        }
        mesh_.generate(0.0, 2.0);
    }

    void run(double t_end, double CFL)
    {
        cout << "2D solver placeholder - will implement after mesh test" << endl;
    }
};

int main()
{
    cout << "2D O2 Injection Test" << endl;
    Solver2D solver(100, 50, 0.020, 0.010);
    solver.run(0.001, 0.3);
    return 0;
}
