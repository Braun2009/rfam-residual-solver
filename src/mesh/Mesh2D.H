/*---------------------------------------------------------------------------*  realFluidAllMach - 2D Structured Mesh
\*---------------------------------------------------------------------------*/

#ifndef Mesh2D_H
#define Mesh2D_H

#include <vector>
#include <cmath>

namespace Foam
{

struct Cell2D
{
    double x, y;      // Cell center
    double dx, dy;    // Cell size
    double volume;    // Cell volume (axisymmetric)
};

class Mesh2D
{
private:
    int Nx_, Ny_;
    double Lx_, Ly_;
    std::vector<std::vector<Cell2D>> cells_;
    
    // Stretching function (hyperbolic tangent)
    double stretch(double xi, double beta) const
    {
        if (beta < 1e-6) return xi;
        return (tanh(beta * (xi - 0.5)) / tanh(0.5 * beta) + 1.0) / 2.0;
    }

public:
    Mesh2D(int Nx, int Ny, double Lx, double Ly)
        : Nx_(Nx), Ny_(Ny), Lx_(Lx), Ly_(Ly)
    {
        cells_.resize(Nx_);
        for (int i = 0; i < Nx_; i++)
            cells_[i].resize(Ny_);
    }

    void generate(double beta_x = 0.0, double beta_y = 2.0)
    {
        // Generate uniform grid in computational space, stretched in physical
        for (int i = 0; i < Nx_; i++)
        {
            for (int j = 0; j < Ny_; j++)
            {
                double xi = (i + 0.5) / Nx_;
                double eta = (j + 0.5) / Ny_;
                
                // Physical coordinates with stretching
                cells_[i][j].x = stretch(xi, beta_x) * Lx_;
                cells_[i][j].y = stretch(eta, beta_y) * Ly_;
                
                // Cell sizes (simplified, uniform in computational space)
                cells_[i][j].dx = Lx_ / Nx_;
                cells_[i][j].dy = Ly_ / Ny_;
                
                // Axisymmetric volume: 2*pi*y*dx*dy
                double y_face = cells_[i][j].y;
                cells_[i][j].volume = 2.0 * M_PI * y_face * cells_[i][j].dx * cells_[i][j].dy;
            }
        }
    }

    Cell2D& operator()(int i, int j) { return cells_[i][j]; }
    const Cell2D& operator()(int i, int j) const { return cells_[i][j]; }
    
    int Nx() const { return Nx_; }
    int Ny() const { return Ny_; }
    double Lx() const { return Lx_; }
    double Ly() const { return Ly_; }
};

} // End namespace Foam

#endif
