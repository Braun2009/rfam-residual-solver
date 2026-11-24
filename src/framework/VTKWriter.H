#ifndef VTKWriter_H
#define VTKWriter_H

/*---------------------------------------------------------------------------*\
    VTKWriter - simple VTK rectilinear grid writer for 1D/2D Domain output.

    Writes ASCII VTK legacy format (*.vtk) consumable by ParaView/Visit.
    Supports scalar fields (density, pressure, temperature, species Y_i)
    and the velocity vector. Works for both 1D (Ny==1) and 2D (Ny>1) domains.

    Usage:
        VTKWriter w("output_0010.vtk", dom);
        w.writeAll(reg);   // writes rho, p, T, U, Y_i for all species in reg
        // or selectively:
        w.writeScalar("rho", [&](int i,int j){ return dom.U(i,j).rho; });
\*---------------------------------------------------------------------------*/

#include "Domain.H"
#include "SpeciesRegistry.H"
#include "TurbulenceOperator.H"
#include <fstream>
#include <string>
#include <functional>
#include <sstream>
#include <iomanip>
#include <cmath>

namespace Foam
{

class VTKWriter
{
private:
    std::ofstream f_;
    const Domain& dom_;
    int Nx_, Ny_;

public:
    VTKWriter(const std::string& filename, const Domain& dom)
        : dom_(dom), Nx_(dom.Nx()), Ny_(dom.Ny())
    {
        f_.open(filename);
        writeHeader();
    }

    ~VTKWriter() { f_.close(); }

    bool good() const { return f_.good(); }

    // Write a named scalar field over all cells
    void writeScalar(const std::string& name,
                     std::function<double(int i, int j)> fn)
    {
        writeScalar3D(name,[&](int i,int j,int){ return fn(i,j); });
    }

    // Write a named scalar field with explicit 3D indexing.
    void writeScalar3D(const std::string& name,
                       std::function<double(int i,int j,int k)> fn)
    {
        f_ << "SCALARS " << name << " float\nLOOKUP_TABLE default\n";
        int Nz = dom_.Nz();
        for (int k = 0; k < Nz; ++k)
            for (int j = 0; j < Ny_; ++j)
                for (int i = 0; i < Nx_; ++i)
                    f_ << fn(i,j,k) << "\n";
    }

    // Write velocity as a VTK VECTORS field
    void writeVelocity()
    {
        f_ << "VECTORS velocity float\n";
        int Nz = dom_.Nz();
        for (int k = 0; k < Nz; ++k)
            for (int j = 0; j < Ny_; ++j)
                for (int i = 0; i < Nx_; ++i)
                    f_ << dom_.W(i,j).U[0] << " "
                       << dom_.W(i,j).U[1] << " "
                       << dom_.W(i,j).U[2] << "\n";
    }

    // Write rho, p, T, velocity, and all species mass fractions
    void writeAll(const SpeciesRegistry& reg)
    {
        writeScalar("density",     [&](int i,int j){ return dom_.U(i,j).rho; });
        writeScalar("pressure",    [&](int i,int j){ return dom_.W(i,j).p;   });
        writeScalar("temperature", [&](int i,int j){ return dom_.W(i,j).T;   });
        writeScalar("soundspeed",  [&](int i,int j){
            return dom_.Q(i,j).acoustic.soundSpeed();
        });
        writeVelocity();
        int nS = dom_.nSpecies();
        for (int g = 0; g < nS && g < reg.size(); ++g)
        {
            std::string sname = "Y_" + reg.name(g);
            writeScalar(sname, [&,g](int i,int j){
                return dom_.W(i,j).Y.size() > (size_t)g
                     ? dom_.W(i,j).Y[g] : 0.0;
            });
        }
    }

    // 1D convenience: write with domain-specific species names
    void writeAll(const std::vector<std::string>& names)
    {
        writeScalar("density",     [&](int i,int j){ return dom_.U(i,j).rho; });
        writeScalar("pressure",    [&](int i,int j){ return dom_.W(i,j).p;   });
        writeScalar("temperature", [&](int i,int j){ return dom_.W(i,j).T;   });
        writeVelocity();
        for (int g = 0; g < (int)names.size(); ++g)
        {
            std::string sname = "Y_" + names[g];
            writeScalar(sname, [&,g](int i,int j){
                return dom_.W(i,j).Y.size() > (size_t)g
                     ? dom_.W(i,j).Y[g] : 0.0;
            });
        }
    }

    // Write turbulence fields k, omega, nu_t from a TurbulenceOperator.
    void writeTurbulence(const TurbulenceOperator& turb)
    {
        writeScalar("k_turb",   [&](int i,int j){ return turb.k(i,j,Nx_);         });
        writeScalar("omega",    [&](int i,int j){ return turb.omega(i,j,Nx_);      });
        writeScalar("nut",      [&](int i,int j){ return turb.nu_t(i,j,Nx_,dom_); });
    }

    // ---- Static helpers for time-series filenames ----

    // Build a zero-padded VTK filename: "prefix_<step>.vtk"
    static std::string filename(const std::string& prefix, int step, int width = 6)
    {
        std::ostringstream s;
        s << prefix << "_" << std::setw(width) << std::setfill('0') << step << ".vtk";
        return s.str();
    }

    // Build a filename encoding physical time: "prefix_t<time_us>.vtk" (µs)
    static std::string filenameTime(const std::string& prefix, double t, int decimals = 2)
    {
        std::ostringstream s;
        s << prefix << "_t" << std::fixed << std::setprecision(decimals)
          << t * 1e6 << "us.vtk";
        return s.str();
    }

private:
    void writeHeader()
    {
        f_ << "# vtk DataFile Version 3.0\n";
        f_ << "RFAM output\n";
        f_ << "ASCII\n";
        f_ << "DATASET RECTILINEAR_GRID\n";

        int Nz = dom_.Nz();
        // Dimensions: VTK uses point count = cell count + 1 in each direction
        f_ << "DIMENSIONS " << Nx_+1 << " " << Ny_+1 << " " << Nz+1 << "\n";

        // X coordinates (cell-face positions)
        f_ << "X_COORDINATES " << Nx_+1 << " float\n";
        double dx = dom_.dx();
        for (int i = 0; i <= Nx_; ++i)
            f_ << i * dx * 1000.0 << " ";   // convert to mm for ParaView
        f_ << "\n";

        // Y coordinates
        double dy = (Ny_ > 1) ? dom_.dy() : dx;
        f_ << "Y_COORDINATES " << Ny_+1 << " float\n";
        for (int j = 0; j <= Ny_; ++j)
            f_ << j * dy * 1000.0 << " ";
        f_ << "\n";

        // Z coordinates
        double dz = (Nz > 1) ? dom_.dz() : dy;
        f_ << "Z_COORDINATES " << Nz+1 << " float\n";
        for (int k = 0; k <= Nz; ++k)
            f_ << k * dz * 1000.0 << " ";
        f_ << "\n";

        f_ << "CELL_DATA " << Nx_ * Ny_ * Nz << "\n";
    }
};

} // End namespace Foam

#endif
