#ifndef MeshReader_H
#define MeshReader_H

/*---------------------------------------------------------------------------*\
    MeshReader - simple mesh I/O for RFAM (reads Vortex chamber profile).

    Level 1: reads a JSON file with {profile: [{x_mm, r_mm}, ...]} format
    (Vortex v5 chamber_profile_3pintle.json) and generates a 2D axisymmetric
    structured grid as Domain cell coordinates + r_wall(x) for BC application.

    Level 2 (future): read general 2D/3D structured grid from file.

    Usage:
        MeshReader reader;
        auto info = reader.loadVortexProfile("config/chamber_profile_3pintle.json");
        Domain dom(info.Nx, info.Ny, info.Lx, info.Ly, nSpecies);
        // info.r_wall[i] gives the wall radius at column i

    Pure std C++, no third-party libraries.
\*---------------------------------------------------------------------------*/

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace Foam
{

struct MeshInfo
{
    int Nx, Ny;
    double Lx, Ly, dx, dy;
    std::vector<double> x_cell;      // (Nx) cell centres [m]
    std::vector<double> r_wall;      // (Nx) wall radius at each column [m]
    double r_throat, r_chamber;      // key dimensions [m]
    int i_throat;                    // column index of throat
};

class MeshReader
{
public:
    // Load Vortex v5 chamber profile JSON. Returns mesh info for Domain creation.
    // Nx = number of axial cells, Ny = number of radial cells.
    // Ly is set to r_max * margin (fits the full chamber radius with headroom).
    static MeshInfo loadVortexProfile(const std::string& path,
                                     int Nx = 80, int Ny = 30,
                                     double margin = 1.2)
    {
        // Parse the JSON for profile entries (minimal parser)
        std::ifstream f(path);
        if (!f.good())
            throw std::runtime_error("MeshReader: cannot open '" + path + "'");

        std::string content((std::istreambuf_iterator<char>(f)),
                             std::istreambuf_iterator<char>());

        std::vector<double> x_prof, r_prof;
        size_t pos = 0;
        while ((pos = content.find("\"x_mm\":", pos)) != std::string::npos)
        {
            pos += 7;
            double x = std::stod(content.substr(pos));
            size_t rp = content.find("\"r_mm\":", pos);
            if (rp == std::string::npos) break;
            double r = std::stod(content.substr(rp + 7));
            x_prof.push_back(x * 1e-3);   // mm → m
            r_prof.push_back(r * 1e-3);
            pos = rp + 10;
        }

        if (x_prof.empty())
            throw std::runtime_error("MeshReader: no profile points found in '" + path + "'");

        MeshInfo info;
        info.Lx = x_prof.back();
        double r_max = *std::max_element(r_prof.begin(), r_prof.end());
        info.Ly = r_max * margin;
        info.Nx = Nx;
        info.Ny = Ny;
        info.dx = info.Lx / Nx;
        info.dy = info.Ly / Ny;

        // Cell centres and interpolated wall radius
        info.x_cell.resize(Nx);
        info.r_wall.resize(Nx);
        info.r_throat = 1e9;
        info.r_chamber = 0.0;
        info.i_throat = 0;

        for (int i = 0; i < Nx; ++i)
        {
            double x = (i + 0.5) * info.dx;
            info.x_cell[i] = x;

            // linear interpolation of r_wall at this x
            double rw = r_prof.back();
            if (x <= x_prof.front()) rw = r_prof.front();
            else if (x >= x_prof.back()) rw = r_prof.back();
            else
            {
                for (size_t k = 1; k < x_prof.size(); ++k)
                {
                    if (x_prof[k] >= x)
                    {
                        double t = (x - x_prof[k-1]) / (x_prof[k] - x_prof[k-1]);
                        rw = r_prof[k-1] + t * (r_prof[k] - r_prof[k-1]);
                        break;
                    }
                }
            }
            info.r_wall[i] = rw;
            if (rw < info.r_throat) { info.r_throat = rw; info.i_throat = i; }
            if (rw > info.r_chamber) info.r_chamber = rw;
        }

        return info;
    }
};

} // End namespace Foam
#endif
