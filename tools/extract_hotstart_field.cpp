// Export a RFAM sector checkpoint as a hot-start field CSV
// (rfam.hotstart-field/v1) that --hotstart-field can re-import.
//
// Usage: extract_hotstart_field <checkpoint.chk> <output.csv>
//        [--grid smoke|pathfinder|intermediate]
//
// This is the self-bootstrap loop: checkpoint -> external field CSV ->
// --hotstart-field import.  Equivalence is validated by comparing the
// exported (p,T,Y,U) field with the imported reconciliation (they must be
// identical up to the thermodynamic backend round trip).

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "framework/Domain.H"
#include "framework/SectorCheckpoint.H"
#include "SectorTurbulenceOperator.H"
#include "mesh/SectorMeshMetrics.H"
#include "GERGResidualTable.H"
#include "thermo/SpeciesRegistry.H"

using namespace Foam;

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        std::cerr << "usage: extract_hotstart_field <checkpoint> <out.csv>"
                  << " [--grid smoke|pathfinder|intermediate]\n";
        return 2;
    }
    const std::string checkpointPath = argv[1];
    const std::string outputPath = argv[2];
    std::string gridName = "smoke";
    for (int a = 3; a < argc; ++a)
    {
        if (std::strcmp(argv[a], "--grid") == 0 && a + 1 < argc)
            gridName = argv[++a];
    }

    const SectorCheckpointInfo info = SectorCheckpoint::inspect(checkpointPath);
    std::cout << "checkpoint: nx=" << info.nx << " ny=" << info.ny
              << " nz=" << info.nz << " nSpecies=" << info.nSpecies
              << " iter=" << info.metadata.iteration
              << " t=" << info.metadata.time << "\n";
    if (info.nx < 4 || info.ny < 4 || info.nz < 4)
    {
        std::cerr << "invalid checkpoint dimensions\n";
        return 1;
    }
    const int Nx = info.nx, Ny = info.ny, Nz = info.nz;

    const std::string profilePath = vortexProfilePath();
    SectorMeshMetrics mesh(profilePath, Nx, Ny, Nz, 120.0);
    Domain dom(Nx, Ny, Nz, mesh.length(), mesh.radiusCenter(Nx / 2, Ny - 2),
               mesh.sectorAngle(), info.nSpecies);
    dom.setSpecies(info.species);
    {
        static GERGResidualTable gergTable;
        const std::string gergPath = std::string(RFAM_SOURCE_DIR)
            + "/tests/reference/gerg_residual_table_v3.csv";
        if (gergTable.load(gergPath))
            dom.thermoModel().setGERGTable(&gergTable);
        else
            std::cerr << "WARNING: GERG table missing; checkpoint recovery "
                      << "of dense LOX cells will fail: " << gergPath << "\n";
    }

    const SectorCheckpointMetadata metadata =
        SectorCheckpoint::read(checkpointPath, dom);
    std::cout << "restored: iter=" << metadata.iteration
              << " t=" << metadata.time << "\n";

    std::ofstream out(outputPath);
    out << "# rfam.hotstart-field/v1\n";
    out << "# time=" << metadata.time << "\n";
    out << "# iteration=" << metadata.iteration << "\n";
    out << "x,y,z,p,T,u,v,w,rho,e";
    for (const std::string& name : info.species) out << "," << name;
    out << "\n";

    int exported = 0;
    // Interior cells only: ghost rows/columns must not be exported,
    // otherwise a re-import maps them onto interior cells and
    // overwrites the physical state (observed: injector-face hot
    // gas replaced by the 700 K background ghost copy).
    for (int i = 1; i < Nx - 1; ++i)
    for (int j = 1; j < Ny - 1; ++j)
    for (int k = 1; k < Nz - 1; ++k)
    // (z is cyclic: k=0 and k=Nz-1 are periodic ghosts)
    {
        const auto& w = dom.W(i, j, k);
        const auto& q = dom.Q(i, j, k);
        const auto center = mesh.cellCenter(i, j, k);
        out << center.x << "," << center.y << "," << center.z << ","
            << w.p << "," << w.T << ","
            << w.U[0] << "," << w.U[1] << "," << w.U[2] << ","
            << q.rho << "," << q.e;
        for (double y : w.Y) out << "," << y;
        out << "\n";
        ++exported;
    }
    std::cout << "exported " << exported << " cells -> " << outputPath << "\n";
    return 0;
}
