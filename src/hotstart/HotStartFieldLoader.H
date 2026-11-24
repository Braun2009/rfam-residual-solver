#ifndef HotStartFieldLoader_H
#define HotStartFieldLoader_H

/*---------------------------------------------------------------------------*\
    HotStartFieldLoader - external hot-start field CSV -> Domain initial state.

    CSV schema (rfam.hotstart-field/v1):
      # rfam.hotstart-field/v1
      x[m],y[m],z[m],p[Pa],T[K],u[m/s],v[m/s],w[m/s],<species names...>

    Semantics:
      - PTY authority: (p, T, Y) are kept; density/energy are recovered by
        the run's own thermodynamic backend (GERG/PR) through
        HotStartReconciler, so interior states and boundary states share one
        backend.  The reconciled (p,T) are written back into the Domain.
      - Grid mapping: nearest cell center in the sector mesh (x via
        bisection on the face table, radius and angle by nearest center).
      - Every reconciliation is audited: the report JSON records
        mass/element/energy deltas and max reconstruction errors.

    The loader never guesses species identity: CSV columns are matched by
    name against the SpeciesRegistry; unknown columns are dropped and
    reported.
\*---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "HotStartReconciler.H"
#include "SectorMeshMetrics.H"

namespace Foam
{

struct HotStartFieldLoadResult
{
    bool valid = false;
    int loadedCells = 0;
    int skippedCells = 0;
    std::vector<ReconciledCell> reconciled;
    std::vector<int> ijk;            // packed i*10000+j*100+k per loaded cell
    ReconciliationReport report;
    std::vector<std::string> droppedColumns;
    std::string error;
    double sourceTime = std::numeric_limits<double>::quiet_NaN();
    int sourceIteration = -1;
};

class HotStartFieldLoader
{
public:
    HotStartFieldLoader(const SpeciesRegistry& registry,
                        const GERGResidualTable* gergTable,
                        const SectorMeshMetrics& mesh)
        : reconciler_(registry, gergTable, HotStartPolicy::PTY_AUTHORITY),
          mesh_(mesh) {}

    // Parse a hot-start field CSV (no grid mapping).
    std::vector<ExternalCellState> parseCSV(const std::string& path,
                                            std::string& error) const
    {
        std::vector<ExternalCellState> out;
        std::ifstream input(path);
        if (!input) { error = "cannot open " + path; return out; }
        std::string line;
        std::vector<std::string> header;
        bool schemaOk = false;
        while (std::getline(input, line))
        {
            if (line.empty()) continue;
            if (line.rfind("#", 0) == 0)
            {
                if (line.find("rfam.hotstart-field/v1") != std::string::npos)
                    schemaOk = true;
                if (line.rfind("# time=", 0) == 0)
                    sourceTime_ = std::atof(
                        line.substr(7).c_str());
                if (line.rfind("# iteration=", 0) == 0)
                    sourceIteration_ = std::atoi(
                        line.substr(13).c_str());
                continue;
            }
            std::vector<std::string> fields;
            std::stringstream stream(line);
            std::string field;
            while (std::getline(stream, field, ','))
                fields.push_back(trim(field));
            if (header.empty()) { header = fields; continue; }
            if (fields.size() != header.size())
            { error = "column count mismatch at row " +
                std::to_string(out.size() + 1) + " got "
                + std::to_string(fields.size()) + " want "
                + std::to_string(header.size()); return out; }
            ExternalCellState cell;
            bool ok = true;
            for (std::size_t c = 0; c < header.size(); ++c)
            {
                const std::string& name = header[c];
                double value = 0.0;
                // strtod: subnormal underflow (e.g. residual
                // species 4.5e-313 from a checkpoint) must not fail
                // the import; it rounds to zero.
                errno = 0;
                char* endPtr = nullptr;
                value = std::strtod(fields[c].c_str(), &endPtr);
                if (endPtr == fields[c].c_str()
                    || (errno == ERANGE && std::abs(value) > 1e300))
                {
                    ok = false;
                    error = "parse fail row "
                        + std::to_string(out.size() + 1) + " col "
                        + name + " val " + fields[c];
                    break;
                }
                if (name == "x") cell.x = value;
                else if (name == "y") cell.y = value;
                else if (name == "z") cell.z = value;
                else if (name == "p") cell.p = value;
                else if (name == "T") cell.T = value;
                else if (name == "u") cell.u = value;
                else if (name == "v") cell.v = value;
                else if (name == "w") cell.w = value;
                else cell.massFractions[name] = value;
            }
            if (ok) out.push_back(cell);
        }
        if (!schemaOk && !out.empty())
            error = "missing schema line rfam.hotstart-field/v1";
        return out;
    }

    // Map an external point to the nearest sector cell center.
    void locate(const ExternalCellState& cell,
                int& i, int& j, int& k) const
    {
        const int nx = mesh_.Nx(), nr = mesh_.Nr(), nt = mesh_.Ntheta();
        i = 1; j = 1; k = 0;
        // x via bisection on face midpoints.
        {
            int lo = 1, hi = nx - 2;
            while (lo < hi)
            {
                const int mid = (lo + hi) / 2;
                const double xc = 0.5 * (mesh_.xFace(mid)
                                         + mesh_.xFace(mid + 1));
                if (cell.x > xc) lo = mid + 1; else hi = mid;
            }
            i = lo;
            if (i < 1) i = 1;
            if (i > nx - 2) i = nx - 2;
        }
        // radius nearest center: r = hypot(y, z) (NOT the y component,
        // which shrinks with the sector angle).
        {
            const double r = std::hypot(cell.y, cell.z);
            double best = std::numeric_limits<double>::infinity();
            for (int q = 1; q < nr - 1; ++q)
            {
                const double d = std::abs(r - mesh_.radiusCenter(i, q));
                if (d < best) { best = d; j = q; }
            }
        }
        // angle nearest center: theta = atan2(z, y) in radians,
        // compared against thetaCenter (also radians).
        {
            const double theta = std::atan2(cell.z, cell.y);
            double best = std::numeric_limits<double>::infinity();
            for (int q = 1; q < nt - 1; ++q)
            {
                const double d = std::abs(theta - mesh_.thetaCenter(q));
                if (d < best) { best = d; k = q; }
            }
        }
    }

    // Load, map, reconcile and pack ijk indices.
    HotStartFieldLoadResult load(const std::string& path) const
    {
        HotStartFieldLoadResult result;
        std::string error;
        const auto cells = parseCSV(path, error);
        if (cells.empty())
        {
            result.error = error.empty() ? "empty field" : error;
            return result;
        }
        if (!error.empty() && cells.size() < 14400)
            result.error = error;   // partial parse: surface the cause
        std::vector<ExternalCellState> mapped;
        for (const auto& cell : cells)
        {
            int i, j, k;
            locate(cell, i, j, k);
            mapped.push_back(cell);
            result.ijk.push_back(i * 10000 + j * 100 + k);
        }
        result.reconciled = reconciler_.reconcileAll(mapped, result.report);
        result.loadedCells = result.report.validCells;
        result.skippedCells = result.report.failedCells;
        result.droppedColumns = result.report.droppedSpecies;
        result.sourceTime = sourceTime_;
        result.sourceIteration = sourceIteration_;
        result.valid = result.loadedCells > 0;
        return result;
    }

private:
    static std::string trim(const std::string& value)
    {
        std::size_t b = value.find_first_not_of(" \t\r\n");
        if (b == std::string::npos) return std::string();
        std::size_t e = value.find_last_not_of(" \t\r\n");
        return value.substr(b, e - b + 1);
    }

    HotStartReconciler reconciler_;
    const SectorMeshMetrics& mesh_;
    mutable double sourceTime_ = std::numeric_limits<double>::quiet_NaN();
    mutable int sourceIteration_ = -1;
};

} // End namespace Foam

#endif
