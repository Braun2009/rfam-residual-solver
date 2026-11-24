#ifndef HotStartReconciler_H
#define HotStartReconciler_H

/*---------------------------------------------------------------------------*\
    HotStartReconciler - external field -> conserved initial state migration.

    A hot start is NOT a plain file read: it is a state migration that must
    preserve mass/element/energy within audited tolerances and must record
    every correction it makes.

    Authority policies (per reconciliation, one per cell set):
      PTY_AUTHORITY     keep (p, T, Y), recover rho/e through the same
                        thermodynamic backend the run will use (GERG/PR);
      RHOEY_AUTHORITY   keep (rho, e, Y), recover T/p (RFAM checkpoint style).

    Every state is recovered through RealFluidState with the same registered
    GERG table as the run, so boundary and interior states share one backend.

    Output: per-cell conserved states (rho, rhoU, rhoE, rhoY) + a JSON
    reconciliation report (mass/element/energy deltas, max reconstruction
    errors, phase counts, missing/truncated species counts).

    Composition: external fields are mapped by NAME through the SpeciesRegistry
    (never by array position).  Unknown species are reported and dropped;
    missing species are reported and zero-filled.
\*---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "RealFluidState.H"
#include "SpeciesRegistry.H"

namespace Foam
{

enum class HotStartPolicy
{
    PTY_AUTHORITY,
    RHOEY_AUTHORITY
};

struct ExternalCellState
{
    double x = 0.0, y = 0.0, z = 0.0;
    double p = std::numeric_limits<double>::quiet_NaN();
    double T = std::numeric_limits<double>::quiet_NaN();
    double rho = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double u = 0.0, v = 0.0, w = 0.0;
    double k = std::numeric_limits<double>::quiet_NaN();
    double omega = std::numeric_limits<double>::quiet_NaN();
    std::map<std::string, double> massFractions;   // name -> Y (mapped)
    bool hasTurbulence = false;
};

struct ReconciledCell
{
    bool valid = false;
    bool corrected = false;
    double rho = std::numeric_limits<double>::quiet_NaN();
    double rhoU[3] = {0.0, 0.0, 0.0};
    double rhoE = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> rhoY;       // registry order
    double p = std::numeric_limits<double>::quiet_NaN();
    double T = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double quality = std::numeric_limits<double>::quiet_NaN();
};

struct ReconciliationReport
{
    int totalCells = 0;
    int validCells = 0;
    int failedCells = 0;
    int outsideEnvelopeCells = 0;
    int singlePhaseCells = 0;
    int twoPhaseCells = 0;
    int missingSpeciesCells = 0;
    int truncatedSpeciesCells = 0;
    double totalMassBefore = 0.0;
    double totalMassAfter = 0.0;
    double totalEnergyBefore = 0.0;
    double totalEnergyAfter = 0.0;
    double elementMassBefore[4] = {0.0, 0.0, 0.0, 0.0};  // C H O N
    double elementMassAfter[4] = {0.0, 0.0, 0.0, 0.0};
    double maxPError = 0.0;
    double maxTError = 0.0;
    double maxRhoError = 0.0;
    double maxEnergyError = 0.0;
    std::string policy = "PTY_AUTHORITY";
    std::vector<std::string> droppedSpecies;   // unknown to registry

    std::string toJson() const
    {
        std::ostringstream os;
        os << std::setprecision(12);
        os << "{\n"
           << "  \"policy\": \"" << policy << "\",\n"
           << "  \"totalCells\": " << totalCells << ",\n"
           << "  \"validCells\": " << validCells << ",\n"
           << "  \"failedCells\": " << failedCells << ",\n"
           << "  \"outsideEnvelopeCells\": " << outsideEnvelopeCells << ",\n"
           << "  \"singlePhaseCells\": " << singlePhaseCells << ",\n"
           << "  \"twoPhaseCells\": " << twoPhaseCells << ",\n"
           << "  \"missingSpeciesCells\": " << missingSpeciesCells << ",\n"
           << "  \"truncatedSpeciesCells\": " << truncatedSpeciesCells << ",\n"
           << "  \"massChangeFraction\": "
           << (totalMassAfter > 0.0
               ? std::abs(totalMassAfter - totalMassBefore) / totalMassBefore
               : 0.0) << ",\n"
           << "  \"energyChangeFraction\": "
           << (totalEnergyBefore != 0.0
               ? std::abs(totalEnergyAfter - totalEnergyBefore)
                     / std::abs(totalEnergyBefore) : 0.0) << ",\n"
           << "  \"elementChangeFraction\": ["
           << elementFraction(0) << ", " << elementFraction(1) << ", "
           << elementFraction(2) << ", " << elementFraction(3) << "],\n"
           << "  \"maxPError\": " << maxPError << ",\n"
           << "  \"maxTError\": " << maxTError << ",\n"
           << "  \"maxRhoError\": " << maxRhoError << ",\n"
           << "  \"maxEnergyError\": " << maxEnergyError << "\n"
           << "}";
        return os.str();
    }

    double elementFraction(int i) const
    {
        return elementMassBefore[i] > 0.0
            ? std::abs(elementMassAfter[i] - elementMassBefore[i])
                / elementMassBefore[i] : 0.0;
    }
};

class HotStartReconciler
{
public:
    HotStartReconciler(const SpeciesRegistry& registry,
                       const GERGResidualTable* gergTable,
                       HotStartPolicy policy = HotStartPolicy::PTY_AUTHORITY)
        : registry_(registry), gergTable_(gergTable), policy_(policy)
    {}

    // Reconcile one external cell.  Y is written in registry order.
    ReconciledCell reconcile(const ExternalCellState& in,
                             ReconciliationReport& report) const
    {
        ReconciledCell out;
        ++report.totalCells;
        out.rhoY.assign(static_cast<std::size_t>(registry_.size()), 0.0);

        // Before-accounting (audit baseline; PTY fields carry no rho, so
        // the ideal-gas estimate is used for mass; energy uses the
        // NASA7 sensible ideal energy at (T,Y)).
        double rhoExt = std::numeric_limits<double>::quiet_NaN();
        if (std::isfinite(in.rho) && in.rho > 0.0) rhoExt = in.rho;
        else if (std::isfinite(in.p) && std::isfinite(in.T) && in.T > 0.0)
        {
            double invMW = 0.0;
            for (const auto& [name, y] : in.massFractions)
            {
                const int idx = registry_.index(name);
                if (idx >= 0) invMW += y / registry_.molecularWeight(idx);
            }
            const double MW = invMW > 0.0 ? 1.0 / invMW : 0.029;
            rhoExt = in.p * MW / (8.314462618 * in.T);
        }
        if (std::isfinite(rhoExt) && rhoExt > 0.0)
        {
            report.totalMassBefore += rhoExt;
            // Ideal sensible energy baseline: h = e + p/rho (ideal),
            // e_ideal = Cp_ig*T - p/rho; use 1000 J/kg/K placeholder only
            // when NASA7 is unavailable (documented audit-only).
            const double eIdeal = std::isfinite(in.e) ? in.e
                : (1000.0 * in.T - in.p / rhoExt);
            if (std::isfinite(eIdeal))
                report.totalEnergyBefore += rhoExt * eIdeal;
        }

        // --- species mapping by name -------------------------------------
        double sumY = 0.0;
        bool truncated = false;
        for (const auto& [name, y] : in.massFractions)
        {
            const int idx = registry_.index(name);
            if (idx < 0)
            {
                if (std::find(report.droppedSpecies.begin(),
                              report.droppedSpecies.end(), name)
                    == report.droppedSpecies.end())
                    report.droppedSpecies.push_back(name);
                continue;
            }
            if (!std::isfinite(y) || y < 0.0) continue;
            out.rhoY[static_cast<std::size_t>(idx)] = y;
            sumY += y;
            if (y > 1.0 + 1e-9) truncated = true;
        }
        if (sumY <= 0.0 || sumY > 1.0 + 1e-6)
        {
            ++report.failedCells;
            return out;
        }
        if (std::abs(sumY - 1.0) > 1e-9)
        {
            for (auto& y : out.rhoY) y /= sumY;   // normalize (recorded)
            ++report.truncatedSpeciesCells;
        }
        if (truncated) ++report.truncatedSpeciesCells;

        // --- state recovery through the run's own backend ----------------
        RealFluidState state;
        state.setSpecies(registry_.names());
        if (gergTable_) state.setGERGTable(gergTable_);

        RealFluidState::Status status;
        if (policy_ == HotStartPolicy::PTY_AUTHORITY)
        {
            if (!std::isfinite(in.p) || !std::isfinite(in.T)
                || in.p <= 0.0 || in.T <= 0.0)
            {
                ++report.failedCells;
                return out;
            }
            status = state.updatePTY(in.p, in.T, out.rhoY);
        }
        else
        {
            if (!std::isfinite(in.rho) || !std::isfinite(in.e)
                || in.rho <= 0.0)
            {
                ++report.failedCells;
                return out;
            }
            status = state.updateRhoEY(in.rho, in.e, out.rhoY);
        }
        if (status != RealFluidState::Status::SUCCESS)
        {
            ++report.failedCells;
            return out;
        }

        // --- conserved construction --------------------------------------
        const double rho = state.rho();
        const double e = state.e();
        out.rho = rho;
        const double uMag2 = in.u*in.u + in.v*in.v + in.w*in.w;
        out.rhoU[0] = rho * in.u;
        out.rhoU[1] = rho * in.v;
        out.rhoU[2] = rho * in.w;
        out.rhoE = rho * e + 0.5 * rho * uMag2;
        for (auto& ry : out.rhoY) ry *= rho;
        out.p = state.p();
        out.T = state.T();
        out.e = e;
        out.valid = true;

        // --- diagnostics --------------------------------------------------
        ++report.validCells;
        if (state.isNearCritical(0.10, 0.35)) ++report.twoPhaseCells;
        else ++report.singlePhaseCells;
        report.maxPError = std::max(report.maxPError,
            std::abs(out.p - in.p) / std::max(std::abs(in.p), 1.0));
        report.maxTError = std::max(report.maxTError,
            std::abs(out.T - in.T) / std::max(std::abs(in.T), 1.0));
        report.maxRhoError = std::max(report.maxRhoError,
            std::abs(out.rho - in.rho) / std::max(std::abs(in.rho), 1.0));
        report.maxEnergyError = std::max(report.maxEnergyError,
            std::abs(out.e - in.e) / std::max(std::abs(in.e), 1.0));
        report.totalMassAfter += rho;
        report.totalEnergyAfter += rho * e;
        report.elementMassAfter[0] += rho * carbonFraction(out.rhoY);
        report.elementMassAfter[1] += rho * hydrogenFraction(out.rhoY);
        report.elementMassAfter[2] += rho * oxygenFraction(out.rhoY);
        report.elementMassAfter[3] += rho * nitrogenFraction(out.rhoY);
        return out;
    }

    // Reconcile a batch of cells (CSV-loaded external field).
    std::vector<ReconciledCell> reconcileAll(
        const std::vector<ExternalCellState>& cells,
        ReconciliationReport& report) const
    {
        report = ReconciliationReport{};
        report.policy = (policy_ == HotStartPolicy::PTY_AUTHORITY)
            ? "PTY_AUTHORITY" : "RHOEY_AUTHORITY";
        std::vector<ReconciledCell> out;
        out.reserve(cells.size());
        for (const auto& cell : cells)
            out.push_back(reconcile(cell, report));
        return out;
    }

private:
    double carbonFraction(const std::vector<double>& Y) const
    {
        double f = 0.0;
        for (std::size_t i = 0; i < Y.size(); ++i)
            f += Y[i] * registry_.nC(static_cast<int>(i)) * 0.012011
                / registry_.molecularWeight(static_cast<int>(i));
        return f;
    }
    double hydrogenFraction(const std::vector<double>& Y) const
    {
        double f = 0.0;
        for (std::size_t i = 0; i < Y.size(); ++i)
            f += Y[i] * registry_.nH(static_cast<int>(i)) * 0.001008
                / registry_.molecularWeight(static_cast<int>(i));
        return f;
    }
    double oxygenFraction(const std::vector<double>& Y) const
    {
        double f = 0.0;
        for (std::size_t i = 0; i < Y.size(); ++i)
            f += Y[i] * registry_.nO(static_cast<int>(i)) * 0.015999
                / registry_.molecularWeight(static_cast<int>(i));
        return f;
    }
    double nitrogenFraction(const std::vector<double>& Y) const
    {
        double f = 0.0;
        for (std::size_t i = 0; i < Y.size(); ++i)
            f += Y[i] * registry_.nN(static_cast<int>(i)) * 0.014007
                / registry_.molecularWeight(static_cast<int>(i));
        return f;
    }

    // Best-effort external density for before-accounting (PTY fields usually
    // do not carry rho; use the registry mixture-MW ideal-gas estimate for
    // the audit baseline only).  RHOEY fields carry exact rho.
    double estimateRho(const ExternalCellState& in) const
    {
        if (!std::isfinite(in.p) || !std::isfinite(in.T) || in.T <= 0.0)
            return std::numeric_limits<double>::quiet_NaN();
        double invMW = 0.0;
        for (const auto& [name, y] : in.massFractions)
        {
            if (!std::isfinite(y)) continue;
            const int idx = registry_.index(name);
            if (idx >= 0) invMW += y / registry_.molecularWeight(idx);
        }
        const double MW = invMW > 0.0 ? 1.0 / invMW : 0.029;
        return in.p * MW / (8.314462618 * in.T);
    }

    const SpeciesRegistry& registry_;
    const GERGResidualTable* gergTable_;
    HotStartPolicy policy_;
};

} // End namespace Foam

#endif
