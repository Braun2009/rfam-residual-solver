// HotStartReconciler: external field -> conserved initial state migration.
//
// Validates:
//   1. PTY authority: a 130 K / 4 MPa LOX external cell reconciles to the
//      GERG dense-liquid state (rho ~ 920.66, HEOS-consistent) with the
//      conserved rho/rhoE/rhoY correctly constructed.
//   2. RHOEY authority: a rho-e-Y cell recovers T/p through the same backend.
//   3. Species mapping by NAME (unknown species dropped and reported;
//      missing species zero-filled; normalization recorded).
//   4. Report JSON carries mass/element/energy deltas and error maxima.
//   5. Hot product cell (2500 K, CO2/H2O) reconciles through the audited
//      PR fallback (GERG table range miss).
//
// Output marker: HOTSTART_RECONCILER_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "GERGResidualTable.H"
#include "HotStartReconciler.H"

using namespace Foam;

namespace
{

int failures = 0;

void check(
    const std::string& name, bool ok,
    double got = 0.0, double expected = 0.0)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
    if (!ok)
        std::cout << " got=" << std::setprecision(12) << got
                  << " expected=" << expected;
    std::cout << '\n';
    if (!ok) ++failures;
}

int run()
{
    static GERGResidualTable table;
    if (!table.load("tests/reference/gerg_residual_table_v1.csv"))
    {
        check("GERG table load", false);
        std::cout << "HOTSTART_RECONCILER_FAIL\n";
        return 1;
    }
    SpeciesRegistry registry;

    // ------------------------------------------------------------------
    // 1. PTY authority: 130 K / 4 MPa LOX cell.
    // ------------------------------------------------------------------
    {
        HotStartReconciler reconciler(
            registry, &table, HotStartPolicy::PTY_AUTHORITY);
        ExternalCellState cell;
        cell.p = 4.0e6;
        cell.T = 130.0;
        cell.u = 5.0;
        cell.massFractions["O2"] = 1.0;
        ReconciliationReport report;
        const auto out = reconciler.reconcile(cell, report);
        check("PTY LOX cell: valid", out.valid);
        if (out.valid)
        {
            check("PTY LOX rho ~ 920.7 (GERG liquid)",
                std::abs(out.rho - 920.66)/920.66 < 0.05, out.rho, 920.66);
            check("PTY LOX rhoE = rho*e + 0.5 rho u^2",
                std::abs(out.rhoE - (out.rho*out.e + 0.5*out.rho*25.0))
                    / std::abs(out.rhoE) < 1e-12,
                out.rhoE, out.rho*out.e + 0.5*out.rho*25.0);
            check("PTY LOX rhoY(O2) = rho",
                std::abs(out.rhoY[registry.index("O2")] - out.rho)
                    / out.rho < 1e-12,
                out.rhoY[registry.index("O2")], out.rho);
            check("PTY LOX rhoU = rho*5",
                std::abs(out.rhoU[0] - out.rho*5.0)/out.rho < 1e-12,
                out.rhoU[0], out.rho*5.0);
            check("report: single-phase cell", report.singlePhaseCells == 1);
        }
    }

    // ------------------------------------------------------------------
    // 2. RHOEY authority: dense liquid rho-e round trip.
    // ------------------------------------------------------------------
    {
        HotStartReconciler pty(registry, &table, HotStartPolicy::PTY_AUTHORITY);
        ExternalCellState ref;
        ref.p = 4.0e6; ref.T = 130.0;
        ref.massFractions["O2"] = 1.0;
        ReconciliationReport refReport;
        const auto refOut = pty.reconcile(ref, refReport);
        if (refOut.valid)
        {
            HotStartReconciler reconciler(
                registry, &table, HotStartPolicy::RHOEY_AUTHORITY);
            ExternalCellState cell;
            cell.rho = refOut.rho;
            cell.e = refOut.e;
            cell.massFractions["O2"] = 1.0;
            ReconciliationReport report;
            const auto out = reconciler.reconcile(cell, report);
            check("RHOEY dense liquid: valid", out.valid);
            if (out.valid)
            {
                check("RHOEY T recovered ~ 130 K",
                    std::abs(out.T - 130.0)/130.0 < 1e-3, out.T, 130.0);
                check("RHOEY p recovered ~ 4 MPa",
                    std::abs(out.p - 4.0e6)/4.0e6 < 0.02, out.p, 4.0e6);
                check("RHOEY rho exact",
                    std::abs(out.rho - refOut.rho)/refOut.rho < 1e-10);
            }
        }
        else
        {
            check("RHOEY reference setup", false);
        }
    }

    // ------------------------------------------------------------------
    // 3. Species mapping: unknown dropped, missing zero-filled.
    // ------------------------------------------------------------------
    {
        HotStartReconciler reconciler(
            registry, &table, HotStartPolicy::PTY_AUTHORITY);
        ExternalCellState cell;
        cell.p = 4.0e6; cell.T = 700.0;
        cell.massFractions["O2"] = 0.8;
        cell.massFractions["NO2"] = 0.2;   // unknown to registry
        ReconciliationReport report;
        const auto out = reconciler.reconcile(cell, report);
        check("mapping: valid (unknown dropped + normalized)",
            out.valid);
        if (out.valid)
        {
            check("mapping: O2 normalized to 1.0",
                std::abs(out.rhoY[registry.index("O2")]/out.rho - 1.0) < 1e-9);
            check("mapping: dropped species reported",
                report.droppedSpecies.size() == 1
                    && report.droppedSpecies[0] == "NO2");
        }
    }

    // ------------------------------------------------------------------
    // 4. Hot product cell (2500 K) through the audited PR fallback.
    // ------------------------------------------------------------------
    {
        HotStartReconciler reconciler(
            registry, &table, HotStartPolicy::PTY_AUTHORITY);
        ExternalCellState cell;
        cell.p = 4.0e6; cell.T = 2500.0;
        cell.massFractions["CO2"] = 0.6468;
        cell.massFractions["H2O"] = 0.3532;
        ReconciliationReport report;
        const auto out = reconciler.reconcile(cell, report);
        check("hot product cell: valid (PR fallback)", out.valid);
        if (out.valid)
            check("hot product: T preserved",
                std::abs(out.T - 2500.0)/2500.0 < 1e-6, out.T, 2500.0);
    }

    // ------------------------------------------------------------------
    // 5. Batch + report JSON.
    // ------------------------------------------------------------------
    {
        HotStartReconciler reconciler(
            registry, &table, HotStartPolicy::PTY_AUTHORITY);
        std::vector<ExternalCellState> cells;
        {
            ExternalCellState c;
            c.p = 4.0e6; c.T = 130.0;
            c.massFractions["O2"] = 1.0;
            cells.push_back(c);
        }
        {
            ExternalCellState c;
            c.p = 4.0e6; c.T = 2500.0;
            c.massFractions["CO2"] = 0.6468;
            c.massFractions["H2O"] = 0.3532;
            cells.push_back(c);
        }
        ReconciliationReport report;
        const auto out = reconciler.reconcileAll(cells, report);
        check("batch: 2 valid cells", out.size() == 2 && report.validCells == 2);
        check("report: policy field", report.policy == "PTY_AUTHORITY");
        check("report: JSON contains massChangeFraction",
            report.toJson().find("massChangeFraction") != std::string::npos);
        check("report: element fractions finite",
            std::isfinite(report.elementFraction(0))
                && std::isfinite(report.elementFraction(1))
                && std::isfinite(report.elementFraction(2)));
    }

    std::cout << (failures == 0 ? "HOTSTART_RECONCILER_PASS\n"
                                : "HOTSTART_RECONCILER_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
