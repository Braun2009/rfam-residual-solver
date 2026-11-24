// Verify the property table against the patched CoolProp backend and use it as
// an independent route for the equilibrium two-phase sound speed.
//  1. off-grid (p,T,Y) states: direct PT flash vs trilinear table interpolation
//  2. (rho,e,Y) states: direct rho-e recovery vs table inversion
//  3. HEM isentropic sound speed: direct CONVERGED_RHOE_ISENTROPE vs the table
#include "CoolPropThermoModel.H"
#include "property_table.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <vector>

using Foam::thermo::EquilibriumSeed;
using Foam::thermo::EquilibriumState;

namespace
{
struct Inversion
{
    bool ok = false;
    double p = 0., T = 0., rho = 0., e = 0.;
};

Inversion invertRhoE(const rfamtable::Table& table, double yO2, double rhoT,
                     double eT, double pSeed, double tSeed)
{
    Inversion result;
    double p = pSeed, t = tSeed;
    rfamtable::Node base;
    const auto residual = [&](double pv, double tv, double& dr, double& de)
    {
        rfamtable::Node n;
        if (!table.interpolate(yO2, tv, pv, n)) return false;
        dr = n.rho - rhoT; de = n.e - eT; return true;
    };
    const auto jacobian = [&](double pv, double tv, double J[2][2])
    {
        const double dp = std::max(1., std::fabs(pv)) * 1e-4;
        const double dt = std::max(1., std::fabs(tv)) * 1e-4;
        double r1, e1, r2, e2;
        if (!residual(pv + dp, tv, r1, e1) || !residual(pv - dp, tv, r2, e2)) return false;
        J[0][0] = (r1 - r2) / (2 * dp); J[1][0] = (e1 - e2) / (2 * dp);
        if (!residual(pv, tv + dt, r1, e1) || !residual(pv, tv - dt, r2, e2)) return false;
        J[0][1] = (r1 - r2) / (2 * dt); J[1][1] = (e1 - e2) / (2 * dt);
        return true;
    };
    for (int iter = 0; iter < 60; ++iter)
    {
        double dr, de;
        if (!residual(p, t, dr, de)) break;
        if (std::fabs(dr) <= 1e-9 * std::max(std::fabs(rhoT), 1.)
            && std::fabs(de) <= 1e-9 * std::max(std::fabs(eT), 1.))
        { result = {true, p, t, rhoT, eT}; return result; }
        double J[2][2];
        if (!jacobian(p, t, J)) break;
        const double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (!std::isfinite(det) || std::fabs(det) < 1e-30) break;
        const double dp = (-dr * J[1][1] + J[0][1] * de) / det;
        const double dt = (-J[0][0] * de + dr * J[1][0]) / det;
        const double norm = std::fabs(dr) / rhoT + std::fabs(de) / std::max(std::fabs(eT), 1.);
        bool advanced = false;
        for (double scale = 1.; scale >= 1. / 1024.; scale *= .5)
        {
            const double pn = p + scale * dp, tn = t + scale * dt;
            if (pn <= 1. || tn <= 20.) continue;
            double r1, e1;
            if (!residual(pn, tn, r1, e1)) continue;
            if (std::fabs(r1) / rhoT + std::fabs(e1) / std::max(std::fabs(eT), 1.) < norm)
            { p = pn; t = tn; advanced = true; break; }
        }
        if (!advanced) break;
    }
    return result;
}

bool tableIsentropePoint(const rfamtable::Table& table, double yO2, double p,
                         double sTarget, double& t, double& rho)
{
    for (int iter = 0; iter < 60; ++iter)
    {
        rfamtable::Node n;
        if (!table.interpolate(yO2, t, p, n)) return false;
        const double err = n.s - sTarget;
        if (std::fabs(err) <= 1e-10 * std::max(std::fabs(sTarget), 1.)) { rho = n.rho; return true; }
        const double dt = std::max(1., t) * 1e-4;
        rfamtable::Node hp, hm;
        if (!table.interpolate(yO2, t + dt, p, hp) || !table.interpolate(yO2, t - dt, p, hm)) return false;
        const double slope = (hp.s - hm.s) / (2 * dt);
        if (!std::isfinite(slope) || slope <= 0.) return false;
        const double step = std::clamp(-err / slope, -.1 * t, .1 * t);
        bool advanced = false;
        for (double scale = 1.; scale >= 1. / 1024.; scale *= .5)
        {
            const double tn = t + scale * step;
            if (tn <= 20.) continue;
            rfamtable::Node trial;
            if (!table.interpolate(yO2, tn, p, trial)) continue;
            if (std::fabs(trial.s - sTarget) < std::fabs(err)) { t = tn; advanced = true; break; }
        }
        if (!advanced) return false;
    }
    return false;
}
}

int main(int argc, char** argv)
{
    std::string tablePath = "work/a406_rootcause/property_table.bin";
    std::string out = "work/a406_rootcause/property_table_verification.json";
    int samples = 2000, recoveries = 300;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--table" && i + 1 < argc) tablePath = argv[++i];
        else if (a == "--out" && i + 1 < argc) out = argv[++i];
        else if (a == "--samples" && i + 1 < argc) samples = std::atoi(argv[++i]);
        else if (a == "--recoveries" && i + 1 < argc) recoveries = std::atoi(argv[++i]);
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    const auto table = rfamtable::Table::read(tablePath);
    auto model = Foam::thermo::makeCoolPropThermoModel();
    std::printf("VERIFY model=%s table=%s nodes=%zu\n",
                model->modelId().c_str(), tablePath.c_str(), table.nodes.size());
    std::mt19937 rng(20260917u);
    std::uniform_real_distribution<double> u01(0., 1.);

    // --- 1. off-grid direct PT vs table interpolation ---------------------
    struct OffGrid
    {
        int tried = 0, misses = 0, directFailures = 0;
        double maxRelRho = 0., maxRelE = 0., maxRelS = 0.;
        double sumRelRho = 0., sumRelE = 0., sumRelS = 0.;
    } off;
    for (int i = 0; i < samples; ++i)
    {
        const double yO2 = table.Y.front() + u01(rng) * (table.Y.back() - table.Y.front());
        const double tv = table.T.front() + u01(rng) * (table.T.back() - table.T.front());
        const double pv = table.p.front() + u01(rng) * (table.p.back() - table.p.front());
        rfamtable::Node n;
        if (!table.interpolate(yO2, tv, pv, n)) { ++off.misses; continue; }
        const auto r = model->initializePTY(pv, tv, {yO2, 1. - yO2});
        if (!r.success() || !r.state) { ++off.directFailures; continue; }
        ++off.tried;
        const double dRho = std::fabs(n.rho - r.state->rho) / r.state->rho;
        const double dE = std::fabs(n.e - r.state->e) / std::max(std::fabs(r.state->e), 1e4);
        const double dS = std::fabs(n.s - r.state->s) / std::max(std::fabs(r.state->s), 1e4);
        off.maxRelRho = std::max(off.maxRelRho, dRho);
        off.maxRelE = std::max(off.maxRelE, dE);
        off.maxRelS = std::max(off.maxRelS, dS);
        off.sumRelRho += dRho; off.sumRelE += dE; off.sumRelS += dS;
    }

    // --- 2. (rho,e) recovery: direct vs table inversion --------------------
    struct Recovery
    {
        int tried = 0, directFailures = 0, tableFailures = 0;
        double maxRelP = 0., maxRelT = 0., maxRelRho = 0., maxRelE = 0.;
    } rec;
    std::uniform_int_distribution<std::size_t> pickY(0, table.Y.size() - 1);
    std::uniform_int_distribution<std::size_t> pickT(1, table.T.size() - 2);
    std::uniform_int_distribution<std::size_t> pickP(1, table.p.size() - 2);
    for (int i = 0; i < recoveries; ++i)
    {
        const std::size_t iy = pickY(rng), it = pickT(rng), ip = pickP(rng);
        const rfamtable::Node& n = table.nodes[table.index(iy, it, ip)];
        if (n.phase < 0) continue;
        const double rhoT = n.rho * (1. + 0.002 * (u01(rng) - .5));
        const double eT = n.e * (1. + 0.002 * (u01(rng) - .5));
        EquilibriumSeed seed;
        seed.pressure = table.p[ip]; seed.temperature = table.T[it];
        // table phase codes are CoolProp iphase_* values (two-phase = 5)
        seed.previousRegion = (n.phase == 5)
            ? Foam::thermo::EquilibriumRegion::TWO_PHASE
            : Foam::thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE;
        const auto direct = model->recoverRhoE(rhoT, eT, {table.Y[iy], 1. - table.Y[iy]}, &seed);
        const auto inv = invertRhoE(table, table.Y[iy], rhoT, eT, table.p[ip], table.T[it]);
        ++rec.tried;
        if (!direct.success()) { ++rec.directFailures; continue; }
        if (!inv.ok) { ++rec.tableFailures; continue; }
        rec.maxRelP = std::max(rec.maxRelP, std::fabs(inv.p - direct.state->p) / direct.state->p);
        rec.maxRelT = std::max(rec.maxRelT, std::fabs(inv.T - direct.state->T) / direct.state->T);
        rec.maxRelRho = std::max(rec.maxRelRho, std::fabs(inv.rho - direct.state->rho) / direct.state->rho);
        rec.maxRelE = std::max(rec.maxRelE, std::fabs(inv.e - direct.state->e) / std::max(std::fabs(direct.state->e), 1e4));
    }

    // --- 3. HEM sound speed: direct vs table isentrope ---------------------
    struct SoundRow
    {
        double yO2, T, p, direct, table, relDiff;
        bool ok;
    };
    // Backend regression references (Y,T,p) + the cell-9272 state.
    const std::vector<std::array<double, 3>> refs = {
        {0.99, 150., 4.0e6}, {0.1, 250., 4.0e6}, {0.9, 170., 4.0e6},
        {0.73, 250., 4.0e6}, {0.9162460633, 189.802, 4.09773e6}};
    std::vector<SoundRow> sound;
    for (const auto& ref : refs)
    {
        SoundRow row{ref[0], ref[1], ref[2], 0., 0., 0., false};
        const auto r = model->initializePTY(ref[2], ref[1], {ref[0], 1. - ref[0]});
        rfamtable::Node n;
        if (r.success() && r.state && table.interpolate(ref[0], ref[1], ref[2], n))
        {
            row.direct = r.state->acoustic.soundSpeed();
            if (n.a2 > 0.)
                std::printf("  TABLE_A2 Y=%.6g T=%.3f p=%.6g tableA=%.6g directA=%.6g rel=%.3e\n",
                            ref[0], ref[1], ref[2], std::sqrt(n.a2), row.direct,
                            std::fabs(std::sqrt(n.a2) - row.direct) / row.direct);
            const double s0 = n.s;
            double a2Sum = 0.;
            bool ok = true;
            for (double eps : {1e-3, 5e-4, 2.5e-4})
            {
                double tl = ref[1], th = ref[1], rl = 0., rh = 0.;
                if (!tableIsentropePoint(table, ref[0], ref[2] * (1. - eps), s0, tl, rl)
                    || !tableIsentropePoint(table, ref[0], ref[2] * (1. + eps), s0, th, rh)
                    || !(rh > rl))
                { ok = false; break; }
                a2Sum += (ref[2] * (1. + eps) - ref[2] * (1. - eps)) / (rh - rl);
            }
            if (ok)
            {
                row.table = std::sqrt(a2Sum / 3.);
                row.relDiff = std::fabs(row.table - row.direct) / row.direct;
                row.ok = true;
            }
        }
        sound.push_back(row);
    }

    // --- 4. independent thermodynamic identity for the sound speed ---------
    // a^2 = (dp/drho)_T + (T/rho^2) * (dp/dT)_rho^2 / Cv, computed entirely
    // from table derivatives (no isentrope integration).
    struct IdentityRow { double yO2, T, p, direct, identity, relIdentity; bool ok; };
    std::vector<IdentityRow> identity;
    for (const auto& ref : refs)
    {
        IdentityRow row{ref[0], ref[1], ref[2], 0., 0., 0., false};
        rfamtable::Node n;
        const auto state = model->initializePTY(ref[2], ref[1], {ref[0], 1. - ref[0]});
        if (!table.interpolate(ref[0], ref[1], ref[2], n)
            || !state.success() || !state.state) { identity.push_back(row); continue; }
        row.direct = state.state->acoustic.soundSpeed();
        const double dp = ref[2] * 1e-3, dt = ref[1] * 1e-3;
        rfamtable::Node pp, pm, tp, tm;
        if (!table.interpolate(ref[0], ref[1], ref[2] + dp, pp)
            || !table.interpolate(ref[0], ref[1], ref[2] - dp, pm)
            || !table.interpolate(ref[0], ref[1] + dt, ref[2], tp)
            || !table.interpolate(ref[0], ref[1] - dt, ref[2], tm))
        { identity.push_back(row); continue; }
        const double drho_dp_T = (pp.rho - pm.rho) / (2 * dp);
        const double drho_dT_p = (tp.rho - tm.rho) / (2 * dt);
        const double de_dT_p = (tp.e - tm.e) / (2 * dt);
        const double de_dp_T = (pp.e - pm.e) / (2 * dp);
        const double dp_drho_T = 1. / drho_dp_T;
        const double dp_dT_rho = -drho_dT_p / drho_dp_T;
        const double cv = de_dT_p - de_dp_T * dp_dT_rho;
        const double a2 = dp_drho_T + (ref[1] / (n.rho * n.rho)) * dp_dT_rho * dp_dT_rho / cv;
        if (std::isfinite(a2) && a2 > 0.)
        {
            row.identity = std::sqrt(a2);
            row.relIdentity = std::fabs(row.identity - row.direct) / row.direct;
            row.ok = true;
        }
        identity.push_back(row);
    }

    std::ofstream j(out);
    j << std::setprecision(10);
    j << "{\n \"model\": \"" << model->modelId() << "\",\n"
      << " \"table\": \"" << tablePath << "\",\n"
      << " \"offGrid\": {\"samples\": " << samples << ", \"interpolated\": " << off.tried
      << ", \"misses\": " << off.misses << ", \"directFlashFailures\": " << off.directFailures
      << ", \"meanRelRho\": " << (off.tried ? off.sumRelRho / off.tried : 0.)
      << ", \"meanRelE\": " << (off.tried ? off.sumRelE / off.tried : 0.)
      << ", \"meanRelS\": " << (off.tried ? off.sumRelS / off.tried : 0.)
      << ", \"maxRelRho\": " << off.maxRelRho << ", \"maxRelE\": " << off.maxRelE
      << ", \"maxRelS\": " << off.maxRelS << "},\n"
      << " \"recovery\": {\"attempted\": " << rec.tried << ", \"directFailures\": " << rec.directFailures
      << ", \"tableFailures\": " << rec.tableFailures
      << ", \"maxRelP\": " << rec.maxRelP << ", \"maxRelT\": " << rec.maxRelT
      << ", \"maxRelRho\": " << rec.maxRelRho << ", \"maxRelE\": " << rec.maxRelE << "},\n"
      << " \"soundSpeed\": [\n";
    for (std::size_t i = 0; i < sound.size(); ++i)
    {
        const auto& s = sound[i];
        j << "  {\"Y_O2\": " << s.yO2 << ", \"T\": " << s.T << ", \"p\": " << s.p
          << ", \"direct\": " << s.direct << ", \"table\": " << s.table
          << ", \"relDiff\": " << s.relDiff << ", \"ok\": " << (s.ok ? "true" : "false") << "}"
          << (i + 1 == sound.size() ? "\n" : ",\n");
    }
    j << " ]\n}\n";
    std::printf("VERIFY_OFFGRID interp=%d/%d misses=%d directFail=%d meanRel rho=%.2e e=%.2e s=%.2e maxRel rho=%.2e\n",
                off.tried, samples, off.misses, off.directFailures,
                off.tried ? off.sumRelRho / off.tried : 0., off.tried ? off.sumRelE / off.tried : 0.,
                off.tried ? off.sumRelS / off.tried : 0., off.maxRelRho);
    std::printf("VERIFY_RECOVERY tried=%d directFail=%d tableFail=%d maxRel p=%.2e T=%.2e rho=%.2e\n",
                rec.tried, rec.directFailures, rec.tableFailures, rec.maxRelP, rec.maxRelT, rec.maxRelRho);
    for (const auto& s : sound)
        std::printf("VERIFY_SOUND Y=%.6g T=%.3f p=%.6g direct=%.6g table=%.6g rel=%.2e ok=%d\n",
                    s.yO2, s.T, s.p, s.direct, s.table, s.relDiff, int(s.ok));
    for (const auto& s : identity)
        std::printf("VERIFY_IDENTITY Y=%.6g T=%.3f p=%.6g direct=%.6g identity=%.6g rel=%.2e ok=%d\n",
                    s.yO2, s.T, s.p, s.direct, s.identity, s.relIdentity, int(s.ok));
    {
        std::ofstream ji(out + ".identity.json");
        ji << std::setprecision(10) << "{\n \"grid\": [";
        for (std::size_t i = 0; i < identity.size(); ++i)
        {
            const auto& s = identity[i];
            ji << (i ? ", " : "") << "{\"Y_O2\": " << s.yO2 << ", \"T\": " << s.T
               << ", \"p\": " << s.p << ", \"direct\": " << s.direct
               << ", \"identity\": " << s.identity << ", \"relDiff\": " << s.relIdentity
               << ", \"ok\": " << (s.ok ? "true" : "false") << "}";
        }
        ji << "]}\n";
    }
    std::printf("VERIFY_JSON %s\n", out.c_str());
    return 0;
}
