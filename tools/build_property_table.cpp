// Build the O2/C3H8 equilibrium property table with the patched native
// CoolProp backend. Verification artifact: it is not read by the solver.
// The PT flash is called directly on the CoolProp state (the wrapper's publish
// would additionally compute the HEM isentrope sound speed for every two-phase
// node, which is far too slow for a table sweep). The wrapper's sensible
// energy offset is calibrated once per pure component so the stored h/e match
// the solver's conservation convention.
#include "CoolPropThermoModel.H"
#include "CoolPropEquilibrium.H"
#include "property_table.H"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using CPState = CoolProp::AbstractState;


// Equilibrium sound speed squared from the thermodynamic identity
//   a^2 = (dp/drho)_T + (T/rho^2) (dp/dT)_rho^2 / Cv
// using equilibrium PT samples; the stencil shrinks until the samples share the
// node's phase. Returns -1 when no reliable stencil is found.
static double identityA2(CoolProp::AbstractState& s, double p, double t, bool twoPhase)
{
    for (double eps : {1e-4, 5e-5, 2.5e-5, 1.25e-5, 6.25e-6})
    {
        try
        {
            const double dp = p * eps, dt = t * eps;
            s.update(CoolProp::PT_INPUTS, p + dp, t);
            const double rp = s.rhomass(), ep = s.umass(); const bool tpp = s.phase() == CoolProp::iphase_twophase;
            s.update(CoolProp::PT_INPUTS, p - dp, t);
            const double rm = s.rhomass(), em = s.umass(); const bool tpm = s.phase() == CoolProp::iphase_twophase;
            s.update(CoolProp::PT_INPUTS, p, t + dt);
            const double rt = s.rhomass(), et = s.umass(); const bool tpt = s.phase() == CoolProp::iphase_twophase;
            s.update(CoolProp::PT_INPUTS, p, t - dt);
            const double rb = s.rhomass(), eb = s.umass(); const bool tpb = s.phase() == CoolProp::iphase_twophase;
            if (twoPhase && !(tpp && tpm && tpt && tpb)) continue;
            const double drho_dp_T = (rp - rm) / (2 * dp);
            const double drho_dT_p = (rt - rb) / (2 * dt);
            const double de_dT_p = (et - eb) / (2 * dt);
            const double de_dp_T = (ep - em) / (2 * dp);
            if (!(drho_dp_T > 0.)) continue;
            const double dp_drho_T = 1. / drho_dp_T;
            const double dp_dT_rho = -drho_dT_p / drho_dp_T;
            const double cv = de_dT_p - de_dp_T * dp_dT_rho;
            if (!std::isfinite(cv) || cv <= 0.) continue;
            const double rho = 0.5 * (rp + rm);
            const double a2 = dp_drho_T + (t / (rho * rho)) * dp_dT_rho * dp_dT_rho / cv;
            if (std::isfinite(a2) && a2 > 0.) return a2;
        }
        catch (const std::exception&) {}
    }
    return -1.;
}

int main(int argc, char** argv)
{
    std::string out = "work/a406_rootcause/property_table.bin";
    double tMin = 80., tMax = 550., pMin = 0.05e6, pMax = 6.5e6;
    int nT = 189, nP = 66;
    bool soundSpeed = false;
    std::vector<double> yNodes = {0., 0.1, 0.25, 0.4, 0.5, 0.6, 0.73, 0.8,
                                  0.9, 0.9162460633, 0.95, 0.99, 1.};
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out = argv[++i];
        else if (a == "--nT" && i + 1 < argc) nT = std::atoi(argv[++i]);
        else if (a == "--nP" && i + 1 < argc) nP = std::atoi(argv[++i]);
        else if (a == "--tmin" && i + 1 < argc) tMin = std::atof(argv[++i]);
        else if (a == "--tmax" && i + 1 < argc) tMax = std::atof(argv[++i]);
        else if (a == "--pmin" && i + 1 < argc) pMin = std::atof(argv[++i]);
        else if (a == "--pmax" && i + 1 < argc) pMax = std::atof(argv[++i]);
        else if (a == "--sound-speed") soundSpeed = true;
        else if (a == "--ylist" && i + 1 < argc)
        {
            yNodes.clear();
            std::string list = argv[++i];
            std::size_t start = 0;
            while (start <= list.size())
            {
                const std::size_t comma = list.find(',', start);
                const std::string item = list.substr(
                    start, comma == std::string::npos ? std::string::npos : comma - start);
                if (!item.empty()) yNodes.push_back(std::atof(item.c_str()));
                if (comma == std::string::npos) break;
                start = comma + 1;
            }
            if (yNodes.size() < 2)
            {
                std::fprintf(stderr, "--ylist needs at least two values\n");
                return 2;
            }
        }
        else
        {
            std::fprintf(stderr, "unknown option %s\n", a.c_str());
            return 2;
        }
    }
    rfamtable::Table table;
    table.Y = yNodes;
    table.T.resize(nT);
    table.p.resize(nP);
    for (int i = 0; i < nT; ++i) table.T[i] = tMin + (tMax - tMin) * i / (nT - 1);
    for (int i = 0; i < nP; ++i) table.p[i] = pMin + (pMax - pMin) * i / (nP - 1);
    table.nodes.assign(table.Y.size() * table.T.size() * table.p.size(), rfamtable::Node{});

    auto model = Foam::thermo::makeCoolPropThermoModel();
    std::printf("TABLE_BUILD model=%s\n", model->modelId().c_str());
    // Calibrate the wrapper's sensible-energy offset per pure component:
    // e_wrapper = u_HEOS - href. Use a single reference state per component.
    double href[2] = {0., 0.};
    {
        const char* names[2] = {"Oxygen", "Propane"};
        for (int comp = 0; comp < 2; ++comp)
        {
            const double pRef = 4.0e6, tRef = 300.;
            std::vector<double> y(2, 0.); y[comp] = 1.;
            const auto r = model->initializePTY(pRef, tRef, y);
            if (!r.success()) throw std::runtime_error("offset calibration failed");
            CPState* s = CPState::factory("HEOS", names[comp]);
            s->update(CoolProp::PT_INPUTS, pRef, tRef);
            href[comp] = s->umass() - r.state->e;
            delete s;
        }
    }
    std::printf("TABLE_OFFSET href_O2=%.10g href_C3H8=%.10g\n", href[0], href[1]);
    const auto start = std::chrono::steady_clock::now();
    long failures = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:failures)
    for (int iy = 0; iy < static_cast<int>(table.Y.size()); ++iy)
    {
        const std::vector<double> y{table.Y[iy], 1. - table.Y[iy]};
        const double offset = y[0] * href[0] + y[1] * href[1];
        CPState* s = (y[0] == 1.) ? CPState::factory("HEOS", "Oxygen")
                    : (y[1] == 1.) ? CPState::factory("HEOS", "Propane")
                                   : CPState::factory("HEOS", "Oxygen&Propane");
        if (y[0] != 1. && y[1] != 1.) s->set_mass_fractions(y);
        for (int it = 0; it < nT; ++it)
            for (int ip = 0; ip < nP; ++ip)
            {
                rfamtable::Node& n = table.nodes[table.index(iy, it, ip)];
                try
                {
                    s->unspecify_phase();
                    s->update(CoolProp::PT_INPUTS, table.p[ip], table.T[it]);
                    const double rho = s->rhomass(), u = s->umass();
                    if (!std::isfinite(rho) || rho <= 0. || !std::isfinite(u)
                        || !std::isfinite(s->smass()))
                    {
                        ++failures;
                        continue;
                    }
                    n.rho = rho; n.e = u - offset; n.h = s->hmass() - offset;
                    n.s = s->smass();
                    // CoolProp's mixture phase() can report iphase_not_imposed
                    // for an equilibrium two-phase point; the quality is the
                    // reliable indicator (Q>=0 => two-phase).
                    n.phase = (s->Q() >= 0.0) ? 5 : static_cast<std::int32_t>(s->phase());
                    // Reject the degenerate trivial split (x=y, rhoL=rhoV).
                    if (s->phase() == CoolProp::iphase_twophase)
                    {
                        const double rhoL = s->saturated_liquid_keyed_output(CoolProp::iDmass);
                        const double rhoV = s->saturated_vapor_keyed_output(CoolProp::iDmass);
                        if (!std::isfinite(rhoL) || !std::isfinite(rhoV)
                            || std::fabs(rhoL - rhoV) < 1e-9 * std::max(rhoL, rhoV))
                            n.phase = -1;
                    }
                    if (soundSpeed && n.phase >= 0)
                    {
                        const double a2 = identityA2(*s, table.p[ip], table.T[it],
                                                     s->phase() == CoolProp::iphase_twophase);
                        n.a2 = a2;
                        if (a2 <= 0.) ++failures;
                    }
                }
                catch (const std::exception&) { ++failures; }
            }
        delete s;
    }
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    table.write(out);
    std::ofstream meta(out + ".json");
    meta << "{\n  \"table\": \"" << out << "\",\n"
         << "  \"model\": \"" << model->modelId() << "\",\n"
         << "  \"yNodes\": " << table.Y.size() << ",\n"
         << "  \"tNodes\": " << nT << ",\n  \"pNodes\": " << nP << ",\n"
         << "  \"tRange_K\": [" << tMin << ", " << tMax << "],\n"
         << "  \"pRange_Pa\": [" << pMin << ", " << pMax << "],\n"
         << "  \"nodes\": " << table.nodes.size() << ",\n"
         << "  \"flashFailures\": " << failures << ",\n"
         << "  \"wallSeconds\": " << seconds << "\n}\n";
    std::printf("TABLE_BUILT nodes=%zu failures=%ld wall=%.1fs file=%s\n",
                table.nodes.size(), failures, seconds, out.c_str());
    return failures == 0 ? 0 : 1;
}
