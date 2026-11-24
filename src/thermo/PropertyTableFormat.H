#ifndef RFAM_PROPERTY_TABLE_FORMAT_H
#define RFAM_PROPERTY_TABLE_FORMAT_H
// Shared binary property-table format (RFAMTAB2) for the O2/C3H8 equilibrium
// states. Axes: Y_O2, T (ascending), p (ascending). Node payload stores the
// RFAM sensible convention (same as CoolPropThermoModel): rho, h, e, s, the
// equilibrium sound speed squared (a2, -1 when not computed) and the CoolProp
// phase code. The table is produced by tools/build_property_table.cpp and
// verified by tools/verify_property_table.cpp; the tetra solver can consume it
// as an explicitly enabled acoustic path (src/thermo/TableAcoustic.H).
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rfamtable
{
struct Node
{
    double rho = 0., h = 0., e = 0., s = 0., a2 = -1.;
    std::int32_t phase = -1;
};
struct Table
{
    std::vector<double> Y, T, p;
    std::vector<Node> nodes; // [iy][it][ip]
    std::size_t index(std::size_t iy, std::size_t it, std::size_t ip) const
    { return (iy * T.size() + it) * p.size() + ip; }
    void write(const std::string& path) const
    {
        std::ofstream out(path, std::ios::binary);
        if (!out) throw std::runtime_error("cannot open table for write: " + path);
        out.write("RFAMTAB2", 8);
        const std::uint64_t ny = Y.size(), nt = T.size(), np = p.size();
        out.write(reinterpret_cast<const char*>(&ny), 8);
        out.write(reinterpret_cast<const char*>(&nt), 8);
        out.write(reinterpret_cast<const char*>(&np), 8);
        out.write(reinterpret_cast<const char*>(Y.data()), 8 * ny);
        out.write(reinterpret_cast<const char*>(T.data()), 8 * nt);
        out.write(reinterpret_cast<const char*>(p.data()), 8 * np);
        out.write(reinterpret_cast<const char*>(nodes.data()), sizeof(Node) * nodes.size());
    }
    static Table read(const std::string& path)
    {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("cannot open table: " + path);
        char magic[8];
        in.read(magic, 8);
        if (std::memcmp(magic, "RFAMTAB2", 8) != 0)
            throw std::runtime_error("bad table magic (need RFAMTAB2)");
        std::uint64_t ny = 0, nt = 0, np = 0;
        in.read(reinterpret_cast<char*>(&ny), 8);
        in.read(reinterpret_cast<char*>(&nt), 8);
        in.read(reinterpret_cast<char*>(&np), 8);
        Table t;
        t.Y.resize(ny); t.T.resize(nt); t.p.resize(np);
        in.read(reinterpret_cast<char*>(t.Y.data()), 8 * ny);
        in.read(reinterpret_cast<char*>(t.T.data()), 8 * nt);
        in.read(reinterpret_cast<char*>(t.p.data()), 8 * np);
        t.nodes.resize(ny * nt * np);
        in.read(reinterpret_cast<char*>(t.nodes.data()), sizeof(Node) * t.nodes.size());
        return t;
    }
    std::size_t bracket(const std::vector<double>& axis, double value) const
    {
        if (value <= axis.front()) return 0;
        if (value >= axis.back()) return axis.size() - 2;
        std::size_t lo = 0, hi = axis.size() - 1;
        while (hi - lo > 1)
        {
            const std::size_t mid = (lo + hi) / 2;
            (axis[mid] <= value ? lo : hi) = mid;
        }
        return lo;
    }
    // Trilinear interpolation; returns false if any corner is missing.
    bool interpolate(double yO2, double Tv, double pv, Node& out) const
    {
        if (Y.size() < 2 || T.size() < 2 || p.size() < 2) return false;
        const std::size_t iy = bracket(Y, yO2), it = bracket(T, Tv), ip = bracket(p, pv);
        const double fy = (yO2 - Y[iy]) / (Y[iy + 1] - Y[iy]);
        const double ft = (Tv - T[it]) / (T[it + 1] - T[it]);
        const double fp = (pv - p[ip]) / (p[ip + 1] - p[ip]);
        double value[6] = {0, 0, 0, 0, 0, 0};
        bool a2Valid = true;
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
                for (int c = 0; c < 2; ++c)
                {
                    const Node& n = nodes[index(iy + a, it + b, ip + c)];
                    if (!(n.rho > 0.) || n.phase < 0) return false;
                    if (!(n.a2 > 0.)) a2Valid = false;
                    const double w = (a ? fy : 1 - fy) * (b ? ft : 1 - ft) * (c ? fp : 1 - fp);
                    value[0] += w * n.rho; value[1] += w * n.h;
                    value[2] += w * n.e; value[3] += w * n.s; value[4] += w * n.a2;
                }
        out.rho = value[0]; out.h = value[1]; out.e = value[2]; out.s = value[3];
        out.a2 = a2Valid ? value[4] : -1.;
        out.phase = nodes[index(iy + (fy > .5), it + (ft > .5), ip + (fp > .5))].phase;
        return true;
    }
};
}
#endif
