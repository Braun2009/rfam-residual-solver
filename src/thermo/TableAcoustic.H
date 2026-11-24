#ifndef RFAM_TABLE_ACOUSTIC_H
#define RFAM_TABLE_ACOUSTIC_H
// Optional table-backed equilibrium sound speed for the tetra HEM path.
// Enabled explicitly with the environment variable RFAM_TETRA_TABLE=<path> to
// an RFAMTAB2 table that carries per-node a2 (built with --sound-speed and
// verified by tools/verify_property_table.cpp). When disabled or out of range
// the caller falls back to the converged EOS isentrope.
#include "PropertyTableFormat.H"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <string>

namespace Foam
{
namespace thermo
{
class TableAcoustic
{
    rfamtable::Table table_;
    std::string path_;
public:
    explicit TableAcoustic(const std::string& path)
      : table_(rfamtable::Table::read(path)), path_(path) {}
    const std::string& path() const { return path_; }
    std::size_t nodeCount() const { return table_.nodes.size(); }
    // True when the 7x7x7 block of table cells around (yO2,T,p) is entirely
    // single phase (and valid). The +-3 node margin keeps the fast path away
    // from the thin near-boundary two-phase band that a 2 K table cannot resolve.
    bool singlePhaseBlock(double yO2, double T, double p) const
    {
        const auto& t = table_;
        if (t.Y.size() < 3 || t.T.size() < 3 || t.p.size() < 3) return false;
        const std::size_t iy = t.bracket(t.Y, yO2), it = t.bracket(t.T, T), ip = t.bracket(t.p, p);
        for (int a = -3; a <= 3; ++a)
        {
            const long long jy = static_cast<long long>(iy) + a;
            if (jy < 0 || jy + 1 >= static_cast<long long>(t.Y.size())) {
                if (a != 0) continue;
                return false;
            }
            for (int b = -3; b <= 3; ++b)
            {
                const long long jt = static_cast<long long>(it) + b;
                if (jt < 0 || jt + 1 >= static_cast<long long>(t.T.size())) {
                    if (b != 0) continue;
                    return false;
                }
                for (int c = -3; c <= 3; ++c)
                {
                    const long long jp = static_cast<long long>(ip) + c;
                    if (jp < 0 || jp + 1 >= static_cast<long long>(t.p.size())) {
                        if (c != 0) continue;
                        return false;
                    }
                    for (int da = 0; da < 2; ++da)
                        for (int db = 0; db < 2; ++db)
                            for (int dc = 0; dc < 2; ++dc)
                            {
                                const rfamtable::Node& n = t.nodes[t.index(
                                    static_cast<std::size_t>(jy + da),
                                    static_cast<std::size_t>(jt + db),
                                    static_cast<std::size_t>(jp + dc))];
                                if (n.phase < 0 || n.phase == 5) return false;  // missing or two-phase
                            }
                }
            }
        }
        return true;
    }
    // Returns the interpolated a^2 when the point is inside the table and all
    // corners carry a computed a2.
    bool lookup(double yO2, double T, double p, double& a2) const
    {
        if (yO2 < table_.Y.front() || yO2 > table_.Y.back()) return false;
        if (T < table_.T.front() || T > table_.T.back()) return false;
        if (p < table_.p.front() || p > table_.p.back()) return false;
        rfamtable::Node n;
        if (!table_.interpolate(yO2, T, p, n)) return false;
        if (!(n.a2 > 0.)) return false;
        a2 = n.a2;
        return true;
    }
    // Equilibrium density from (p, h, Y) -- used by the interface-consistent
    // flux. The flux needs only a *consistent* density for the interface state
    // (the cell states themselves stay exact through the recovery), so Y and p
    // snap to the nearest table node and h is inverted along the temperature
    // grid (h is monotone in T at fixed Y,p); rho is interpolated linearly.
    // Returns false when the point is outside the table or the h bracket --
    // the caller then keeps the unmodified AUSM density (counted separately).
    bool stateFromPHY(double yO2, double p, double h, double& rho) const
    {
        const auto& Y = table_.Y;
        const auto& T = table_.T;
        const auto& P = table_.p;
        if (Y.empty() || T.size() < 2 || P.size() < 2) return false;
        const auto nearest = [](const std::vector<double>& g, double v, std::size_t& idx)
        {
            const double lo = g.front(), hi = g.back();
            const double dlo = g.size() > 1 ? g[1] - g[0] : 1.0;
            const double dhi = g.size() > 1 ? g.back() - g[g.size() - 2] : 1.0;
            if (v < lo - 0.5 * dlo || v > hi + 0.5 * dhi) return false;
            const auto it = std::lower_bound(g.begin(), g.end(), v);
            std::size_t i = static_cast<std::size_t>(it - g.begin());
            if (i >= g.size()) i = g.size() - 1;
            if (i > 0 && std::abs(g[i - 1] - v) < std::abs(g[i] - v)) --i;
            idx = i;
            return true;
        };
        std::size_t iy = 0, ip = 0;
        if (!nearest(Y, yO2, iy) || !nearest(P, p, ip)) return false;
        const auto hAt = [&](std::size_t it)
        { return table_.nodes[table_.index(iy, it, ip)].h; };
        if (!(hAt(0) <= h && h <= hAt(T.size() - 1))) return false;
        std::size_t lo = 0, hi = T.size() - 1;
        while (hi - lo > 1)
        {
            const std::size_t mid = (lo + hi) / 2;
            if (hAt(mid) <= h) lo = mid; else hi = mid;
        }
        const double h0 = hAt(lo), h1 = hAt(hi);
        const double w = h1 > h0 ? (h - h0) / (h1 - h0) : 0.0;
        const auto& n0 = table_.nodes[table_.index(iy, lo, ip)];
        const auto& n1 = table_.nodes[table_.index(iy, hi, ip)];
        rho = n0.rho + w * (n1.rho - n0.rho);
        return rho > 0. && std::isfinite(rho);
    }
    // Acoustic fallback variant: the last-resort sound speed must be available
    // even just outside the tabulated box (A-route cells 698987/785575 sit at
    // p=1.578 MPa where the local derivative routes all fail). The pressure
    // sensitivity of a^2 is weak compared with temperature, so clamping p into
    // the box is a documented approximation of this fallback branch only (the
    // fast path keeps strict out-of-range semantics).
    bool lookupClamped(double yO2, double T, double p, double& a2, bool& clamped) const
    {
        clamped=false;
        const double yc=std::clamp(yO2,table_.Y.front(),table_.Y.back());
        const double tc=std::clamp(T,table_.T.front(),table_.T.back());
        const double pc=std::clamp(p,table_.p.front(),table_.p.back());
        clamped=(yc!=yO2||tc!=T||pc!=p);
        if(lookup(yc,tc,pc,a2))return true;
        // Some nodes carry no a2 (the table build fails ~1 % of them, mostly in
        // the near-critical band). Walk outward node by node and take the first
        // valid interpolated value: the field-state check showed a2 varies
        // slowly (0.11 % mean, 0.49 % worst over the operating box), so a
        // nearest-valid-node value is acceptable for this fallback branch.
        const double dY=table_.Y.size()>1?table_.Y[1]-table_.Y[0]:0.;
        const double dT=table_.T.size()>1?table_.T[1]-table_.T[0]:0.;
        const double dp=table_.p.size()>1?table_.p[1]-table_.p[0]:0.;
        for(int k=1;k<=6;++k)
            for(int iy=-k;iy<=k;++iy)
                for(int it=-k;it<=k;++it)
                    for(int ip=-k;ip<=k;++ip)
                    {
                        if(std::max({std::abs(iy),std::abs(it),std::abs(ip)})!=k)continue;
                        const double yt=std::clamp(yc+iy*dY,table_.Y.front(),table_.Y.back());
                        const double tt=std::clamp(tc+it*dT,table_.T.front(),table_.T.back());
                        const double pt=std::clamp(pc+ip*dp,table_.p.front(),table_.p.back());
                        if(lookup(yt,tt,pt,a2)){clamped=true;return true;}
                    }
        return false;
    }
};
// Function-local singleton loaded once from RFAM_TETRA_TABLE. Returns nullptr
// when the variable is unset or the table cannot be read.
inline const TableAcoustic* tableAcoustic()
{
    static std::unique_ptr<TableAcoustic> instance;
    static std::once_flag once;
    std::call_once(once, []()
    {
        const char* path = std::getenv("RFAM_TETRA_TABLE");
        if (!path || !*path) return;
        try { instance = std::make_unique<TableAcoustic>(path); }
        catch (...) { instance.reset(); }
    });
    return instance.get();
}
}
}
#endif
