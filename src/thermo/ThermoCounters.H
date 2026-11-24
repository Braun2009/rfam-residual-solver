#ifndef RFAM_THERMO_COUNTERS_H
#define RFAM_THERMO_COUNTERS_H
// Lightweight opt-in counters for the equilibrium EOS call paths. Enabled with
// RFAM_THERMO_COUNTERS=1. Increments are relaxed atomics: negligible next to an
// EOS flash, and zero work when disabled.
#include <atomic>
#include <cmath>
#include <cstdlib>
namespace Foam
{
namespace thermo
{
struct ThermoCounters
{
    std::atomic<unsigned long long> pt{0};
    std::atomic<unsigned long long> isentrope{0};
    std::atomic<unsigned long long> identity{0};
    std::atomic<unsigned long long> recover{0};
    std::atomic<unsigned long long> recoverNewtonFail{0};
    std::atomic<unsigned long long> recoverScanRestart{0};
    std::atomic<unsigned long long> twoPhasePublish{0};
    std::atomic<unsigned long long> tableAcousticHit{0};
    std::atomic<unsigned long long> speciesRepairEvents{0};
    std::atomic<unsigned long long> pureFallbackEvents{0};
    std::atomic<double> maxPureSnapFraction{0.};
    std::atomic<double> maxSpeciesRepairRel{0.};
    std::atomic<double> maxRelRho{0.};
    std::atomic<double> maxRelEnergy{0.};
};
inline ThermoCounters& thermoCounters()
{
    static ThermoCounters counters;
    return counters;
}
inline bool thermoCountersEnabled()
{
    static const bool enabled = std::getenv("RFAM_THERMO_COUNTERS") != nullptr;
    return enabled;
}
inline void thermoCount(std::atomic<unsigned long long>& counter)
{
    if (thermoCountersEnabled()) counter.fetch_add(1, std::memory_order_relaxed);
}
inline void thermoMax(std::atomic<double>& target, double value)
{
    if (!thermoCountersEnabled() || !std::isfinite(value)) return;
    double current = target.load(std::memory_order_relaxed);
    while (value > current
           && !target.compare_exchange_weak(current, value,
                                            std::memory_order_relaxed))
    {}
}
}
}
#endif
