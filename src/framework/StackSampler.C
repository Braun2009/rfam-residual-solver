#include "StackSampler.H"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

namespace Foam
{
namespace telemetry
{
namespace
{

constexpr int kMaxSites = 128;

// A sample names the instrumented loop that was executing. `label` is a string
// literal chosen by the owner of that loop, so pointer identity is a complete
// key: no hashing, no allocation, and no comparison beyond a short scan that
// only runs when a new label appears.
struct Site
{
    const char* label = nullptr;
    unsigned long long samples = 0;
    unsigned long long firstMs = 0;
    unsigned long long lastMs = 0;
};

Site g_sites[kMaxSites];
std::atomic<int> g_siteCount(0);
std::atomic<unsigned long long> g_totalSamples(0);
std::atomic<std::uint64_t> g_lastCaptureNs(0);
std::atomic<int> g_lastSite(-1);
std::atomic<unsigned long long> g_lastStampMs(0);

unsigned long long monotonicMs()
{
    return static_cast<unsigned long long>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count());
}

std::uint64_t monotonicNs()
{
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count());
}

int findOrRegister(const char* label)
{
    const char* key = label ? label : "unnamed";
    const int count = g_siteCount.load(std::memory_order_acquire);
    for (int i = 0; i < count; ++i)
        if (g_sites[i].label == key) return i;
    const int slot = g_siteCount.fetch_add(1, std::memory_order_acq_rel);
    if (slot >= kMaxSites) return -1;
    g_sites[slot].label = key;
    g_sites[slot].firstMs = monotonicMs();
    return slot;
}

struct Entry
{
    Site site;
};

} // End anonymous namespace

struct StackSampler::Impl
{
    std::string runDir;
    int hz = 5;
    std::uint64_t periodNs = 200000000ULL;
    std::atomic<bool> stopped{false};
};

void StackSampler::start(int hz, const std::string& runDir)
{
    if (running_) return;
    if (!impl_) impl_ = new Impl();
    impl_->runDir = runDir;
    impl_->hz = std::max(0, std::min(hz, 200));
    impl_->periodNs = impl_->hz > 0 ? static_cast<std::uint64_t>(1.0e9 / impl_->hz) : 0ULL;
    impl_->stopped.store(false);
    running_ = impl_->hz > 0;
}

void StackSampler::stop()
{
    if (!impl_) return;
    impl_->stopped.store(true);
    running_ = false;
}

void StackSampler::captureIfDue(const char* label) noexcept
{
    if (!running_ || !impl_ || impl_->periodNs == 0) return;
    if (impl_->stopped.load(std::memory_order_relaxed)) return;
    const std::uint64_t now = monotonicNs();
    const std::uint64_t last = g_lastCaptureNs.load(std::memory_order_relaxed);
    if (now - last < impl_->periodNs) return;
    // Every OpenMP worker can reach this at once; only one may record. Losing a
    // sample is cheaper than serialising the progress path.
    std::uint64_t expected = last;
    if (!g_lastCaptureNs.compare_exchange_strong(expected, now)) return;

    const int slot = findOrRegister(label);
    if (slot < 0) return;
    Site& site = g_sites[slot];
    site.samples += 1;
    site.lastMs = monotonicMs();
    g_lastSite.store(slot, std::memory_order_relaxed);
    g_lastStampMs.store(site.lastMs, std::memory_order_relaxed);
    g_totalSamples.fetch_add(1, std::memory_order_relaxed);
}

json::Value StackSampler::report() const
{
    json::Value root = json::Value::object();
    root.set("enabled", json::Value(running_));
    root.set("samples", json::Value(static_cast<long long>(g_totalSamples.load())));
    if (!running_ || !impl_) return root;
    root.set("hz", json::Value(impl_->hz));
    root.set("method", json::Value("call-site sampling (no signals, no unwinding)"));

    const int count = std::min(g_siteCount.load(std::memory_order_acquire), kMaxSites);
    std::vector<Entry> entries;
    for (int i = 0; i < count; ++i)
        if (g_sites[i].samples > 0) entries.push_back(Entry{g_sites[i]});
    std::sort(entries.begin(), entries.end(),
              [](const Entry& a, const Entry& b) { return a.site.samples > b.site.samples; });
    root.set("unique", json::Value(static_cast<long long>(entries.size())));

    const unsigned long long total = std::max<unsigned long long>(g_totalSamples.load(), 1ULL);
    json::Value top = json::Value::array();
    const std::size_t limit = std::min<std::size_t>(entries.size(), 8);
    for (std::size_t i = 0; i < limit; ++i)
    {
        json::Value item = json::Value::object();
        item.set("label", json::Value(std::string(entries[i].site.label)));
        item.set("n", json::Value(static_cast<long long>(entries[i].site.samples)));
        item.set("share", json::Value(static_cast<double>(entries[i].site.samples)
                                      / static_cast<double>(total)));
        top.push(item);
    }
    root.set("top", top);

    const int last = g_lastSite.load(std::memory_order_relaxed);
    if (last >= 0 && last < kMaxSites && g_sites[last].label)
    {
        json::Value current = json::Value::object();
        current.set("label", json::Value(std::string(g_sites[last].label)));
        current.set("age_ms", json::Value(static_cast<long long>(
            monotonicMs() - g_lastStampMs.load(std::memory_order_relaxed))));
        root.set("current", current);
    }
    return root;
}

std::string StackSampler::histogramText() const
{
    std::ostringstream out;
    out << "# RFAM call-site histogram\n";
    out << "# Each sample names the instrumented loop that was executing.\n";
    out << "# The call chain is reconstructed as: phase -> sub-phase -> site.\n";
    out << "# samples=" << g_totalSamples.load() << "\n";
    const int count = std::min(g_siteCount.load(std::memory_order_acquire), kMaxSites);
    std::vector<Entry> entries;
    for (int i = 0; i < count; ++i)
        if (g_sites[i].samples > 0) entries.push_back(Entry{g_sites[i]});
    std::sort(entries.begin(), entries.end(),
              [](const Entry& a, const Entry& b) { return a.site.samples > b.site.samples; });
    out << "# unique_sites=" << entries.size() << "\n";
    const unsigned long long total = std::max<unsigned long long>(g_totalSamples.load(), 1ULL);
    for (const Entry& entry : entries)
    {
        out << "count=" << entry.site.samples
            << " share=" << (100.0 * static_cast<double>(entry.site.samples) / static_cast<double>(total))
            << "%  " << entry.site.label << "\n";
    }
    return out.str();
}

} // End namespace telemetry
} // End namespace Foam
