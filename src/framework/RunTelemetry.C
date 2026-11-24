#include "RunTelemetry.H"
#include "StackSampler.H"
#include "ThermoDataProvenance.H"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace Foam
{
namespace telemetry
{

// ===========================================================================
// Static tables
// ===========================================================================

const char* severityName(Severity s)
{
    switch (s)
    {
        case Severity::Requested:     return "requested";
        case Severity::Approximation: return "approximation";
        case Severity::Downgrade:     return "downgrade";
        case Severity::FailClosed:    return "fail_closed";
    }
    return "unknown";
}

namespace
{

const char* const kMetricNames[] = {
    "residual", "residual_ratio", "dt_pseudo", "dt_physical", "physical_time",
    "chamber_pressure", "mach_max",
    "mdot_lox", "mdot_fuel", "mdot_torch", "mdot_out",
    "mass_kg", "mass_error", "cfl_limit",
    "temperature_min", "temperature_max", "pressure_min", "pressure_max",
};
static_assert(sizeof(kMetricNames) / sizeof(kMetricNames[0]) ==
              static_cast<std::size_t>(Metric::Count_),
              "metric table and enum are out of sync");

const char* const kCountNames[] = {
    "frozen", "limited", "stiffened", "two_phase_cells", "rejected", "retries",
    "thermo_pt", "thermo_isentrope", "thermo_recover", "thermo_recover_fail",
    "species_repair", "pure_fallback", "boundary_fallback", "interface_fallback",
    "outlet_shortcut_hits", "outlet_exact_faces",
    "chem_accepted", "chem_rejected", "chem_skipped_cold", "chem_skipped_two_phase",
};
static_assert(sizeof(kCountNames) / sizeof(kCountNames[0]) ==
              static_cast<std::size_t>(Count::Count_),
              "count table and enum are out of sync");

constexpr std::size_t kMetricCount = static_cast<std::size_t>(Metric::Count_);
constexpr std::size_t kCountCount = static_cast<std::size_t>(Count::Count_);

// ---------------------------------------------------------------------------
// Phase registry. Static storage, one mutex used only at registration time.
// ---------------------------------------------------------------------------

constexpr int kMaxPhases = 192;

struct PhaseSlot
{
    const char* name;
    const char* file;
    int line;
};

// Function-local statics: a driver's static-init phase registration must not
// depend on this translation unit's globals having been constructed first.
struct PhaseRegistry
{
    PhaseSlot slots[kMaxPhases];
    std::atomic<int> count{0};
    std::mutex mutex;
};

PhaseRegistry& registry()
{
    static PhaseRegistry instance;
    return instance;
}

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

std::string nowIso()
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    struct tm tm;
    gmtime_r(&ts.tv_sec, &tm);
    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S", &tm);
    char full[48];
    std::snprintf(full, sizeof(full), "%s.%03ldZ", buffer, ts.tv_nsec / 1000000L);
    return std::string(full);
}


std::uint64_t nowNs()
{
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count());
}

bool makeDirs(const std::string& path)
{
    if (path.empty()) return false;
    std::string partial;
    for (std::size_t i = 0; i < path.size(); ++i)
    {
        partial.push_back(path[i]);
        if (path[i] == '/' || i + 1 == path.size())
        {
            if (partial == "/") continue;
            if (::mkdir(partial.c_str(), 0775) != 0 && errno != EEXIST) return false;
        }
    }
    return true;
}

bool writeAtomic(const std::string& path, const std::string& bytes, bool doFsync)
{
    const std::string temp = path + ".tmp." + std::to_string(static_cast<long long>(::getpid()));
    const int fd = ::open(temp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
    if (fd < 0) return false;
    const char* data = bytes.data();
    std::size_t left = bytes.size();
    while (left > 0)
    {
        const ssize_t n = ::write(fd, data, left);
        if (n <= 0) { ::close(fd); ::unlink(temp.c_str()); return false; }
        data += n;
        left -= static_cast<std::size_t>(n);
    }
    if (doFsync) ::fdatasync(fd);
    ::close(fd);
    if (::rename(temp.c_str(), path.c_str()) != 0)
    {
        ::unlink(temp.c_str());
        return false;
    }
    return true;
}

void appendFile(const std::string& path, const std::string& bytes)
{
    const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0664);
    if (fd < 0) return;
    const char* data = bytes.data();
    std::size_t left = bytes.size();
    while (left > 0)
    {
        const ssize_t n = ::write(fd, data, left);
        if (n <= 0) break;
        data += n;
        left -= static_cast<std::size_t>(n);
    }
    ::close(fd);
}

std::string readFile(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) return std::string();
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::string exePath()
{
    char buffer[4096];
    const ssize_t n = ::readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (n <= 0) return std::string();
    buffer[n] = 0;
    return std::string(buffer);
}

std::string hostName()
{
    char buffer[256];
    if (::gethostname(buffer, sizeof(buffer)) != 0) return "unknown";
    buffer[sizeof(buffer) - 1] = 0;
    return std::string(buffer);
}

std::string userName()
{
    const char* user = std::getenv("USER");
    if (!user || !*user) user = std::getenv("LOGNAME");
    return user ? std::string(user) : std::string("unknown");
}

std::string cwdPath()
{
    char buffer[4096];
    if (!::getcwd(buffer, sizeof(buffer))) return "?";
    return std::string(buffer);
}

// ---------------------------------------------------------------------------
// /proc sampling. Failures produce NaN, never a plausible-looking guess.
// ---------------------------------------------------------------------------

struct ProcSample
{
    double cpuPercent = std::nan("");
    double rssMb = std::nan("");
    long threads = -1;
    char state = '?';
    unsigned long long readBytes = 0;
    unsigned long long writeBytes = 0;
    unsigned long long cpuTicks = 0;
};

struct CpuTicks
{
    unsigned long long total = 0;
    unsigned long long idle = 0;
};

CpuTicks readCpuTicks()
{
    CpuTicks ticks;
    std::ifstream input("/proc/stat");
    std::string line;
    if (!std::getline(input, line)) return ticks;
    if (line.compare(0, 4, "cpu ") != 0) return ticks;
    std::istringstream stream(line.substr(4));
    unsigned long long values[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 10; ++i) if (!(stream >> values[i])) break;
    for (int i = 0; i < 10; ++i) ticks.total += values[i];
    ticks.idle = values[3] + values[4];
    return ticks;
}

ProcSample readProcSample(int pid)
{
    ProcSample sample;
    const std::string stat = readFile("/proc/" + std::to_string(pid) + "/stat");
    const std::size_t close = stat.rfind(')');
    if (close == std::string::npos) return sample;
    std::istringstream stream(stat.substr(close + 1));
    std::vector<std::string> fields;
    std::string token;
    while (stream >> token) fields.push_back(token);
    // fields[0] == field 3 of proc(5): pid and comm were consumed above.
    if (fields.size() > 21)
    {
        sample.state = fields[0].empty() ? '?' : fields[0][0];
        sample.cpuTicks = std::strtoull(fields[11].c_str(), nullptr, 10)
                        + std::strtoull(fields[12].c_str(), nullptr, 10);
        sample.threads = std::atol(fields[17].c_str());
        const long pages = std::atol(fields[21].c_str());
        sample.rssMb = static_cast<double>(pages) * static_cast<double>(::sysconf(_SC_PAGESIZE))
                     / (1024.0 * 1024.0);
    }
    std::ifstream io("/proc/" + std::to_string(pid) + "/io");
    std::string line;
    while (std::getline(io, line))
    {
        if (line.compare(0, 6, "rchar:") == 0) sample.readBytes = std::strtoull(line.c_str() + 6, nullptr, 10);
        else if (line.compare(0, 7, "wchar:") == 0) sample.writeBytes = std::strtoull(line.c_str() + 7, nullptr, 10);
    }
    return sample;
}

struct MemorySample
{
    double totalGb = std::nan("");
    double usedGb = std::nan("");
};

MemorySample readMemory()
{
    MemorySample sample;
    std::ifstream input("/proc/meminfo");
    std::string line;
    double totalKb = 0.0, availableKb = 0.0;
    bool haveTotal = false, haveAvailable = false;
    while (std::getline(input, line))
    {
        if (line.compare(0, 9, "MemTotal:") == 0)
        {
            totalKb = std::strtod(line.c_str() + 9, nullptr);
            haveTotal = true;
        }
        else if (line.compare(0, 13, "MemAvailable:") == 0)
        {
            availableKb = std::strtod(line.c_str() + 13, nullptr);
            haveAvailable = true;
        }
    }
    if (haveTotal)
    {
        sample.totalGb = totalKb / (1024.0 * 1024.0);
        if (haveAvailable) sample.usedGb = (totalKb - availableKb) / (1024.0 * 1024.0);
    }
    return sample;
}

double readLoad1()
{
    std::ifstream input("/proc/loadavg");
    double load1 = std::nan("");
    if (!(input >> load1)) return std::nan("");
    return load1;
}

// ---------------------------------------------------------------------------
// Live snapshot. Only atomics: the heartbeat reads it while the solver writes,
// and neither side may ever block the other.
// ---------------------------------------------------------------------------

struct Snapshot
{
    std::atomic<int> phaseId{-1};
    std::atomic<int> phaseActive{0};
    std::atomic<std::uint64_t> phaseStartNs{0};
    std::atomic<std::uint64_t> phaseEndNs{0};
    std::atomic<const char*> subName{nullptr};
    std::atomic<const char*> subFile{nullptr};
    std::atomic<int> subLine{0};
    std::atomic<std::uint64_t> subStartNs{0};

    std::atomic<std::uint64_t> progDone{0};
    std::atomic<std::uint64_t> progTotal{0};
    std::atomic<const char*> progUnit{"?"};
    std::atomic<std::uint64_t> progStartNs{0};
    std::atomic<std::uint64_t> lastProgressNs{0};

    std::atomic<long long> stepIndex{-1};
    std::atomic<long long> stepTotal{-1};
    std::atomic<long long> stepAccepted{0};

    std::atomic<double> metrics[kMetricCount];
    std::atomic<unsigned long long> counts[kCountCount];

    Snapshot()
    {
        for (std::size_t i = 0; i < kMetricCount; ++i)
            metrics[i].store(std::nan(""), std::memory_order_relaxed);
        for (std::size_t i = 0; i < kCountCount; ++i)
            counts[i].store(0ULL, std::memory_order_relaxed);
    }
};

// ---------------------------------------------------------------------------
// Degradation ledger. Counts are atomics; the detail string of a slot is
// immutable once `published` is set with release ordering, which is what lets
// the heartbeat read it without a lock.
// ---------------------------------------------------------------------------

constexpr int kMaxLedger = 64;

struct LedgerSlot
{
    std::atomic<int> published{0};
    std::atomic<unsigned long long> occurrences{0};
    Severity severity = Severity::Approximation;
    std::string code;
    std::string detail;
};

struct Ledger
{
    std::mutex mutex;
    std::atomic<int> used{0};
    LedgerSlot slots[kMaxLedger];

    // Returns true the first time a code is seen. Callers use that to decide
    // whether to emit an event line: a downgrade that repeats every monitor
    // step should raise its count, not flood the event stream.
    bool note(const char* code, Severity severity, const std::string& detail)
    {
        std::lock_guard<std::mutex> lock(mutex);
        const int count = used.load(std::memory_order_relaxed);
        for (int i = 0; i < count; ++i)
        {
            if (slots[i].code == code)
            {
                slots[i].occurrences.fetch_add(1, std::memory_order_relaxed);
                if (static_cast<int>(severity) > static_cast<int>(slots[i].severity))
                    slots[i].severity = severity;
                return false;
            }
        }
        if (count >= kMaxLedger) return false;
        slots[count].code = code;
        slots[count].detail = detail;
        slots[count].severity = severity;
        slots[count].occurrences.store(1ULL, std::memory_order_relaxed);
        used.store(count + 1, std::memory_order_relaxed);
        slots[count].published.store(1, std::memory_order_release);
        return true;
    }

    // Exact-level counts: a fail-closed event is not also counted as a
    // downgrade in the summary, so the three numbers a reader sees always add
    // up to the number of ledger entries.
    unsigned long long countExactly(Severity level) const
    {
        unsigned long long total = 0;
        const int count = used.load(std::memory_order_acquire);
        for (int i = 0; i < count; ++i)
        {
            if (slots[i].published.load(std::memory_order_acquire) == 0) continue;
            if (slots[i].severity != level) continue;
            total += slots[i].occurrences.load(std::memory_order_relaxed);
        }
        return total;
    }

    Severity worst() const
    {
        Severity level = Severity::Requested;
        const int count = used.load(std::memory_order_acquire);
        for (int i = 0; i < count; ++i)
        {
            if (slots[i].published.load(std::memory_order_acquire) == 0) continue;
            if (static_cast<int>(slots[i].severity) > static_cast<int>(level))
                level = slots[i].severity;
        }
        return level;
    }

    json::Value dump() const
    {
        json::Value items = json::Value::array();
        const int count = used.load(std::memory_order_acquire);
        for (int i = 0; i < count; ++i)
        {
            if (slots[i].published.load(std::memory_order_acquire) == 0) continue;
            json::Value item = json::Value::object();
            item.set("code", json::Value(slots[i].code));
            item.set("severity", json::Value(severityName(slots[i].severity)));
            item.set("n", json::Value(slots[i].occurrences.load(std::memory_order_relaxed)));
            item.set("detail", json::Value(slots[i].detail));
            items.push(item);
        }
        return items;
    }
};

} // End anonymous namespace

// ===========================================================================
// phase()
// ===========================================================================

Phase phase(const char* name, const char* file, int line)
{
    PhaseRegistry& table = registry();
    std::lock_guard<std::mutex> lock(table.mutex);
    const int count = table.count.load(std::memory_order_relaxed);
    for (int i = 0; i < count; ++i)
        if (std::strcmp(table.slots[i].name, name) == 0)
            return Phase{i, table.slots[i].name, table.slots[i].file, table.slots[i].line};
    if (count >= kMaxPhases) return Phase{-1, name, file, line};
    table.slots[count] = PhaseSlot{name, file, line};
    table.count.store(count + 1, std::memory_order_release);
    return Phase{count, name, file, line};
}

namespace
{

const PhaseSlot* phaseSlot(int id)
{
    const PhaseRegistry& table = registry();
    const int count = table.count.load(std::memory_order_acquire);
    if (id < 0 || id >= count) return nullptr;
    return &table.slots[id];
}

// ---------------------------------------------------------------------------
// Signal plumbing. The handler does one async-signal-safe thing: it writes a
// one-byte code to a pre-opened descriptor and sets an atomic. Everything else
// happens on the heartbeat thread, which is never blocked by solver work.
// ---------------------------------------------------------------------------

std::atomic<int> g_signalNumber(0);
std::atomic<int> g_signalFd(-1);
std::atomic<int> g_signalCount(0);

extern "C" void rfamTelemetrySignalHandler(int signo)
{
    const int fd = g_signalFd.load(std::memory_order_relaxed);
    if (fd >= 0)
    {
        unsigned char code = static_cast<unsigned char>(signo & 0xFF);
        ssize_t ignored = ::write(fd, &code, 1);
        (void)ignored;
    }
    g_signalNumber.store(signo, std::memory_order_relaxed);
    // Duplicate signals are normal (`timeout`, Ctrl-C twice, a service manager
    // that retries). Exiting on the second one would throw away the summary --
    // the record is the reason the handler exists, so only a fourth delivery
    // abandons it.
    if (g_signalCount.fetch_add(1, std::memory_order_relaxed) + 1 > 3)
        ::_exit(128 + signo);
}

} // End anonymous namespace

// ===========================================================================
// Session::Impl
// ===========================================================================

struct Session::Impl
{
    Options options;
    bool active = false;
    std::string runId;
    std::string runDir;
    std::string rootDir;

    json::Value manifest = json::Value::object();
    json::Value extras = json::Value::object();

    Snapshot snap;
    Ledger ledger;

    // argv-derived state
    std::vector<std::string> argv;
    std::vector<std::string> argFlags;   // "--flag" spellings that were present

    // manifest field sources
    std::mutex manifestMutex;

    // files
    std::string pathStatus, pathEvents, pathHeartbeat, pathSessionLog,
                pathManifest, pathSummary, pathStacks, pathSignal;

    std::mutex eventMutex;      // solver-thread events only
    std::mutex summaryMutex;
    std::atomic<bool> summaryWritten{false};
    std::atomic<bool> finishing{false};

    // heartbeat bookkeeping
    std::thread heartbeatThread;
    std::mutex heartbeatMutex;
    std::condition_variable heartbeatWake;
    std::atomic<bool> heartbeatStop{false};
    std::int64_t heartbeatPeriodNs = 1000000000LL;
    std::uint64_t tickCount = 0;
    std::uint64_t missedTicks = 0;
    double maxJitterMs = 0.0;
    double lastJitterMs = 0.0;
    double lastWriteMs = 0.0;
    std::uint64_t startSteadyNs = 0;

    // proc deltas
    unsigned long long lastCpuTicks = 0;
    std::uint64_t lastCpuSampleNs = 0;
    CpuTicks lastSystemCpu;
    bool haveSystemCpu = false;

    // stall bookkeeping
    bool inStall = false;
    std::uint64_t stallSinceNs = 0;
    double stallSeconds = 0.0;

    StackSampler sampler;

    std::string joinPath(const char* name) const { return runDir + "/" + name; }

    bool flagPresent(const std::string& flag) const
    {
        return std::find(argFlags.begin(), argFlags.end(), flag) != argFlags.end();
    }

    void emitEvent(const json::Value& fields, const char* type);
    void note(const char* type, const std::string& human, const json::Value& fields);
    void runHeartbeat();
    bool writeStatus();
    void writeSummary(const char* verdict, int exitCode, const json::Value& extra);
    void flushManifest();
    json::Value snapshotJson(bool includeSystem);
    void onTermination();
    void pruneOldLink();
};

// ---------------------------------------------------------------------------

void Session::Impl::emitEvent(const json::Value& fields, const char* type)
{
    json::Value record = json::Value::object();
    record.set("t", json::Value(nowIso()));
    record.set("type", json::Value(type));
    record.set("pid", json::Value(static_cast<long long>(::getpid())));
    for (const auto& member : fields.members()) record.set(member.first, member.second);

    std::lock_guard<std::mutex> lock(eventMutex);
    appendFile(pathEvents, record.dump(false) + "\n");
    if (options.mirrorSessionLog)
    {
        std::string line = nowIso();
        line += "  ";
        line += type;
        for (const auto& member : fields.members())
        {
            line += "  ";
            line += member.first;
            line += "=";
            const json::Value& value = member.second;
            if (value.isString()) line += value.asString();
            else line += value.dump(false);
        }
        line += "\n";
        appendFile(pathSessionLog, line);
    }
}

void Session::Impl::note(const char* type, const std::string& human, const json::Value& fields)
{
    (void)human;
    emitEvent(fields, type);
}
void Session::Impl::onTermination()
{
    const int signo = g_signalNumber.load(std::memory_order_relaxed);
    if (signo == 0) return;
    if (finishing.exchange(true)) return;
    const PhaseSlot* slot = phaseSlot(snap.phaseId.load(std::memory_order_relaxed));
    {
        json::Value record = json::Value::object();
        record.set("t", json::Value(nowIso()));
        record.set("type", json::Value("SIGNAL"));
        record.set("signal", json::Value(static_cast<long long>(signo)));
        record.set("phase", json::Value(slot ? slot->name : "unknown"));
        record.set("step", json::Value(snap.stepIndex.load(std::memory_order_relaxed)));
        appendFile(pathHeartbeat, record.dump(false) + "\n");
    }
    json::Value extra = json::Value::object();
    extra.set("signal", json::Value(static_cast<long long>(signo)));
    writeSummary("terminated", 128 + signo, extra);
    ::_exit(128 + signo);
}

void Session::Impl::runHeartbeat()
{
    // Exact-deadline scheduling: tick n happens at start + n*period. `sleep_for`
    // accumulates the cost of each write; on a 9p mount over a Windows drive
    // that drift is visible within a few minutes, and a drifting heartbeat is
    // indistinguishable from a stalling solver.
    const std::int64_t period = heartbeatPeriodNs;
    const std::int64_t startedNs = static_cast<std::int64_t>(startSteadyNs);

    // First snapshot immediately: a viewer that starts before the first tick
    // must never see a missing file and guess.
    {
        const std::uint64_t writeStart = nowNs();
        writeStatus();
        lastWriteMs = static_cast<double>(nowNs() - writeStart) / 1.0e6;
        ++tickCount;
    }

    std::unique_lock<std::mutex> lock(heartbeatMutex);
    for (std::uint64_t tick = 0; !heartbeatStop.load(std::memory_order_relaxed); ++tick)
    {
        const std::int64_t deadline = startedNs + static_cast<std::int64_t>(tick + 1) * period;
        const std::int64_t nowSigned = static_cast<std::int64_t>(nowNs());
        if (nowSigned + 200000LL < deadline)
        {
            heartbeatWake.wait_for(lock, std::chrono::nanoseconds(deadline - nowSigned));
            if (heartbeatStop.load(std::memory_order_relaxed)) break;
        }
        const std::int64_t actual = static_cast<std::int64_t>(nowNs());
        const double jitterMs = static_cast<double>(actual - deadline) / 1.0e6;
        lastJitterMs = jitterMs;
        if (jitterMs > maxJitterMs) maxJitterMs = jitterMs;
        if (jitterMs > static_cast<double>(period) / 1.0e6) ++missedTicks;
        ++tickCount;

        const std::uint64_t writeStart = nowNs();
        const bool wrote = writeStatus();
        lastWriteMs = static_cast<double>(nowNs() - writeStart) / 1.0e6;
        if (!wrote)
        {
            json::Value record = json::Value::object();
            record.set("t", json::Value(nowIso()));
            record.set("type", json::Value("HEARTBEAT_WRITE_FAILED"));
            record.set("path", json::Value(pathStatus));
            record.set("errno", json::Value(static_cast<long long>(errno)));
            appendFile(pathHeartbeat, record.dump(false) + "\n");
        }

        // Stall detection: no progress published for longer than the threshold.
        const std::uint64_t lastProgress = snap.lastProgressNs.load(std::memory_order_relaxed);
        const std::uint64_t nowValue = nowNs();
        const double idleSeconds = lastProgress ? static_cast<double>(nowValue - lastProgress) / 1.0e9 : 0.0;
        const bool stalledNow = idleSeconds > static_cast<double>(options.stallSeconds);
        if (stalledNow && !inStall)
        {
            inStall = true;
            stallSinceNs = nowValue;
            const PhaseSlot* slot = phaseSlot(snap.phaseId.load(std::memory_order_relaxed));
            json::Value record = json::Value::object();
            record.set("t", json::Value(nowIso()));
            record.set("type", json::Value("STALL"));
            record.set("phase", json::Value(slot ? slot->name : "unknown"));
            record.set("file", json::Value(slot ? slot->file : ""));
            record.set("line", json::Value(slot ? slot->line : 0));
            record.set("idle_s", json::Value(idleSeconds));
            record.set("step", json::Value(snap.stepIndex.load(std::memory_order_relaxed)));
            appendFile(pathHeartbeat, record.dump(false) + "\n");

        }
        else if (!stalledNow && inStall)
        {
            inStall = false;
            stallSeconds = static_cast<double>(nowValue - stallSinceNs) / 1.0e9;
            json::Value record = json::Value::object();
            record.set("t", json::Value(nowIso()));
            record.set("type", json::Value("STALL_CLEARED"));
            record.set("stalled_s", json::Value(stallSeconds));
            appendFile(pathHeartbeat, record.dump(false) + "\n");

        }

        if (g_signalNumber.load(std::memory_order_relaxed) != 0)
        {
            lock.unlock();
            onTermination();
            return;
        }
    }
}

// ---------------------------------------------------------------------------

json::Value Session::Impl::snapshotJson(bool includeSystem)
{
    json::Value root = json::Value::object();
    const std::uint64_t nowValue = nowNs();
    root.set("schema", json::Value("rfam.run.status/1"));
    root.set("run_id", json::Value(runId));
    root.set("pid", json::Value(static_cast<long long>(::getpid())));
    root.set("seq", json::Value(static_cast<long long>(tickCount)));
    root.set("wall", json::Value(nowIso()));
    root.set("uptime_s", json::Value(static_cast<double>(nowValue - startSteadyNs) / 1.0e9));

    json::Value heartbeat = json::Value::object();
    heartbeat.set("period_ms", json::Value(static_cast<long long>(heartbeatPeriodNs / 1000000LL)));
    heartbeat.set("ticks", json::Value(static_cast<long long>(tickCount)));
    heartbeat.set("missed", json::Value(static_cast<long long>(missedTicks)));
    heartbeat.set("jitter_ms", json::Value(lastJitterMs));
    heartbeat.set("max_jitter_ms", json::Value(maxJitterMs));
    heartbeat.set("write_ms", json::Value(lastWriteMs));
    root.set("heartbeat", heartbeat);

    const int phaseIdValue = snap.phaseId.load(std::memory_order_relaxed);
    const PhaseSlot* slot = phaseSlot(phaseIdValue);
    json::Value phaseInfo = json::Value::object();
    phaseInfo.set("name", json::Value(slot ? slot->name : "startup"));
    phaseInfo.set("file", json::Value(slot ? slot->file : ""));
    phaseInfo.set("line", json::Value(slot ? slot->line : 0));
    const std::uint64_t phaseStart = snap.phaseStartNs.load(std::memory_order_relaxed);
    const bool phaseActive = snap.phaseActive.load(std::memory_order_relaxed) != 0;
    phaseInfo.set("active", json::Value(phaseActive));
    if (phaseActive)
        phaseInfo.set("age_s", json::Value(phaseStart ? static_cast<double>(nowValue - phaseStart) / 1.0e9 : 0.0));
    else
    {
        const std::uint64_t ended = snap.phaseEndNs.load(std::memory_order_relaxed);
        phaseInfo.set("age_s", json::Value(phaseStart ? static_cast<double>(ended - phaseStart) / 1.0e9 : 0.0));
        phaseInfo.set("since_end_s", json::Value(ended ? static_cast<double>(nowValue - ended) / 1.0e9 : 0.0));
    }
    const char* subName = snap.subName.load(std::memory_order_relaxed);
    if (subName)
    {
        json::Value sub = json::Value::object();
        sub.set("name", json::Value(subName));
        sub.set("file", json::Value(snap.subFile.load(std::memory_order_relaxed) ? snap.subFile.load(std::memory_order_relaxed) : ""));
        sub.set("line", json::Value(snap.subLine.load(std::memory_order_relaxed)));
        const std::uint64_t subStart = snap.subStartNs.load(std::memory_order_relaxed);
        sub.set("age_s", json::Value(subStart ? static_cast<double>(nowValue - subStart) / 1.0e9 : 0.0));
        phaseInfo.set("sub", sub);
    }
    root.set("phase", phaseInfo);

    const std::uint64_t progDone = snap.progDone.load(std::memory_order_relaxed);
    const std::uint64_t progTotal = snap.progTotal.load(std::memory_order_relaxed);
    json::Value progressInfo = json::Value::object();
    progressInfo.set("done", json::Value(static_cast<long long>(progDone)));
    progressInfo.set("total", json::Value(static_cast<long long>(progTotal)));
    progressInfo.set("unit", json::Value(snap.progUnit.load(std::memory_order_relaxed)));
    const std::uint64_t progStart = snap.progStartNs.load(std::memory_order_relaxed);
    double rate = std::nan("");
    if (progStart && nowValue > progStart && progDone > 0)
        rate = static_cast<double>(progDone) / (static_cast<double>(nowValue - progStart) / 1.0e9);
    progressInfo.set("rate_per_s", json::Value(rate));
    const std::uint64_t lastProgress = snap.lastProgressNs.load(std::memory_order_relaxed);
    progressInfo.set("idle_s", json::Value(lastProgress ? static_cast<double>(nowValue - lastProgress) / 1.0e9 : 0.0));
    root.set("progress", progressInfo);

    json::Value step = json::Value::object();
    step.set("index", json::Value(snap.stepIndex.load(std::memory_order_relaxed)));
    step.set("of", json::Value(snap.stepTotal.load(std::memory_order_relaxed)));
    step.set("accepted", json::Value(snap.stepAccepted.load(std::memory_order_relaxed)));
    root.set("step", step);

    json::Value metrics = json::Value::object();
    for (std::size_t i = 0; i < kMetricCount; ++i)
    {
        const double value = snap.metrics[i].load(std::memory_order_relaxed);
        if (std::isnan(value)) continue;
        metrics.set(kMetricNames[i], json::Value(value));
    }
    root.set("metrics", metrics);

    json::Value counters = json::Value::object();
    for (std::size_t i = 0; i < kCountCount; ++i)
    {
        const unsigned long long value = snap.counts[i].load(std::memory_order_relaxed);
        if (value == 0ULL) continue;
        counters.set(kCountNames[i], json::Value(static_cast<long long>(value)));
    }
    root.set("counters", counters);

    json::Value ledgerInfo = json::Value::object();
    ledgerInfo.set("downgrades", json::Value(static_cast<long long>(ledger.countExactly(Severity::Downgrade))));
    ledgerInfo.set("approximations", json::Value(static_cast<long long>(ledger.countExactly(Severity::Approximation))));
    ledgerInfo.set("fail_closed", json::Value(static_cast<long long>(ledger.countExactly(Severity::FailClosed))));
    ledgerInfo.set("strict", json::Value(options.strictNoDowngrade));
    ledgerInfo.set("worst", json::Value(severityName(ledger.worst())));
    ledgerInfo.set("items", ledger.dump());
    root.set("ledger", ledgerInfo);

    root.set("stalled", json::Value(inStall));
    if (inStall) root.set("stalled_s", json::Value(static_cast<double>(nowValue - stallSinceNs) / 1.0e9));

    if (options.sampleStacks) root.set("sampler", sampler.report());

    // ---- process ----
    ProcSample proc = readProcSample(static_cast<int>(::getpid()));
    if (lastCpuSampleNs && nowValue > lastCpuSampleNs)
    {
        const double wallSeconds = static_cast<double>(nowValue - lastCpuSampleNs) / 1.0e9;
        const double tickSeconds = static_cast<double>(proc.cpuTicks - lastCpuTicks) / 100.0;
        if (wallSeconds > 0.0) proc.cpuPercent = 100.0 * tickSeconds / wallSeconds;
    }
    lastCpuTicks = proc.cpuTicks;
    lastCpuSampleNs = nowValue;

    json::Value procInfo = json::Value::object();
    procInfo.set("cpu_percent", json::Value(proc.cpuPercent));
    procInfo.set("rss_mb", json::Value(proc.rssMb));
    procInfo.set("threads", json::Value(static_cast<long long>(proc.threads)));
    procInfo.set("state", json::Value(std::string(1, proc.state)));
    procInfo.set("read_bytes", json::Value(static_cast<long long>(proc.readBytes)));
    procInfo.set("write_bytes", json::Value(static_cast<long long>(proc.writeBytes)));
    root.set("proc", procInfo);

    if (includeSystem)
    {
        const CpuTicks ticks = readCpuTicks();
        json::Value system = json::Value::object();
        if (haveSystemCpu && ticks.total > lastSystemCpu.total)
        {
            const double dTotal = static_cast<double>(ticks.total - lastSystemCpu.total);
            const double dIdle = static_cast<double>(ticks.idle - lastSystemCpu.idle);
            system.set("cpu_percent", json::Value(100.0 * (1.0 - dIdle / std::max(dTotal, 1.0))));
        }
        lastSystemCpu = ticks;
        haveSystemCpu = true;
        system.set("cores", json::Value(static_cast<long long>(::sysconf(_SC_NPROCESSORS_ONLN))));
        const MemorySample memory = readMemory();
        system.set("mem_total_gb", json::Value(memory.totalGb));
        system.set("mem_used_gb", json::Value(memory.usedGb));
        system.set("load1", json::Value(readLoad1()));
        root.set("system", system);
    }
    return root;
}

bool Session::Impl::writeStatus()
{
    const json::Value status = snapshotJson(true);
    return writeAtomic(pathStatus, status.dump(false), options.fsyncStatus);
}

// ---------------------------------------------------------------------------

void Session::Impl::writeSummary(const char* verdict, int exitCode, const json::Value& extra)
{
    if (summaryWritten.exchange(true)) return;
    std::lock_guard<std::mutex> lock(summaryMutex);

    json::Value summary = json::Value::object();
    summary.set("schema", json::Value("rfam.run.summary/1"));
    summary.set("run_id", json::Value(runId));
    summary.set("verdict", json::Value(verdict));
    summary.set("exit_code", json::Value(static_cast<long long>(exitCode)));
    summary.set("finished", json::Value(nowIso()));
    summary.set("uptime_s", json::Value(static_cast<double>(nowNs() - startSteadyNs) / 1.0e9));

    const int phaseIdValue = snap.phaseId.load(std::memory_order_relaxed);
    const PhaseSlot* slot = phaseSlot(phaseIdValue);
    summary.set("last_phase", json::Value(slot ? slot->name : "unknown"));
    summary.set("last_phase_file", json::Value(slot ? slot->file : ""));
    summary.set("last_phase_line", json::Value(slot ? slot->line : 0));
    summary.set("step_index", json::Value(snap.stepIndex.load(std::memory_order_relaxed)));
    summary.set("step_accepted", json::Value(snap.stepAccepted.load(std::memory_order_relaxed)));

    json::Value metrics = json::Value::object();
    for (std::size_t i = 0; i < kMetricCount; ++i)
    {
        const double value = snap.metrics[i].load(std::memory_order_relaxed);
        if (std::isnan(value)) continue;
        metrics.set(kMetricNames[i], json::Value(value));
    }
    summary.set("metrics", metrics);

    json::Value counters = json::Value::object();
    for (std::size_t i = 0; i < kCountCount; ++i)
        counters.set(kCountNames[i], json::Value(static_cast<long long>(snap.counts[i].load(std::memory_order_relaxed))));
    summary.set("counters", counters);

    json::Value ledgerInfo = json::Value::object();
    ledgerInfo.set("downgrades", json::Value(static_cast<long long>(ledger.countExactly(Severity::Downgrade))));
    ledgerInfo.set("approximations", json::Value(static_cast<long long>(ledger.countExactly(Severity::Approximation))));
    ledgerInfo.set("fail_closed", json::Value(static_cast<long long>(ledger.countExactly(Severity::FailClosed))));
    ledgerInfo.set("strict", json::Value(options.strictNoDowngrade));
    ledgerInfo.set("items", ledger.dump());
    summary.set("ledger", ledgerInfo);

    json::Value heartbeat = json::Value::object();
    heartbeat.set("ticks", json::Value(static_cast<long long>(tickCount)));
    heartbeat.set("missed", json::Value(static_cast<long long>(missedTicks)));
    heartbeat.set("max_jitter_ms", json::Value(maxJitterMs));
    heartbeat.set("max_write_ms", json::Value(lastWriteMs));
    summary.set("heartbeat", heartbeat);

    if (options.sampleStacks) summary.set("sampler", sampler.report());
    for (const auto& member : extra.members()) summary.set(member.first, member.second);

    writeAtomic(pathSummary, summary.dump(true), options.fsyncStatus);
    if (options.sampleStacks) writeAtomic(pathStacks, sampler.histogramText(), options.fsyncStatus);

    if (options.mirrorSessionLog)
    {
        std::string line = nowIso();
        line += "  SUMMARY  verdict=";
        line += verdict;
        line += "  exit_code=";
        line += std::to_string(exitCode);
        line += "  downgrades=";
        line += std::to_string(ledger.countExactly(Severity::Downgrade));
        line += "\n";
        appendFile(pathSessionLog, line);
    }
}

void Session::Impl::flushManifest()
{
    std::lock_guard<std::mutex> lock(manifestMutex);
    writeAtomic(pathManifest, manifest.dump(true), options.fsyncStatus);
}

void Session::Impl::pruneOldLink()
{
    // `<root>/latest` is a convenience pointer; recreating it is cheap and its
    // staleness (a link into a deleted directory) is harmless.
    const std::string link = rootDir + "/latest";
    const std::string temp = link + ".tmp." + std::to_string(static_cast<long long>(::getpid()));
    ::unlink(temp.c_str());
    if (::symlink(runId.c_str(), temp.c_str()) == 0)
    {
        if (::rename(temp.c_str(), link.c_str()) != 0) ::unlink(temp.c_str());
    }
}

// ===========================================================================
// Session
// ===========================================================================

Session& Session::instance()
{
    static Session singleton;
    return singleton;
}

Session::Session() : impl_(new Impl()) {}

Session::~Session()
{
    if (!impl_) return;
    if (!impl_->summaryWritten.load() && impl_->active)
    {
        impl_->writeSummary("abandoned", 0, json::Value::object());
    }
    impl_->heartbeatStop.store(true);
    impl_->heartbeatWake.notify_all();
    if (impl_->heartbeatThread.joinable()) impl_->heartbeatThread.join();
    if (g_signalFd.load() >= 0) ::close(g_signalFd.load());
    // impl_ is intentionally leaked: atexit handlers and static destructors may
    // still touch it after this destructor runs.
}

namespace
{
Session* g_session = nullptr;
}

void Session::start(int argc, char** argv, const Options& options)
{
    active_ = false;
    Impl& impl = *impl_;
    impl.options = options;
    g_session = this;

    const char* envMode = std::getenv("RFAM_TELEMETRY");
    if (envMode && (std::strcmp(envMode, "0") == 0 || std::strcmp(envMode, "off") == 0
                    || std::strcmp(envMode, "false") == 0))
        impl.options.enabled = false;

    if (!impl.options.enabled) return;

    // Stack sampling can be switched off independently: it is the only part of
    // this layer that runs code in a signal handler, so it must be separable
    // when a target build behaves differently under it.
    {
        const char* envSample = std::getenv("RFAM_TELEMETRY_SAMPLE");
        if (envSample && (std::strcmp(envSample, "0") == 0 || std::strcmp(envSample, "off") == 0))
            impl.options.sampleStacks = false;
    }
    {
        const char* envSampleHz = std::getenv("RFAM_TELEMETRY_SAMPLE_HZ");
        if (envSampleHz && *envSampleHz)
            impl.options.sampleHz = std::atoi(envSampleHz);
    }
    const char* envHz = std::getenv("RFAM_TELEMETRY_HZ");
    if (envHz && *envHz)
    {
        const double hz = std::strtod(envHz, nullptr);
        if (hz >= 0.1 && hz <= 100.0)
            impl.heartbeatPeriodNs = static_cast<std::int64_t>(1.0e9 / hz);
    }
    if (impl.options.heartbeatMs > 0)
        impl.heartbeatPeriodNs = static_cast<std::int64_t>(impl.options.heartbeatMs) * 1000000LL;

    // ---- run identity ----
    char stamp[32];
    {
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        struct tm tm;
        gmtime_r(&ts.tv_sec, &tm);
        std::strftime(stamp, sizeof(stamp), "%Y%m%dT%H%M%S", &tm);
    }
    impl.runId = impl.options.runName.empty()
        ? std::string(stamp) + "_" + hostName() + "_" + std::to_string(static_cast<long long>(::getpid()))
        : impl.options.runName;

    impl.rootDir = impl.options.rootDir;
    if (impl.rootDir.empty())
    {
        const char* envDir = std::getenv("RFAM_TELEMETRY_DIR");
        impl.rootDir = (envDir && *envDir) ? std::string(envDir) : std::string("runs/telemetry");
    }
    impl.runDir = impl.rootDir + "/" + impl.runId;
    runId_ = impl.runId;
    runDir_ = impl.runDir;

    if (!makeDirs(impl.runDir)) return;   // no writable directory: stay silent

    impl.pathStatus = impl.joinPath("status.json");
    impl.pathEvents = impl.joinPath("events.jsonl");
    impl.pathHeartbeat = impl.joinPath("heartbeat.jsonl");
    impl.pathSessionLog = impl.joinPath("session.log");
    impl.pathManifest = impl.joinPath("manifest.json");
    impl.pathSummary = impl.joinPath("summary.json");
    impl.pathStacks = impl.joinPath("stacks.txt");
    impl.pathSignal = impl.joinPath("signals.log");

    impl.startSteadyNs = nowNs();
    impl.snap.lastProgressNs.store(impl.startSteadyNs);
    impl.snap.phaseStartNs.store(impl.startSteadyNs);

    // ---- argv ----
    for (int i = 0; i < argc; ++i)
    {
        impl.argv.push_back(argv[i] ? argv[i] : "");
        if (argv[i] && argv[i][0] == '-' && argv[i][1] == '-')
        {
            std::string flag(argv[i]);
            const std::size_t equals = flag.find('=');
            if (equals != std::string::npos) flag = flag.substr(0, equals);
            impl.argFlags.push_back(flag);
        }
    }

    // ---- manifest (provisional: everything known before the driver resolves
    //      mesh, thermo and chemistry) ----
    impl.manifest = json::Value::object();
    impl.manifest.set("schema", json::Value("rfam.run.manifest/1"));
    impl.manifest.set("run_id", json::Value(impl.runId));
    impl.manifest.set("started", json::Value(nowIso()));
    impl.manifest.set("host", json::Value(hostName()));
    impl.manifest.set("user", json::Value(userName()));
    impl.manifest.set("pid", json::Value(static_cast<long long>(::getpid())));
    impl.manifest.set("cwd", json::Value(cwdPath()));

    json::Value binary = json::Value::object();
    const std::string exe = exePath();
    binary.set("path", json::Value(exe));
    if (!exe.empty())
    {
        struct stat info;
        if (::stat(exe.c_str(), &info) == 0)
        {
            binary.set("bytes", json::Value(static_cast<long long>(info.st_size)));
            binary.set("mtime", json::Value(static_cast<long long>(info.st_mtime)));
            const std::string bytes = readFile(exe);
            if (!bytes.empty())
                binary.set("sha256", json::Value(thermoProvenance::sha256Hex(bytes)));
        }
    }
#ifdef RFAM_GIT_COMMIT
    binary.set("git_commit", json::Value(RFAM_GIT_COMMIT));
#endif
    binary.set("compiled", json::Value(std::string(__DATE__) + " " + __TIME__));
    binary.set("compiler", json::Value(std::string(__VERSION__)));
    impl.manifest.set("binary", binary);

    json::Value argvValue = json::Value::array();
    for (const auto& item : impl.argv) argvValue.push(json::Value(item));
    impl.manifest.set("argv", argvValue);

    json::Value env = json::Value::object();
    const char* keys[] = {"OMP_NUM_THREADS", "RFAM_TETRA_TABLE", "RFAM_STATE_ACCEPTANCE_REL",
                          "RFAM_TELEMETRY_DIR", "RFAM_GIT_COMMIT", "SUNDIALS_ROOT", ""};
    for (int i = 0; keys[i][0]; ++i)
    {
        const char* value = std::getenv(keys[i]);
        if (value) env.set(keys[i], json::Value(std::string(value)));
    }
    impl.manifest.set("env", env);
    impl.manifest.set("telemetry", json::Value::object());
    impl.manifest.set("effective", json::Value::object());
    impl.manifest.set("files", json::Value::object());
    impl.flushManifest();
    impl.pruneOldLink();

    // ---- signal handlers ----
    const std::string signalPath = impl.pathSignal;
    const int signalFd = ::open(signalPath.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0664);
    g_signalFd.store(signalFd);
    struct sigaction action;
    std::memset(&action, 0, sizeof(action));
    action.sa_handler = rfamTelemetrySignalHandler;
    sigemptyset(&action.sa_mask);
    action.sa_flags = SA_RESTART;
    for (int signo : {SIGINT, SIGTERM, SIGHUP, SIGQUIT})
        ::sigaction(signo, &action, nullptr);
    ::signal(SIGPIPE, SIG_IGN);

    // ---- stack sampler ----
    if (impl.options.sampleStacks)
        impl.sampler.start(impl.options.sampleHz, impl.runDir);

    // ---- heartbeat ----
    active_ = true;
    impl.active = true;
    {
        json::Value startEvent = json::Value::object();
        startEvent.set("run_dir", json::Value(impl.runDir));
        startEvent.set("argv", argvValue);
        impl.emitEvent(startEvent, "START");
    }
    impl.heartbeatThread = std::thread([&impl]() { impl.runHeartbeat(); });
    pthread_setname_np(impl.heartbeatThread.native_handle(), "rfam-heartbeat");
}

json::Value& Session::describe()
{
    return impl_->manifest;
}

void Session::defineField(const std::string& key, const std::string& value,
                          const char* flag, const char* forceSource)
{
    if (!impl_ || !impl_->active) return;
    std::lock_guard<std::mutex> lock(impl_->manifestMutex);
    std::string spelling = flag ? std::string(flag) : ("--" + key);
    if (!flag)
        for (char& c : spelling) if (c == '_') c = '-';
    const char* source = forceSource;
    if (!source) source = impl_->flagPresent(spelling) ? "cli" : "default";

    json::Value item = json::Value::object();
    item.set("value", json::Value(value));
    item.set("source", json::Value(source));
    impl_->manifest["effective"][key] = item;
}

void Session::flushManifest()
{
    if (!impl_ || !impl_->active) return;
    impl_->flushManifest();
}

// ---------------------------------------------------------------------------
// Hot path
// ---------------------------------------------------------------------------

void Session::enter(const Phase& p)
{
    if (!active_) return;
    const std::uint64_t now = nowNs();
    impl_->snap.phaseId.store(p.id, std::memory_order_relaxed);
    impl_->snap.phaseActive.store(1, std::memory_order_relaxed);
    impl_->snap.phaseStartNs.store(now, std::memory_order_relaxed);
    impl_->snap.subName.store(nullptr, std::memory_order_relaxed);
    impl_->snap.subFile.store(nullptr, std::memory_order_relaxed);
    impl_->snap.lastProgressNs.store(now, std::memory_order_relaxed);

    json::Value fields = json::Value::object();
    fields.set("phase", json::Value(p.name));
    fields.set("file", json::Value(p.file));
    fields.set("line", json::Value(p.line));
    fields.set("step", json::Value(impl_->snap.stepIndex.load(std::memory_order_relaxed)));
    impl_->emitEvent(fields, "PHASE_BEGIN");
}

void Session::leave(const Phase& p)
{
    if (!active_) return;
    const std::uint64_t now = nowNs();
    const std::uint64_t started = impl_->snap.phaseStartNs.load(std::memory_order_relaxed);
    // A phase that has ended but was not followed by another one must be
    // labelled as such: reporting "thermo_setup, 30 s" for a run that finished
    // thermo setup in 0.7 s is exactly the kind of true-looking lie this layer
    // exists to prevent.
    impl_->snap.phaseActive.store(0, std::memory_order_relaxed);
    impl_->snap.phaseEndNs.store(now, std::memory_order_relaxed);
    json::Value fields = json::Value::object();
    fields.set("phase", json::Value(p.name));
    fields.set("seconds", json::Value(static_cast<double>(now - started) / 1.0e9));
    impl_->emitEvent(fields, "PHASE_END");
    impl_->snap.lastProgressNs.store(now, std::memory_order_relaxed);
}

void Session::subPhase(const char* name, const char* file, int line)
{
    if (!active_) return;
    impl_->snap.subName.store(name, std::memory_order_relaxed);
    impl_->snap.subFile.store(file, std::memory_order_relaxed);
    impl_->snap.subLine.store(line, std::memory_order_relaxed);
    impl_->snap.subStartNs.store(nowNs(), std::memory_order_relaxed);
}

void Session::endSubPhase()
{
    if (!active_) return;
    impl_->snap.subName.store(nullptr, std::memory_order_relaxed);
}

void Session::setStep(long long index, long long total, long long accepted)
{
    if (!active_) return;
    impl_->snap.stepIndex.store(index, std::memory_order_relaxed);
    impl_->snap.stepTotal.store(total, std::memory_order_relaxed);
    impl_->snap.stepAccepted.store(accepted, std::memory_order_relaxed);
    impl_->snap.lastProgressNs.store(nowNs(), std::memory_order_relaxed);
}

void Session::metric(Metric m, double value)
{
    if (!active_) return;
    const std::size_t index = static_cast<std::size_t>(m);
    if (index >= kMetricCount) return;
    impl_->snap.metrics[index].store(value, std::memory_order_relaxed);
}

void Session::count(Count c, unsigned long long value)
{
    if (!active_) return;
    const std::size_t index = static_cast<std::size_t>(c);
    if (index >= kCountCount) return;
    impl_->snap.counts[index].store(value, std::memory_order_relaxed);
}

void Session::progress(unsigned long long done, unsigned long long total, const char* unit)
{
    if (!active_) return;
    if (done == 0)
    {
        impl_->snap.progStartNs.store(nowNs(), std::memory_order_relaxed);
        impl_->snap.progDone.store(0, std::memory_order_relaxed);
    }
    impl_->snap.progTotal.store(total, std::memory_order_relaxed);
    impl_->snap.progUnit.store(unit, std::memory_order_relaxed);
    impl_->snap.progDone.store(done, std::memory_order_relaxed);
    impl_->snap.lastProgressNs.store(nowNs(), std::memory_order_relaxed);
    impl_->sampler.captureIfDue(unit);
}

void Session::bumpProgress(unsigned long long n)
{
    if (!active_) return;
    impl_->snap.progDone.fetch_add(n, std::memory_order_relaxed);
    impl_->snap.lastProgressNs.store(nowNs(), std::memory_order_relaxed);
    impl_->sampler.captureIfDue(nullptr);
}

void Session::noteProgress()
{
    if (!active_) return;
    impl_->snap.lastProgressNs.store(nowNs(), std::memory_order_relaxed);
    impl_->sampler.captureIfDue(nullptr);
}

// ---------------------------------------------------------------------------
// Events
// ---------------------------------------------------------------------------

void Session::event(const char* type, const json::Value& fields)
{
    if (!active_) return;
    impl_->emitEvent(fields, type);
}

void Session::degrade(const char* code, Severity severity, const std::string& detail)
{
    if (!active_) return;
    const bool firstTime = impl_->ledger.note(code, severity, detail);
    if (firstTime)
    {
        json::Value fields = json::Value::object();
        fields.set("code", json::Value(code));
        fields.set("severity", json::Value(severityName(severity)));
        fields.set("detail", json::Value(detail));
        impl_->emitEvent(fields, "DOWNGRADE");
    }

    // --strict-no-downgrade: a run that had to approximate is not a run whose
    // result may be quoted. FailClosed entries are recorded only; their caller
    // is already aborting, and turning a clean exception into _exit() would
    // throw away the driver's own failure diagnostics.
    if (impl_->options.strictNoDowngrade
        && static_cast<int>(severity) >= static_cast<int>(Severity::Downgrade))
    {
        std::fprintf(stdout, "RFAM_STRICT_DOWNGRADE_FAILED code=%s detail=%s\n",
                     code, detail.c_str());
        std::fflush(stdout);
        impl_->writeSummary("strict_downgrade", 3, json::Value::object());
        ::_exit(3);
    }
}

void Session::warning(const std::string& text)
{
    if (!active_) return;
    json::Value fields = json::Value::object();
    fields.set("text", json::Value(text));
    impl_->emitEvent(fields, "WARNING");
}

void Session::checkpoint(const std::string& path, const std::string& kind)
{
    if (!active_) return;
    json::Value fields = json::Value::object();
    fields.set("path", json::Value(path));
    fields.set("kind", json::Value(kind));
    impl_->emitEvent(fields, "CHECKPOINT");
}

void Session::finish(const char* verdict, int exitCode, const json::Value& extra)
{
    if (!impl_ || !impl_->active) return;
    if (impl_->finishing.exchange(true)) return;
    impl_->writeSummary(verdict, exitCode, extra);
    impl_->heartbeatStop.store(true);
    impl_->heartbeatWake.notify_all();
    if (impl_->heartbeatThread.joinable()
        && impl_->heartbeatThread.get_id() != std::this_thread::get_id())
        impl_->heartbeatThread.join();
    impl_->sampler.stop();
    impl_->active = false;
    active_ = false;
}

bool Session::terminationRequested() const
{
    return g_signalNumber.load(std::memory_order_relaxed) != 0;
}

int Session::terminationSignal() const
{
    return g_signalNumber.load(std::memory_order_relaxed);
}

json::Value Session::samplerReport() const
{
    if (!impl_ || !impl_->active) return json::Value::object();
    return impl_->sampler.report();
}

} // End namespace telemetry
} // End namespace Foam
