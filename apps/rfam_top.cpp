// rfam_top -- a btop-style dashboard for a running RFAM simulation.
//
// It reads the telemetry directory written by Foam::telemetry (status.json,
// manifest.json, events.jsonl, summary.json) and samples the host's own
// counters. Two rules shape everything here:
//
//   1. The dashboard never invents a number. A field the solver did not write
//      renders as "--"; a stale status file renders as STALE with its real age;
//      a dead PID renders as PROCESS GONE. A green bar over a dead run would be
//      worse than no dashboard at all.
//   2. The layout never jumps. Column widths are fixed, numeric precision is
//      fixed, and the sparkline keeps a constant length, so the same run looks
//      the same at every refresh -- which is what makes a change in the picture
//      meaningful.

#include "CursesCompat.H"
#include "JsonLite.H"
#include "RunTelemetry.H"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <locale.h>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <dirent.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace
{

// ===========================================================================
// Options
// ===========================================================================

struct Options
{
    std::string runDir;
    std::string rootDir = "runs/telemetry";
    int refreshMs = 250;
    int dataMs = 1000;
    bool ascii = false;
    bool plain = false;
    bool once = false;
    bool noGpu = false;
};

Options g_options;
std::atomic<bool> g_quit(false);

void onSignal(int) { g_quit.store(true); }

// ===========================================================================
// Data model
// ===========================================================================

struct LedgerItem
{
    std::string code;
    std::string severity;
    std::string detail;
    long long n = 0;
};

struct Live
{
    bool valid = false;
    std::string error;
    Foam::json::Value root;

    std::string runId;
    long long pid = 0;
    long long seq = 0;
    std::string wall;
    double uptime = 0.0;
    double fileAge = 1e9;      // seconds since status.json was last replaced

    std::string phase = "--";
    std::string phaseFile;
    long long phaseLine = 0;
    double phaseAge = 0.0;
    bool phaseActive = false;
    double phaseSinceEnd = 0.0;
    std::string sub;
    std::string subFile;
    long long subLine = 0;
    double subAge = 0.0;

    long long done = 0, total = 0;
    std::string unit = "--";
    double rate = NAN;
    double idle = 0.0;

    long long step = -1, stepOf = -1, accepted = 0;

    std::map<std::string, double> metrics;
    std::map<std::string, long long> counters;

    long long downgrades = 0, approximations = 0, failClosed = 0;
    std::string worst = "--";
    std::vector<LedgerItem> items;

    bool stalled = false;
    double stalledS = 0.0;

    long long hbTicks = 0, hbMissed = 0;
    double hbJitter = NAN, hbMaxJitter = NAN, hbWriteMs = NAN;
    long long hbPeriodMs = 1000;

    double cpu = NAN, rss = NAN;
    long long threads = -1;
    std::string procState = "?";

    std::string stackNow;
    long long stackSamples = 0;

    double sysCpu = NAN, memUsed = NAN, memTotal = NAN, load1 = NAN;
    long long cores = 0;

    bool processAlive = false;
    bool hasSummary = false;
    std::string verdict;
};

struct EventLine
{
    std::string time;
    std::string type;
    std::string text;
};

struct StepPoint
{
    double ratio = NAN;
    double residual = NAN;
    long long step = -1;
};

struct HostState
{
    std::vector<double> corePercent;
    double cpuPercent = NAN;
    double memUsedGb = NAN;
    double memTotalGb = NAN;
    double load1 = NAN;
    double procCpu = NAN;
    double procRssMb = NAN;
    long long procThreads = -1;
    char procState = '?';
    bool gpuAvailable = false;
    std::string gpuName;
    double gpuUtil = NAN;
    double gpuMemUsed = NAN;
    double gpuMemTotal = NAN;
    double gpuTemp = NAN;
    std::string gpuError;
};

HostState g_host;
std::mutex g_hostMutex;
std::atomic<bool> g_hostStop(false);
std::atomic<long long> g_targetPid(0);

// ===========================================================================
// Small utilities
// ===========================================================================

double nowSeconds()
{
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

bool readWholeFile(const std::string& path, std::string& out)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) return false;
    std::ostringstream buffer;
    buffer << input.rdbuf();
    out = buffer.str();
    return true;
}

bool fileAgeSeconds(const std::string& path, double& age)
{
    struct stat info;
    if (::stat(path.c_str(), &info) != 0) return false;
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    age = static_cast<double>(now.tv_sec - info.st_mtim.tv_sec)
        + static_cast<double>(now.tv_nsec - info.st_mtim.tv_nsec) / 1e9;
    if (age < 0.0) age = 0.0;
    return true;
}

bool processExists(long long pid)
{
    if (pid <= 0) return false;
    struct stat info;
    return ::stat(("/proc/" + std::to_string(pid)).c_str(), &info) == 0;
}

std::string formatDuration(double seconds)
{
    if (!std::isfinite(seconds) || seconds < 0) return "--";
    const long long total = static_cast<long long>(seconds);
    const long long hours = total / 3600;
    const long long minutes = (total % 3600) / 60;
    const long long secs = total % 60;
    char buffer[32];
    if (hours > 0) std::snprintf(buffer, sizeof(buffer), "%lldh%02lldm%02llds", hours, minutes, secs);
    else if (minutes > 0) std::snprintf(buffer, sizeof(buffer), "%lldm%02llds", minutes, secs);
    else std::snprintf(buffer, sizeof(buffer), "%llds", secs);
    return std::string(buffer);
}

std::string formatNumber(double value, int digits = 3)
{
    if (!std::isfinite(value)) return "--";
    char buffer[48];
    const double magnitude = std::fabs(value);
    if (magnitude != 0.0 && (magnitude < 1e-3 || magnitude >= 1e6))
        std::snprintf(buffer, sizeof(buffer), "%.*e", std::max(1, digits), value);
    else
        std::snprintf(buffer, sizeof(buffer), "%.*f", digits, value);
    return std::string(buffer);
}

std::string formatCount(double value)
{
    if (!std::isfinite(value)) return "--";
    char buffer[48];
    if (std::fabs(value) >= 1e9) std::snprintf(buffer, sizeof(buffer), "%.2fG", value / 1e9);
    else if (std::fabs(value) >= 1e6) std::snprintf(buffer, sizeof(buffer), "%.2fM", value / 1e6);
    else if (std::fabs(value) >= 1e3) std::snprintf(buffer, sizeof(buffer), "%.1fk", value / 1e3);
    else std::snprintf(buffer, sizeof(buffer), "%.0f", value);
    return std::string(buffer);
}

std::string padRight(const std::string& text, std::size_t width)
{
    if (text.size() >= width) return text.substr(0, width);
    return text + std::string(width - text.size(), ' ');
}

std::string padLeft(const std::string& text, std::size_t width)
{
    if (text.size() >= width) return text.substr(0, width);
    return std::string(width - text.size(), ' ') + text;
}

std::string elide(const std::string& text, std::size_t width)
{
    if (text.size() <= width) return text;
    if (width <= 3) return text.substr(0, width);
    return text.substr(0, width - 3) + "...";
}

std::string shortName(const std::string& path)
{
    // Show the repository-relative tail of a source path: the full path is in
    // the manifest and the events; the dashboard only needs the file identity.
    const std::size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return path;
    const std::string file = path.substr(slash + 1);
    const std::size_t previous = slash == 0 ? std::string::npos : path.find_last_of('/', slash - 1);
    if (previous == std::string::npos || previous == 0) return file;
    return path.substr(previous + 1);
}

// ===========================================================================
// Reading telemetry files
// ===========================================================================

Live loadLive(const std::string& runDir)
{
    Live live;
    const std::string path = runDir + "/status.json";
    std::string text;
    if (!readWholeFile(path, text))
    {
        live.error = "status.json not found";
        return live;
    }
    std::string error;
    if (!Foam::json::parse(text, live.root, error))
    {
        live.error = "status.json unreadable: " + error;
        return live;
    }
    double age = 1e9;
    fileAgeSeconds(path, age);
    live.fileAge = age;

    const Foam::json::Value& root = live.root;
    live.valid = true;
    live.runId = root["run_id"].asString();
    live.pid = root["pid"].asInt();
    live.seq = root["seq"].asInt();
    live.wall = root["wall"].asString();
    live.uptime = root["uptime_s"].asDouble();

    live.phase = root["phase"]["name"].asString("--");
    live.phaseFile = root["phase"]["file"].asString();
    live.phaseLine = root["phase"]["line"].asInt();
    live.phaseAge = root["phase"]["age_s"].asDouble();
    live.phaseActive = root["phase"]["active"].asBool(false);
    live.phaseSinceEnd = root["phase"]["since_end_s"].asDouble(0.0);
    if (root["phase"]["sub"].isObject())
    {
        live.sub = root["phase"]["sub"]["name"].asString();
        live.subFile = root["phase"]["sub"]["file"].asString();
        live.subLine = root["phase"]["sub"]["line"].asInt();
        live.subAge = root["phase"]["sub"]["age_s"].asDouble();
    }

    live.done = root["progress"]["done"].asInt();
    live.total = root["progress"]["total"].asInt();
    live.unit = root["progress"]["unit"].asString("--");
    live.rate = root["progress"]["rate_per_s"].asDouble(NAN);
    live.idle = root["progress"]["idle_s"].asDouble();

    live.step = root["step"]["index"].asInt(-1);
    live.stepOf = root["step"]["of"].asInt(-1);
    live.accepted = root["step"]["accepted"].asInt();

    for (const auto& member : root["metrics"].members())
        live.metrics[member.first] = member.second.asDouble(NAN);
    for (const auto& member : root["counters"].members())
        live.counters[member.first] = member.second.asInt();

    live.downgrades = root["ledger"]["downgrades"].asInt();
    live.approximations = root["ledger"]["approximations"].asInt();
    live.failClosed = root["ledger"]["fail_closed"].asInt();
    live.worst = root["ledger"]["worst"].asString("--");
    for (const auto& item : root["ledger"]["items"].elements())
    {
        LedgerItem entry;
        entry.code = item["code"].asString();
        entry.severity = item["severity"].asString();
        entry.detail = item["detail"].asString();
        entry.n = item["n"].asInt();
        live.items.push_back(entry);
    }

    live.stalled = root["stalled"].asBool();
    live.stalledS = root["stalled_s"].asDouble();

    live.hbTicks = root["heartbeat"]["ticks"].asInt();
    live.hbMissed = root["heartbeat"]["missed"].asInt();
    live.hbPeriodMs = root["heartbeat"]["period_ms"].asInt(1000);
    live.hbJitter = root["heartbeat"]["jitter_ms"].asDouble(NAN);
    live.hbMaxJitter = root["heartbeat"]["max_jitter_ms"].asDouble(NAN);
    live.hbWriteMs = root["heartbeat"]["write_ms"].asDouble(NAN);

    live.cpu = root["proc"]["cpu_percent"].asDouble(NAN);
    live.rss = root["proc"]["rss_mb"].asDouble(NAN);
    live.threads = root["proc"]["threads"].asInt(-1);
    live.procState = root["proc"]["state"].asString("?");

    live.stackSamples = root["sampler"]["samples"].asInt();
    if (root["sampler"]["current"].isObject())
        live.stackNow = root["sampler"]["current"]["stack"].asString();

    live.sysCpu = root["system"]["cpu_percent"].asDouble(NAN);
    live.memUsed = root["system"]["mem_used_gb"].asDouble(NAN);
    live.memTotal = root["system"]["mem_total_gb"].asDouble(NAN);
    live.load1 = root["system"]["load1"].asDouble(NAN);
    live.cores = root["system"]["cores"].asInt();

    live.processAlive = processExists(live.pid);
    return live;
}

void readEvents(const std::string& runDir, std::vector<EventLine>& events,
                std::vector<StepPoint>& steps)
{
    std::string text;
    if (!readWholeFile(runDir + "/events.jsonl", text)) return;
    std::istringstream stream(text);
    std::string line;
    events.clear();
    steps.clear();
    while (std::getline(stream, line))
    {
        if (line.empty()) continue;
        Foam::json::Value record;
        std::string error;
        if (!Foam::json::parse(line, record, error)) continue;
        EventLine event;
        event.time = record["t"].asString();
        event.type = record["type"].asString();
        std::ostringstream detail;
        for (const auto& member : record.members())
        {
            if (member.first == "t" || member.first == "type" || member.first == "pid") continue;
            if (!detail.str().empty()) detail << ' ';
            detail << member.first << '=';
            if (member.second.isString()) detail << member.second.asString();
            else detail << member.second.dump(false);
        }
        event.text = detail.str();
        if (event.type == "STEP")
        {
            StepPoint point;
            point.step = record["step"].asInt(-1);
            point.residual = record["residual"].asDouble(NAN);
            point.ratio = record["ratio"].asDouble(NAN);
            steps.push_back(point);
        }
        events.push_back(event);
    }
}

std::string loadSummaryVerdict(const std::string& runDir, bool& found)
{
    std::string text;
    found = false;
    if (!readWholeFile(runDir + "/summary.json", text)) return std::string();
    Foam::json::Value summary;
    std::string error;
    if (!Foam::json::parse(text, summary, error)) return std::string();
    found = true;
    return summary["verdict"].asString("unknown");
}

// ===========================================================================
// Host sampling (independent thread: never blocks the UI)
// ===========================================================================

struct CpuLine
{
    unsigned long long total = 0;
    unsigned long long idle = 0;
};

std::vector<CpuLine> readCpuLines()
{
    std::vector<CpuLine> lines;
    std::ifstream input("/proc/stat");
    std::string line;
    while (std::getline(input, line))
    {
        if (line.compare(0, 3, "cpu") != 0) break;
        CpuLine entry;
        std::istringstream stream(line.substr(3));
        unsigned long long values[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        for (int i = 0; i < 10; ++i) if (!(stream >> values[i])) break;
        for (int i = 0; i < 10; ++i) entry.total += values[i];
        entry.idle = values[3] + values[4];
        lines.push_back(entry);
    }
    return lines;
}

void readProcHost(long long pid, HostState& state)
{
    if (pid <= 0) return;
    std::string text;
    if (!readWholeFile("/proc/" + std::to_string(pid) + "/stat", text)) return;
    const std::size_t close = text.rfind(')');
    if (close == std::string::npos) return;
    std::istringstream stream(text.substr(close + 1));
    std::vector<std::string> fields;
    std::string token;
    while (stream >> token) fields.push_back(token);
    if (fields.size() > 21)
    {
        state.procState = fields[0].empty() ? '?' : fields[0][0];
        state.procThreads = std::atol(fields[17].c_str());
        const long pages = std::atol(fields[21].c_str());
        state.procRssMb = static_cast<double>(pages) * static_cast<double>(::sysconf(_SC_PAGESIZE))
                        / (1024.0 * 1024.0);
    }
}

void readGpu(HostState& state)
{
    FILE* pipe = ::popen(
        "nvidia-smi --query-gpu=name,utilization.gpu,memory.used,memory.total,temperature.gpu"
        " --format=csv,noheader,nounits 2>/dev/null", "r");
    if (!pipe)
    {
        state.gpuError = "nvidia-smi not found";
        return;
    }
    char buffer[512];
    std::string first;
    if (std::fgets(buffer, sizeof(buffer), pipe)) first = buffer;
    const int status = ::pclose(pipe);
    if (status != 0 || first.empty())
    {
        state.gpuError = "nvidia-smi unavailable";
        state.gpuAvailable = false;
        return;
    }
    // name may contain commas; the five numeric fields are last.
    std::vector<std::string> parts;
    std::istringstream stream(first);
    std::string part;
    while (std::getline(stream, part, ',')) parts.push_back(part);
    if (parts.size() < 5)
    {
        state.gpuError = "unexpected nvidia-smi output";
        state.gpuAvailable = false;
        return;
    }
    state.gpuAvailable = true;
    state.gpuName = parts[0];
    state.gpuUtil = std::atof(parts[parts.size() - 4].c_str());
    state.gpuMemUsed = std::atof(parts[parts.size() - 3].c_str());
    state.gpuMemTotal = std::atof(parts[parts.size() - 2].c_str());
    state.gpuTemp = std::atof(parts[parts.size() - 1].c_str());
    state.gpuError.clear();
}

void hostSamplerLoop()
{
    std::vector<CpuLine> previous = readCpuLines();
    int tick = 0;
    while (!g_hostStop.load(std::memory_order_relaxed))
    {
        const std::vector<CpuLine> current = readCpuLines();
        HostState state;
        {
            std::lock_guard<std::mutex> lock(g_hostMutex);
            state = g_host;
        }
        state.corePercent.assign(current.size() > 1 ? current.size() - 1 : 0, 0.0);
        double totalDelta = 0.0, idleDelta = 0.0;
        for (std::size_t i = 1; i < current.size() && i < previous.size(); ++i)
        {
            const double dTotal = static_cast<double>(current[i].total - previous[i].total);
            const double dIdle = static_cast<double>(current[i].idle - previous[i].idle);
            if (dTotal > 0.0) state.corePercent[i - 1] = 100.0 * (1.0 - dIdle / dTotal);
            totalDelta += dTotal;
            idleDelta += dIdle;
        }
        if (totalDelta > 0.0) state.cpuPercent = 100.0 * (1.0 - idleDelta / totalDelta);

        std::ifstream meminfo("/proc/meminfo");
        std::string line;
        double totalKb = 0, availableKb = 0;
        while (std::getline(meminfo, line))
        {
            if (line.compare(0, 9, "MemTotal:") == 0) totalKb = std::strtod(line.c_str() + 9, nullptr);
            else if (line.compare(0, 13, "MemAvailable:") == 0) availableKb = std::strtod(line.c_str() + 13, nullptr);
        }
        if (totalKb > 0)
        {
            state.memTotalGb = totalKb / (1024.0 * 1024.0);
            state.memUsedGb = (totalKb - availableKb) / (1024.0 * 1024.0);
        }
        std::ifstream load("/proc/loadavg");
        double load1 = 0;
        if (load >> load1) state.load1 = load1;

        readProcHost(g_targetPid.load(), state);
        if (!g_options.noGpu && (tick % 4) == 0) readGpu(state);

        {
            std::lock_guard<std::mutex> lock(g_hostMutex);
            g_host = state;
        }
        previous = current;
        ++tick;
        for (int slept = 0; slept < 250 && !g_hostStop.load(); slept += 25)
            std::this_thread::sleep_for(std::chrono::milliseconds(25));
    }
}

// ===========================================================================
// Rendering primitives
// ===========================================================================

struct Rect
{
    int y = 0, x = 0, h = 0, w = 0;
    bool empty() const { return h <= 0 || w <= 0; }
    Rect inner() const { return Rect{y + 1, x + 1, h - 2, w - 2}; }
};

int g_attrNormal = 0;
int g_attrTitle = 0;
int g_attrOk = 0;
int g_attrWarn = 0;
int g_attrBad = 0;
int g_attrDim = 0;
int g_attrAccent = 0;
int g_attrBar = 0;
int g_attrBarHot = 0;

bool g_color = false;

void put(int y, int x, const std::string& text, int attr = 0)
{
    if (y < 0 || x < 0) return;
    if (y >= getmaxy(stdscr) || x >= getmaxx(stdscr)) return;
    if (attr) wattron(stdscr, attr);
    wmove(stdscr, y, x);
    waddnstr(stdscr, text.c_str(), static_cast<int>(text.size()));
    if (attr) wattroff(stdscr, attr);
}

void fill(int y, int x, int h, int w, char c = ' ')
{
    if (h <= 0 || w <= 0) return;
    const std::string blank(static_cast<std::size_t>(w), c);
    for (int row = 0; row < h; ++row) put(y + row, x, blank);
}

const char* boxChars(int index)
{
    static const char* utf8[] = {"\u2502", "\u2500", "\u250c", "\u2510", "\u2514", "\u2518", "\u251c", "\u2524"};
    static const char* ascii[] = {"|", "-", "+", "+", "+", "+", "+", "+"};
    return g_options.ascii ? ascii[index] : utf8[index];
}

void drawBox(const Rect& rect, const std::string& title, int titleAttr)
{
    if (rect.h < 2 || rect.w < 4) return;
    fill(rect.y, rect.x, rect.h, rect.w, ' ');
    // horizontal edges
    std::string top = boxChars(2);
    top += std::string(static_cast<std::size_t>(rect.w - 2), '-');
    top += boxChars(3);
    std::string bottom = boxChars(4);
    bottom += std::string(static_cast<std::size_t>(rect.w - 2), '-');
    bottom += boxChars(5);
    put(rect.y, rect.x, top);
    put(rect.y + rect.h - 1, rect.x, bottom);
    for (int row = 1; row < rect.h - 1; ++row)
    {
        put(rect.y + row, rect.x, boxChars(0));
        put(rect.y + row, rect.x + rect.w - 1, boxChars(0));
    }
    if (!title.empty() && rect.w > 6)
    {
        put(rect.y, rect.x + 2, " " + elide(title, static_cast<std::size_t>(rect.w - 5)) + " ", titleAttr);
    }
}

std::string bar(double fraction, int width, bool ascii)
{
    fraction = std::max(0.0, std::min(1.0, fraction));
    const int filled = static_cast<int>(std::lround(fraction * width));
    if (ascii)
        return std::string(static_cast<std::size_t>(filled), '#')
             + std::string(static_cast<std::size_t>(width - filled), '.');
    static const char* blocks[] = {"\u2588", "\u2589", "\u258a", "\u258b", "\u258c", "\u258d", "\u258e", "\u258f"};
    (void)blocks;
    return std::string(static_cast<std::size_t>(filled), '#')
         + std::string(static_cast<std::size_t>(width - filled), '.');
}

std::string sparkline(const std::vector<double>& values, int width, bool ascii)
{
    if (values.empty() || width <= 0) return std::string();
    const std::size_t count = std::min<std::size_t>(values.size(), static_cast<std::size_t>(width));
    const std::size_t first = values.size() - count;
    double low = values[first], high = values[first];
    for (std::size_t i = first; i < values.size(); ++i)
    {
        low = std::min(low, values[i]);
        high = std::max(high, values[i]);
    }
    const double span = (high - low) > 0 ? (high - low) : 1.0;
    static const char* utf8[] = {"\u2581", "\u2582", "\u2583", "\u2584", "\u2585", "\u2586", "\u2587", "\u2588"};
    static const char* asciiLevels[] = {"_", ".", "-", "=", "+", "*", "#", "@"};
    std::string out;
    for (std::size_t i = first; i < values.size(); ++i)
    {
        const double normalised = (values[i] - low) / span;
        int level = static_cast<int>(normalised * 7.0);
        level = std::max(0, std::min(7, level));
        out += ascii ? asciiLevels[level] : utf8[level];
    }
    return out;
}

// ===========================================================================
// Dashboard state
// ===========================================================================

struct Dashboard
{
    std::string runDir;
    Live live;
    std::vector<EventLine> events;
    std::vector<StepPoint> steps;
    std::vector<std::string> eventView;
    Foam::json::Value manifest;
    bool manifestLoaded = false;
    double manifestAge = -1.0;

    enum class Overlay { None, Manifest, Events, Summary, Help };
    Overlay overlay = Overlay::None;
    int overlayScroll = 0;
    std::vector<std::string> overlayLines;

    bool paused = false;
    double lastDataLoad = 0.0;
    double startedAt = nowSeconds();
};

std::vector<double> residentSeries(const std::vector<StepPoint>& steps, int maximum)
{
    std::vector<double> values;
    const std::size_t count = std::min<std::size_t>(steps.size(), static_cast<std::size_t>(maximum));
    for (std::size_t i = steps.size() - count; i < steps.size(); ++i)
    {
        double value = steps[i].ratio;
        if (!std::isfinite(value)) value = steps[i].residual;
        if (std::isfinite(value)) values.push_back(value);
    }
    return values;
}

void reloadData(Dashboard& board)
{
    const Live previous = board.live;
    board.live = loadLive(board.runDir);
    if (board.live.valid && board.live.pid > 0) g_targetPid.store(board.live.pid);
    readEvents(board.runDir, board.events, board.steps);

    bool hasSummary = false;
    const std::string verdict = loadSummaryVerdict(board.runDir, hasSummary);
    board.live.hasSummary = hasSummary;
    board.live.verdict = verdict;

    double age = -1.0;
    if (fileAgeSeconds(board.runDir + "/manifest.json", age)) board.manifestAge = age;
    if (!board.manifestLoaded)
    {
        std::string text;
        if (readWholeFile(board.runDir + "/manifest.json", text))
        {
            std::string error;
            if (Foam::json::parse(text, board.manifest, error)) board.manifestLoaded = true;
        }
    }
    (void)previous;
}

// ===========================================================================
// Panels
// ===========================================================================

void drawHeader(int width, const Dashboard& board)
{
    const Live& live = board.live;
    std::string state = "STARTING";
    int attr = g_attrDim;
    if (live.valid)
    {
        if (!live.processAlive) { state = "PROCESS GONE"; attr = g_attrBad; }
        else if (live.hasSummary) { state = "FINISHED " + live.verdict; attr = g_attrOk; }
        else if (live.fileAge > 3.0 * (live.hbPeriodMs / 1000.0)) { state = "TELEMETRY STALE"; attr = g_attrBad; }
        else if (live.stalled) { state = "STALLED"; attr = g_attrBad; }
        else { state = "RUNNING"; attr = g_attrOk; }
    }
    std::string left = " RFAM TOP  " + (live.runId.empty() ? std::string("(no run)") : live.runId);
    if (live.pid > 0) left += "   pid " + std::to_string(live.pid);
    left = padRight(elide(left, static_cast<std::size_t>(std::max(10, width - 34))), static_cast<std::size_t>(std::max(10, width - 34)));
    std::string right = "[" + state + "] " + formatDuration(live.uptime) + " ";
    right = padLeft(right, 33);
    put(0, 0, left, g_attrTitle | (g_color ? COLOR_PAIR(6) : 0));
    put(0, static_cast<int>(left.size()), right, attr | A_BOLD);
}

void drawCpuPanel(const Rect& rect, const Dashboard& board)
{
    drawBox(rect, "HARDWARE", g_attrTitle);
    const Rect inner = rect.inner();
    if (inner.empty()) return;

    HostState host;
    {
        std::lock_guard<std::mutex> lock(g_hostMutex);
        host = g_host;
    }

    int row = inner.y;
    const int cores = static_cast<int>(host.corePercent.size());
    char buffer[256];
    std::snprintf(buffer, sizeof(buffer), "CPU %s%%  %lld cores   load %s",
                  formatNumber(host.cpuPercent, 1).c_str(),
                  static_cast<long long>(board.live.cores > 0 ? board.live.cores : cores),
                  formatNumber(host.load1, 2).c_str());
    put(row++, inner.x, padRight(buffer, static_cast<std::size_t>(inner.w)), g_attrAccent);
    if (row >= inner.y + inner.h) return;

    // Core bars: as many per row as fit, at least one column, capped at 4.
    const int barWidth = 10;
    const int cellWidth = barWidth + 9;              // "c00 " + bar + " 100%"
    const int columns = std::max(1, std::min(4, inner.w / cellWidth));
    const int coreRows = (cores + columns - 1) / std::max(columns, 1);
    const int available = std::max(0, inner.y + inner.h - row - 4);  // reserve MEM/GPU/PROC
    for (int index = 0; index < cores; ++index)
    {
        const int column = index % columns;
        const int line = index / columns;
        if (line >= available || line >= coreRows) break;
        const double percent = host.corePercent[static_cast<std::size_t>(index)];
        char text[64];
        std::snprintf(text, sizeof(text), "c%02d %s %3.0f%%", index,
                      bar(percent / 100.0, barWidth, g_options.ascii).c_str(), percent);
        const int attr = percent > 90.0 ? g_attrBarHot : g_attrBar;
        put(row + line, inner.x + column * cellWidth, text, attr);
    }
    row += std::min(available, coreRows);
    if (row >= inner.y + inner.h) return;

    // Memory
    const double used = std::isfinite(host.memUsedGb) ? host.memUsedGb : board.live.memUsed;
    const double total = std::isfinite(host.memTotalGb) ? host.memTotalGb : board.live.memTotal;
    std::snprintf(buffer, sizeof(buffer), "MEM  %s %5.1f/%.1f GB",
                  bar(total > 0 ? used / total : 0.0, barWidth, g_options.ascii).c_str(), used, total);
    put(row++, inner.x, padRight(buffer, static_cast<std::size_t>(inner.w)), g_attrDim);
    if (row >= inner.y + inner.h) return;

    // GPU
    if (host.gpuAvailable)
        std::snprintf(buffer, sizeof(buffer), "GPU  %s %3.0f%% %5.1f/%.1f GB %.0fC",
                      bar(host.gpuUtil / 100.0, barWidth, g_options.ascii).c_str(),
                      host.gpuUtil, host.gpuMemUsed, host.gpuMemTotal, host.gpuTemp);
    else
        std::snprintf(buffer, sizeof(buffer), "GPU  -- %s", host.gpuError.empty() ? "not present" : host.gpuError.c_str());
    put(row++, inner.x, padRight(buffer, static_cast<std::size_t>(inner.w)), g_attrDim);
    if (row >= inner.y + inner.h) return;

    // Process
    std::snprintf(buffer, sizeof(buffer), "SIM  cpu %5.0f%%  rss %6.0f MB  threads %lld  state %c",
                  board.live.cpu, board.live.rss, board.live.threads,
                  board.live.procState.empty() ? '?' : board.live.procState[0]);
    put(row++, inner.x, padRight(buffer, static_cast<std::size_t>(inner.w)), g_attrDim);
}

void drawRunPanel(const Rect& rect, const Dashboard& board)
{
    drawBox(rect, "RUN", g_attrTitle);
    const Rect inner = rect.inner();
    if (inner.empty()) return;
    const Live& live = board.live;
    int row = inner.y;
    char buffer[512];

    const int width = inner.w;
    auto line = [&](const std::string& text, int attr) {
        if (row >= inner.y + inner.h) return;
        put(row++, inner.x, padRight(elide(text, static_cast<std::size_t>(width)), static_cast<std::size_t>(width)), attr);
    };

    if (!board.manifestLoaded)
    {
        line("manifest: not available", g_attrWarn);
    }
    else
    {
        const Foam::json::Value& manifest = board.manifest;
        line("binary  " + shortName(manifest["binary"]["path"].asString()), g_attrNormal);
        std::string commit = manifest["binary"]["git_commit"].asString();
        if (commit.empty()) commit = "?";
        line("commit  " + commit + "   compiled " + manifest["binary"]["compiled"].asString(), g_attrDim);
        const Foam::json::Value& effective = manifest["effective"];
        for (const char* key : {"thermo_backend", "mechanism", "turbulence", "inlet_mode"})
        {
            const Foam::json::Value* entry = effective.find(key);
            if (!entry) continue;
            const std::string value = (*entry)["value"].asString();
            const std::string source = (*entry)["source"].asString();
            char text[256];
            std::snprintf(text, sizeof(text), "%-14s %-18s %s", key, value.c_str(),
                          source == "default" ? "(default)" : "");
            line(text, source == "default" ? g_attrWarn : g_attrNormal);
        }
    }
    line("status age " + formatDuration(live.fileAge) + "   seq " + std::to_string(live.seq), g_attrDim);
    (void)buffer;
}

void drawPhasePanel(const Rect& rect, const Dashboard& board)
{
    drawBox(rect, "PHASE", g_attrTitle);
    const Rect inner = rect.inner();
    if (inner.empty()) return;
    const Live& live = board.live;
    int row = inner.y;
    const int width = inner.w;
    auto line = [&](const std::string& text, int attr) {
        if (row >= inner.y + inner.h) return;
        put(row++, inner.x, padRight(elide(text, static_cast<std::size_t>(width)), static_cast<std::size_t>(width)), attr);
    };

    line(live.phase + "   " + formatDuration(live.phaseAge), g_attrAccent | A_BOLD);
    if (!live.sub.empty())
        line("  \u2192 " + live.sub + "   " + formatDuration(live.subAge), g_attrNormal);
    const std::string file = live.subFile.empty() ? live.phaseFile : live.subFile;
    const long long fileLine = live.subFile.empty() ? live.phaseLine : live.subLine;
    if (!file.empty())
        line("  " + shortName(file) + ":" + std::to_string(fileLine), g_attrDim);

    if (live.total > 0)
    {
        const double fraction = static_cast<double>(live.done) / std::max<double>(1.0, static_cast<double>(live.total));
        char text[256];
        std::snprintf(text, sizeof(text), "%s %s/%s %s",
                      bar(fraction, 20, g_options.ascii).c_str(),
                      formatCount(static_cast<double>(live.done)).c_str(),
                      formatCount(static_cast<double>(live.total)).c_str(),
                      live.unit.c_str());
        line(text, g_attrNormal);
        std::snprintf(text, sizeof(text), "  %s/s   idle %s",
                      formatCount(live.rate).c_str(), formatDuration(live.idle).c_str());
        line(text, g_attrDim);
    }
    else
    {
        line("  progress not reported for this phase", g_attrDim);
        line("  idle " + formatDuration(live.idle), g_attrDim);
    }

    if (!live.stackNow.empty())
        line("  stack " + elide(live.stackNow, static_cast<std::size_t>(std::max(10, width - 8))), g_attrDim);
}

void drawStepPanel(const Rect& rect, const Dashboard& board)
{
    drawBox(rect, "STEP", g_attrTitle);
    const Rect inner = rect.inner();
    if (inner.empty()) return;
    const Live& live = board.live;
    int row = inner.y;
    char buffer[512];

    std::snprintf(buffer, sizeof(buffer), "step %lld / %lld   accepted %lld",
                  live.step, live.stepOf, live.accepted);
    put(row++, inner.x, padRight(buffer, static_cast<std::size_t>(inner.w)), g_attrAccent);

    auto metric = [&](const char* key) -> double {
        auto found = live.metrics.find(key);
        return found == live.metrics.end() ? NAN : found->second;
    };
    std::snprintf(buffer, sizeof(buffer), "res %s   ratio %s",
                  formatNumber(metric("residual")).c_str(),
                  formatNumber(metric("residual_ratio"), 2).c_str());
    put(row++, inner.x, padRight(buffer, static_cast<std::size_t>(inner.w)), g_attrNormal);

    std::snprintf(buffer, sizeof(buffer), "dt %s s   t_phys %s s",
                  formatNumber(metric("dt_physical")).c_str(),
                  formatNumber(metric("physical_time")).c_str());
    put(row++, inner.x, padRight(buffer, static_cast<std::size_t>(inner.w)), g_attrNormal);

    std::snprintf(buffer, sizeof(buffer), "p_ch %s MPa  Mach %s  T[%s, %s]",
                  formatNumber(metric("chamber_pressure") / 1e6, 4).c_str(),
                  formatNumber(metric("mach_max"), 2).c_str(),
                  formatNumber(metric("temperature_min"), 0).c_str(),
                  formatNumber(metric("temperature_max"), 0).c_str());
    put(row++, inner.x, padRight(buffer, static_cast<std::size_t>(inner.w)), g_attrNormal);

    std::snprintf(buffer, sizeof(buffer), "mdot in %s/%s/%s  out %s  dm %s",
                  formatNumber(metric("mdot_lox"), 3).c_str(),
                  formatNumber(metric("mdot_fuel"), 3).c_str(),
                  formatNumber(metric("mdot_torch"), 3).c_str(),
                  formatNumber(metric("mdot_out")).c_str(),
                  formatNumber(metric("mass_error")).c_str());
    put(row++, inner.x, padRight(buffer, static_cast<std::size_t>(inner.w)), g_attrNormal);
}

void drawLedgerPanel(const Rect& rect, const Dashboard& board)
{
    const Live& live = board.live;
    std::string title = "LEDGER";
    if (live.downgrades == 0 && live.approximations == 0 && live.failClosed == 0)
        title += "   no downgrade recorded";
    else
        title += "   downgrade " + std::to_string(live.downgrades)
               + "  approx " + std::to_string(live.approximations)
               + "  fail-closed " + std::to_string(live.failClosed);
    drawBox(rect, title, live.downgrades > 0 || live.failClosed > 0 ? g_attrBad : g_attrTitle);
    const Rect inner = rect.inner();
    if (inner.empty()) return;
    int row = inner.y;
    if (live.items.empty())
    {
        put(row, inner.x, padRight("nothing was approximated, substituted or skipped",
                                   static_cast<std::size_t>(inner.w)), g_attrOk);
        return;
    }
    for (const auto& item : live.items)
    {
        if (row >= inner.y + inner.h) break;
        const bool severe = item.severity == "downgrade" || item.severity == "fail_closed";
        const std::string marker = severe ? "\u26a0" : "~";
        const std::string text = " " + marker + " " + item.severity + "  " + item.code
                               + "  n=" + std::to_string(item.n) + "  " + item.detail;
        put(row++, inner.x, padRight(elide(text, static_cast<std::size_t>(inner.w)),
                                     static_cast<std::size_t>(inner.w)),
            severe ? g_attrBad : g_attrWarn);
    }
}

void drawEventsPanel(const Rect& rect, const Dashboard& board)
{
    drawBox(rect, "EVENTS (tail)", g_attrTitle);
    const Rect inner = rect.inner();
    if (inner.empty()) return;
    const std::size_t maximum = static_cast<std::size_t>(std::max(0, inner.h));
    const std::size_t count = std::min(maximum, board.events.size());
    int row = inner.y;
    for (std::size_t i = board.events.size() - count; i < board.events.size(); ++i)
    {
        const EventLine& event = board.events[i];
        int attr = g_attrDim;
        if (event.type == "DOWNGRADE") attr = g_attrBad;
        else if (event.type == "PHASE_BEGIN" || event.type == "PHASE_END") attr = g_attrNormal;
        else if (event.type == "WARNING" || event.type == "CHECKPOINT") attr = g_attrWarn;
        std::string time = event.time;
        const std::size_t t = time.find('T');
        if (t != std::string::npos) time = time.substr(t + 1, 8);
        const std::string text = " " + time + "  " + padRight(event.type, 12) + " " + event.text;
        put(row++, inner.x, padRight(elide(text, static_cast<std::size_t>(inner.w)),
                                     static_cast<std::size_t>(inner.w)), attr);
    }
}

void drawSparkPanel(const Rect& rect, const Dashboard& board)
{
    drawBox(rect, "RESIDUAL", g_attrTitle);
    const Rect inner = rect.inner();
    if (inner.empty()) return;
    const std::vector<double> series = residentSeries(board.steps, std::max(1, inner.w - 14));
    if (series.empty())
    {
        put(inner.y, inner.x, " no STEP events yet", g_attrDim);
        return;
    }
    const std::string spark = sparkline(series, std::max(1, inner.w - 14), g_options.ascii);
    const double last = series.back();
    char text[128];
    std::snprintf(text, sizeof(text), " %s  %s", spark.c_str(), formatNumber(last, 2).c_str());
    put(inner.y, inner.x, padRight(text, static_cast<std::size_t>(inner.w)), g_attrAccent);
    if (inner.h > 1)
    {
        char range[128];
        double low = series[0], high = series[0];
        for (double value : series) { low = std::min(low, value); high = std::max(high, value); }
        std::snprintf(range, sizeof(range), " min %s   max %s   n=%zu",
                      formatNumber(low, 2).c_str(), formatNumber(high, 2).c_str(), series.size());
        put(inner.y + 1, inner.x, padRight(range, static_cast<std::size_t>(inner.w)), g_attrDim);
    }
}

void drawKeys(int y, int width, const Dashboard& board)
{
    std::string text;
    if (board.overlay == Dashboard::Overlay::None)
        text = " q quit  space pause  m manifest  e events  s summary  h help  r refresh";
    else
        text = " up/down scroll  pgup/pgdn  esc close";
    put(y, 0, padRight(elide(text, static_cast<std::size_t>(width)), static_cast<std::size_t>(width)),
        g_attrDim | A_REVERSE);
}

void drawDashboard(const Dashboard& board)
{
    const int height = getmaxy(stdscr);
    const int width = getmaxx(stdscr);
    werase(stdscr);
    if (height < 20 || width < 76)
    {
        put(0, 0, "terminal too small: need at least 76x20", g_attrBad);
        wnoutrefresh(stdscr);
        doupdate();
        return;
    }

    drawHeader(width, board);

    const int top = 1;
    const int bottom = height - 2;
    const int bodyHeight = bottom - top + 1;
    const int leftWidth = std::max(34, std::min(60, width * 45 / 100));
    const int rightWidth = width - leftWidth + 1;

    const int ledgerHeight = std::max(3, std::min(7, static_cast<int>(board.live.items.size()) + 2));
    const int eventsHeight = std::max(3, std::min(9, bodyHeight / 3));
    const int sparkHeight = 4;
    int columnHeight = bodyHeight - ledgerHeight - eventsHeight - sparkHeight;
    if (columnHeight < 12) columnHeight = bodyHeight - ledgerHeight - 3;
    if (columnHeight < 8) columnHeight = bodyHeight;

    // Left column: hardware.
    Rect hardware{top, 0, columnHeight, leftWidth};
    drawCpuPanel(hardware, board);

    // Right column: run, phase, step.
    const int runHeight = 6;
    const int phaseHeight = 7;
    Rect run{top, leftWidth, runHeight, rightWidth};
    Rect phasePanel{top + runHeight, leftWidth, phaseHeight, rightWidth};
    Rect stepPanel{top + runHeight + phaseHeight, leftWidth,
                   std::max(0, columnHeight - runHeight - phaseHeight), rightWidth};
    drawRunPanel(run, board);
    drawPhasePanel(phasePanel, board);
    if (stepPanel.h >= 2) drawStepPanel(stepPanel, board);

    int cursor = top + columnHeight;
    Rect ledger{cursor, 0, std::min(ledgerHeight, bottom - cursor + 1), width};
    drawLedgerPanel(ledger, board);
    cursor += ledger.h;

    Rect spark{cursor, 0, std::min(sparkHeight, bottom - cursor + 1), width};
    drawSparkPanel(spark, board);
    cursor += spark.h;

    Rect events{cursor, 0, std::max(0, bottom - cursor + 1), width};
    if (events.h >= 2) drawEventsPanel(events, board);

    drawKeys(height - 1, width, board);
    wnoutrefresh(stdscr);
    doupdate();
}

// ===========================================================================
// Overlays
// ===========================================================================

void buildOverlayLines(Dashboard& board)
{
    board.overlayLines.clear();
    board.overlayScroll = 0;
    if (board.overlay == Dashboard::Overlay::Help)
    {
        board.overlayLines = {
            "RFAM TOP -- what this shows",
            "",
            "HARDWARE  host CPU per core, memory, GPU (nvidia-smi), and the solver's",
            "          own RSS/thread count, sampled independently of the solver.",
            "RUN       binary path and git commit, plus the effective values of the",
            "          settings that change the answer, and where each came from.",
            "PHASE     the function the solver is in right now, its source location,",
            "          how long it has been there, and sub-phase progress.",
            "STEP      residual, dt, chamber pressure, mass flows for the last step.",
            "RESIDUAL  history of residual ratios read from events.jsonl.",
            "LEDGER    every approximation or downgrade the solver recorded.",
            "EVENTS    tail of events.jsonl.",
            "",
            "Absent data renders as '--'. A stale status.json renders as TELEMETRY",
            "STALE with its true age. A dead PID renders as PROCESS GONE.",
        };
        return;
    }
    if (board.overlay == Dashboard::Overlay::Manifest)
    {
        std::string text;
        if (readWholeFile(board.runDir + "/manifest.json", text))
        {
            std::istringstream stream(text);
            std::string line;
            while (std::getline(stream, line)) board.overlayLines.push_back(line);
        }
        else board.overlayLines.push_back("manifest.json not available");
        return;
    }
    if (board.overlay == Dashboard::Overlay::Events)
    {
        for (const auto& event : board.events)
            board.overlayLines.push_back(event.time + "  " + padRight(event.type, 14) + " " + event.text);
        if (board.overlayLines.empty()) board.overlayLines.push_back("no events recorded");
        board.overlayScroll = static_cast<int>(board.overlayLines.size());
        return;
    }
    if (board.overlay == Dashboard::Overlay::Summary)
    {
        std::string text;
        if (readWholeFile(board.runDir + "/summary.json", text))
        {
            std::istringstream stream(text);
            std::string line;
            while (std::getline(stream, line)) board.overlayLines.push_back(line);
        }
        else board.overlayLines.push_back("summary.json not written yet (run still in progress)");
        return;
    }
}

void drawOverlay(const Dashboard& board)
{
    const int height = getmaxy(stdscr);
    const int width = getmaxx(stdscr);
    const char* title = "HELP";
    if (board.overlay == Dashboard::Overlay::Manifest) title = "MANIFEST (effective configuration)";
    else if (board.overlay == Dashboard::Overlay::Events) title = "EVENTS (full)";
    else if (board.overlay == Dashboard::Overlay::Summary) title = "SUMMARY";
    drawBox(Rect{0, 0, height - 1, width}, title, g_attrTitle);
    const int visible = height - 3;
    int start = board.overlayScroll;
    if (board.overlay == Dashboard::Overlay::Events)
        start = static_cast<int>(board.overlayLines.size()) - visible;
    start = std::max(0, std::min(start, static_cast<int>(board.overlayLines.size()) - visible));
    if (start < 0) start = 0;
    for (int row = 0; row < visible; ++row)
    {
        const std::size_t index = static_cast<std::size_t>(start + row);
        if (index >= board.overlayLines.size()) break;
        put(1 + row, 1, padRight(elide(board.overlayLines[index], static_cast<std::size_t>(width - 2)),
                                 static_cast<std::size_t>(width - 2)), g_attrNormal);
    }
    drawKeys(height - 1, width, board);
    wnoutrefresh(stdscr);
    doupdate();
}

// ===========================================================================
// Plain (non-interactive) mode
// ===========================================================================

int plainMode(const std::string& runDir, bool once)
{
    for (;;)
    {
        const Live live = loadLive(runDir);
        if (!live.valid)
        {
            std::printf("RFAM_TOP_ERROR %s\n", live.error.c_str());
            return 2;
        }
        std::printf(
            "run=%s pid=%lld seq=%lld state=%s status_age=%.1fs "
            "phase=%s phase_active=%d phase_age=%.1fs since_phase_end=%.1fs "
            "sub=%s sub_age=%.1fs progress=%lld/%lld unit=%s idle=%.1fs "
            "step=%lld/%lld res=%.6g ratio=%.6g "
            "downgrades=%lld approximations=%lld fail_closed=%lld "
            "hb_ticks=%lld hb_missed=%lld hb_jitter_ms=%.3f\n",
            live.runId.c_str(), live.pid, live.seq,
            live.processAlive ? (live.hasSummary ? "finished" : "running") : "gone",
            live.fileAge,
            live.phase.c_str(), live.phaseActive ? 1 : 0, live.phaseAge, live.phaseSinceEnd,
            live.sub.empty() ? "-" : live.sub.c_str(), live.subAge,
            live.done, live.total, live.unit.c_str(), live.idle,
            live.step, live.stepOf,
            live.metrics.count("residual") ? live.metrics.at("residual") : NAN,
            live.metrics.count("residual_ratio") ? live.metrics.at("residual_ratio") : NAN,
            live.downgrades, live.approximations, live.failClosed,
            live.hbTicks, live.hbMissed, live.hbJitter);
        std::fflush(stdout);
        // --once reports whether the status could be read, not whether the
        // process is still alive: a finished run is a perfectly good answer.
        if (once) return 0;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
}

// ===========================================================================
// Run discovery
// ===========================================================================

std::string discoverRun(const std::string& root)
{
    DIR* dir = ::opendir(root.c_str());
    if (!dir) return std::string();
    std::string best;
    double bestScore = -1e18;
    while (struct dirent* entry = ::readdir(dir))
    {
        if (entry->d_name[0] == '.') continue;
        const std::string candidate = root + "/" + entry->d_name;
        double age = 0;
        if (!fileAgeSeconds(candidate + "/status.json", age)) continue;
        // Prefer a run whose telemetry is still being written; a finished run
        // still wins over nothing at all, but never over a live one.
        double score = -age;
        if (age < 30.0) score += 1e6;
        if (score > bestScore) { bestScore = score; best = candidate; }
    }
    ::closedir(dir);
    return best;
}

int runInteractive(const std::string& runDir)
{
    Dashboard board;
    board.runDir = runDir;

    if (!initscr())
    {
        std::fprintf(stderr, "rfam_top: cannot initialise the terminal\n");
        return 2;
    }
    cbreak();
    noecho();
    nonl();
    keypad(stdscr, true);
    wtimeout(stdscr, 100);
    curs_set(0);
    if (has_colors())
    {
        start_color();
        assume_default_colors(-1, -1);
        init_pair(1, COLOR_CYAN, -1);
        init_pair(2, COLOR_GREEN, -1);
        init_pair(3, COLOR_YELLOW, -1);
        init_pair(4, COLOR_RED, -1);
        init_pair(5, COLOR_WHITE, -1);
        init_pair(6, COLOR_MAGENTA, -1);
        init_pair(7, COLOR_BLUE, -1);
        g_color = true;
    }
    g_attrNormal = g_color ? COLOR_PAIR(5) : 0;
    g_attrTitle = (g_color ? COLOR_PAIR(1) : 0) | A_BOLD;
    g_attrOk = g_color ? COLOR_PAIR(2) : 0;
    g_attrWarn = g_color ? COLOR_PAIR(3) : 0;
    g_attrBad = (g_color ? COLOR_PAIR(4) : 0) | A_BOLD;
    g_attrDim = A_DIM;
    g_attrAccent = g_color ? COLOR_PAIR(1) : 0;
    g_attrBar = g_color ? COLOR_PAIR(7) : 0;
    g_attrBarHot = g_color ? COLOR_PAIR(4) : 0;

    g_hostStop.store(false);
    std::thread hostThread(hostSamplerLoop);

    reloadData(board);
    board.lastDataLoad = nowSeconds();

    bool running = true;
    while (running && !g_quit.load())
    {
        const int key = wgetch(stdscr);
        if (key != ERR)
        {
            if (board.overlay != Dashboard::Overlay::None)
            {
                if (key == 27 || key == 'q' || key == 'Q') { board.overlay = Dashboard::Overlay::None; }
                else if (key == KEY_UP) { board.overlayScroll = std::max(0, board.overlayScroll - 1); }
                else if (key == KEY_DOWN) { board.overlayScroll += 1; }
                else if (key == KEY_PPAGE) { board.overlayScroll = std::max(0, board.overlayScroll - 20); }
                else if (key == KEY_NPAGE) { board.overlayScroll += 20; }
                else if (key == KEY_RESIZE) { /* ncurses handled it */ }
            }
            else
            {
                switch (key)
                {
                    case 'q': case 'Q': running = false; break;
                    case ' ': board.paused = !board.paused; break;
                    case 'r': reloadData(board); board.lastDataLoad = nowSeconds(); break;
                    case 'm': board.overlay = Dashboard::Overlay::Manifest; buildOverlayLines(board); break;
                    case 'e': board.overlay = Dashboard::Overlay::Events; buildOverlayLines(board); break;
                    case 's': board.overlay = Dashboard::Overlay::Summary; buildOverlayLines(board); break;
                    case 'h': case '?': board.overlay = Dashboard::Overlay::Help; buildOverlayLines(board); break;
                    case KEY_RESIZE: break;
                    default: break;
                }
            }
        }

        const double now = nowSeconds();
        if (!board.paused && now - board.lastDataLoad >= g_options.dataMs / 1000.0)
        {
            reloadData(board);
            board.lastDataLoad = now;
        }
        if (board.overlay != Dashboard::Overlay::None) drawOverlay(board);
        else drawDashboard(board);
        std::this_thread::sleep_for(std::chrono::milliseconds(g_options.refreshMs));
    }

    g_hostStop.store(true);
    if (hostThread.joinable()) hostThread.join();
    endwin();
    return 0;
}

void usage()
{
    std::printf(
        "usage: rfam_top [options]\n"
        "  --run DIR        telemetry run directory to watch\n"
        "  --root DIR       search here for the newest run (default runs/telemetry)\n"
        "  --refresh MS     screen refresh period (default 250)\n"
        "  --data MS        telemetry reload period (default 1000)\n"
        "  --ascii          plain ASCII box drawing and sparklines\n"
        "  --no-gpu         do not poll nvidia-smi\n"
        "  --plain          non-interactive: one status line per second\n"
        "  --once           non-interactive: print one status line and exit\n"
        "  --help           this text\n");
}

} // End anonymous namespace

int main(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        auto next = [&](const char* name) -> std::string {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "%s requires a value\n", name);
                std::exit(2);
            }
            return argv[++i];
        };
        if (arg == "--run") g_options.runDir = next("--run");
        else if (arg == "--root") g_options.rootDir = next("--root");
        else if (arg == "--refresh") g_options.refreshMs = std::atoi(next("--refresh").c_str());
        else if (arg == "--data") g_options.dataMs = std::atoi(next("--data").c_str());
        else if (arg == "--ascii") g_options.ascii = true;
        else if (arg == "--no-gpu") g_options.noGpu = true;
        else if (arg == "--plain") g_options.plain = true;
        else if (arg == "--once") { g_options.plain = true; g_options.once = true; }
        else if (arg == "--help" || arg == "-h") { usage(); return 0; }
        else { std::fprintf(stderr, "unknown option %s\n", arg.c_str()); usage(); return 2; }
    }

    if (g_options.refreshMs < 50) g_options.refreshMs = 50;
    if (g_options.dataMs < 100) g_options.dataMs = 100;

    const std::string envRoot = std::getenv("RFAM_TELEMETRY_DIR")
        ? std::string(std::getenv("RFAM_TELEMETRY_DIR")) : std::string();
    if (g_options.rootDir == "runs/telemetry" && !envRoot.empty()) g_options.rootDir = envRoot;

    std::string runDir = g_options.runDir;
    if (runDir.empty()) runDir = discoverRun(g_options.rootDir);
    if (runDir.empty())
    {
        std::fprintf(stderr, "rfam_top: no telemetry run found under %s\n", g_options.rootDir.c_str());
        return 2;
    }
    // Accept a root/latest symlink, a run directory, or a directory that holds
    // the telemetry files directly.
    if (::access((runDir + "/status.json").c_str(), R_OK) != 0)
    {
        if (::access((runDir + "/latest/status.json").c_str(), R_OK) == 0) runDir += "/latest";
    }

    if (g_options.plain || !::isatty(STDOUT_FILENO)) return plainMode(runDir, g_options.once);

    setlocale(LC_ALL, "");
    struct sigaction action;
    std::memset(&action, 0, sizeof(action));
    action.sa_handler = onSignal;
    sigemptyset(&action.sa_mask);
    ::sigaction(SIGINT, &action, nullptr);
    ::sigaction(SIGTERM, &action, nullptr);

    return runInteractive(runDir);
}
