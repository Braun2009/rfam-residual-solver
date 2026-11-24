#ifndef RunTelemetry_H
#define RunTelemetry_H

#include "JsonLite.H"

#include <atomic>
#include <cstdint>
#include <string>

namespace Foam
{
namespace telemetry
{

// ===========================================================================
// Run telemetry: what the solver is actually doing, recorded continuously and
// independently of the stdout log.
//
// The design rule this file exists to enforce: a long run must never be able to
// fail, hang or silently downgrade without leaving a machine-readable trace.
// Three consequences drive every decision below.
//
//  1. The heartbeat thread never touches solver state that a lock protects, and
//     never allocates. It reads atomics and writes its own files. If the solver
//     is wedged inside a 40-second residual loop, the heartbeat still ticks --
//     and the tick shows the wedge.
//  2. There are two writers and they never share a file: the solver thread owns
//     `events.jsonl`, the heartbeat thread owns `status.json` and
//     `heartbeat.jsonl`. No cross-thread file lock can turn a stall into a
//     deadlock.
//  3. Every downgrade is recorded with a code and a severity, not just printed.
//     `--strict-no-downgrade` turns the ledger into a gate.
//
// Files written into one run directory:
//
//   manifest.json    effective configuration, every field tagged with the
//                    source that produced it (cli / env / default / resolved)
//   events.jsonl     solver-thread events, one JSON object per line
//   status.json      latest snapshot, replaced atomically at the heartbeat rate
//   heartbeat.jsonl  heartbeat-thread events (stall, io error, signal)
//   session.log      human-readable mirror of both streams
//   summary.json     final verdict, degradation ledger, counters
//   stacks.txt       sampled call-stack histogram (P2)
// ===========================================================================

enum class Severity { Requested, Approximation, Downgrade, FailClosed };

const char* severityName(Severity);

// Scalar quantities shown per step. Fixed enum so the writer and every reader
// agree on spelling without a schema negotiation.
enum class Metric
{
    Residual, ResidualRatio, DtPseudo, DtPhysical, PhysicalTime,
    ChamberPressure, MachMax,
    MdotLox, MdotFuel, MdotTorch, MdotOut,
    MassKg, MassError, CflLimit,
    TemperatureMin, TemperatureMax, PressureMin, PressureMax,
    Count_
};

const char* metricName(Metric);

// Monotonic counters. Absolute values: callers publish the current total.
enum class Count
{
    Frozen, Limited, Stiffened, TwoPhaseCells, Rejected, Retries,
    ThermoPt, ThermoIsentrope, ThermoRecover, ThermoRecoverFail,
    SpeciesRepair, PureFallback, BoundaryFallback, InterfaceFallback,
    OutletShortcutHits, OutletExactFaces,
    ChemAccepted, ChemRejected, ChemSkippedCold, ChemSkippedTwoPhase,
    Count_
};

const char* countName(Count);

// A named execution phase. Registration is idempotent; `__FILE__`/`__LINE__`
// give the reader the concrete source location without a debugger.
struct Phase
{
    int id = -1;
    const char* name = "";
    const char* file = "";
    int line = 0;
};

Phase phase(const char* name, const char* file, int line);

struct Options
{
    bool enabled = true;
    std::string rootDir;              // empty -> env or <cwd>/runs/telemetry
    std::string runName;              // empty -> timestamped id
    int heartbeatMs = 1000;
    // Call-stack capture rate, in captures per second, taken at the progress
    // sites (never from a signal handler). 0 disables the histogram.
    bool sampleStacks = true;
    int sampleHz = 5;
    int stallSeconds = 5;             // no progress for this long -> STALL event
    bool strictNoDowngrade = false;
    bool fsyncStatus = false;         // 9p / network filesystems: usually off
    bool mirrorSessionLog = true;
};

class Session
{
public:
    static Session& instance();

    // Called once from main() before anything else of interest happens.
    // Safe to call when disabled: every later call becomes a no-op.
    void start(int argc, char** argv, const Options& options);
    bool active() const { return active_; }
    const std::string& runDir() const { return runDir_; }
    const std::string& runId() const { return runId_; }

    // Fill the manifest once the driver knows the mesh, the thermo backend and
    // the resolved configuration. Rewrites manifest.json atomically.
    json::Value& describe();

    // Effective-configuration table. `flag` is the command-line spelling; pass
    // nullptr to derive "--key-with-dashes" from `key`. The source tag is
    // "cli" when that flag appears in argv, otherwise "default".
    void defineField(const std::string& key, const std::string& value,
                     const char* flag = nullptr,
                     const char* forceSource = nullptr);
    void flushManifest();

    // Hot path. All noexcept, none of them allocate.
    void enter(const Phase&);
    void leave(const Phase&);
    void subPhase(const char* name, const char* file, int line);
    void endSubPhase();
    void setStep(long long index, long long total, long long accepted);
    void metric(Metric, double value);
    void count(Count, unsigned long long value);
    void progress(unsigned long long done, unsigned long long total,
                  const char* unit);
    void bumpProgress(unsigned long long n);
    void noteProgress();

    // Structured events. Safe from any thread; serialised by one mutex that is
    // never taken by the heartbeat.
    void event(const char* type, const json::Value& fields);
    void degrade(const char* code, Severity, const std::string& detail);
    void warning(const std::string& text);
    void checkpoint(const std::string& path, const std::string& kind);

    // Final verdict. Idempotent; also runs from atexit and from the signal path.
    void finish(const char* verdict, int exitCode, const json::Value& extra);

    bool terminationRequested() const;
    int terminationSignal() const;

    // Sampler diagnostics (P2).
    json::Value samplerReport() const;

    ~Session();

private:
    Session();
    struct Impl;
    Impl* impl_ = nullptr;
    bool active_ = false;
    std::string runDir_, runId_;
};

// Free-function shorthands: the call sites read better without `.instance()`.
inline bool active() { return Session::instance().active(); }
inline void enter(const Phase& p) { Session::instance().enter(p); }
inline void leave(const Phase& p) { Session::instance().leave(p); }
inline void subPhase(const char* n, const char* f, int l) { Session::instance().subPhase(n, f, l); }
inline void endSubPhase() { Session::instance().endSubPhase(); }
inline void setStep(long long i, long long n, long long a) { Session::instance().setStep(i, n, a); }
inline void metric(Metric m, double v) { Session::instance().metric(m, v); }
inline void count(Count c, unsigned long long v) { Session::instance().count(c, v); }
inline void progress(unsigned long long d, unsigned long long t, const char* u)
{ Session::instance().progress(d, t, u); }
inline void bumpProgress(unsigned long long n) { Session::instance().bumpProgress(n); }
inline void noteProgress() { Session::instance().noteProgress(); }
inline void warning(const std::string& text) { Session::instance().warning(text); }
inline void degrade(const char* code, Severity s, const std::string& d)
{ Session::instance().degrade(code, s, d); }
inline void checkpoint(const std::string& path, const std::string& kind)
{ Session::instance().checkpoint(path, kind); }
inline bool terminationRequested() { return Session::instance().terminationRequested(); }

// Convenience macro so the file/line of a phase is never hand-typed.
#define RFAM_PHASE(name) \
    static const ::Foam::telemetry::Phase rfamPhase_##name = \
        ::Foam::telemetry::phase(#name, __FILE__, __LINE__)

#define RFAM_ENTER(name) ::Foam::telemetry::enter(rfamPhase_##name)
#define RFAM_LEAVE(name) ::Foam::telemetry::leave(rfamPhase_##name)

} // End namespace telemetry
} // End namespace Foam

#endif
