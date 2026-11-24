#ifndef StackSampler_H
#define StackSampler_H

#include "JsonLite.H"

#include <string>

namespace Foam
{
namespace telemetry
{

// ---------------------------------------------------------------------------
// Call-site sampling.
//
// Two earlier designs were rejected here, both after measurement:
//
//   1. SIGPROF + backtrace() in the handler. It produced good-looking stacks and
//      it crashed the simulation: once under steady sampling of the main thread,
//      reliably within 75 s once it sampled OpenMP workers. glibc's unwinder is
//      not async-signal-safe in every context.
//   2. backtrace() from ordinary code at the progress sites. Safe, but on this
//      toolchain (gcc 9, -O2, OpenMP outlined regions) it returns the leaf frame,
//      the thread-start frame, and garbage in between -- even with
//      -fno-omit-frame-pointer. A histogram of garbage is worse than none.
//
// What survives is the part that cannot lie: the sampler records the LABEL of
// the instrumented loop that called it -- "faces", "cells_recovered",
// "outlet_faces". The label is a string literal chosen by the code that owns
// the loop, combined with the heartbeat's phase and sub-phase it reconstructs
// the chain that matters:
//
//     steady_step -> lusgs -> cells_recovered
//
// The result is a flat profile by call site: "this run spent 812 of 1000
// samples in the recovery loop". No signals, no unwinding, no possibility of
// corrupting the run, and no possibility of reporting a frame that is not real.
// ---------------------------------------------------------------------------

class StackSampler
{
public:
    void start(int hz, const std::string& runDir);
    void stop();

    // Called from the progress sites with their own __FILE__/__LINE__.
    void captureIfDue(const char* label) noexcept;

    bool running() const { return running_; }

    json::Value report() const;
    std::string histogramText() const;

private:
    struct Impl;
    Impl* impl_ = nullptr;
    bool running_ = false;
};

} // End namespace telemetry
} // End namespace Foam

#endif
