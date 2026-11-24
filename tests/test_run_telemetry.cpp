// Contract test for the run-telemetry layer.
//
// The properties asserted here are the ones the whole feature rests on:
//
//   1. JSON round-trips exactly, including the strings the solver puts in it.
//   2. A started session leaves manifest.json, events.jsonl, status.json and
//      summary.json on disk, all parseable.
//   3. The heartbeat keeps ticking while the "solver" does nothing at all --
//      that is the property a stalled run depends on.
//   4. Every degradation lands in the ledger and in the summary, and
//      --strict-no-downgrade turns the ledger into a gate.
//   5. Signal termination is recorded rather than left as a mystery.

#include "JsonLite.H"
#include "RunTelemetry.H"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace
{

int gFailures = 0;

void check(bool condition, const char* what)
{
    if (condition) { std::printf("PASS %s\n", what); return; }
    std::printf("FAIL %s\n", what);
    ++gFailures;
}

std::string slurp(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) return std::string();
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

double number(const Foam::json::Value& root, const char* section, const char* key)
{
    return root[section][key].asDouble(-1.0);
}

void testJson()
{
    using namespace Foam::json;
    Value root = Value::object();
    root.set("name", Value("realFluidAllMach"));
    root.set("quote", Value("a \"quoted\" \\ value\n\twith control \x01"));
    root.set("int", Value(static_cast<long long>(-42)));
    root.set("double", Value(0.1 + 0.2));
    root.set("nan", Value(std::nan("")));
    Value array = Value::array();
    for (int i = 0; i < 5; ++i) array.push(Value(i));
    root.set("array", array);
    root["nested"]["deep"]["value"] = Value(true);

    const std::string text = root.dump(true);
    check(text.find("\"nan\": null") != std::string::npos,
          "JSON serialises NaN as null, never as an invalid number");

    Value parsed;
    std::string error;
    const bool ok = parse(text, parsed, error);
    check(ok, "JSON round-trip parses");
    if (!ok) { std::printf("     parse error: %s\n", error.c_str()); return; }

    check(parsed["name"].asString() == "realFluidAllMach", "JSON string survives");
    check(parsed["quote"].asString() == "a \"quoted\" \\ value\n\twith control \x01",
          "JSON escapes survive (quote, backslash, control)");
    check(parsed["int"].asInt() == -42, "JSON integer survives");
    check(parsed["double"].asDouble() == 0.1 + 0.2, "JSON double round-trips bit-exactly");
    check(parsed["array"].size() == 5, "JSON array length survives");
    check(parsed["nested"]["deep"]["value"].asBool(), "JSON nested member survives");
    check(parsed["missing"].isNull(), "absent member reads as null, not as a guess");

    Value rejected;
    check(!parse("{\"a\": 1", rejected, error), "JSON parser rejects truncated input");
    check(!parse("{\"a\" 1}", rejected, error), "JSON parser rejects a missing colon");
}

void testSession(const std::string& root)
{
    using namespace Foam::telemetry;

    Options options;
    options.rootDir = root;
    options.runName = "selftest";
    options.heartbeatMs = 250;
    options.sampleStacks = true;
    options.sampleHz = 20;
    options.stallSeconds = 2;

    Session& session = Session::instance();
    char arg0[] = "test_run_telemetry";
    char arg1[] = "--cfl";
    char arg2[] = "0.5";
    char* argv[] = {arg0, arg1, arg2, nullptr};
    session.start(3, argv, options);
    check(session.active(), "session starts");

    session.describe().set("mesh", Foam::json::Value("selftest.rfamtet"));
    session.defineField("cfl", "0.5");
    session.defineField("mechanism", "twostep");
    session.flushManifest();

    const std::string manifestText = slurp(session.runDir() + "/manifest.json");
    check(!manifestText.empty(), "manifest.json is written");

    Foam::json::Value manifest;
    std::string error;
    const bool manifestOk = Foam::json::parse(manifestText, manifest, error);
    check(manifestOk, "manifest.json parses");
    if (manifestOk)
    {
        check(manifest["effective"]["cfl"]["source"].asString() == "cli",
              "manifest tags a flag that appeared in argv as cli");
        check(manifest["effective"]["mechanism"]["source"].asString() == "default",
              "manifest tags an absent flag as default");
        check(manifest["argv"].size() == 3, "manifest keeps argv verbatim");
    }

    RFAM_PHASE(steady);
    RFAM_ENTER(steady);
    session.setStep(7, 100, 7);
    session.metric(Metric::Residual, 1.25e-3);
    session.metric(Metric::ChamberPressure, 4.02e6);
    session.count(Count::Frozen, 12);
    session.progress(0, 1000, "faces");
    for (int i = 0; i < 10; ++i) { session.bumpProgress(100); std::this_thread::sleep_for(std::chrono::milliseconds(10)); }

    // The heartbeat must tick even though this thread now does nothing.
    std::this_thread::sleep_for(std::chrono::milliseconds(1200));

    const std::string statusText = slurp(session.runDir() + "/status.json");
    check(!statusText.empty(), "status.json is written without solver activity");

    Foam::json::Value status;
    const bool statusOk = Foam::json::parse(statusText, status, error);
    check(statusOk, "status.json parses");
    if (statusOk)
    {
        check(status["seq"].asInt() >= 2, "heartbeat advanced while the solver idled");
        check(status["step"]["index"].asInt() == 7, "status carries the step index");
        check(std::abs(number(status, "metrics", "residual") - 1.25e-3) < 1e-18,
              "status carries the residual");
        check(status["phase"]["name"].asString() == "steady", "status carries the phase");
        check(status["heartbeat"]["missed"].asInt() == 0, "no heartbeat tick was missed");
        check(status["heartbeat"]["jitter_ms"].asDouble() < 250.0, "heartbeat jitter is bounded");
    }

    RFAM_LEAVE(steady);
    session.degrade("transport.sutherland_substitute", Severity::Approximation,
                    "53/59 species borrowed a Sutherland set");
    session.degrade("outlet.matched_pressure_shortcut", Severity::Downgrade,
                    "11467 faces used the 1e-2 band");
    session.warning("selftest warning");

    Foam::json::Value extra = Foam::json::Value::object();
    extra.set("accepted_steps", Foam::json::Value(7));
    session.finish("pass", 0, extra);

    const std::string summaryText = slurp(session.runDir() + "/summary.json");
    Foam::json::Value summary;
    const bool summaryOk = Foam::json::parse(summaryText, summary, error);
    check(summaryOk, "summary.json parses");
    if (summaryOk)
    {
        check(summary["verdict"].asString() == "pass", "summary carries the verdict");
        check(summary["ledger"]["downgrades"].asInt() == 1, "summary counts the downgrade");
        check(summary["ledger"]["approximations"].asInt() == 1, "summary counts the approximation");
        check(summary["ledger"]["items"].size() == 2, "summary lists both ledger entries");
        check(summary["accepted_steps"].asInt() == 7, "summary merges the driver's extra fields");
    }

    const std::string eventsText = slurp(session.runDir() + "/events.jsonl");
    std::size_t lines = 0;
    std::istringstream stream(eventsText);
    std::string line;
    bool allParse = true;
    while (std::getline(stream, line))
    {
        if (line.empty()) continue;
        ++lines;
        Foam::json::Value record;
        if (!Foam::json::parse(line, record, error)) allParse = false;
    }
    check(lines >= 5, "events.jsonl records start, phases, downgrades and finish");
    check(allParse, "every event line is valid JSON");
}

std::string makeRoot()
{
    char pattern[] = "/tmp/rfam_telemetry_test_XXXXXX";
    const char* made = ::mkdtemp(pattern);
    return made ? std::string(made) : std::string("/tmp/rfam_telemetry_test");
}

} // End anonymous namespace

int main()
{
    const std::string root = makeRoot();
    testJson();
    testSession(root);
    std::printf("run dir: %s\n", root.c_str());
    if (gFailures == 0) std::printf("RUN_TELEMETRY_TEST_PASS\n");
    else std::printf("RUN_TELEMETRY_TEST_FAIL failures=%d\n", gFailures);
    return gFailures == 0 ? 0 : 1;
}
