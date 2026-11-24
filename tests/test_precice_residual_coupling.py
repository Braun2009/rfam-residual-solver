#!/usr/bin/env python3
"""Drive two real preCICE participants and audit the accepted windows."""
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


def run(binary, config, floor):
    outputs = {}
    with tempfile.TemporaryDirectory(prefix="rfam-precice-") as tmp:
        work = Path(tmp)
        env = dict(os.environ, OMP_NUM_THREADS="1")
        processes = []
        try:
            for role in ("Fluid", "Solid"):
                cmd = [str(binary), role, str(config)]
                if role == "Fluid":
                    cmd.append(str(floor))
                with (work / f"{role}.log").open("w") as log:
                    processes.append(subprocess.Popen(cmd, cwd=work, env=env,
                                                      stdout=log, stderr=subprocess.STDOUT))
            end = time.monotonic() + 35
            while True:
                states = [p.poll() for p in processes]
                if any(code not in (0, None) for code in states):
                    raise RuntimeError(f"participant failed: {states}")
                if all(code == 0 for code in states):
                    break
                if time.monotonic() > end:
                    raise RuntimeError("preCICE coupling timed out")
                time.sleep(0.05)
        finally:
            for p in processes:
                if p.poll() is None:
                    p.kill()
                p.wait()
            for role in ("Fluid", "Solid"):
                log = work / f"{role}.log"
                if log.exists():
                    outputs[role] = log.read_text(errors="replace")
                    print(outputs[role], end="", flush=True)

    data = {}
    for role, log in outputs.items():
        prefix = f"RFAM_PRECICE {role} "
        line, = [row for row in log.splitlines() if row.startswith(prefix)]
        values = dict(item.split("=", 1) for item in line[len(prefix):].split())
        data[role] = {key: float(value) for key, value in values.items()}
    return data


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def audit(data):
    fluid, solid = data["Fluid"], data["Solid"]
    for values in data.values():
        check(all(math.isfinite(v) for v in values.values()), "non-finite diagnostics")
        check(values["windows"] == 3 and values["restores"] >= 3, "missing accepted windows or rollbacks")
        check(abs(values["physical_time_s"] - 3e-6) < 1e-15, "rejected trial advanced physical time")
        check(values["max_balance_error_J"] < 1e-8, "participant energy imbalance")
    check(fluid["solver_accepted_steps"] == fluid["windows"],
          "discarded fluid trials leaked into the solver's accepted-step counter")
    check(fluid["heat_J"] > 0, "fluid did not heat solid")
    check(abs(fluid["heat_J"] - solid["heat_J"]) < 1e-10, "mapped interface power does not balance")
    check(0 < fluid["max_inner_defect"] <= 1e-12, "physical-step residual did not converge")
    check(fluid["max_inner_iterations"] > 0, "fluid did not use inner iterations")


if __name__ == "__main__":
    binary, config = (Path(arg).resolve() for arg in sys.argv[1:])
    on = run(binary, config, 0.2)
    off = run(binary, config, 1.0)
    audit(on)
    audit(off)
    check(abs(on["Fluid"]["heat_J"] - off["Fluid"]["heat_J"]) < 1e-9,
          "preconditioning changed the accepted heat result")
    print("RFAM_PRECICE_COUPLING_PASS")
