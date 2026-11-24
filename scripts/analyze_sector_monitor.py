#!/usr/bin/env python3
"""Audit time convergence from rfam_v20_sector_smoke MONITOR lines."""

import argparse
import math
import re
import sys


MONITOR_RE = re.compile(
    r"MONITOR iter=(?P<iter>\d+) t=(?P<time>\S+) mass=(?P<mass>\S+) "
    r"mdot_in/out=(?P<mdot_in>[^/]+)/(?P<mdot_out>\S+) "
    r"(?:mdot_axis/wall=(?P<axis_out>[^/]+)/(?P<wall_out>\S+) )?"
    r"(?:mdot_patch=(?P<background_in>[^/]+)/(?P<torch_in>[^/]+)/"
    r"(?P<lox_in>[^/]+)/(?P<fuel_in>\S+) )?"
    r"T=(?P<tmin>[^/]+)/(?P<tmax>\S+) Pc=(?P<pc>\S+) "
    r"throat=(?P<throat>\S+) cstar=(?P<cstar>\S+) "
    r"qmax/qavg=(?P<qmax>[^/]+)/(?P<qavg>\S+)"
    r"(?: mdot_target=(?P<mdot_target>\S+) "
    r"throat/target=(?P<throat_target_ratio>\S+))?"
)

AUDIT_METRICS = ("mass", "throat", "cstar", "qmax", "qavg", "tmax")


def parse_monitors(lines):
    rows = []
    for line in lines:
        match = MONITOR_RE.search(line)
        if not match:
            continue
        row = {"iter": int(match.group("iter"))}
        row.update(
            (name, float(value))
            for name, value in match.groupdict().items()
            if name != "iter" and value is not None
        )
        rows.append(row)
    return rows


def linear_slope(rows, metric):
    x_mean = sum(row["iter"] for row in rows) / len(rows)
    y_mean = sum(row[metric] for row in rows) / len(rows)
    denominator = sum((row["iter"] - x_mean) ** 2 for row in rows)
    if denominator == 0.0:
        return 0.0
    return sum(
        (row["iter"] - x_mean) * (row[metric] - y_mean) for row in rows
    ) / denominator


def audit(rows, span_limit, trend_limit):
    results = {}
    for metric in AUDIT_METRICS:
        first = rows[0][metric]
        last = rows[-1][metric]
        scale = max(abs(last), 1.0e-300)
        slope = linear_slope(rows, metric)
        span_relative = abs(last - first) / scale
        trend_relative = abs(slope * 10000.0) / scale
        results[metric] = {
            "first": first,
            "last": last,
            "span_relative": span_relative,
            "slope_per_1000": slope * 1000.0,
            "trend_relative_per_10000": trend_relative,
            "pass": span_relative <= span_limit and trend_relative <= trend_limit,
        }
    return results


def mass_closure(rows):
    required=("axis_out","wall_out")
    if any(any(name not in row for name in required) for row in rows):
        return None
    intervals=[]
    for first,last in zip(rows,rows[1:]):
        delta_time=last["time"]-first["time"]
        if delta_time<=0.0:
            continue
        inventory_rate=(last["mass"]-first["mass"])/delta_time
        first_net=(first["mdot_in"]-first["mdot_out"]
                   -first["axis_out"]-first["wall_out"])
        last_net=(last["mdot_in"]-last["mdot_out"]
                  -last["axis_out"]-last["wall_out"])
        boundary_rate=0.5*(first_net+last_net)
        flow_scale=max(0.5*(abs(first["mdot_in"])+abs(last["mdot_in"])),
                       0.5*(abs(first["mdot_out"])+abs(last["mdot_out"])),
                       1.0e-300)
        intervals.append({
            "first_iter":first["iter"],
            "last_iter":last["iter"],
            "inventory_rate":inventory_rate,
            "boundary_rate":boundary_rate,
            "relative":abs(inventory_rate-boundary_rate)/flow_scale,
        })
    if not intervals:
        return None
    worst=max(intervals,key=lambda interval:interval["relative"])
    return {
        "maximum":worst["relative"],
        "rms":math.sqrt(sum(interval["relative"]**2 for interval in intervals)
                         /len(intervals)),
        "intervals":len(intervals),
        "worst":worst,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", help="solver log path, or - for standard input")
    parser.add_argument("--window", type=int, default=20, help="latest monitor points")
    parser.add_argument("--span-limit", type=float, default=0.01,
                        help="maximum endpoint relative change")
    parser.add_argument("--trend-limit", type=float, default=0.005,
                        help="maximum fitted relative drift per 10000 iterations")
    parser.add_argument("--wall-flow-limit", type=float, default=1.0e-10,
                        help="maximum absolute outward wall mass flow in kg/s")
    parser.add_argument("--background-flow-limit", type=float, default=1.0e-10,
                        help="maximum absolute non-injector headwall mass flow in kg/s")
    parser.add_argument("--injection-flow-limit", type=float, default=0.02,
                        help="maximum injector numerical/target mass-flow mismatch")
    parser.add_argument("--torch-flow", type=float, default=0.1348 / 3.0,
                        help="full-ramp prescribed torch mass flow in kg/s")
    parser.add_argument("--lox-flow", type=float, default=0.9333,
                        help="full-ramp prescribed LOX mass flow in kg/s")
    parser.add_argument("--fuel-flow", type=float, default=0.3500,
                        help="full-ramp prescribed fuel mass flow in kg/s")
    parser.add_argument("--ramp-initial", type=float, default=0.1,
                        help="initial fraction of prescribed injector flow")
    parser.add_argument("--ramp-duration", type=float, default=2.0e-4,
                        help="linear injector ramp duration in seconds")
    parser.add_argument("--mass-closure-limit", type=float, default=0.02,
                        help="maximum normalized inventory/flux residual")
    args = parser.parse_args()

    if args.window < 2:
        parser.error("--window must be at least 2")
    if args.ramp_duration <= 0.0:
        parser.error("--ramp-duration must be positive")
    stream = sys.stdin if args.log == "-" else open(args.log, encoding="utf-8")
    try:
        rows = parse_monitors(stream)
    finally:
        if stream is not sys.stdin:
            stream.close()
    if len(rows) < args.window:
        print(f"insufficient monitor points: found {len(rows)}, need {args.window}",
              file=sys.stderr)
        return 2

    rows = rows[-args.window:]
    results = audit(rows, args.span_limit, args.trend_limit)
    print(f"window iter={rows[0]['iter']}..{rows[-1]['iter']} points={len(rows)}")
    print("metric          first          last    span[%]   slope/1k  trend/10k[%] status")
    for metric in AUDIT_METRICS:
        result = results[metric]
        status = "PASS" if result["pass"] else "DRIFT"
        print(f"{metric:<8} {result['first']:>13.6g} {result['last']:>13.6g} "
              f"{100.0 * result['span_relative']:>10.4f} "
              f"{result['slope_per_1000']:>11.6g} "
              f"{100.0 * result['trend_relative_per_10000']:>14.4f} {status}")
    wall_rows = [row for row in rows if "wall_out" in row]
    wall_pass = False
    if len(wall_rows) == len(rows):
        maximum_wall_flow=max(abs(row["wall_out"]) for row in wall_rows)
        wall_pass=maximum_wall_flow<=args.wall_flow_limit
        print(f"max_abs_wall_mdot={maximum_wall_flow:.6g} "
              f"{'PASS' if wall_pass else 'LEAK'}")
    else:
        print("max_abs_wall_mdot=unavailable")
    background_rows=[row for row in rows if "background_in" in row]
    background_pass=False
    if len(background_rows)==len(rows):
        maximum_background_flow=max(abs(row["background_in"]) for row in background_rows)
        background_pass=maximum_background_flow<=args.background_flow_limit
        print(f"max_abs_background_mdot={maximum_background_flow:.6g} "
              f"{'PASS' if background_pass else 'LEAK'}")
    else:
        print("max_abs_background_mdot=unavailable")
    injection_rows=[row for row in rows
                    if all(name in row for name in
                           ("torch_in","lox_in","fuel_in","mdot_target"))]
    injection_pass=False
    if len(injection_rows)==len(rows):
        mismatches=[]
        patch_mismatches={"torch_in":[],"lox_in":[],"fuel_in":[]}
        prescribed={"torch_in":args.torch_flow,"lox_in":args.lox_flow,
                    "fuel_in":args.fuel_flow}
        for row in injection_rows:
            ramp=args.ramp_initial+(1.0-args.ramp_initial)*min(
                row["time"]/args.ramp_duration,1.0)
            actual=row["torch_in"]+row["lox_in"]+row["fuel_in"]
            ramped_target=ramp*row["mdot_target"]
            mismatches.append(abs(actual-ramped_target)
                              /max(abs(ramped_target),1.0e-300))
            for name,target in prescribed.items():
                patch_target=ramp*target
                patch_mismatches[name].append(
                    abs(row[name]-patch_target)/max(abs(patch_target),1.0e-300))
        maximum_patch_mismatch=max(max(values)
                                   for values in patch_mismatches.values())
        injection_pass=(max(mismatches)<=args.injection_flow_limit
                        and maximum_patch_mismatch<=args.injection_flow_limit)
        latest=injection_rows[-1]
        latest_ramp=args.ramp_initial+(1.0-args.ramp_initial)*min(
            latest["time"]/args.ramp_duration,1.0)
        latest_actual=latest["torch_in"]+latest["lox_in"]+latest["fuel_in"]
        print(f"injector_mdot_actual/target={latest_actual:.6g}/"
              f"{latest_ramp*latest['mdot_target']:.6g} "
              f"max_total_mismatch={100.0*max(mismatches):.4f}% "
              f"max_patch_mismatch={100.0*maximum_patch_mismatch:.4f}% "
              f"{'PASS' if injection_pass else 'DRIFT'}")
    else:
        print("injector_mdot_actual/target=unavailable")
    closure=mass_closure(rows)
    closure_pass=False
    if closure is None:
        print("mass_closure=unavailable")
    else:
        closure_pass=closure["maximum"]<=args.mass_closure_limit
        print(f"mass_closure_intervals={closure['intervals']} "
              f"max={100.0*closure['maximum']:.4f}% "
              f"rms={100.0*closure['rms']:.4f}% "
              f"{'PASS' if closure_pass else 'FAIL'}")
        worst=closure["worst"]
        print(f"mass_closure_worst={worst['first_iter']}..{worst['last_iter']} "
              f"inventory_rate={worst['inventory_rate']:.6g} "
              f"boundary_rate={worst['boundary_rate']:.6g} kg/s")
    converged = (all(result["pass"] for result in results.values())
                 and wall_pass and background_pass and injection_pass and closure_pass)
    print("SECTOR_TIME_CONVERGED" if converged else "SECTOR_TIME_NOT_CONVERGED")
    return 0 if converged else 1


if __name__ == "__main__":
    sys.exit(main())
