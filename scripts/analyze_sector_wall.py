#!/usr/bin/env python3
"""Independently audit an RFAM sector wall CSV."""

import argparse
import csv
import math
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path")
    parser.add_argument("--power-tolerance",type=float,default=1.0e-9,
                        help="relative tolerance for q*area versus power column")
    args=parser.parse_args()

    with open(args.csv_path,newline="",encoding="utf-8") as stream:
        rows=list(csv.DictReader(stream))
    if not rows:
        print("wall CSV has no data rows",file=sys.stderr)
        return 2
    required=("x_m","theta_rad","Tgas_K","qwall_W_m2","area_m2",
              "wall_power_W","y_plus")
    missing=[name for name in required if name not in rows[0]]
    if missing:
        print("wall CSV missing columns: "+", ".join(missing),file=sys.stderr)
        return 2

    numeric=[]
    maximum_power_error=0.0
    for row in rows:
        item={name:float(row[name]) for name in required}
        if not all(math.isfinite(value) for value in item.values()):
            print("wall CSV contains non-finite data",file=sys.stderr)
            return 1
        expected_power=item["qwall_W_m2"]*item["area_m2"]
        maximum_power_error=max(maximum_power_error,
            abs(item["wall_power_W"]-expected_power)
            /max(abs(expected_power),1.0e-300))
        numeric.append(item)

    total_area=sum(item["area_m2"] for item in numeric)
    total_power=sum(item["wall_power_W"] for item in numeric)
    if total_area<=0.0:
        print("wall CSV total area is not positive",file=sys.stderr)
        return 1
    positive_area=sum(item["area_m2"] for item in numeric if item["qwall_W_m2"]>0.0)
    negative_area=sum(item["area_m2"] for item in numeric if item["qwall_W_m2"]<0.0)
    maximum=max(numeric,key=lambda item:item["qwall_W_m2"])
    minimum=min(numeric,key=lambda item:item["qwall_W_m2"])

    print(f"rows={len(numeric)} area_m2={total_area:.12g}")
    print(f"power_W={total_power:.12g} qavg_W_m2={total_power/total_area:.12g}")
    print(f"positive_area_fraction={positive_area/total_area:.6g} "
          f"negative_area_fraction={negative_area/total_area:.6g}")
    print(f"qmax_W_m2={maximum['qwall_W_m2']:.12g} x_m={maximum['x_m']:.12g} "
          f"theta_rad={maximum['theta_rad']:.12g} y_plus={maximum['y_plus']:.12g}")
    print(f"qmin_W_m2={minimum['qwall_W_m2']:.12g} x_m={minimum['x_m']:.12g} "
          f"theta_rad={minimum['theta_rad']:.12g} y_plus={minimum['y_plus']:.12g}")
    print(f"Tgas_min/max_K={min(item['Tgas_K'] for item in numeric):.12g}/"
          f"{max(item['Tgas_K'] for item in numeric):.12g}")
    print(f"max_power_identity_relative_error={maximum_power_error:.6g}")
    if maximum_power_error>args.power_tolerance:
        print("SECTOR_WALL_AUDIT_FAIL")
        return 1
    print("SECTOR_WALL_AUDIT_PASS")
    return 0


if __name__=="__main__":
    sys.exit(main())
