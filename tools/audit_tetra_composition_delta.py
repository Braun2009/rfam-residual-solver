#!/usr/bin/env python3
"""Compare species fields without inferring physical time or outlet flux.

The report is a diagnostic of cell-centred composition changes. A threshold
front is not a propagation speed, and axial slices are not injector exits.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from tetra_field_to_vtu import read_field, read_mesh


def parse_field(spec: str) -> tuple[str, str]:
    label, sep, path = spec.partition("=")
    if not sep or not label or not path:
        raise argparse.ArgumentTypeError("field must be LABEL=PATH")
    return label, path


def field_data(path: str, cell_count: int):
    centroids, fields, composition, species_count = read_field(path)
    if len(fields) != cell_count or len(composition) != cell_count:
        raise ValueError(f"field/mesh cell mismatch: {path}")
    if species_count < 2:
        raise ValueError(f"at least two species are required: {path}")
    y = composition[:, 0].astype(np.float64)
    if not np.isfinite(y).all() or np.any(y < -1e-6) or np.any(y > 1.000001):
        raise ValueError(f"invalid O2 mass fractions: {path}")
    return centroids.astype(np.float64), fields.astype(np.float64), y


def range_for(mask: np.ndarray, z_mm: np.ndarray):
    if not np.any(mask):
        return None
    return {
        "min_mm": float(z_mm[mask].min()),
        "max_mm": float(z_mm[mask].max()),
        "count": int(np.count_nonzero(mask)),
    }


def summarize(label: str, path: str, centroids: np.ndarray,
              fields: np.ndarray, y: np.ndarray, reference: np.ndarray,
              volumes: np.ndarray, bin_mm: float):
    delta = np.abs(y - reference)
    z_mm = centroids[:, 2] * 1e3
    max_cell = int(np.argmax(delta))
    thresholds = {}
    for threshold in (1e-7, 1e-5, 1e-3, 1e-2):
        thresholds[str(threshold)] = range_for(delta > threshold, z_mm)

    edges = np.arange(
        np.floor(float(z_mm.min()) / bin_mm) * bin_mm,
        np.ceil(float(z_mm.max()) / bin_mm) * bin_mm + bin_mm,
        bin_mm,
    )
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (z_mm >= lo) & (z_mm < hi)
        if not np.any(mask):
            continue
        bins.append({
            "z_min_mm": float(lo),
            "z_max_mm": float(hi),
            "cells": int(np.count_nonzero(mask)),
            "max_abs_delta": float(delta[mask].max()),
            "mean_abs_delta": float(delta[mask].mean()),
            "changed_gt_1e-5": int(np.count_nonzero(delta[mask] > 1e-5)),
        })

    rho = fields[:, 5]
    volume = volumes.astype(np.float64)
    total_mass = float(np.sum(rho * volume))
    return {
        "label": label,
        "source": path,
        "max_abs_delta": float(delta.max()),
        "max_delta_cell": max_cell,
        "max_delta_xyz_mm": (centroids[max_cell] * 1e3).tolist(),
        "reference_Y_O2_at_max": float(reference[max_cell]),
        "field_Y_O2_at_max": float(y[max_cell]),
        "mean_abs_delta": float(delta.mean()),
        "volume_weighted_abs_delta": float(np.sum(volume * delta) / np.sum(volume)),
        "mass_kg": total_mass,
        "O2_mass_kg": float(np.sum(rho * y * volume)),
        "Y_O2_min": float(y.min()),
        "Y_O2_max": float(y.max()),
        "cells_Y_O2_gt_0p9": int(np.count_nonzero(y > 0.9)),
        "furthest_Y_O2_gt_0p9_mm": (
            float(z_mm[y > 0.9].max()) if np.any(y > 0.9) else None
        ),
        "delta_ranges_by_threshold": thresholds,
        "z_bins": bins,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--reference", required=True,
                        help="reference field used for composition deltas")
    parser.add_argument("--field", action="append", required=True,
                        type=parse_field, help="LABEL=FIELD.rfamtetf")
    parser.add_argument("--bin-mm", type=float, default=1.0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    if not np.isfinite(args.bin_mm) or args.bin_mm <= 0:
        parser.error("--bin-mm must be positive")

    _, cells, volumes = read_mesh(args.mesh)
    reference_centroids, _, reference_y = field_data(args.reference, len(cells))
    report = {
        "scope": "cell-centred composition delta only; not physical time, outlet flux or ignition evidence",
        "mesh": args.mesh,
        "reference": args.reference,
        "bin_mm": args.bin_mm,
        "fields": [],
    }
    for label, path in args.field:
        centroids, fields, y = field_data(path, len(cells))
        if not np.allclose(centroids, reference_centroids, rtol=0.0, atol=2e-7):
            raise ValueError(f"field centroid mismatch against reference: {path}")
        report["fields"].append(
            summarize(label, path, centroids, fields, y, reference_y, volumes, args.bin_mm)
        )

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for item in report["fields"]:
        print(item["label"],
              "max_abs_delta=", item["max_abs_delta"],
              "changed_gt_1e-5=", item["delta_ranges_by_threshold"]["1e-05"],
              "max_z_mm=", item["delta_ranges_by_threshold"]["1e-05"]["max_mm"]
              if item["delta_ranges_by_threshold"]["1e-05"] else None)
    print("COMPOSITION_DELTA_AUDIT_WRITTEN", output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
