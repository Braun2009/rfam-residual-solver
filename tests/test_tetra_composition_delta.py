#!/usr/bin/env python3
"""Regression for the composition-delta audit contract."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from audit_tetra_composition_delta import summarize  # noqa: E402


def main() -> None:
    centroids = np.array(
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.001], [0.0, 0.0, 0.003]],
        dtype=float,
    )
    fields = np.zeros((3, 8), dtype=float)
    fields[:, 5] = 2.0
    reference = np.array([0.784, 0.784, 0.784], dtype=float)
    current = np.array([0.784, 0.7, 0.784], dtype=float)
    result = summarize(
        "synthetic", "synthetic.rfamtetf", centroids, fields, current,
        reference, np.ones(3, dtype=float), 1.0,
    )
    assert np.isclose(result["max_abs_delta"], 0.084)
    assert result["delta_ranges_by_threshold"]["0.01"]["count"] == 1
    assert np.isclose(result["delta_ranges_by_threshold"]["0.01"]["max_mm"], 1.0)
    assert result["z_bins"][1]["changed_gt_1e-5"] == 1
    print("TETRA_COMPOSITION_DELTA_PASS")


if __name__ == "__main__":
    main()
