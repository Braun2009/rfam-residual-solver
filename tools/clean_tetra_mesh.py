#!/usr/bin/env python3
"""Improve tetrahedral mesh quality without touching the boundary surface.

Boundary vertices are held fixed, so the fluid surface, the inlet patches and
every boundary zone classification are preserved exactly. Interior vertices are
relaxed toward the average of their incident cell centroids; any relaxation that
would invert a cell is rejected. Remaining bad cells are repaired by collapsing
their shortest edge when every incident cell stays positive.

Quality is reported with the same mean-ratio metric the mesh audits use.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

TET_FACES = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]])


def cell_volumes(points, cells):
    x = points[cells]
    return np.einsum('ij,ij->i', x[:, 1] - x[:, 0],
                     np.cross(x[:, 2] - x[:, 0], x[:, 3] - x[:, 0])) / 6.0


def mean_ratio(points, cells):
    x = points[cells]
    volume = np.abs(cell_volumes(points, cells))
    edge2 = ((x[:, [0, 0, 0, 1, 1, 2]] - x[:, [1, 2, 3, 2, 3, 3]]) ** 2).sum(axis=(1, 2))
    return 12.0 * np.power(3.0 * volume, 2.0 / 3.0) / edge2


def boundary_edges(faces):
    e = np.sort(faces[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
    _, first, counts = np.unique(e, axis=0, return_index=True, return_counts=True)
    return e[first[counts == 1]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--iterations', type=int, default=30)
    parser.add_argument('--relax', type=float, default=0.4)
    parser.add_argument('--collapse-rounds', type=int, default=6)
    parser.add_argument('--quality-floor', type=float, default=1.0e-3,
                        help='mean-ratio below this is repaired')
    parser.add_argument('--report', default=None)
    args = parser.parse_args()

    data = np.load(args.input)
    points = np.array(data['points'], dtype=np.float64)
    cells = np.array(data['cells'], dtype=np.int64)
    boundary = np.array(data['boundary'], dtype=np.int64)
    zones = np.array(data['boundary_zone'], dtype=np.int32)
    boundary_owner = np.array(data['boundary_owner'], dtype=np.int64)
    inlet_patch = np.array(data['inlet_patch'], dtype=np.int32)

    report = {'input': args.input, 'input_cells': int(len(cells)),
              'input_min_volume_mm3': float(cell_volumes(points, cells).min())}
    q0 = mean_ratio(points, cells)
    report['input_quality'] = {'min': float(q0.min()),
                               'p001': float(np.percentile(q0, 0.1)),
                               'median': float(np.median(q0)),
                               'below_floor': int((q0 < args.quality_floor).sum())}

    # Boundary-adjacent vertices are fixed: the surface must not move or the
    # inlet patches and wall geometry would change.
    boundary_vertices = np.unique(boundary)
    interior = np.ones(len(points), dtype=bool)
    interior[boundary_vertices] = False
    report['interior_vertices'] = int(interior.sum())

    cells_of_vertex = [[] for _ in range(len(points))]
    for index, tet in enumerate(cells):
        for vertex in tet:
            cells_of_vertex[vertex].append(index)

    for iteration in range(args.iterations):
        centroids = points[cells].mean(axis=1)
        accumulated = np.zeros_like(points)
        counts = np.zeros(len(points))
        for vertex in np.flatnonzero(interior):
            incident = cells_of_vertex[vertex]
            if not incident:
                continue
            accumulated[vertex] = centroids[incident].sum(axis=0)
            counts[vertex] = len(incident)
        move = np.zeros_like(points)
        active = interior & (counts > 0)
        move[active] = args.relax * (accumulated[active] / counts[active, None]
                                     - points[active])
        trial = points + move
        volume = cell_volumes(trial, cells)
        bad = volume <= 0.0
        if bad.any():
            # Revert the vertices of every inverted cell, then re-check.
            revert = np.zeros(len(points), dtype=bool)
            revert[np.unique(cells[bad])] = True
            trial[revert] = points[revert]
            volume = cell_volumes(trial, cells)
            if (volume <= 0.0).any():
                # Fall back to a halved move.
                trial = points + 0.5 * move
                volume = cell_volumes(trial, cells)
                if (volume <= 0.0).any():
                    report['smoothing_stopped_at'] = iteration
                    break
        points = trial
    q1 = mean_ratio(points, cells)
    report['after_smoothing'] = {'min': float(q1.min()),
                                 'p001': float(np.percentile(q1, 0.1)),
                                 'median': float(np.median(q1)),
                                 'below_floor': int((q1 < args.quality_floor).sum())}

    # Edge collapse for the remaining bad cells (interior vertices only).
    edges = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
    resolved = 0
    for _ in range(args.collapse_rounds):
        q = mean_ratio(points, cells)
        bad = np.flatnonzero(q < args.quality_floor)
        if len(bad) == 0:
            break
        progress = 0
        for cell in bad:
            tet = cells[cell]
            lengths = []
            for a, b in edges:
                lengths.append((np.linalg.norm(points[tet[a]] - points[tet[b]]),
                                tet[a], tet[b]))
            lengths.sort()
            done = False
            for _, keep, remove in lengths:
                if not interior[keep] or not interior[remove]:
                    continue
                trial = points.copy()
                trial[keep] = 0.5 * (points[keep] + points[remove])
                candidate = cells.copy()
                candidate[candidate == remove] = keep
                key = np.sort(candidate, axis=1)
                if np.any(np.diff(key, axis=1) == 0):
                    continue  # would create a degenerate cell
                volume = cell_volumes(trial, candidate)
                if (volume <= 0.0).any():
                    continue
                points, cells = trial, candidate
                interior[remove] = False
                # Rebuild the vertex->cell map for the moved vertex.
                cells_of_vertex[keep] = [i for i, t in enumerate(cells) if keep in t]
                progress += 1
                resolved += 1
                done = True
                break
            if not done:
                continue
        if progress == 0:
            break
    q2 = mean_ratio(points, cells)
    volumes = cell_volumes(points, cells)
    report['after_collapse'] = {'collapsed_vertices': resolved,
                                'min': float(q2.min()),
                                'p001': float(np.percentile(q2, 0.1)),
                                'median': float(np.median(q2)),
                                'below_floor': int((q2 < args.quality_floor).sum()),
                                'min_volume': float(volumes.min()),
                                'negative_volumes': int((volumes <= 0).sum())}
    report['boundary_unchanged'] = True

    used = np.unique(cells)
    remap = np.full(len(points), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    points = points[used]
    cells = remap[cells]
    if np.any(cells < 0):
        raise SystemExit('cell referenced a removed vertex')
    # Boundary faces must reference surviving vertices (they were fixed).
    if np.any(remap[boundary] < 0):
        raise SystemExit('boundary face referenced a removed vertex')

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, points=points, cells=cells,
                        boundary=remap[boundary], boundary_zone=zones,
                        boundary_owner=boundary_owner, inlet_patch=inlet_patch)
    report['output'] = str(out)
    report['output_cells'] = int(len(cells))
    report['output_points'] = int(len(points))
    if args.report:
        Path(args.report).write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
