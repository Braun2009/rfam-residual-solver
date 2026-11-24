#!/usr/bin/env python3
"""Restore CAD-derived inlet contours by conforming tetrahedral edge operations.

Project rim vertices onto reference polylines, optionally contract an explicitly
identified folded rim edge, and bisect every incident tetrahedron along long
rim edges. New rim nodes lie on the reference polyline at half arclength.
Neither isolated bad cells nor boundary triangles are deleted to hide failures.
Every surviving tetrahedron must keep positive volume; full audit is separate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment

TET_FACES = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]])


def perimeter(faces):
    edges, count = np.unique(np.sort(faces[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1), axis=0, return_counts=True)
    if np.any(count > 2):
        raise ValueError('Nonmanifold inlet cap')
    return edges[count == 1]


def cap_area(points, faces):
    x = points[faces]
    return float(.5 * np.linalg.norm(np.cross(x[:, 1]-x[:, 0], x[:, 2]-x[:, 0]), axis=1).sum())


def volumes(points, cells):
    x = points[cells]
    return np.einsum('ij,ij->i', x[:, 1]-x[:, 0], np.cross(x[:, 2]-x[:, 0], x[:, 3]-x[:, 0])) / 6


def polyline_loops(points, edges):
    neighbours = {}
    for u, v in edges:
        neighbours.setdefault(int(u), []).append(int(v))
        neighbours.setdefault(int(v), []).append(int(u))
    if any(len(n) != 2 for n in neighbours.values()):
        raise ValueError('Reference perimeter is not a set of closed loops')
    unused = set(neighbours)
    loops = []
    while unused:
        first = min(unused); path = [first]; previous = None; current = first
        while True:
            nxt = next(v for v in neighbours[current] if v != previous)
            if nxt == first:
                break
            if nxt in path:
                raise ValueError('Self-connected reference perimeter')
            path.append(nxt); previous, current = current, nxt
        unused.difference_update(path)
        xy = points[path, :2]
        length = np.linalg.norm(np.roll(xy, -1, axis=0)-xy, axis=1)
        if np.any(length <= 0):
            raise ValueError('Zero reference perimeter segment')
        loops.append((xy, np.r_[0., np.cumsum(length)]))
    return loops


def on_loop(loop, distance):
    xy, arc = loop; distance %= arc[-1]
    j = min(int(np.searchsorted(arc, distance, side='right')-1), len(xy)-1)
    t = (distance-arc[j])/(arc[j+1]-arc[j])
    return (1-t)*xy[j] + t*xy[(j+1) % len(xy)]


def project(loop, xy):
    vertices, arc = loop; end = np.roll(vertices, -1, axis=0); vec = end-vertices
    t = np.clip(np.einsum('ij,ij->i', xy-vertices, vec)/(vec*vec).sum(1), 0, 1)
    closest = vertices+t[:, None]*vec
    sq = ((closest-xy)**2).sum(1); j = int(np.argmin(sq))
    return closest[j], float(arc[j]+t[j]*(arc[j+1]-arc[j])), float(np.sqrt(sq[j]))


class LocalElements:
    """Mutable incident elements; the millions of untouched cells stay in NumPy."""
    def __init__(self, elements):
        self.elements = {}; self.adjacency = {}; self.next_id = 0
        for el in elements:
            self.add(el)

    def add(self, el):
        idx = self.next_id; self.next_id += 1; el = tuple(int(x) for x in el)
        self.elements[idx] = el
        for v in el:
            self.adjacency.setdefault(v, set()).add(idx)

    def split_edge(self, u, v, midpoint):
        ids = list(self.adjacency.get(u, set()) & self.adjacency.get(v, set()))
        for idx in ids:
            el = self.elements.pop(idx)
            for n in el:
                self.adjacency[n].remove(idx)
            self.add([midpoint if n == v else n for n in el])
            self.add([midpoint if n == u else n for n in el])
        return len(ids)

    def incident(self, u, v):
        return [self.elements[i] for i in self.adjacency.get(u, set()) & self.adjacency.get(v, set())]

    def array(self, width):
        return np.array(list(self.elements.values()), dtype=np.int64).reshape(-1, width)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True); parser.add_argument('--reference', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--max-movement-mm', type=float, default=.04)
    parser.add_argument('--edge-mm', type=float, default=.1)
    parser.add_argument('--area-relative-tolerance', type=float, default=.003)
    parser.add_argument('--collapse-edge', type=int, nargs=2, metavar=('REMOVE', 'KEEP'))
    parser.add_argument('--expected-input-sha256', help='Guard node-index repair against a different input mesh')
    args = parser.parse_args()
    if not 0 < args.edge_mm <= .3:
        raise ValueError('Inlet edge target must be in (0, 0.3] mm')
    root = Path(args.mesh); ref = Path(args.reference)
    input_hash = hashlib.sha256((root/'fluid_volume_mm.npz').read_bytes()).hexdigest()
    if args.expected_input_sha256 and input_hash != args.expected_input_sha256:
        raise ValueError('Input mesh hash differs from the node-index repair record')
    m = np.load(root/'fluid_volume_mm.npz'); r = np.load(ref/'fluid_volume_mm.npz')
    points = m['points'].copy(); cells = m['cells'].copy()
    boundary = m['boundary'].copy(); labels = m['inlet_patch'].copy()
    input_volume = float(volumes(points, cells).sum())
    old = json.loads((root/'mesh_audit.json').read_text())['inlets']
    new = json.loads((ref/'mesh_audit.json').read_text())['inlets']
    dist = np.linalg.norm(np.array([g['centroid_cad_mm'] for g in old])[:, None, :]-np.array([g['centroid_cad_mm'] for g in new])[None, :, :], axis=2)
    rows, cols = linear_sum_assignment(dist)
    if len(rows) != 64 or np.any(dist[rows, cols] > .05):
        raise ValueError('Inlet matching failed')
    loops = []; rim_location = {}; projection = []; groups = []
    for i, j in zip(rows, cols):
        faces = boundary[labels == old[i]['id']]
        ref_faces = r['boundary'][r['inlet_patch'] == new[j]['id']]
        ids = np.unique(perimeter(faces)); loop_ids = []
        for loop in polyline_loops(r['points'], perimeter(ref_faces)):
            loop_ids.append(len(loops)); loops.append(loop)
        motions = []
        for vertex in ids:
            candidates = [(k, *project(loops[k], points[vertex, :2])) for k in loop_ids]
            k, xy, distance, motion = min(candidates, key=lambda x: x[-1])
            if motion > args.max_movement_mm:
                raise ValueError(f'Inlet {old[i]["id"]}: motion {motion} exceeds {args.max_movement_mm} mm')
            points[vertex, :2] = xy; rim_location[int(vertex)] = (k, distance); motions.append(motion)
        projection.append({'id': old[i]['id'], 'max_node_movement_mm': max(motions)})
        groups.append((old[i], float(cap_area(m['points'], faces)), float(cap_area(r['points'], ref_faces))))
    collapse = None
    if args.collapse_edge:
        remove, keep = args.collapse_edge
        if remove not in rim_location or keep not in rim_location or rim_location[remove][0] != rim_location[keep][0]:
            raise ValueError('Only an edge on one matched inlet perimeter may be contracted')
        if not np.any(np.all(perimeter(boundary[labels > 0]) == sorted([remove, keep]), axis=1)):
            raise ValueError('Specified contraction is not a perimeter edge')
        length = float(np.linalg.norm(m['points'][remove]-m['points'][keep]))
        if length > .35:
            raise ValueError('Folded rim edge is larger than a nominal inlet surface cell')
        cells[cells == remove] = keep; boundary[boundary == remove] = keep
        collapsed_cells = np.any(np.diff(np.sort(cells, axis=1), axis=1) == 0, axis=1)
        collapsed_faces = np.any(np.diff(np.sort(boundary, axis=1), axis=1) == 0, axis=1)
        cells = cells[~collapsed_cells]; boundary = boundary[~collapsed_faces]; labels = labels[~collapsed_faces]
        collapse = {'remove_node': remove, 'keep_node': keep, 'input_edge_length_mm': length,
                    'collapsed_degenerate_cells': int(collapsed_cells.sum()), 'collapsed_degenerate_faces': int(collapsed_faces.sum())}
        del rim_location[remove]
    rim = np.array(list(rim_location))
    incident = np.any(np.isin(cells, rim), axis=1)
    if np.any(volumes(points, cells[incident]) <= 0):
        raise ValueError('Projection/contraction would invert a surviving incident tetrahedron')
    changed_faces = np.any(np.isin(boundary, rim), axis=1)
    local_cells = LocalElements(cells[incident]); local_faces = LocalElements(boundary[changed_faces])
    static_cells = cells[~incident]; static_faces = boundary[~changed_faces]
    cap_edges = perimeter(boundary[labels > 0]); queue = [tuple(int(x) for x in e) for e in cap_edges]
    added = []; npoints = len(points); splits = 0
    def coordinate(v):
        return points[v] if v < npoints else added[v-npoints]
    while queue:
        u, v = queue.pop(); ku, su = rim_location[u]; kv, sv = rim_location[v]
        if ku != kv:
            raise ValueError('Inlet edge connects different reference loops')
        total = loops[ku][1][-1]; delta = (sv-su+total/2) % total-total/2
        if abs(delta) <= args.edge_mm:
            continue
        new_arc = (su+delta/2) % total
        xy = on_loop(loops[ku], new_arc)
        midpoint = np.array([*xy, coordinate(u)[2]])
        star = local_cells.incident(u, v)
        if not star:
            raise ValueError('Perimeter edge lost all incident cells')
        new_id = npoints+len(added)
        # Validate both children against the unmodified parent ordering.
        for tet in star:
            child = np.array([coordinate(n) for n in tet])
            for replace in [u, v]:
                x = child.copy(); x[tet.index(replace)] = midpoint
                signed = np.dot(x[1]-x[0], np.cross(x[2]-x[0], x[3]-x[0]))/6
                if signed <= 0:
                    raise ValueError(f'CAD arc edge split {u},{v} would invert an incident child')
        added.append(midpoint); rim_location[new_id] = (ku, new_arc)
        local_cells.split_edge(u, v, new_id); local_faces.split_edge(u, v, new_id)
        queue.extend([(u, new_id), (new_id, v)]); splits += 1
        if splits > 20000:
            raise ValueError('Unexpected inlet-refinement growth')
    points = np.vstack([points, np.array(added).reshape(-1, 3)])
    cells = np.vstack([static_cells, local_cells.array(4)])
    boundary = np.vstack([static_faces, local_faces.array(3)])
    if np.any(volumes(points, cells) <= 0):
        raise ValueError('Repaired mesh has inverted tetrahedra')
    # Recompute cap connectivity/areas, then match by centroid for acceptance.
    inlet_z = float(old[0]['centroid_cad_mm'][2])
    cap = boundary[np.max(np.abs(points[boundary, 2]-inlet_z), axis=1) < 1e-7]
    import trimesh
    used = np.unique(cap); ids = np.full(len(points), -1); ids[used] = np.arange(len(used))
    caps = trimesh.Trimesh(points[used], ids[cap], process=False).split(only_watertight=False)
    if len(caps) != 64:
        raise ValueError(f'Inlet cap connectivity changed to {len(caps)}')
    cent = np.array([x.centroid for x in caps])
    cost = np.linalg.norm(cent[:, None, :]-np.array([x[0]['centroid_cad_mm'] for x in groups])[None, :, :], axis=2)
    ii, jj = linear_sum_assignment(cost); areas = []
    for ci, gi in zip(ii, jj):
        group, before, reference = groups[gi]; after = float(caps[ci].area)
        areas.append({'id': group['id'], 'family': group['design_family'], 'before_mm2': before,
                      'after_mm2': after, 'reference_mm2': reference, 'relative_area_error': after/reference-1})
    maximum_error = max(abs(x['relative_area_error']) for x in areas)
    report = {'input_sha256': input_hash,
              'max_movement_limit_mm': args.max_movement_mm, 'projection': projection, 'edge_contraction': collapse,
              'reference_arc_max_length_mm': args.edge_mm, 'conforming_edge_splits': splits,
              'inverted_surviving_cells': 0, 'inlet_patches': len(caps), 'areas': areas,
              'max_relative_area_error': maximum_error, 'volume_change_mm3': float(volumes(points, cells).sum()-input_volume),
              'full_audit_required': True}
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.with_suffix('.projection.json').write_text(json.dumps(report, indent=2))
    if maximum_error > args.area_relative_tolerance:
        raise ValueError(f'Inlet area relative error {maximum_error} exceeds {args.area_relative_tolerance}')
    np.savez_compressed(out, points=points, cells=cells)
    used = np.unique(boundary); ids = np.full(len(points), -1); ids[used] = np.arange(len(used))
    np.savez_compressed(out.with_suffix('.boundary.npz'), vertices=points[used], faces=ids[boundary])
    print('INLET_CONTOURS_RESTORED', 'splits', splits, 'max_area_error', maximum_error, 'cells', len(cells), flush=True)

if __name__ == '__main__':
    main()
