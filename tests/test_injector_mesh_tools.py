"""Geometry/topology regressions for the nozzle mesh tools (Python 3.12)."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import trimesh

TOOLS = Path(__file__).resolve().parents[1] / 'tools'
sys.path.insert(0, str(TOOLS))
from rebuild_injector_fluid import load_welded
from restore_injector_inlet_contours import (
    LocalElements, TET_FACES, on_loop, perimeter, polyline_loops, project, volumes,
)


def boundary_of(cells):
    faces = cells[:, TET_FACES].reshape(-1, 3)
    _, first, count = np.unique(np.sort(faces, axis=1), axis=0,
                                return_index=True, return_counts=True)
    if np.any(count > 2):
        raise AssertionError('Nonmanifold test mesh')
    return faces[first[count == 1]]


class InjectorMeshToolsTest(unittest.TestCase):
    def test_stl_triangle_soup_is_welded_before_connectivity_check(self):
        cube = trimesh.creation.box()
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'cube.stl'
            cube.export(p)
            welded = load_welded(p)
            self.assertTrue(welded.is_watertight)
            self.assertEqual(len(welded.vertices), 8)
            self.assertAlmostEqual(welded.volume, 1.0)
            cube.update_faces(np.arange(len(cube.faces)) != 0)
            cube.export(p)
            with self.assertRaises(ValueError):
                load_welded(p)

    def test_edge_split_preserves_volume_and_shared_faces(self):
        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                           [0, 0, 1], [0, 0, -1], [.5, 0, 0]], float)
        cells = np.array([[0, 1, 2, 3], [0, 2, 1, 4]])
        before_boundary = boundary_of(cells)
        star = LocalElements(cells)
        skin = LocalElements(before_boundary)
        self.assertEqual(star.split_edge(0, 1, 5), 2)
        self.assertEqual(skin.split_edge(0, 1, 5), 2)
        after = star.array(4)
        self.assertEqual(len(after), 4)
        self.assertTrue(np.all(volumes(points, after) > 0))
        self.assertAlmostEqual(volumes(points, after).sum(), 1/3)
        actual = {tuple(f) for f in np.sort(boundary_of(after), axis=1)}
        expected = {tuple(f) for f in np.sort(skin.array(3), axis=1)}
        self.assertEqual(actual, expected)
        self.assertEqual(len(actual), 8)

    def test_annular_cap_keeps_inner_and_outer_loops(self):
        p = np.array([[-2, -2, 0], [2, -2, 0], [2, 2, 0], [-2, 2, 0],
                      [-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]], float)
        faces = []
        for i in range(4):
            j = (i+1) % 4
            faces.extend([[i, j, j+4], [i, j+4, i+4]])
        loops = polyline_loops(p, perimeter(np.array(faces)))
        self.assertEqual(len(loops), 2)
        self.assertEqual(sorted(float(arc[-1]) for _, arc in loops), [8., 16.])
        with self.assertRaises(ValueError):
            polyline_loops(p, np.array([[0, 1], [1, 2]]))

    def test_reference_arc_wrap_and_projection(self):
        loop = (np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float), np.arange(5.))
        np.testing.assert_allclose(on_loop(loop, 4.5), [.5, 0])
        xy, arc, error = project(loop, np.array([1.2, .25]))
        np.testing.assert_allclose(xy, [1., .25])
        self.assertAlmostEqual(arc, 1.25)
        self.assertAlmostEqual(error, .2)

    def test_optimizer_multiset_preserves_unchanged_optimized_cell(self):
        # The optimizer retains the original two cells and emits a second
        # mesh: one unchanged cell plus two children of the other cell.
        p = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                      [3, 0, 0], [4, 0, 0], [3, 1, 0], [3, 0, 1], [3.5, 0, 0]], float)
        old = np.array([[0, 1, 2, 3], [4, 5, 6, 7]])
        combined = np.vstack([old, old[:1], [[4, 8, 6, 7], [8, 5, 6, 7]]])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            np.savez(root/'old.npz', points=p[:8], cells=old)
            np.savez(root/'new.npz', points=p, cells=combined)
            result = subprocess.run([sys.executable, str(TOOLS/'separate_optimized_tetra.py'),
                                     '--original', str(root/'old.npz'), '--combined', str(root/'new.npz'),
                                     '--out', str(root/'result.npz')], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with np.load(root/'result.npz') as data:
                self.assertEqual(len(data['cells']), 3)
                self.assertAlmostEqual(volumes(data['points'], data['cells']).sum(), 1/3)
            report = json.loads((root/'result.separation.json').read_text())
            self.assertEqual(report['original_cells_removed_exactly_once'], 2)


if __name__ == '__main__':
    unittest.main()
