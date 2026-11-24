#!/usr/bin/env python3
"""Remove an unchanged input-mesh multiset from a combined optimizer output.

Gmsh/Netgen on this discrete-volume input retained original cells alongside
optimized cells. This tool only operates when *every* original node and cell
is present exactly. It removes one copy of each original cell, preserving any
unchanged cell that also occurs in the optimized mesh. The remaining mesh must
be audited independently by mesh_injector_tetra.py.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',required=True);p.add_argument('--combined',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();old=np.load(a.original);new=np.load(a.combined)
    v=old['points'];c=old['cells'];points=new['points'];cells=new['cells']
    distance,node_map=cKDTree(points).query(v)
    if not np.all(distance==0) or len(np.unique(node_map))!=len(v):
        raise ValueError('Optimizer did not retain all original vertices exactly; cannot separate')
    orig=np.ascontiguousarray(np.sort(node_map[c],axis=1),dtype=np.int64)
    combined=np.ascontiguousarray(np.sort(cells,axis=1),dtype=np.int64)
    dtype=np.dtype([(str(i),np.int64) for i in range(4)])
    orig_keys=orig.view(dtype).ravel();new_keys=combined.view(dtype).ravel()
    orig_unique,orig_counts=np.unique(orig_keys,return_counts=True)
    new_unique,first,counts=np.unique(new_keys,return_index=True,return_counts=True)
    pos=np.searchsorted(new_unique,orig_unique)
    if np.any(pos>=len(new_unique)) or not np.array_equal(new_unique[pos],orig_unique):
        raise ValueError('Input tetrahedral cells are not a complete subset of combined output')
    if np.any(orig_counts!=1):raise ValueError('Original mesh has duplicate cells')
    # Exactly one input copy is removed. Equal optimized cells remain.
    keep=np.ones(len(cells),dtype=bool);keep[first[pos]]=False
    remaining=cells[keep]
    keys=np.sort(remaining,axis=1)
    if len(np.unique(keys,axis=0))!=len(remaining):raise ValueError('Optimized remainder contains duplicate cells')
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out,points=points,cells=remaining)
    report={'input_vertices_matched_exactly':True,'all_input_cells_present':True,
            'original_cells_removed_exactly_once':len(c),'combined_cells':len(cells),
            'remaining_cells':len(remaining),'remaining_duplicate_cells':0,
            'remaining_requires_independent_audit':True}
    out.with_suffix('.separation.json').write_text(json.dumps(report,indent=2))
    print('OPTIMIZED_MESH_SEPARATED',report,flush=True)

if __name__=='__main__':main()
