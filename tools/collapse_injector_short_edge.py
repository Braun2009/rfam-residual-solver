#!/usr/bin/env python3
"""Contract one sub-tolerance edge, retaining geometric and topological evidence.

This is a mesh operation, not removal of a selected bad cell. Every incident
cell is remapped; only topologically degenerate tetrahedra caused by the edge
contraction are removed. A separate full audit is mandatory after this step.
"""
import argparse,json
from pathlib import Path
import numpy as np


def signed_volume(points,cells):
    x=points[cells]
    return np.einsum('ij,ij->i',x[:,1]-x[:,0],np.cross(x[:,2]-x[:,0],x[:,3]-x[:,0]))/6


def quality(points,cells):
    x=points[cells];v=np.abs(signed_volume(points,cells))
    edge2=((x[:,[0,0,0,1,1,2]]-x[:,[1,2,3,2,3,3]])**2).sum(axis=(1,2))
    return 12*(3*v)**(2/3)/edge2


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',required=True);p.add_argument('--out',required=True)
    p.add_argument('--tolerance-mm',type=float,default=.001)
    a=p.parse_args();d=np.load(a.mesh);points=d['points'].copy();cells=d['cells'].copy()
    q=quality(points,cells);bad=np.flatnonzero(q<.01)
    if len(bad)!=1:raise ValueError('Local repair expects exactly one identified sliver')
    target=cells[bad[0]];pairs=np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
    lengths=np.linalg.norm(points[target[pairs[:,0]]]-points[target[pairs[:,1]]],axis=1)
    edge=target[pairs[np.argmin(lengths)]];length=float(lengths.min())
    if not 0<length<a.tolerance_mm:raise ValueError('Short edge exceeds the recorded repair tolerance')
    old_volume=float(signed_volume(points,cells).sum())
    incident=np.any(np.isin(cells,edge),axis=1)
    old_incident_sign=signed_volume(points,cells[incident])
    points[edge[0]]=(points[edge[0]]+points[edge[1]])/2
    cells[cells==edge[1]]=edge[0]
    sorted_cells=np.sort(cells,axis=1)
    degenerate=np.any(np.diff(sorted_cells,axis=1)==0,axis=1)
    keep=~degenerate
    # The operation must not invert any nondegenerate incident tetrahedron.
    new_sign=signed_volume(points,cells[incident])
    surviving=~degenerate[incident]
    if np.any(old_incident_sign[surviving]*new_sign[surviving]<=0):
        raise ValueError('Edge contraction would invert an incident tetrahedron')
    cells=cells[keep]
    if len(np.unique(np.sort(cells,axis=1),axis=0))!=len(cells):raise ValueError('Contraction created duplicate cells')
    after=quality(points,cells)
    if after.min()<=q.min() or np.sum(after<.01)>0:raise ValueError('Contraction did not repair minimum quality')
    used=np.unique(cells);ids=np.full(len(points),-1);ids[used]=np.arange(len(used))
    cells=ids[cells];points=points[used]
    new_volume=float(signed_volume(points,cells).sum())
    report={'edge_node_ids_before':edge.tolist(),'edge_length_mm':length,
            'max_vertex_movement_mm':length/2,'tolerance_mm':a.tolerance_mm,
            'incident_cells':int(incident.sum()),'collapsed_degenerate_cells':int(degenerate.sum()),
            'surviving_incident_inversions':0,'duplicate_remaining_cells':0,
            'volume_change_mm3':new_volume-old_volume,
            'quality_min_before':float(q.min()),'quality_min_after':float(after.min()),
            'quality_below_0p01_after':int((after<.01).sum()),'cells_after':len(cells),
            'full_audit_required':True}
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(out,points=points,cells=cells)
    out.with_suffix('.collapse.json').write_text(json.dumps(report,indent=2))
    print('LOCAL_EDGE_CONTRACTION',report,flush=True)

if __name__=='__main__':main()
