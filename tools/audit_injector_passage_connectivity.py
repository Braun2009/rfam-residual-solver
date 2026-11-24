#!/usr/bin/env python3
"""Audit passage-only connected components upstream of outlet mergers.

Selection stops BEFORE the observed fuel/LOX/torch passages join the exterior.
The output is a geometry mask, not a physical field or a proven exit surface.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import struct
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from tetra_field_to_vtu import read_mesh

FAMILIES={'lox':([10,11,12],280.9,3),
          'fuel':([20,21,22],276.8,60),
          'torch':([30],283.3,1)}


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as file:
        for piece in iter(lambda:file.read(1024*1024),b''):
            digest.update(piece)
    return digest.hexdigest()


def selection(mesh_path,include_torch=True):
    points,cells,volume=read_mesh(mesh_path)
    raw=Path(mesh_path).read_bytes()
    npoints,ncells,nfaces=struct.unpack_from('<QQQ',raw,8)
    if len(cells)!=ncells or len(points)!=npoints:raise ValueError('mesh header mismatch')
    faces=np.frombuffer(raw,'<i4',nfaces*6,32+npoints*24+ncells*24).reshape(-1,6)
    last_z=points[cells][:,:,2].max(axis=1)
    tag=np.zeros(ncells,dtype=np.uint8)
    details={}
    for number,(family,(zones,cutoff_mm,required)) in enumerate(FAMILIES.items(),start=1):
        if family=='torch' and not include_torch:continue
        inside=last_z<cutoff_mm*1e-3
        remap=np.full(ncells,-1,dtype=np.int32)
        remap[inside]=np.arange(np.count_nonzero(inside),dtype=np.int32)
        internal=faces[faces[:,4]>=0][:,3:5]
        pairs=internal[inside[internal[:,0]]&inside[internal[:,1]]]
        edge=remap[pairs]
        graph=coo_matrix((np.ones(len(edge),dtype=np.uint8),(edge[:,0],edge[:,1])),
                         shape=(np.count_nonzero(inside),)*2).tocsr()
        ncomponents,labels=connected_components(graph,directed=False)
        inlet=faces[(faces[:,4]<0)&np.isin(faces[:,5],zones)]
        owners=inlet[:,3]
        if not inside[owners].all():raise ValueError('inlet lies outside passage cutoff')
        roots=np.unique(labels[remap[owners]])
        if len(roots)!=required:
            raise ValueError(f'{family} passage count changed: {len(roots)} vs {required}')
        selected=np.flatnonzero(inside)[np.isin(labels,roots)]
        if np.any(tag[selected]):raise ValueError('passage masks overlap')
        tag[selected]=number
        details[family]={'cutoff_z_mm':cutoff_mm,'components':len(roots),
            'cells':len(selected),'volume_cm3':float(np.sum(volume[selected])*1e6),
            'highest_selected_vertex_z_mm':float(np.max(last_z[selected])*1e3),
            'zone_ids':zones,'graph_component_count':int(ncomponents)}
    return tag,details


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',required=True)
    ap.add_argument('--out',required=True)
    ap.add_argument('--mask-out',help='optional compressed uint8 cell mask')
    args=ap.parse_args()
    tags,details=selection(args.mesh)
    result={'scope':'geometry connectivity only, not injector exit or initialized CFD field',
        'mesh':args.mesh,'mesh_sha256':sha256(args.mesh),
        'family_ids':{'lox':1,'fuel':2,'torch':3},
        'selected_cells':int(np.count_nonzero(tags)),
        'families':details}
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if args.mask_out:
        np.savez_compressed(args.mask_out,tag=tags,mesh_sha256=result['mesh_sha256'])
    print('INJECTOR_PASSAGE_CONNECTIVITY_PASS',result['selected_cells'],
          [(family,data['components'],data['cells']) for family,data in details.items()])


if __name__=='__main__':main()
