#!/usr/bin/env python3
"""Export raw Gmsh/Netgen optimization output for mandatory mesh auditing.

Some discrete-volume runs retain the original mesh alongside the optimized
mesh. Such output must be separated with separate_optimized_tetra.py, then
independently audited; a successful optimizer call is not a valid mesh claim.
"""
import argparse
import time
from pathlib import Path
import numpy as np
import gmsh


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--method',default='Netgen',choices=['Netgen','Relocate3D','default'])
    p.add_argument('--write-msh',action='store_true')
    a=p.parse_args();d=np.load(a.mesh);points=d['points'];cells=d['cells'];faces=d['boundary']
    boundary_nodes=np.unique(faces);inside=np.setdiff1d(np.arange(len(points)),boundary_nodes)
    start=time.time();gmsh.initialize();gmsh.option.setNumber('General.Terminal',1)
    try:
        gmsh.model.add('nozzle_volume')
        s=gmsh.model.addDiscreteEntity(2)
        v=gmsh.model.addDiscreteEntity(3,boundary=[s])
        gmsh.model.mesh.addNodes(2,s,boundary_nodes+1,points[boundary_nodes].ravel())
        gmsh.model.mesh.addNodes(3,v,inside+1,points[inside].ravel())
        gmsh.model.mesh.addElementsByType(s,2,[],(faces+1).ravel())
        gmsh.model.mesh.addElementsByType(v,4,[],(cells+1).ravel())
        gmsh.model.mesh.optimize('' if a.method=='default' else a.method,force=True,niter=3)
        tags,coords,_=gmsh.model.mesh.getNodes();coords=np.asarray(coords).reshape(-1,3)
        ids=np.empty(int(max(tags))+1,dtype=np.int64);ids[tags]=np.arange(len(tags))
        entities=gmsh.model.getEntities(3)
        print('VOLUME_ENTITIES',entities,flush=True)
        for dim,tag in entities:
            etypes,etags,_=gmsh.model.mesh.getElements(dim,tag)
            print('ENTITY_ELEMENTS',dim,tag,[(int(t),len(e)) for t,e in zip(etypes,etags)],flush=True)
        _,conn=gmsh.model.mesh.getElementsByType(4)
        c=ids[np.asarray(conn).reshape(-1,4)]
        out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
        if a.write_msh:
            gmsh.write(str(out.with_suffix('.msh')))
        np.savez_compressed(out,points=coords,cells=c)
        print('OPTIMIZER_RAW_REQUIRES_AUDIT',len(c),'seconds',time.time()-start,flush=True)
    finally:gmsh.finalize()

if __name__=='__main__':main()
