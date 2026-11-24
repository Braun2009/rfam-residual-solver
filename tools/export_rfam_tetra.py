#!/usr/bin/env python3
"""Convert audited tetra NPZ to RFAMTET1 binary mesh for the C++ loader."""
import argparse,struct
from pathlib import Path
import numpy as np

LOCAL=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]],dtype=np.int32)

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args()
 d=np.load(a.input)
 # RFAM runtime is SI: CAD review mesh is mm, convert coordinates to m.
 points=np.asarray(d['points'],dtype='<f8')*1e-3;cells=np.asarray(d['cells'],dtype='<i4')
 # Build one oriented face record per tetrahedron face.
 raw=cells[:,LOCAL].reshape(-1,3); keys=np.sort(raw,axis=1)
 uniq,first,inv,count=np.unique(keys,axis=0,return_index=True,return_inverse=True,return_counts=True)
 if np.any(count>2):raise ValueError('nonmanifold face')
 face=[];owner=[];neigh=[]
 order=np.argsort(inv,kind='stable')
 starts=np.r_[0,np.cumsum(count)[:-1]]
 for group,(idx,n,start) in enumerate(zip(first,count,starts)):
  rec=raw[idx];o=idx//4
  neigh.append(int(order[start+1]//4) if n==2 else -1)
  face.append(rec);owner.append(o)
 # boundary zone map from audited boundary faces.
 zone=np.zeros(len(face),dtype='<i4');
 bfaces=np.sort(np.asarray(d['boundary'],dtype=np.int32),axis=1);bzone=np.asarray(d['boundary_zone'],dtype=np.int32)
 lookup={tuple(x):int(z) for x,z in zip(bfaces,bzone)}
 for i,x in enumerate(face):
  if neigh[i]<0: zone[i]=lookup.get(tuple(sorted(x)),0)
 # Positive volumes are already audited; compute fresh in SI (points are m).
 x=points[cells];vol=np.abs(np.einsum('ij,ij->i',x[:,1]-x[:,0],np.cross(x[:,2]-x[:,0],x[:,3]-x[:,0])))/6
 if np.any(~np.isfinite(vol)) or np.any(vol<=0):raise ValueError('invalid tetra volume')
 out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
 with out.open('wb') as f:
  f.write(b'RFAMTET1');f.write(struct.pack('<QQQ',len(points),len(cells),len(face)))
  f.write(points.tobytes());f.write(cells.tobytes());f.write(np.asarray(vol,dtype='<f8').tobytes())
  for n,o,nn,z in zip(face,owner,neigh,zone):
   f.write(struct.pack('<3i3i',*n,int(o),int(nn),int(z)))
 print('RFAM_TETRA_EXPORTED',out,'points',len(points),'cells',len(cells),'faces',len(face),'boundary',sum(x<0 for x in neigh),'volume_m3',float(vol.sum()))
if __name__=='__main__':main()
