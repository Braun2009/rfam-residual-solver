import sys, struct, numpy as np
def zones(path):
    with open(path,'rb') as f:
        m=f.read(8)
        if m!=b'RFAMTET1': return None
        np_,nc,nf=struct.unpack('<QQQ',f.read(24))
        f.seek(np_*24+nc*16+nc*8,1)
        raw=f.read(nf*28)
        faces=np.frombuffer(raw,dtype=np.dtype([('n','<i4',3),('own','<i4'),('nb','<i4'),('zone','<i4')]))
        pts=np.frombuffer(open(path,'rb').read(),dtype=np.uint8)  # placeholder
    return faces
for p in sys.argv[1:]:
    try:
        with open(p,'rb') as f:
            m=f.read(8)
            if m!=b'RFAMTET1': print(p,'BADMAGIC'); continue
            np_,nc,nf=struct.unpack('<QQQ',f.read(24))
            pts=np.frombuffer(f.read(np_*24),dtype='<f8').reshape(np_,3)
            f.seek(nc*16+nc*8,1)
            faces=np.frombuffer(f.read(nf*28),dtype=np.dtype([('n','<i4',3),('own','<i4'),('nb','<i4'),('zone','<i4')]))
        b=faces[faces['nb']<0]; p3=pts[b['n']]
        a=0.5*np.linalg.norm(np.cross(p3[:,1]-p3[:,0],p3[:,2]-p3[:,0]),axis=1)
        lo={}; 
        for z in np.unique(b['zone']):
            m=b['zone']==z; lo[int(z)]=a[m].sum()*1e6
        lox=sum(v for k,v in lo.items() if 10<=k<=12)
        fuel=sum(v for k,v in lo.items() if 20<=k<=22)
        print('%-58s cells=%-8d LOX=%7.2f mm2 FUEL=%8.2f' % (p.split('/')[-2]+'/'+p.split('/')[-1], nc, lox, fuel))
    except Exception as e:
        print(p,'ERR',e)
