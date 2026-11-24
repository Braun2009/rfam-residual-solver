#!/usr/bin/env python3
"""Package an audited nozzle mesh, its CAD sources, and replayable inlet repair.

No solver deployment. Sources/outputs are copied and hashed; existing packages
are left intact. All paths inside the package metadata are relative.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',required=True)
    parser.add_argument('--out',required=True)
    args=parser.parse_args()
    run=Path(args.run).resolve();out=Path(args.out).resolve()
    repo=Path(__file__).resolve().parents[1]
    if out.exists():
        raise ValueError('Package already exists; use a new directory or archive it first')
    audit=json.loads((run/'inlet_restored/mesh_audit.json').read_text())
    projection=json.loads((run/'inlet_restored_tetra_mm.projection.json').read_text())
    if audit['fluid_components_by_shared_faces']!=1 or audit['inlet_connected_patches']!=64:
        raise ValueError('Mesh connectivity is not the accepted review candidate')
    if audit['quality_mean_ratio_below_0p01'] or projection['max_relative_area_error']>.003:
        raise ValueError('Mesh fails the documented diagnostic quality/area checks')
    out.mkdir(parents=True)
    def copy(src,relative):
        dst=out/relative;dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dst)
        if digest(src)!=digest(dst):raise IOError(f'Copy hash mismatch: {relative}')
    for name in ['fluid_volume_mm.npz','fluid_volume_mm.vtu','boundary_patches_mm.vtu']:
        copy(run/'inlet_restored'/name,'mesh/'+name)
    for name in ['mesh_audit.json','geometry_audit.json']:
        copy(run/'inlet_restored'/name,'reports/'+name)
    for src,name in [('nozzle_only/geometry_manifest.json','geometry_manifest.json'),
                     ('inlet_restored_tetra_mm.projection.json','inlet_restoration.json'),
                     ('candidate_tetra_mm.collapse.json','short_edge_contraction.json'),
                     ('optimized_separated_mm.separation.json','optimizer_separation.json'),
                     ('surface_remesh/remesh_report.json','surface_remesh.json')]:
        copy(run/src,'reports/'+name)
    for p in (run/'inlet_restored/review').iterdir():
        if p.is_file():copy(p,'review/'+p.name)
    for name in ['injector_v3.STL','injector_top.STL']:
        copy(run/'source'/name,'source/'+name)
    copy(repo/'config/chamber_profile_3pintle.json','source/chamber_profile_3pintle.json')
    for name in ['geometry_manifest.json','fluid_boundary_mm.npz','assembled_head_mm.npz']:
        copy(run/'nozzle_only'/name,'geometry/'+name)
    for seed,folder in [('candidate','candidate_before_inlet'),('nozzle_tetra','reference_inlets')]:
        for name in ['fluid_volume_mm.npz','mesh_audit.json']:
            copy(run/seed/name,f'seeds/{folder}/'+name)
    for name in ['rebuild_injector_fluid.py','mesh_injector_tetra.py','inspect_injector_volume.py',
                 'audit_injector_geometry.py','remesh_injector_surface.py','prepare_injector_boundary.py',
                 'optimize_injector_volume.py','separate_optimized_tetra.py','collapse_injector_short_edge.py',
                 'restore_injector_inlet_contours.py','package_injector_review.py','requirements_mesh_review.txt']:
        copy(repo/'tools'/name,'tools/'+name)
    copy(repo/'tests/test_injector_mesh_tools.py','tests/test_injector_mesh_tools.py')
    copy(repo/'docs/mesh_rebuild_status.md','README.md')
    seed_hash=digest(out/'seeds/candidate_before_inlet/fluid_volume_mm.npz')
    if seed_hash!=projection['input_sha256']:
        raise ValueError('Saved inlet-repair seed differs from the documented input')
    commands={
        'tests':['python','tests/test_injector_mesh_tools.py'],
        'audit_saved_mesh':['python','tools/mesh_injector_tetra.py','--geometry','geometry','--out','checks/topology',
                            '--boundary','mesh/fluid_volume_mm.npz','--input-tets','mesh/fluid_volume_mm.npz'],
        'geometry_audit':['python','tools/audit_injector_geometry.py','--mesh','checks/topology',
                          '--geometry','geometry','--reference-mesh','seeds/reference_inlets',
                          '--out','checks/geometry_audit.json'],
        'rebuild_geometry':['python','tools/rebuild_injector_fluid.py','--body','source/injector_v3.STL',
                            '--top','source/injector_top.STL','--profile','source/chamber_profile_3pintle.json',
                            '--out','replay/geometry','--board-z','270.9527893066406','--inlet-z','269.2'],
        'replay_inlet_repair':['python','tools/restore_injector_inlet_contours.py','--mesh','seeds/candidate_before_inlet',
                               '--reference','seeds/reference_inlets','--out','replay/inlet_restored_tetra_mm.npz',
                               '--expected-input-sha256',seed_hash,'--max-movement-mm','0.1',
                               '--collapse-edge','107168','107172','--edge-mm','0.1'],
        'audit_replay':['python','tools/mesh_injector_tetra.py','--geometry','geometry','--out','replay/audited',
                        '--boundary','replay/inlet_restored_tetra_mm.boundary.npz',
                        '--input-tets','replay/inlet_restored_tetra_mm.npz'],
        'render_saved_mesh':['python','tools/inspect_injector_volume.py','--mesh','checks/topology',
                             '--geometry','geometry','--out','checks/review']}
    lineage=[
        {'stage':'nozzle_only','tool':'rebuild_injector_fluid.py','description':'Authoritative STL solids, two mating planes, Boolean subtract, outlet-connected component'},
        {'stage':'nozzle_tetra','tool':'mesh_injector_tetra.py','description':'Original boundary, no quality refinement; establishes unsmoothed inlet reference'},
        {'stage':'refined','tool':'mesh_injector_tetra.py','description':'Recovered tetrahedral boundary, nobisect=True, radius-edge 1.6, maxvolume=15 mm3; raw logs and arrays retained in original run'},
        {'stage':'surface_remesh','tool':'remesh_injector_surface.py','description':'Base 1.5 mm, injector 0.3 mm, throat 0.5 mm surface targets; projection and deviations recorded'},
        {'stage':'remeshed_tetra','tool':'mesh_injector_tetra.py','description':'Prepared surface boundary, quality=True, radius-edge=1.6, maxvolume=15 mm3, Steiner limit 1000000'},
        {'stage':'optimized_verified','tool':'optimize_injector_volume.py + separate_optimized_tetra.py','description':'Netgen optimization; remove exactly one proven unchanged input-mesh copy, audit remainder'},
        {'stage':'candidate','tool':'collapse_injector_short_edge.py','description':'Recorded 0.000116 mm short-edge contraction, full audit; saved as replay seed'},
        {'stage':'inlet_restored','tool':'restore_injector_inlet_contours.py','description':'CAD inlet projection, folded-rim repair and 0.1 mm arc subdivision; exact input hash required'}]
    versions={n:importlib.metadata.version(n) for n in ['numpy','scipy','trimesh','manifold3d','tetgen','meshio','vtk','matplotlib','pymeshlab','gmsh','rtree','shapely']}
    metadata={'status':'GEOMETRY_TOPOLOGY_REVIEW_ONLY','python':sys.version,'versions':versions,
              'source_commit_before_changes':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),
              'packaged_tools_are_worktree_snapshot':True,'full_pipeline_has_not_been_replayed_end_to_end':True,
              'commands_from_package_root':commands,'lineage':lineage,
              'local_node_repairs_are_guarded_by_seed_hash':True,
              'production_cfd_accepted':False,'rfam_unstructured_import_implemented':False,
              'missing':['prism wall layers','production geometry tolerance/mesh convergence','RFAM unstructured operators and ghosts']}
    (out/'reproduction.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
    quick='''# 从本目录复核（Python 3.12）\n\n运行参数以 reproduction.json 的结构化命令为准。\n\n```bash\npython -m pip install -r tools/requirements_mesh_review.txt\npython tests/test_injector_mesh_tools.py\npython tools/mesh_injector_tetra.py --geometry geometry --out checks/topology --boundary mesh/fluid_volume_mm.npz --input-tets mesh/fluid_volume_mm.npz\npython tools/audit_injector_geometry.py --mesh checks/topology --geometry geometry --reference-mesh seeds/reference_inlets --out checks/geometry_audit.json\n```\n\nVTU/NPZ 为真实体网格；review 里的 STL 是剖切体单元的表面，仅供查看。\n该检查包没有流场，也不能被现有 RFAM 结构化求解器直接读取。\n'''
    (out/'VERIFY.md').write_text(quick,encoding='utf-8')
    sums={str(p.relative_to(out)).replace('\\','/'): {'bytes':p.stat().st_size,'sha256':digest(p)}
          for p in sorted(out.rglob('*')) if p.is_file()}
    (out/'SHA256SUMS.json').write_text(json.dumps(sums,indent=2),encoding='utf-8')
    print('REVIEW_PACKAGE_WRITTEN',out,'files',len(sums),'bytes',sum(x['bytes'] for x in sums.values()),flush=True)

if __name__=='__main__':main()
