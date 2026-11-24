"""Build a ParaView state that opens the nozzle mesh with a ready-made half cut."""
from paraview.simple import (
    OpenDataFile, Show, Hide, GetActiveViewOrCreate, Render, ResetCamera,
    ColorBy, SaveState, Clip, SetActiveSource, GetActiveSource, Interact,
    GetDisplayProperties,
)

mesh = r'D:/realFluidAllMach/work/mesh_rebuild/delivery/mesh/fluid_volume_mm.npz'
vtu = r'D:/realFluidAllMach/work/mesh_rebuild/delivery/mesh/fluid_volume_mm.vtu'
state = r'D:/realFluidAllMach/work/mesh_rebuild/delivery/view_stl/paraview_nozzle_state.pvsm'

reader = OpenDataFile(vtu)
reader.UpdatePipeline()
bounds = reader.GetDataInformation().GetBounds()
print('bounds', bounds)

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1600, 1000]
display = Show(reader, view)
display.Representation = 'Surface With Edges'
display.EdgeColor = [0.25, 0.32, 0.38]
cell_arrays = [name for name in reader.CellData.keys()] if hasattr(reader, 'CellData') else []
print('cell arrays:', cell_arrays)
if cell_arrays:
    ColorBy(display, ('CELLS', cell_arrays[0]))
else:
    display.DiffuseColor = [0.72, 0.78, 0.84]

# Half cut along the transverse axis through the pintle axis, so the interior
# passages and the metal between them are visible.
centre_y = bounds[2] + 0.5 * (bounds[3] - bounds[2])
clip = Clip(Input=reader)
clip.ClipType = 'Plane'
clip.ClipType.Origin = [0.0, centre_y, 0.0]
clip.ClipType.Normal = [0.0, 1.0, 0.0]
Show(clip, view)
clip_display = GetDisplayProperties(clip, view)
clip_display.Representation = 'Surface With Edges'
clip_display.DiffuseColor = [0.62, 0.72, 0.80]
Hide(reader, view)

SetActiveSource(clip)
Render()
ResetCamera()
camera = view.GetActiveCamera()
camera.Azimuth(35); camera.Elevation(18); camera.Zoom(1.05)
Render()
SaveState(state)
print('STATE_SAVED', state)
