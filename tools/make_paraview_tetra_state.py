"""Build a ParaView state for a tetra cold-flow VTU (Mach + composition views)."""
import sys
from paraview.simple import (
    OpenDataFile, Show, Hide, GetActiveViewOrCreate, Render, ResetCamera,
    ColorBy, SaveState, Clip, Slice, SaveScreenshot, GetDisplayProperties,
    SetActiveSource, GetColorTransferFunction, GetScalarBar, Text,
)

vtu = sys.argv[1] if len(sys.argv) > 1 else \
    r'D:/realFluidAllMach/work/mesh_rebuild/tetra_coldflow_gas.vtu'
state = sys.argv[2] if len(sys.argv) > 2 else \
    r'D:/realFluidAllMach/work/mesh_rebuild/tetra_coldflow_gas.pvsm'
png = sys.argv[3] if len(sys.argv) > 3 else \
    r'D:/realFluidAllMach/work/mesh_rebuild/tetra_coldflow_mach.png'

label = sys.argv[4] if len(sys.argv)>4 else ''
reader = OpenDataFile(vtu)
reader.UpdatePipeline()
bounds = reader.GetDataInformation().GetBounds()
print('bounds', bounds)

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1600, 1000]
display = Show(reader, view)
display.Representation = 'Surface'
ColorBy(display, ('CELLS', 'Mach'))
lut = GetColorTransferFunction('Mach')
lut.RescaleTransferFunction(0.0, 2.0) if hasattr(lut, 'RescaleTransferFunction') else None
bar = GetScalarBar(lut, view)
bar.Title = 'Mach'
bar.ComponentTitle = ''

# Half cut through the pintle axis so the interior mixing is visible.
centre_y = bounds[2] + 0.5 * (bounds[3] - bounds[2])
clip = Clip(Input=reader)
clip.ClipType = 'Plane'
clip.ClipType.Origin = [0.0, centre_y, 0.0]
clip.ClipType.Normal = [0.0, 1.0, 0.0]
clip_display = Show(clip, view)
clip_display.Representation = 'Surface'
ColorBy(clip_display, ('CELLS', 'Mach'))
Hide(reader, view)

if label:
    annotation=Text()
    annotation.Text=label
    annotation_display=Show(annotation,view)
    annotation_display.WindowLocation='Upper Left Corner'
    annotation_display.FontSize=16
    annotation_display.Color=[1.0,0.75,0.25]

SetActiveSource(clip)
Render()
ResetCamera()
camera = view.GetActiveCamera()
camera.Azimuth(30); camera.Elevation(16); camera.Zoom(1.05)
Render()
SaveState(state)
SaveScreenshot(png, view, ImageResolution=[1600, 1000])
print('STATE_SAVED', state)
print('SCREENSHOT_SAVED', png)
