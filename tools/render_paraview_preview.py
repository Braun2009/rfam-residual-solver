"""Open the saved ParaView state and write a PNG preview (headless)."""
from paraview.simple import LoadState, Render, GetActiveViewOrCreate, SaveScreenshot

state = r'D:/realFluidAllMach/work/mesh_rebuild/delivery/view_stl/paraview_nozzle_state.pvsm'
png = r'D:/realFluidAllMach/work/mesh_rebuild/delivery/view_stl/paraview_preview.png'

LoadState(state)
view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1600, 1000]
Render()
SaveScreenshot(png, view, ImageResolution=[1600, 1000])
print('SCREENSHOT_SAVED', png)
