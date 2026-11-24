@echo off
rem Open the core files of the A-route injector simulation (Notepad4, one instance)
set NP=D:\Notepad4\build\bin\Release\x64\Notepad4.exe
set R=D:\realFluidAllMach
start "" "%NP%" "%R%\solver_3d_tetra.cpp" "%R%\src\mesh\TetraFiniteVolume.C" "%R%\src\mesh\TetraFiniteVolume.H" "%R%\src\mesh\TetraLUSGS.H" "%R%\src\numerics\flux\AUSMPlusUpFlux.C" "%R%\src\thermo\CoolPropThermoModel.C" "%R%\src\thermo\CoolPropEquilibrium.H" "%R%\src\thermo\StateAcceptance.H" "%R%\src\thermo\TableAcoustic.H" "%R%\src\variables\CompleteCellState.C" "%R%\src\mesh\TetraStateRecovery.H" "%R%\src\mesh\TetraGradients.H" "%R%\src\mesh\TetraViscousFlux.H"
