/*---------------------------------------------------------------------------*\
    test_vtk_output - VTKWriter smoke test: writes a 2D field and verifies
    the output file is well-formed (non-empty, correct header, cell count).
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/VTKWriter.H"
#include "src/turbulence/TurbulenceOperator.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/SolverBase.H"

using namespace Foam;
using namespace std;

static int g_pass=0, g_fail=0;
static void check(const string& n, bool ok)
{
    cout << "  [" << (ok?"PASS":"FAIL") << "] " << n << "\n";
    ok ? ++g_pass : ++g_fail;
}

int main()
{
    cout << "=== VTKWriter output smoke test ===\n\n";

    // 1. Filename helpers
    string f0 = VTKWriter::filename("rfam", 42, 6);
    string ft = VTKWriter::filenameTime("rfam", 1.5e-4, 2);
    cout << "  filename(42)   = " << f0 << "\n";
    cout << "  filenameTime() = " << ft << "\n\n";
    check("filename ends with .vtk", f0.find(".vtk") != string::npos);
    check("filename zero-padded step", f0.find("000042") != string::npos);
    check("filenameTime contains 'us'", ft.find("us") != string::npos);

    // 2. Build a small 2D domain and write a VTK file
    const int Nx=8, Ny=4;
    Domain dom(Nx, Ny, 0.01, 0.005, 2);
    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
    {
        PrimitiveVars& w=dom.stagedW(i,j);
        w.p=1e6; w.T=300.0+j*100;
        w.U[0]=50.0; w.U[1]=w.U[2]=0.0;
        w.Y[0]=1.0; w.Y[1]=0.0;
    }
    dom.primitiveToConserved();

    const string tmpvtk = "/tmp/rfam_vtk_test.vtk";
    {
        VTKWriter vtk(tmpvtk, dom);
        check("VTKWriter opened OK", vtk.good());
        vtk.writeScalar("temperature", [&](int i,int j){ return dom.W(i,j).T; });
        vtk.writeScalar("pressure",    [&](int i,int j){ return dom.W(i,j).p; });
        vtk.writeVelocity();
    }

    // 3. Check file content
    ifstream f(tmpvtk);
    check("VTK file exists and opens", f.good());
    string line; int lines=0;
    bool hasHeader=false, hasDataset=false, hasCellData=false;
    while (getline(f, line))
    {
        ++lines;
        if (line.find("vtk DataFile") != string::npos) hasHeader=true;
        if (line.find("RECTILINEAR_GRID") != string::npos) hasDataset=true;
        if (line.find("CELL_DATA") != string::npos)
        {
            hasCellData=true;
            // Verify cell count matches Nx*Ny
            size_t pos = line.find("CELL_DATA ");
            if (pos != string::npos)
            {
                int n = stoi(line.substr(pos+10));
                check("CELL_DATA count == Nx*Ny", n == Nx*Ny);
            }
        }
    }
    check("VTK header present",  hasHeader);
    check("RECTILINEAR_GRID set", hasDataset);
    check("CELL_DATA line found", hasCellData);
    check("File has >10 lines",   lines > 10);

    // 4. Turbulence field output (compile-check only, values not tested)
    {
        TurbulenceOperator turb(1.5e-5);
        SolverBase solver(dom);
        solver.addOperator(unique_ptr<Operator>(
            new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
        solver.addOperator(unique_ptr<Operator>(new TurbulenceOperator(1.5e-5)));
        // Don't run — just ensure compile + write path compiles cleanly.
        const string tmpturb = "/tmp/rfam_vtk_turb.vtk";
        VTKWriter vtk2(tmpturb, dom);
        vtk2.writeTurbulence(turb);
        check("turbulence VTK write compiles+runs", vtk2.good());
    }

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "VTK_OUTPUT_PASS" : "VTK_OUTPUT_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
