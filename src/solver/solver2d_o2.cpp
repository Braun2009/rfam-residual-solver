#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include "../mesh/Mesh2D.H"
#include "../variables/ConservedVars.H"
#include "../numerics/flux/AUSMPlusUpFlux.H"
#include "../thermo/RealFluidState.H"
#include "../thermo/EquationOfState/PengRobinsonEOS.H"
#include "../thermo/MixingRules/SpeciesCriticalDB.H"
#include "../thermo/MixingRules/BinaryInteractionDB.H"

using namespace Foam;
using namespace std;

class Solver2D {
private:
    Mesh2D mesh_;
    vector<vector<ConservedVars>> U_, W_old_;
    vector<vector<PrimitiveVars>> W_;
    AUSMPlusUpFlux fluxScheme_;
    SpeciesCriticalDB criticalDB_;
    BinaryInteractionDB kijDB_;
    PengRobinsonEOS eos_;
    double t_;
    int iter_;

public:
    Solver2D(int Nx, int Ny, double Lx, double Ly)
        : mesh_(Nx, Ny, Lx, Ly), 
          criticalDB_(SpeciesCriticalDB::createDefault()),
          kijDB_(BinaryInteractionDB::createDefault()),
          eos_(criticalDB_, kijDB_),
          t_(0.0), iter_(0)
    {
        U_.resize(Nx); W_.resize(Nx);
        for (int i=0; i<Nx; i++) {
            U_[i].resize(Ny);
            W_[i].resize(Ny);
        }
        mesh_.generate(0.0, 2.0);
    }

    void initialize() {
        cout << "Init 2D domain: " << mesh_.Nx() << "x" << mesh_.Ny() << endl;
        double r_noz = 0.0005;
        for (int i=0; i<mesh_.Nx(); i++) {
            for (int j=0; j<mesh_.Ny(); j++) {
                double x = mesh_(i,j).x, y = mesh_(i,j).y;
                if (x<0.001 && y<r_noz) {
                    W_[i][j].p=6e6; W_[i][j].T=100; W_[i][j].U[0]=10; W_[i][j].Y[0]=1; W_[i][j].Y[1]=0;
                } else {
                    W_[i][j].p=0.5e6; W_[i][j].T=300; W_[i][j].U[0]=1; W_[i][j].Y[0]=0; W_[i][j].Y[1]=1;
                }
                W_[i][j].U[1]=W_[i][j].U[2]=0;
                StateConversion::primitiveToConserved(W_[i][j], U_[i][j]);
            }
        }
    }

    double computeTimestep(double CFL) {
        double dt=1e10;
        for (int i=0; i<mesh_.Nx(); i++) {
            for (int j=0; j<mesh_.Ny(); j++) {
                RealFluidState s; s.updatePTY(W_[i][j].p, W_[i][j].T, W_[i][j].Y);
                double a=s.a(), u=sqrt(W_[i][j].U[0]*W_[i][j].U[0]+W_[i][j].U[1]*W_[i][j].U[1]);
                double ds=min(mesh_(i,j).dx, mesh_(i,j).dy);
                dt=min(dt, CFL*ds/(u+a+1e-10));
            }
        }
        return dt;
    }

    void step(double dt) {
        auto Un = U_;
        vector<double> nx(3,0), ny(3,0); nx[0]=1; ny[1]=1;
        for (int i=1; i<mesh_.Nx()-1; i++) {
            for (int j=1; j<mesh_.Ny()-1; j++) {
                ConservedVars fxL,fxR,fyL,fyR;
                fluxScheme_.computeFlux(W_[i-1][j],W_[i][j],nx,fxL);
                fluxScheme_.computeFlux(W_[i][j],W_[i+1][j],nx,fxR);
                fluxScheme_.computeFlux(W_[i][j-1],W_[i][j],ny,fyL);
                fluxScheme_.computeFlux(W_[i][j],W_[i][j+1],ny,fyR);
                double dx=mesh_(i,j).dx, dy=mesh_(i,j).dy, v=mesh_(i,j).volume, y=mesh_(i,j).y;
                double Ax=2*M_PI*y*dy, Ayl=2*M_PI*(y-0.5*dy)*dx, Ayh=2*M_PI*(y+0.5*dy)*dx;
                Un[i][j].rho = U_[i][j].rho - dt/v*((fxR.rho-fxL.rho)*Ax+(fyR.rho*Ayh-fyL.rho*Ayl));
                for (int k=0;k<3;k++) Un[i][j].rhoU[k]=U_[i][j].rhoU[k]-dt/v*((fxR.rhoU[k]-fxL.rhoU[k])*Ax+(fyR.rhoU[k]*Ayh-fyL.rhoU[k]*Ayl));
                Un[i][j].rhoE = U_[i][j].rhoE - dt/v*((fxR.rhoE-fxL.rhoE)*Ax+(fyR.rhoE*Ayh-fyL.rhoE*Ayl));
                for (size_t k=0;k<U_[i][j].rhoY.size();k++) Un[i][j].rhoY[k]=U_[i][j].rhoY[k]-dt/v*((fxR.rhoY[k]-fxL.rhoY[k])*Ax+(fyR.rhoY[k]*Ayh-fyL.rhoY[k]*Ayl));
            }
        }
        applyBC(Un);
        U_=Un;
        for (int i=0;i<mesh_.Nx();i++) for (int j=0;j<mesh_.Ny();j++) StateConversion::conservedToPrimitive(U_[i][j],W_[i][j]);
        t_+=dt; iter_++;
    }

    void applyBC(vector<vector<ConservedVars>>& U) {
        double r=0.0005;
        for (int j=0;j<mesh_.Ny();j++) {
            PrimitiveVars b;
            if (mesh_(0,j).y<r) { b.p=6e6;b.T=100;b.U[0]=10;b.Y[0]=1;b.Y[1]=0; }
            else { b.p=0.5e6;b.T=300;b.U[0]=1;b.Y[0]=0;b.Y[1]=1; }
            b.U[1]=b.U[2]=0; StateConversion::primitiveToConserved(b,U[0][j]);
        }
        for (int j=0;j<mesh_.Ny();j++) U[mesh_.Nx()-1][j]=U[mesh_.Nx()-2][j];
        for (int i=0;i<mesh_.Nx();i++) { U[i][0]=U[i][1]; U[i][0].rhoU[1]=-U[i][1].rhoU[1]; }
        for (int i=0;i<mesh_.Nx();i++) {
            PrimitiveVars f; f.p=0.5e6;f.T=300;f.U[0]=1;f.Y[0]=0;f.Y[1]=1;f.U[1]=f.U[2]=0;
            StateConversion::primitiveToConserved(f,U[i][mesh_.Ny()-1]);
        }
    }

    void writeVTK(const string& fn) {
        ofstream f(fn);
        f<<"# vtk DataFile Version 3.0\n2D O2\nASCII\nDATASET RECTILINEAR_GRID\n";
        f<<"DIMENSIONS "<<mesh_.Nx()<<" "<<mesh_.Ny()<<" 1\n";
        f<<"X_COORDINATES "<<mesh_.Nx()<<" float\n";
        for (int i=0;i<mesh_.Nx();i++) f<<mesh_(i,0).x*1000<<" ";
        f<<"\nY_COORDINATES "<<mesh_.Ny()<<" float\n";
        for (int j=0;j<mesh_.Ny();j++) f<<mesh_(0,j).y*1000<<" ";
        f<<"\nZ_COORDINATES 1 float\n0\n";
        f<<"CELL_DATA "<<mesh_.Nx()*mesh_.Ny()<<"\n";
        f<<"SCALARS density float\nLOOKUP_TABLE default\n";
        for (int j=0;j<mesh_.Ny();j++) { for (int i=0;i<mesh_.Nx();i++) {
            RealFluidState s;s.updatePTY(W_[i][j].p,W_[i][j].T,W_[i][j].Y); f<<s.rho()<<"\n";
        }}
        f<<"SCALARS pressure float\nLOOKUP_TABLE default\n";
        for (int j=0;j<mesh_.Ny();j++) for (int i=0;i<mesh_.Nx();i++) f<<W_[i][j].p<<"\n";
        f<<"SCALARS temperature float\nLOOKUP_TABLE default\n";
        for (int j=0;j<mesh_.Ny();j++) for (int i=0;i<mesh_.Nx();i++) f<<W_[i][j].T<<"\n";
        f<<"SCALARS Y_O2 float\nLOOKUP_TABLE default\n";
        for (int j=0;j<mesh_.Ny();j++) for (int i=0;i<mesh_.Nx();i++) f<<W_[i][j].Y[0]<<"\n";
        f<<"VECTORS velocity float\n";
        for (int j=0;j<mesh_.Ny();j++) for (int i=0;i<mesh_.Nx();i++) f<<W_[i][j].U[0]<<" "<<W_[i][j].U[1]<<" 0\n";
        f.close();
    }

    void run(double te, double CFL, int wi) {
        cout<<"Run: t_end="<<te<<"s, CFL="<<CFL<<"\n";
        while (t_<te) {
            double dt=computeTimestep(CFL); step(dt);
            if (iter_%100==0) cout<<"  iter="<<iter_<<", t="<<t_<<"s\n";
            if (iter_%wi==0) { ostringstream o; o<<"out_"<<iter_<<".vtk"; writeVTK(o.str()); }
        }
        cout<<"Done: "<<iter_<<" iters\n"; writeVTK("final.vtk");
    }
};

int main() {
    cout<<"2D Transcritical O2 Injection\n";
    Solver2D s(100,50,0.020,0.010);
    s.initialize();
    s.run(0.005, 0.3, 500);
    return 0;
}
