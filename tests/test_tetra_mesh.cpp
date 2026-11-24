#include "TetraMesh.H"
#include <fstream>
#include <iostream>
#include <stdexcept>
int main(int argc,char**argv)
{
    if(argc!=2){std::cerr<<"usage: test_tetra_mesh mesh.rfamtet\n";return 2;}
    try{Foam::TetraMesh m;m.read(argv[1]);if(m.cellCount()==0||m.faceCount()==0) return 1;std::size_t boundary=0;for(std::size_t i=0;i<m.faceCount();++i)boundary+=m.face(i).neighbour<0;if(boundary==0) return 1;std::cout<<"TETRA_MESH_PASS points="<<m.pointCount()<<" cells="<<m.cellCount()<<" faces="<<m.faceCount()<<" boundary="<<boundary<<"\n";return 0;}
    catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
