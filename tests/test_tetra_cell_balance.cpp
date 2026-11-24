// Face-table invariants and per-cell area balance for the tetra loader.
//
// Checks that cellFaces(c) really references faces owned by or neighbouring c,
// then verifies sum(A_outward)=0 for every cell. A failure here means the face
// table or the orientation rule is wrong, not that the physics is hard.
#include "TetraMesh.H"
#include <cmath>
#include <iostream>
#include <vector>

using namespace Foam;

int main(int argc,char**argv)
{
    if(argc!=2){std::cerr<<"usage: test_tetra_face_table <mesh.rfamtet>\n";return 2;}
    try
    {
        TetraMesh mesh;mesh.read(argv[1]);
        std::size_t badLink=0,badLocal=0;
        for(std::size_t c=0;c<mesh.cellCount();++c)
        {
            const auto& faces=mesh.cellFaces(c);
            if(faces.size()!=4){++badLocal;continue;}
            for(std::int32_t fi:faces)
            {
                const auto& f=mesh.face(fi);
                if(f.owner!=(std::int32_t)c&&f.neighbour!=(std::int32_t)c)++badLink;
            }
        }
        std::cout<<"table: cells="<<mesh.cellCount()<<" badLocal="<<badLocal
                 <<" facesNotIncidentToCell="<<badLink<<"\n";
        double worstRel=0.0;std::size_t worstCell=0;
        for(std::size_t c=0;c<mesh.cellCount();++c)
        {
            TetraPoint sum{0,0,0};double total=0.0;
            for(std::int32_t fi:mesh.cellFaces(c))
            {
                // faceAreaVector is oriented outward from the OWNER, so a cell
                // that owns the face on its neighbour side must negate it.
                auto a=mesh.faceAreaVector(fi);
                if(mesh.face(fi).owner!=(std::int32_t)c){a.x=-a.x;a.y=-a.y;a.z=-a.z;}
                sum.x+=a.x;sum.y+=a.y;sum.z+=a.z;
                total+=std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
            }
            const double imbalance=std::sqrt(sum.x*sum.x+sum.y*sum.y+sum.z*sum.z);
            const double rel=imbalance/std::max(total,1e-30);
            if(rel>worstRel){worstRel=rel;worstCell=c;}
        }
        std::cout<<"worst relative area imbalance="<<worstRel<<" at cell "<<worstCell<<"\n";
        // Detail for the worst cell: nodes, face areas and orientation dots.
        for(std::int32_t fi:mesh.cellFaces(worstCell))
        {
            const auto& f=mesh.face(fi);
            auto a=mesh.faceAreaVector(fi);
            if(f.owner!=(std::int32_t)worstCell){a.x=-a.x;a.y=-a.y;a.z=-a.z;}
            const auto fc=mesh.faceCentroid(fi);
            const auto oc=mesh.cellCentroid(worstCell);
            std::cout<<"  face "<<fi<<" nodes "<<f.nodes[0]<<","<<f.nodes[1]<<","<<f.nodes[2]
                     <<" owner "<<f.owner<<" neigh "<<f.neighbour<<" zone "<<f.boundaryZone
                     <<" |A_outward|="<<std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z)
                     <<" dot="<<(a.x*(fc.x-oc.x)+a.y*(fc.y-oc.y)+a.z*(fc.z-oc.z))<<"\n";
        }
        if(badLink||badLocal){std::cerr<<"face table is inconsistent\n";return 1;}
        if(!(worstRel<1e-9)){std::cerr<<"cell faces are not consistently closed\n";return 1;}
        std::cout<<"TETRA_FACE_TABLE_PASS\n";
        return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
