#include "TetraMesh.H"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>

namespace Foam
{
namespace
{
using Key=std::array<std::int32_t,3>;
}
void TetraMesh::read(const std::string& path)
{
 std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("cannot open tetra mesh: "+path);
 char magic[8];f.read(magic,8);if(std::string(magic,8)!="RFAMTET1")throw std::runtime_error("invalid RFAM tetra mesh magic");
 std::uint64_t np=0,nc=0,nf=0;f.read((char*)&np,8);f.read((char*)&nc,8);f.read((char*)&nf,8);
 points_.resize(np);cells_.resize(nc);volumes_.resize(nc);faces_.resize(nf);cellFaces_.assign(nc,{});
 f.read((char*)points_.data(),np*sizeof(TetraPoint));f.read((char*)cells_.data(),nc*sizeof(cells_[0]));
 for(auto& x:volumes_)f.read((char*)&x,sizeof(x));
 f.read((char*)faces_.data(),nf*sizeof(TetraFace));
 if(!f)throw std::runtime_error("truncated tetra mesh");
 for(std::size_t i=0;i<nf;++i){if(faces_[i].owner<0||faces_[i].owner>=(int)nc)throw std::runtime_error("invalid tetra face owner");cellFaces_[faces_[i].owner].push_back((int)i);if(faces_[i].neighbour>=0)cellFaces_[faces_[i].neighbour].push_back((int)i);}
 validate();
}
TetraPoint TetraMesh::cellCentroid(std::size_t i) const
{
 const auto& c=cells_.at(i); TetraPoint r{};
 for(int n:c){r.x+=points_[n].x;r.y+=points_[n].y;r.z+=points_[n].z;}
 return {r.x/4.0,r.y/4.0,r.z/4.0};
}
TetraPoint TetraMesh::faceCentroid(std::size_t i) const
{
 const auto& f=faces_.at(i);TetraPoint r{};
 for(int n:f.nodes){r.x+=points_[n].x;r.y+=points_[n].y;r.z+=points_[n].z;}
 return {r.x/3.0,r.y/3.0,r.z/3.0};
}
TetraPoint TetraMesh::faceAreaVector(std::size_t i) const
{
 const auto& f=faces_.at(i);const auto&a=points_[f.nodes[0]];const auto&b=points_[f.nodes[1]];const auto&c=points_[f.nodes[2]];
 TetraPoint v{(b.y-a.y)*(c.z-a.z)-(b.z-a.z)*(c.y-a.y),(b.z-a.z)*(c.x-a.x)-(b.x-a.x)*(c.z-a.z),(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x)};
 const auto fc=faceCentroid(i);const auto oc=cellCentroid(f.owner);if((v.x*(fc.x-oc.x)+v.y*(fc.y-oc.y)+v.z*(fc.z-oc.z))<0){v.x=-v.x;v.y=-v.y;v.z=-v.z;}return {0.5*v.x,0.5*v.y,0.5*v.z};
}
void TetraMesh::validate() const
{
 if(points_.empty()||cells_.empty()||faces_.empty())throw std::runtime_error("empty tetra mesh");
 for(std::size_t i=0;i<cells_.size();++i){const auto&c=cells_[i];for(int n:c)if(n<0||(std::size_t)n>=points_.size())throw std::runtime_error("tetra node out of range");if(!std::isfinite(volumes_[i])||!(volumes_[i]>0))throw std::runtime_error("nonpositive tetra volume");if(cellFaces_[i].size()!=4)throw std::runtime_error("tetra does not have four faces");}
 for(const auto& x:faces_){for(auto n:x.nodes)if(n<0||(std::size_t)n>=points_.size())throw std::runtime_error("face node out of range");if(x.neighbour>=0&&x.neighbour==x.owner)throw std::runtime_error("self-neighbour face");}
}
}
