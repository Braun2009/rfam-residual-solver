#ifndef TetraMesh_H
#define TetraMesh_H

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace Foam
{

struct TetraPoint { double x=0.0,y=0.0,z=0.0; };
struct TetraFace
{
    std::array<std::int32_t,3> nodes{};
    std::int32_t owner=-1;
    std::int32_t neighbour=-1;
    std::int32_t boundaryZone=0;
};

// Unstructured linear-tetrahedron mesh used by the future RFAM finite-volume
// path. Boundary faces have neighbour=-1; internal faces have two owners.
// The loader is deliberately independent of STL/VTU viewers.
class TetraMesh
{
    std::vector<TetraPoint> points_;
    std::vector<std::array<std::int32_t,4>> cells_;
    std::vector<double> volumes_;
    std::vector<TetraFace> faces_;
    std::vector<std::vector<std::int32_t>> cellFaces_;

public:
    void read(const std::string& path);
    void validate() const;

    std::size_t pointCount() const { return points_.size(); }
    std::size_t cellCount() const { return cells_.size(); }
    std::size_t faceCount() const { return faces_.size(); }
    const TetraPoint& point(std::size_t i) const { return points_.at(i); }
    const std::array<std::int32_t,4>& cell(std::size_t i) const { return cells_.at(i); }
    double volume(std::size_t i) const { return volumes_.at(i); }
    const TetraFace& face(std::size_t i) const { return faces_.at(i); }
    const std::vector<std::int32_t>& cellFaces(std::size_t i) const { return cellFaces_.at(i); }

    TetraPoint cellCentroid(std::size_t i) const;
    TetraPoint faceCentroid(std::size_t i) const;
    TetraPoint faceAreaVector(std::size_t i) const;

};

} // namespace Foam
#endif
