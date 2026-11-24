#include "TetraMesh.H"
#include <cmath>
#include <iostream>

// Total volume and face-area sanity for the audited nozzle mesh.
// The expected value is the SI volume from the mesh audit (mm^3 -> m^3).
int main(int argc, char** argv)
{
    if (argc != 2) return 2;
    try
    {
        Foam::TetraMesh m; m.read(argv[1]);
        double sum = 0.0;
        for (std::size_t i = 0; i < m.cellCount(); ++i)
        {
            if (!(m.volume(i) > 0)) return 1;
            sum += m.volume(i);
        }
        double area = 0.0;
        for (std::size_t i = 0; i < m.faceCount(); ++i)
        {
            const auto a = m.faceAreaVector(i);
            const double q = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
            if (!(q > 0 && std::isfinite(q))) return 1;
            area += q;
        }
        const double expected = 2.958603075278e-3;
        const double diff = std::abs(sum - expected);
        std::cout << "total_volume_m3=" << sum << " diff=" << diff << std::endl;
        if (diff > 1e-10) return 1;
        std::cout << "TETRA_GEOMETRY_PASS volume=" << sum
                  << " area_vector_sum=" << area << std::endl;
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
