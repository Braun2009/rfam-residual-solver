#ifndef SectorMeshMetrics_H
#define SectorMeshMetrics_H

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Foam
{

inline std::string vortexProfilePath()
{
    if(const char* configured=std::getenv("VORTEX_PROFILE")) return configured;
    std::ifstream local("config/chamber_profile_3pintle.json");
    if(local.good()) return "config/chamber_profile_3pintle.json";
#ifdef RFAM_SOURCE_DIR
    const std::string sourceProfile=std::string(RFAM_SOURCE_DIR)
        +"/config/chamber_profile_3pintle.json";
    std::ifstream source(sourceProfile);
    if(source.good()) return sourceProfile;
#endif
    return "/mnt/d/Vortex/config/chamber_profile_3pintle.json";
}

struct FaceAreaVector
{
    double x, y, z;
};

class SectorMeshMetrics
{
    int nx_, nr_, nt_;
    double length_, sectorAngle_;
    std::vector<double> xFaces_, wallFaces_, etaFaces_, thetaFaces_;

    static void readProfile(const std::string& path,
                            std::vector<double>& x,
                            std::vector<double>& r)
    {
        std::ifstream file(path);
        if (!file) throw std::runtime_error("cannot open Vortex profile: " + path);
        const std::string content((std::istreambuf_iterator<char>(file)),
                                  std::istreambuf_iterator<char>());
        size_t pos = 0;
        while ((pos = content.find("\"x_mm\":", pos)) != std::string::npos)
        {
            const double xv = std::stod(content.substr(pos + 7)) * 1e-3;
            const size_t rp = content.find("\"r_mm\":", pos);
            if (rp == std::string::npos) break;
            const double rv = std::stod(content.substr(rp + 7)) * 1e-3;
            x.push_back(xv);
            r.push_back(rv);
            pos = rp + 7;
        }
        if (x.size() < 2) throw std::runtime_error("Vortex profile has no usable points");
    }

    static double interpolate(const std::vector<double>& x,
                              const std::vector<double>& y,
                              double value)
    {
        if (value <= x.front()) return y.front();
        if (value >= x.back()) return y.back();
        const auto upper = std::lower_bound(x.begin(), x.end(), value);
        const size_t hi = static_cast<size_t>(upper - x.begin());
        const size_t lo = hi - 1;
        const double fraction = (value - x[lo]) / (x[hi] - x[lo]);
        return y[lo] + fraction * (y[hi] - y[lo]);
    }

public:
    SectorMeshMetrics(const std::string& profilePath,
                      int nx, int nr, int nt, double sectorDegrees)
        : nx_(nx), nr_(nr), nt_(nt), sectorAngle_(sectorDegrees * M_PI / 180.0)
    {
        if (nx < 4 || nr < 4 || nt < 4) throw std::invalid_argument("sector mesh dimensions too small");
        std::vector<double> profileX, profileR;
        readProfile(profilePath, profileX, profileR);
        length_ = profileX.back();

        xFaces_.resize(nx_ + 1);
        wallFaces_.resize(nx_ + 1);
        etaFaces_.resize(nr_ + 1);
        thetaFaces_.resize(nt_ + 1);

        // Structured solver baseline. CAD nozzle geometry is reconstructed
        // separately; this mesh has no embedded solids or local refinement.
        // Physical faces are 1..N-1. The two outer faces extend one full
        // ghost cell beyond each boundary (signed eta mirrors the axis).
        // Do not clamp ghost faces onto physical boundaries: wall/axis
        // operators need distinct ghost centres and nonzero spacing.
        for (int i = 0; i <= nx_; ++i)
        {
            xFaces_[i] = length_ * (i-1) / (nx_-2);
            wallFaces_[i] = interpolate(profileX, profileR, xFaces_[i]);
        }
        for (int j = 0; j <= nr_; ++j)
            etaFaces_[j] = static_cast<double>(j-1) / (nr_-2);
        // k=0 and k=nt-1 are periodic halo cells. Interior faces 1 and
        // nt-1 therefore delimit the complete physical sector.
        for (int k = 0; k <= nt_; ++k)
            thetaFaces_[k] = sectorAngle_ * (k-1) / (nt_-2);
    }

    int Nx() const { return nx_; }
    int Nr() const { return nr_; }
    int Ntheta() const { return nt_; }
    double length() const { return length_; }
    double sectorAngle() const { return sectorAngle_; }
    double xFace(int i) const { return xFaces_[i]; }
    double wallFace(int i) const { return wallFaces_[i]; }
    double etaFace(int j) const { return etaFaces_[j]; }
    double thetaFace(int k) const { return thetaFaces_[k]; }
    double thetaCenter(int k) const
    { return 0.5 * (thetaFaces_[k] + thetaFaces_[k+1]); }

    double wallCenter(int i) const { return 0.5 * (wallFaces_[i] + wallFaces_[i+1]); }
    double radiusCenter(int i, int j) const
    {
        return wallCenter(i) * 0.5 * (etaFaces_[j] + etaFaces_[j+1]);
    }

    FaceAreaVector cellCenter(int i, int j, int k) const
    {
        const double theta = thetaCenter(k);
        const double radius = radiusCenter(i,j);
        return {0.5*(xFaces_[i]+xFaces_[i+1]),
                radius*std::cos(theta), radius*std::sin(theta)};
    }

    double volume(int i, int j, int k) const
    {
        (void)k;
        const double dx = xFaces_[i+1] - xFaces_[i];
        const double dtheta = thetaFaces_[1] - thetaFaces_[0];
        const double rw = wallCenter(i);
        const double r0 = rw * etaFaces_[j];
        const double r1 = rw * etaFaces_[j+1];
        return 0.5 * (r1*r1 - r0*r0) * dtheta * dx;
    }

    double axialArea(int faceI, int j, int k) const
    {
        (void)k;
        const double dtheta = thetaFaces_[1] - thetaFaces_[0];
        const double r0 = wallFaces_[faceI] * etaFaces_[j];
        const double r1 = wallFaces_[faceI] * etaFaces_[j+1];
        return 0.5 * (r1*r1 - r0*r0) * dtheta;
    }

    double radialArea(int i, int faceJ, int k) const
    {
        (void)k;
        const double dx = xFaces_[i+1] - xFaces_[i];
        const double dtheta = thetaFaces_[1] - thetaFaces_[0];
        return dx * wallCenter(i) * etaFaces_[faceJ] * dtheta;
    }

    double azimuthalArea(int i, int j, int faceK) const
    {
        (void)faceK;
        const double dx = xFaces_[i+1] - xFaces_[i];
        const double rw = wallCenter(i);
        return dx * rw * (etaFaces_[j+1] - etaFaces_[j]);
    }

    double totalVolume() const
    {
        double total = 0.0;
        for (int i = 1; i < nx_-1; ++i)
        for (int j = 1; j < nr_-1; ++j)
        for (int k = 1; k < nt_-1; ++k)
            total += volume(i,j,k);
        return total;
    }

    FaceAreaVector axialAreaVector(int faceI, int j, int k, double sign) const
    {
        const double theta0 = thetaFaces_[k];
        const double theta1 = thetaFaces_[k+1];
        const double rw = wallFaces_[faceI];
        const double area = 0.5 * rw * rw
            * (etaFaces_[j+1]*etaFaces_[j+1] - etaFaces_[j]*etaFaces_[j])
            * (theta1 - theta0);
        return {sign * area, 0.0, 0.0};
    }

    // Area vector for eta=constant surface, positive sign points toward
    // increasing eta. The axial component exactly accounts for wall slope.
    FaceAreaVector radialAreaVector(int i, int faceJ, int k, double sign) const
    {
        const double eta = etaFaces_[faceJ];
        const double theta0 = thetaFaces_[k];
        const double theta1 = thetaFaces_[k+1];
        const double thetaMid = 0.5 * (theta0 + theta1);
        const double rLeft = wallFaces_[i];
        const double rRight = wallFaces_[i+1];
        const double dx = xFaces_[i+1] - xFaces_[i];
        const double axial = -0.5 * eta * eta * (rRight*rRight-rLeft*rLeft)
                             * (theta1-theta0);
        const double radial = eta * 0.5 * (rLeft+rRight) * dx
                              * 2.0 * std::sin(0.5*(theta1-theta0));
        return {sign*axial, sign*radial*std::cos(thetaMid),
                sign*radial*std::sin(thetaMid)};
    }

    FaceAreaVector azimuthalAreaVector(int i, int j, int faceK, double sign) const
    {
        const double theta = thetaFaces_[faceK];
        const double dx = xFaces_[i+1] - xFaces_[i];
        const double radialWidth = 0.5 * (wallFaces_[i]+wallFaces_[i+1])
                                 * (etaFaces_[j+1]-etaFaces_[j]);
        const double area = dx * radialWidth;
        return {0.0, sign*(-std::sin(theta))*area,
                sign*std::cos(theta)*area};
    }
};

} // End namespace Foam

#endif
