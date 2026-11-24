#include "FPVTable.H"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

namespace Foam
{

FPVTable::FPVTable()
:
    nZ_(0),
    nC_(0),
    nh_(0)
{
}

bool FPVTable::loadFromCSV(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "ERROR: Cannot open FPV table: " << filename << std::endl;
        return false;
    }
    
    std::cout << "Loading FPV table from: " << filename << std::endl;
    
    // Read header
    std::string header;
    std::getline(file, header);
    
    // Parse species names from header
    std::stringstream ss(header);
    std::string token;
    species_names_.clear();
    
    int col = 0;
    while (std::getline(ss, token, ','))
    {
        if (col >= 5 && col < 18)  // Y_species columns (5-17)
        {
            species_names_.push_back(token);
        }
        col++;
    }
    
    std::cout << "  Species: " << species_names_.size() << std::endl;
    
    // Read data
    std::vector<double> Z_raw, C_raw, h_raw;
    std::vector<double> T_raw, rho_raw, omega_raw;
    std::vector<std::vector<double>> Y_raw(species_names_.size());
    
    std::string line;
    int lineCount = 0;
    
    while (std::getline(file, line))
    {
        std::stringstream ls(line);
        std::string val;
        
        double Z, C, h, T, rho;
        std::vector<double> Y(species_names_.size());
        double omega;
        
        // Parse: Z, C, h, T, rho, Y[0-12], omega_dot
        std::getline(ls, val, ','); Z = std::stod(val);
        std::getline(ls, val, ','); C = std::stod(val);
        std::getline(ls, val, ','); h = std::stod(val);
        std::getline(ls, val, ','); T = std::stod(val);
        std::getline(ls, val, ','); rho = std::stod(val);
        
        for (size_t i = 0; i < species_names_.size(); i++)
        {
            std::getline(ls, val, ',');
            Y[i] = std::stod(val);
        }
        
        std::getline(ls, val, ',');
        omega = std::stod(val);
        
        Z_raw.push_back(Z);
        C_raw.push_back(C);
        h_raw.push_back(h);
        T_raw.push_back(T);
        rho_raw.push_back(rho);
        omega_raw.push_back(omega);
        
        for (size_t i = 0; i < species_names_.size(); i++)
        {
            Y_raw[i].push_back(Y[i]);
        }
        
        lineCount++;
    }
    
    file.close();
    
    std::cout << "  Read " << lineCount << " data points" << std::endl;
    
    // Extract unique grid values
    auto uniqueSorted = [](const std::vector<double>& v) {
        std::vector<double> u = v;
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), u.end()), u.end());
        return u;
    };
    
    Z_grid_ = uniqueSorted(Z_raw);
    C_grid_ = uniqueSorted(C_raw);
    h_grid_ = uniqueSorted(h_raw);
    
    nZ_ = Z_grid_.size();
    nC_ = C_grid_.size();
    nh_ = h_grid_.size();
    
    std::cout << "  Grid: Z=" << nZ_ << ", C=" << nC_ << ", h=" << nh_ << std::endl;
    
    // Reorganize data into 3D grid
    int totalPoints = nZ_ * nC_ * nh_;
    T_.resize(totalPoints, 0.0);
    rho_.resize(totalPoints, 0.0);
    omega_dot_.resize(totalPoints, 0.0);
    Y_.resize(species_names_.size());
    for (auto& Yi : Y_) Yi.resize(totalPoints, 0.0);
    
    // Fill grid (assumes data is ordered)
    for (int idx = 0; idx < lineCount; idx++)
    {
        T_[idx] = T_raw[idx];
        rho_[idx] = rho_raw[idx];
        omega_dot_[idx] = omega_raw[idx];
        for (size_t i = 0; i < species_names_.size(); i++)
        {
            Y_[i][idx] = Y_raw[i][idx];
        }
    }
    
    std::cout << "FPV table loaded successfully!" << std::endl;
    return true;
}

void FPVTable::findBracket(
    const std::vector<double>& grid,
    double value,
    int& i0, int& i1,
    double& weight
) const
{
    // Clamp to bounds
    if (value <= grid.front())
    {
        i0 = 0;
        i1 = 0;
        weight = 0.0;
        return;
    }
    if (value >= grid.back())
    {
        i0 = grid.size() - 1;
        i1 = i0;
        weight = 0.0;
        return;
    }
    
    // Binary search
    auto it = std::lower_bound(grid.begin(), grid.end(), value);
    i1 = std::distance(grid.begin(), it);
    i0 = i1 - 1;
    
    weight = (value - grid[i0]) / (grid[i1] - grid[i0]);
}

double FPVTable::interpolate3D(
    const std::vector<double>& data,
    double Z, double C, double h
) const
{
    int iZ0, iZ1, iC0, iC1, ih0, ih1;
    double wZ, wC, wh;
    
    findBracket(Z_grid_, Z, iZ0, iZ1, wZ);
    findBracket(C_grid_, C, iC0, iC1, wC);
    findBracket(h_grid_, h, ih0, ih1, wh);
    
    // Trilinear interpolation
    double c000 = data[index(iZ0, iC0, ih0)];
    double c001 = data[index(iZ0, iC0, ih1)];
    double c010 = data[index(iZ0, iC1, ih0)];
    double c011 = data[index(iZ0, iC1, ih1)];
    double c100 = data[index(iZ1, iC0, ih0)];
    double c101 = data[index(iZ1, iC0, ih1)];
    double c110 = data[index(iZ1, iC1, ih0)];
    double c111 = data[index(iZ1, iC1, ih1)];
    
    double c00 = c000 * (1-wh) + c001 * wh;
    double c01 = c010 * (1-wh) + c011 * wh;
    double c10 = c100 * (1-wh) + c101 * wh;
    double c11 = c110 * (1-wh) + c111 * wh;
    
    double c0 = c00 * (1-wC) + c01 * wC;
    double c1 = c10 * (1-wC) + c11 * wC;
    
    return c0 * (1-wZ) + c1 * wZ;
}

double FPVTable::getOmegaDot(double Z, double C, double h) const
{
    return interpolate3D(omega_dot_, Z, C, h);
}

void FPVTable::print() const
{
    std::cout << "FPV Table Summary:" << std::endl;
    std::cout << "  Z: [" << Z_min() << ", " << Z_max() << "], " << nZ_ << " points" << std::endl;
    std::cout << "  C: [" << C_min() << ", " << C_max() << "], " << nC_ << " points" << std::endl;
    std::cout << "  h: [" << h_min() << ", " << h_max() << "], " << nh_ << " points" << std::endl;
    std::cout << "  Species: " << species_names_.size() << std::endl;
}

} // End namespace Foam
