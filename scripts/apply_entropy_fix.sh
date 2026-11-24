#!/bin/bash
cd /root/rfam-tetra
cp apps/solver_3d_tetra.cpp.backup apps/solver_3d_tetra.cpp

# 1. Add the tolerance parameter to Config struct (after supplyEntropyStrict)
sed -i '/bool supplyEntropyStrict=false;/a\    double supplyEntropyTolerance=1e-7; // relative entropy margin tolerance (allow losses)' apps/solver_3d_tetra.cpp

# 2. Add command-line argument parsing (after --supply-entropy-strict)
sed -i '/else if(a=="--supply-entropy-strict")c.supplyEntropyStrict=true;/a\        else if(a=="--supply-entropy-tolerance")c.supplyEntropyTolerance=std::atof(need("--supply-entropy-tolerance"));' apps/solver_3d_tetra.cpp

# 3. Replace hardcoded 1e-7 with cfg.supplyEntropyTolerance (three occurrences)
sed -i 's/margin<-1e-7\*std::max(std::abs(zone.supplyEntropy),1\.)/margin<-cfg.supplyEntropyTolerance*std::max(std::abs(zone.supplyEntropy),1.)/g' apps/solver_3d_tetra.cpp
sed -i 's/zone\.entropyMargin<-1e-7\*std::max(std::abs(zone.supplyEntropy),1\.)/zone.entropyMargin<-cfg.supplyEntropyTolerance*std::max(std::abs(zone.supplyEntropy),1.)/g' apps/solver_3d_tetra.cpp

# 4. Add help text
sed -i '/the supply entropy (default: report only)/a\        "  --supply-entropy-tolerance R  relative entropy margin tolerance (default 1e-7);\\n"\n        "                         allows s_face - s_supply >= -R*max(|s_supply|, 1) for losses\\n"' apps/solver_3d_tetra.cpp

echo "Patch applied successfully"
