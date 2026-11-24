#!/bin/bash
#------------------------------------------------------------------------------
# Build script for realFluidAllMach
#------------------------------------------------------------------------------

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  realFluidAllMach Build Script                             ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Check for build type
BUILD_TYPE=${1:-Release}
echo -e "${YELLOW}Build type: ${BUILD_TYPE}${NC}"

# Create build directory
if [ ! -d "build" ]; then
    mkdir build
    echo "Created build/ directory"
fi

cd build

# Configure
echo ""
echo -e "${YELLOW}Configuring CMake...${NC}"
cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ..

if [ $? -ne 0 ]; then
    echo -e "${RED}✗ CMake configuration failed${NC}"
    exit 1
fi

# Build
echo ""
echo -e "${YELLOW}Building (using all cores)...${NC}"
make -j$(nproc)

if [ $? -ne 0 ]; then
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}✓ Build successful${NC}"

# Check what was built
echo ""
echo "Built artifacts:"
echo "  Solver: $(ls realFluidAllMachFoam 2>/dev/null && echo '✓' || echo '✗')"
echo "  Libraries:"
echo "    - librealFluidThermo.a: $(ls librealFluidThermo.a 2>/dev/null && echo '✓' || echo '✗')"
echo "    - librealFluidNumerics.a: $(ls librealFluidNumerics.a 2>/dev/null && echo '✓' || echo '✗')"
echo "    - librealFluidCombustion.a: $(ls librealFluidCombustion.a 2>/dev/null && echo '✓' || echo '✗')"
echo "    - librealFluidCHT.a: $(ls librealFluidCHT.a 2>/dev/null && echo '✓' || echo '✗')"

# Check for tests
if [ -f "test_pr_roots" ]; then
    echo "  Tests: ✓"
else
    echo "  Tests: ✗ (GTest not found)"
fi

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  Build complete. Next steps:                               ║"
echo "║                                                            ║"
echo "║  1. Run tests:     ./scripts/run_tests.sh                 ║"
echo "║  2. Run solver:    ./build/realFluidAllMachFoam            ║"
echo "║  3. Documentation: make doc                                ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
