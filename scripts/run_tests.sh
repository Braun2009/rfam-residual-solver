#!/bin/bash
#------------------------------------------------------------------------------
# Run tests for realFluidAllMach
#------------------------------------------------------------------------------

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  realFluidAllMach Test Suite                               ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Check if build directory exists
if [ ! -d "build" ]; then
    echo -e "${RED}✗ Build directory not found. Run ./scripts/build.sh first.${NC}"
    exit 1
fi

cd build

# Check if tests were built
if [ ! -f "test_pr_roots" ]; then
    echo -e "${YELLOW}⚠ Tests not built (GTest not found). Skipping.${NC}"
    exit 0
fi

# Run all tests with CTest
echo -e "${YELLOW}Running test suite...${NC}"
echo ""

ctest --output-on-failure --verbose

RESULT=$?

echo ""
if [ $RESULT -eq 0 ]; then
    echo -e "${GREEN}╔════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║  ✓ ALL TESTS PASSED                                        ║${NC}"
    echo -e "${GREEN}║                                                            ║${NC}"
    echo -e "${GREEN}║  Thermo kernel is thermodynamically consistent.            ║${NC}"
    echo -e "${GREEN}║  Ready to proceed to Stage 1 (conservative solver).        ║${NC}"
    echo -e "${GREEN}╚════════════════════════════════════════════════════════════╝${NC}"
    exit 0
else
    echo -e "${RED}╔════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║  ✗ TESTS FAILED                                            ║${NC}"
    echo -e "${RED}║                                                            ║${NC}"
    echo -e "${RED}║  DO NOT PROCEED. Fix failures first.                       ║${NC}"
    echo -e "${RED}╚════════════════════════════════════════════════════════════╝${NC}"
    exit 1
fi
