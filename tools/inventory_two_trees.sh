#!/usr/bin/env bash
# Capability inventory between the two RFAM trees that must become one.
OLD=/root/rfam-tetra
NEW=/root/realFluidAllMach

echo "===== 1. driver diff (old solver_3d_tetra_react.cpp vs new rfam_work copy)"
if [ -f "$OLD/solver_3d_tetra_react.cpp" ] && [ -f "$NEW/rfam_work/solver_3d_tetra_react.cpp" ]; then
  a="$OLD/solver_3d_tetra_react.cpp"; b="$NEW/rfam_work/solver_3d_tetra_react.cpp"
  echo "old_lines=$(wc -l < "$a") new_lines=$(wc -l < "$b") differing=$(diff "$a" "$b" | grep -c '^[<>]')"
  echo "--- hunk headers touching supply / GERG / PSY / outlet"
  diff "$a" "$b" | grep -E '^[0-9].*(initializePSY|gerg|GERG|supply|outlet)' | head -15
  echo "--- flags only in old:"
  comm -23 <(grep -oE 'a=="--[A-Za-z-]+"' "$a" | sed 's/a=="//;s/"//' | sort -u) \
           <(grep -oE 'a=="--[A-Za-z-]+"' "$b" | sed 's/a=="//;s/"//' | sort -u) | tr '\n' ' '
  echo
  echo "--- flags only in new:"
  comm -13 <(grep -oE 'a=="--[A-Za-z-]+"' "$a" | sed 's/a=="//;s/"//' | sort -u) \
           <(grep -oE 'a=="--[A-Za-z-]+"' "$b" | sed 's/a=="//;s/"//' | sort -u) | tr '\n' ' '
  echo
else
  echo "driver missing in one tree"
fi

echo
echo "===== 2. source files only in the OLD tree (src/, max 25)"
comm -23 <(cd "$OLD/src" 2>/dev/null && find . -type f \( -name '*.H' -o -name '*.C' -o -name '*.h' -o -name '*.cpp' \) | sort) \
         <(cd "$NEW/src" 2>/dev/null && find . -type f \( -name '*.H' -o -name '*.C' -o -name '*.h' -o -name '*.cpp' \) | sort) | head -25

echo
echo "===== 3. source files only in the NEW tree (src/, max 25)"
comm -13 <(cd "$OLD/src" 2>/dev/null && find . -type f \( -name '*.H' -o -name '*.C' -o -name '*.h' -o -name '*.cpp' \) | sort) \
         <(cd "$NEW/src" 2>/dev/null && find . -type f \( -name '*.H' -o -name '*.C' -o -name '*.h' -o -name '*.cpp' \) | sort) | head -25

echo
echo "===== 4. drivers present in each tree"
echo "--- old tree executables with 'tetra' in the name:"
find "$OLD" -maxdepth 2 -type f -executable -name '*tetra*' 2>/dev/null | head -10
echo "--- new tree executables with 'tetra' in the name:"
find "$NEW" -maxdepth 2 -type f -executable -name '*tetra*' 2>/dev/null | head -10

echo
echo "===== 5. mesh types available"
echo "--- tetra readers:"; grep -rl "RFAMTET1" "$NEW/src" 2>/dev/null | head -5
echo "--- hex/structured mesh readers:"; grep -rl "struct.*Mesh\|class.*Mesh" "$NEW/src/framework" 2>/dev/null | head -8
