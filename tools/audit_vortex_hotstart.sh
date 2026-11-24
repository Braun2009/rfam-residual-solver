#!/bin/bash
# ============================================================================
# RFAM Vortex hot-start 成果核验脚本 (2026-08-03)
#
# 用法（WSL 内）:
#   bash /mnt/f/realFluidAllMach/tools/audit_vortex_hotstart.sh [--full]
#
# 核验内容:
#   L1 参考数据完整性（SHA256）
#   L2 核心测试（GERG 热力学链）
#   L3 输运/湍流/restart 组件测试
#   L4 hot-start 生命周期（checkpoint -> restart -> 续算）
#   L5 全量 ctest（--full）
# ============================================================================
set -u
cd /mnt/f/realFluidAllMach || exit 1

FULL=0
[ "${1:-}" = "--full" ] && FULL=1

PASS=0; FAIL=0; SKIP=0
ok()   { PASS=$((PASS+1)); echo "  [PASS] $1"; }
bad()  { FAIL=$((FAIL+1)); echo "  [FAIL] $1"; }
skip() { SKIP=$((SKIP+1)); echo "  [SKIP] $1"; }

echo "============================================================"
echo " RFAM Vortex hot-start 核验  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo " git: $(git log --oneline -1 2>/dev/null)"
echo "============================================================"

# ---------------------------------------------------------------- L1 数据
echo ""
echo "--- L1: 参考数据 SHA256 ---"
GERG_SHA=$(sha256sum tests/reference/gerg_residual_table_v3.csv | awk '{print $1}')
if [ "$GERG_SHA" = "2481ea116e362ff01ec503bfe47b821664d5b747571809f4ac4b32eec3014f5b" ]; then
    ok "GERG residual 表 SHA256 ($GERG_SHA)"
else
    bad "GERG residual 表 SHA256 (got $GERG_SHA)"
fi
CRYO_SHA=$(sha256sum tests/reference/cryogenic_transport_runtime_coolprop72_v1.csv | awk '{print $1}')
if [ "$CRYO_SHA" = "2bf00a17b4c35eb242e28548e1d4673e8d9c3f37412f47e690c9054cab8c3f5c" ]; then
    ok "低温输运表 SHA256 ($CRYO_SHA)"
else
    bad "低温输运表 SHA256 (got $CRYO_SHA)"
fi

# ---------------------------------------------------------------- L2 核心测试
echo ""
echo "--- L2: GERG 热力学链 ---"
CORE_TESTS=(test_gerg_residual test_realfluid_gerg test_hybrid_helmholtz \
            test_hybrid_bubble test_pr_thermo_identities test_pr_tp_flash \
            test_pr_ph_flash test_pr_rho_e_flash test_pr_phase_caloric)
mkdir -p build-audit
# 复用已配置 SUNDIALS 的构建树（本机 /mnt/f/opt/sundials）
SUNDIALS_FLAG=
[ -d /mnt/f/opt/sundials ] && SUNDIALS_FLAG="-DSUNDIALS_ROOT=/mnt/f/opt/sundials"
cmake -S . -B build-audit -DCMAKE_BUILD_TYPE=Release $SUNDIALS_FLAG > /dev/null 2>&1
cmake --build build-audit --parallel 16 --target "${CORE_TESTS[@]}" > /dev/null 2>&1
for t in "${CORE_TESTS[@]}"; do
    if [ -x "build-audit/$t" ]; then
        if (cd build-audit && ctest -R "^$t$" > /dev/null 2>&1); then
            ok "$t"
        else
            bad "$t"
        fi
    else
        skip "$t (构建失败)"
    fi
done

# ---------------------------------------------------------------- L3 组件
echo ""
echo "--- L3: 输运/湍流/restart 组件 ---"
COMP_TESTS=(test_cryogenic_transport_table test_cryogenic_backend \
            test_droplet_evaporation test_wale test_sector_turbulence \
            test_sector_checkpoint test_sector_restart_equivalence_warm \
            test_sector_injection test_step_transaction test_restart_publication)
cmake --build build-audit --parallel 16 --target "${COMP_TESTS[@]}" > /dev/null 2>&1
for t in "${COMP_TESTS[@]}"; do
    if [ -x "build-audit/$t" ]; then
        if (cd build-audit && ctest -R "^$t$" > /dev/null 2>&1); then
            ok "$t"
        else
            bad "$t"
        fi
    else
        skip "$t (构建失败)"
    fi
done

# ---------------------------------------------------------------- L4 hot-start
echo ""
echo "--- L4: hot-start 生命周期 (3D 扇区诊断) ---"
cmake --build build-audit --parallel 16 --target rfam_v20_sector_smoke > /dev/null 2>&1
if [ -x build-audit/rfam_v20_sector_smoke ]; then
    rm -f /tmp/hs_audit.chk
    OUT1=$(timeout 300 ./build-audit/rfam_v20_sector_smoke \
        --reacting --allow-single-phase-reacting-diagnostic \
        --steps 60 --checkpoint-every 30 --checkpoint /tmp/hs_audit.chk 2>&1)
    if echo "$OUT1" | grep -q "checkpoint /tmp/hs_audit.chk iter=60"; then
        ok "阶段1: 60 步 + checkpoint (iter=60)"
    else
        bad "阶段1: 60 步 + checkpoint"
        echo "$OUT1" | grep -E "ERROR|error|what" | head -3
    fi
    OUT2=$(timeout 300 ./build-audit/rfam_v20_sector_smoke \
        --reacting --allow-single-phase-reacting-diagnostic \
        --steps 90 --restart /tmp/hs_audit.chk 2>&1)
    if echo "$OUT2" | grep -q "Restarted /tmp/hs_audit.chk at iter=60"; then
        ok "阶段2: restart 恢复 (iter=60)"
    else
        bad "阶段2: restart 恢复"
        echo "$OUT2" | grep -E "what|ERROR" | head -3
    fi
    if echo "$OUT2" | grep -q "SECTOR_3D_STEP_PASS"; then
        ok "阶段3: 续算至 90 步 SECTOR_3D_STEP_PASS"
    else
        bad "阶段3: 续算"
    fi
else
    skip "rfam_v20_sector_smoke (构建失败)"
fi

# ---------------------------------------------------------------- L5 全量
if [ "$FULL" = "1" ]; then
    echo ""
    echo "--- L5: 全量 ctest ---"
    cmake --build build-audit --parallel 16 > /dev/null 2>&1
    (cd build-audit && ctest 2>&1 | tail -3)
fi

echo ""
echo "============================================================"
echo " 核验结果: PASS=$PASS FAIL=$FAIL SKIP=$SKIP"
echo "============================================================"
[ "$FAIL" = "0" ] || exit 1
