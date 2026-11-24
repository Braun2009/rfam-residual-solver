#!/bin/bash
set -e

echo "=== Thread 1: RFAM 编译（WSL） ==="
cd /d/realFluidAllMach
mkdir -p build_wsl && cd build_wsl

# 配置
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-O3 -march=native" \
  > cmake.log 2>&1

echo "CMake 配置完成，开始编译..."
make -j8 > make.log 2>&1 &
MAKE_PID=$!
echo "编译进程 PID: $MAKE_PID"
echo "监控命令: tail -f /d/realFluidAllMach/build_wsl/make.log"

echo ""
echo "=== Thread 2: v4 STL 快速检查（Python） ==="
cd /d/Vortex/cad/injector
timeout 300 python3 << 'EOFPY' > v4_quick_check.log 2>&1 &
import struct, os
with open('injector_v4.STL', 'rb') as f:
    header = f.read(80)
    print(f"Header: {header[:20]}")
    nfaces_bytes = f.read(4)
    nfaces = struct.unpack('<I', nfaces_bytes)[0]
    print(f"Declared faces: {nfaces:,}")
    expected_size = 84 + nfaces * 50
    actual_size = os.path.getsize('injector_v4.STL')
    print(f"Expected size: {expected_size:,} bytes")
    print(f"Actual size: {actual_size:,} bytes")
    print(f"Match: {expected_size == actual_size}")
    
    # 尝试读取前10个面
    print("\nFirst 10 faces:")
    for i in range(min(10, nfaces)):
        f.seek(84 + i * 50)
        face_data = f.read(50)
        if len(face_data) < 50:
            print(f"Face {i}: TRUNCATED")
            break
        normal = struct.unpack('<fff', face_data[0:12])
        v1 = struct.unpack('<fff', face_data[12:24])
        v2 = struct.unpack('<fff', face_data[24:36])
        v3 = struct.unpack('<fff', face_data[36:48])
        print(f"Face {i}: normal={normal[0]:.3f},{normal[1]:.3f},{normal[2]:.3f}")
EOFPY
V4_PID=$!
echo "v4 检查进程 PID: $V4_PID"

echo ""
echo "=== Thread 3: A406 Fluent 环境搜索 ==="
sshpass -p '159357' ssh -o StrictHostKeyChecking=no -p 60683 root@yun.1617k.com '
  echo "搜索 Fluent/ANSYS 安装..."
  find /opt /usr/local /ansys* /ANSYS* -name "fluent" -type f 2>/dev/null | head -5 || echo "未找到 Fluent"
  echo ""
  echo "检查是否有 Docker Fluent 镜像..."
  docker images 2>/dev/null | grep -i fluent || echo "无 Docker Fluent"
' > /d/realFluidAllMach/a406_fluent_search.log 2>&1 &
A406_PID=$!
echo "A406 搜索进程 PID: $A406_PID"

echo ""
echo "=== 所有后台任务已启动 ==="
echo "1. RFAM 编译：tail -f /d/realFluidAllMach/build_wsl/make.log"
echo "2. v4 检查：cat /d/Vortex/cad/injector/v4_quick_check.log"
echo "3. A406 搜索：cat /d/realFluidAllMach/a406_fluent_search.log"
