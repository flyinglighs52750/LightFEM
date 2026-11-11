# LightFEM

LightFEM 是面向二维小变形线弹性问题的轻量级有限元求解器。当前版本支持平面应力/平面应变 Q4 静力分析，从文本输入、模型检查、CSR 装配和边界消元到 Jacobi 预条件 CG 求解与 VTK 结果写出均可直接运行。

## 已实现

- CMake 工程、核心静态库、命令行程序和 CTest 测试；
- 节点、材料、Q4 单元、位移边界、节点载荷、体力和边界面载荷；
- 按节点 ID 排序的确定性自由度编号和 ID 查询索引；
- 分节式 `.fem` 文本输入及带行号的语法错误；
- ID、引用关系、材料参数、Q4 Jacobian、网格连通性、刚体模态约束、载荷和输出设置检查；
- 四节点双线性等参单元刚度矩阵（2×2 Gauss 积分），支持平面应力和应变；
- 基于网格拓扑预分配的 CSR 全局刚度矩阵与节点载荷累加；
- RCM 节点/自由度重排序，自动报告重排前后 CSR `nnz`、带宽和 profile；
- 预计算“单元局部矩阵项 → CSR values 下标”映射，数值装配不再搜索 CSR 行；
- 对称 Dirichlet 消元、Jacobi 预条件共轭梯度（PCG）与约束反力恢复；
- 恒定体力、恒定边界面力和压力的等效节点载荷；`fixed/pinned`、滚动及对称边界简写；
- OpenMP 单元刚度计算和 CSR SpMV；支持原子、单元着色、线程局部缓冲三种并行装配策略；
- Legacy ASCII VTK 输出：节点位移、反力、位移模量、单元中心应变/应力及面内 von Mises 应力；
- 节点位移与高斯点应变/应力 CSV 导出；VTK 变形放大系数 `deformation_scale`；
- 合法悬臂模型和顺时针错误单元样例。

## 构建

在 Visual Studio 2022 Developer PowerShell 或具有 CMake 的终端中执行：

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## 检查或求解模型

```powershell
.\build\Release\lightfem.exe check .\examples\cantilever.fem
```

合法模型返回退出码 `0`。语法错误、模型错误和程序错误分别返回退出码 `2`、`3` 和 `4`。

检查错误样例：

```powershell
.\build\Release\lightfem.exe check .\examples\invalid_clockwise.fem
```

该样例应返回 `E_NONPOSITIVE_JACOBIAN`。详细输入规范见 [docs/input_format.md](docs/input_format.md)。

求解合法算例，并按 `[settings]` 中的 `output` 写出 VTK：

```powershell
.\build\Release\lightfem.exe solve .\examples\cantilever.fem
```

OpenMP 着色装配示例：

```powershell
.\build\Release\lightfem.exe solve .\examples\cantilever_parallel.fem
```

体力、边界面力及压力示例：

```powershell
.\build\Release\lightfem.exe solve .\examples\cantilever_extended_loads.fem
```

也可以在命令行覆盖输出路径：

```powershell
.\build\Release\lightfem.exe solve .\examples\cantilever.fem .\cantilever.vtk
```

生成的 `.vtk` 可用 ParaView 打开。命令行会报告自由度数、CSR 非零项数、PCG 迭代次数及初末残差。求解同时会写出 `<stem>_nodes.csv`（节点位移、反力、单元中心应力）与 `<stem>_gauss_points.csv`（各高斯点应变/应力）；在 `[settings]` 中设置 `deformation_scale`（默认 `1.0`，`0` 表示不变形）可控制 VTK 中的位移放大比例，`csv_prefix` 可覆盖 CSV 文件名前缀。

在 `[settings]` 中设置 `reordering rcm` 可启用 Reverse Cuthill-McKee 重排序；`reordering none` 为默认值。`PROFILE` 行会同时给出重排前后矩阵的 `nnz`、带宽和 profile，理论上 RCM 不改变 `nnz`。

若构建时检测到 OpenMP，可用 `assembly atomic|coloring|thread_local` 选择并行装配，并以 `threads 4` 指定线程数（`0` 使用 OpenMP 运行时默认值）。单元刚度计算和 CSR SpMV 会并行执行；`PARALLEL` 行会输出实际线程容量、颜色数及刚度/装配耗时。`atomic` 内存占用低但可能有写竞争，`coloring` 对同色单元无写冲突，`thread_local` 无竞争但额外占用约 `线程数 × nnz × 8` 字节。未检测到 OpenMP 时仍可使用 `assembly serial`。

## 当前分析范围

仅支持二维、线弹性、小变形、单材料或多材料 Q4 网格的线性静力问题。单元应力和应变为单元中心值；当前载荷为恒定节点力、体力、边界面力或压力，不包含随时间/空间变化的载荷、热载荷、非线性和动力学。

## 目录结构

```text
LightFEM/
├─ CMakeLists.txt
├─ include/lightfem/    公共接口
├─ src/                 实现与命令行入口
├─ tests/               无第三方依赖的自动测试
├─ examples/            合法和错误输入样例
├─ tools/               算例生成脚本
└─ docs/                输入格式说明
```

## 复现论文实验

学位论文第 4 章的性能实验全部由命令行版本测得，不经过图形界面。规模扫描用的结构化网格由 `tools/gen_mesh.py` 生成，无需手工建模：

```powershell
# 200×200 结构化网格：与论文 4.4 节同规模
python tools/gen_mesh.py --n 200 -o mesh_200x200.fem
.\build\Release\lightfem.exe solve mesh_200x200.fem mesh_200x200.vtk
```

该算例的输出与论文 4.4 节一致：`nodes=40401`、`elements=40000`、`dofs=80802`、`nnz=1444804`，`PROFILE` 行给出带宽 405→803、profile 32 522 001→42 949 601（表 4-3）。

装配策略对比（表 4-4）通过命令行参数切换：

```powershell
python tools/gen_mesh.py --n 200 --assembly serial      --threads 8 -o a_serial.fem
python tools/gen_mesh.py --n 200 --assembly atomic      --threads 8 -o a_atomic.fem
python tools/gen_mesh.py --n 200 --assembly coloring    --threads 8 -o a_coloring.fem
python tools/gen_mesh.py --n 200 --assembly thread_local --threads 8 -o a_thread_local.fem
```

逐个 `solve` 后读取 `TIMING` 行中的 `element_stiffness_ms`、`assembly_ms` 与 `pcg_ms`。论文中的数字是同批次 7 轮取中位数，绝对值随机器状态波动，只有同批次内可比。

**注意**：不同批次之间不要比较绝对耗时；单元刚度计算、装配、PCG 三个阶段需要在同一进程批次内采集。

## 图形界面

`LightFEM Studio` 是独立的前后处理界面，负责几何草绘、材料/网格/边界/载荷设置与结果云图显示，本身不参与矩阵存储、装配与求解。界面只生成 `.fem` 文件并调用 `lightfem.exe`。Windows 免安装包见 [dist/LightFEM-Studio-0.8.0-win64.zip](dist/LightFEM-Studio-0.8.0-win64.zip)，解压后双击 `LightFEM Studio.exe` 即可运行，无需安装。

## 许可证

MIT License，见 [LICENSE](LICENSE)。

## 引用

若在论文或报告中引用本工作，请标注：

> 钟诺熙. 面向二维线弹性问题的轻量级有限元求解器设计——稀疏矩阵存储与装配优化[D]. 广州: 华南理工大学, 2026.
