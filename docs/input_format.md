# LightFEM 输入文件格式

输入文件采用 UTF-8 纯文本，推荐扩展名为 `.fem`。空行会被忽略，`#` 后面的内容为注释。字段之间使用一个或多个空格分隔。

## settings

```text
[settings]
plane_mode plane_stress
solver pcg
tolerance 1.0e-10
max_iterations 5000
reordering rcm
assembly coloring
threads 4
output results.vtk
deformation_scale 100.0
csv_prefix results
```

- `plane_mode`：`plane_stress` 或 `plane_strain`；
- `solver`：当前只接受 `pcg`；
- `tolerance`：满足 `0 < tolerance < 1`；
- `max_iterations`：正整数；
- `reordering`：`none`（默认）或 `rcm`。后者按节点图执行 Reverse Cuthill-McKee 重排序；
- `assembly`：`serial`（默认）、`atomic`、`coloring` 或 `thread_local`；后三者要求构建时检测到 OpenMP；
- `threads`：非负整数；`0`（默认）使用 OpenMP 运行时默认线程数，`1` 可用于串行对照；
- `output`：结果文件名，建议使用 `.vtk`。
- `deformation_scale`：VTK 中位移放大比例，非负数，`0` 表示不变形，默认 `1.0`；
- `csv_prefix`：CSV 导出文件名前缀（输出 `<prefix>_nodes.csv` 与 `<prefix>_gauss_points.csv`），为空时由 VTK 输出文件主名推导。

所有设置都有默认值，但同一设置不允许重复出现。

## materials

```text
[materials]
# id E(Pa) nu thickness(m)
1 2.10e11 0.30 0.01
```

材料编号必须为正整数且不得重复；`E > 0`，`-1 < nu < 0.5`，厚度必须为正数。

## nodes

```text
[nodes]
# id x y
1 0.0 0.0
2 1.0 0.0
```

节点编号必须为正整数且不得重复，坐标必须为有限数值。当前范围不允许两个节点具有完全相同的坐标。

## elements

```text
[elements]
# id n1 n2 n3 n4 material_id
1 1 2 3 4 1
```

Q4 节点必须沿单元边界逆时针排列，四个节点不得重复。程序在单元中心和四个 2×2 高斯点检查 Jacobian，非正或接近零的 Jacobian 会被拒绝。

## boundary

```text
[boundary]
# node_id ux uy
1 0.0 0.0
4 0.0 free
```

`ux`、`uy` 可以是数值，也可以使用 `free`、`none` 或 `-` 表示不约束。程序检查重复约束，并根据二维刚体位移矩阵判断约束是否消除了两个平动和一个转动模态。

也可使用边界类型简写：

```text
1 fixed        # ux=0, uy=0；pinned 同义
4 roller_x     # ux=0；symmetry_x 同义
4 roller_y     # uy=0；symmetry_y 同义
```

## loads

```text
[loads]
# node_id fx fy
3 0.0 -500.0
```

同一节点允许出现多条载荷记录；求解器会累加这些载荷，检查程序会给出提示。载荷施加在已约束自由度上时也会产生警告。

## body_forces

```text
[body_forces]
# element_id bx(N/m^3) by(N/m^3)
1 0.0 -78500.0
```

每条记录为一个 Q4 单元内恒定体力密度。程序通过 2×2 Gauss 积分形成一致等效节点载荷，并自动乘以材料厚度。`body_loads` 可作为节名别名。

## edge_loads

```text
[edge_loads]
# element_id local_edge traction tx(Pa) ty(Pa)
2 2 traction 0.0 -2.0e6
# element_id local_edge pressure p(Pa); p>0 沿单元外法线反向（向内）
2 2 pressure 1.0e6
```

局部边号依次为 `1:n1-n2`、`2:n2-n3`、`3:n3-n4`、`4:n4-n1`。`traction` 是恒定边界面力（力/面积），按线段长度和材料厚度转换为等效节点力；`pressure` 是法向压力，正值向单元内部作用。网格须保持逆时针 Q4 连通性。`surface_loads`、`tractions` 和 `pressures` 可作为节名别名，但记录格式不变。

## 求解与 VTK 输出

```powershell
lightfem solve model.fem
lightfem solve model.fem result.vtk
```

`solve` 会先执行全部模型检查；检查通过后，采用 2×2 Gauss 积分形成 Q4 单元刚度，按网格拓扑预分配 CSR 全局矩阵，并以对称消元施加数值位移边界。线性方程组由 Jacobi 预条件共轭梯度法求解。第二种调用形式会覆盖 `settings` 中的 `output`。

CSR 模式创建后会预计算每个单元 8×8 局部刚度项在 `values` 数组中的直接下标；数值装配阶段不再对 CSR 行做查找。启用 `reordering rcm` 时，程序会输出重排前后的 `nnz`、带宽和 profile；重排序不改变物理节点 ID、连通性、载荷或数值解。

若 OpenMP 可用，所有策略均并行计算单元刚度，CSR 的矩阵向量乘法也按行并行。`atomic` 直接原子累加全局 `values`；`coloring` 将共享节点的单元分到不同颜色并逐色并行，消除同色写冲突；`thread_local` 为每个线程建立完整 CSR 数值缓冲，最终归并。后者额外内存约为 `threads × nnz × sizeof(double)`。求解输出中的 `PARALLEL` 行给出所选策略、线程数、颜色数和两个阶段耗时。

VTK 采用 ASCII Legacy Unstructured Grid 格式，包含：

- 节点数据：`displacement`、`reaction`、`displacement_magnitude`；
- 单元数据（单元中心）：`strain_xx`、`strain_yy`、`strain_xy`、`stress_xx`、`stress_yy`、`stress_xy`、`von_mises_in_plane`。

## 当前模型边界

当前阶段仅支持一个连通的二维 Q4 网格、恒定节点集中力/体力/边界面力/压力和位移边界。随空间或时间变化的分布载荷、梁壳单元、接触、非线性和动力分析均不在当前范围内。
