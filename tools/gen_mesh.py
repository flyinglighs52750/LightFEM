#!/usr/bin/env python3
"""生成 N×N 结构化 Q4 网格的 LightFEM 输入文件。

用于复现论文第 4.4 节的规模扫描：N=200 时得到 40 000 个单元、
80 802 个自由度，与论文一致。

用法：
    python tools/gen_mesh.py --n 200 -o mesh_200x200.fem
    python tools/gen_mesh.py --n 200 --assembly coloring --threads 8 --reordering rcm -o a.fem
"""

from __future__ import annotations

import argparse
import sys


def build(nx: int, ny: int, length: float, height: float, thickness: float,
          young: float, poisson: float, traction: float, assembly: str,
          threads: int, reordering: str, output: str) -> str:
    """按行优先（先 x 后 y）编号生成 .fem 文本。"""
    dx = length / nx
    dy = height / ny
    nx1 = nx + 1

    def nid(i: int, j: int) -> int:
        # i: 列号 0..nx, j: 行号 0..ny
        return j * nx1 + i + 1

    lines: list[str] = []
    lines.append("# 由 tools/gen_mesh.py 生成的结构化 Q4 网格算例")
    lines.append(f"# 网格 {nx}x{ny}，单元 {nx * ny}，节点 {nx1 * (ny + 1)}")
    lines.append("")
    lines.append("[settings]")
    lines.append("plane_mode plane_stress")
    lines.append("solver pcg")
    lines.append("tolerance 1.0e-10")
    lines.append("max_iterations 5000")
    lines.append(f"reordering {reordering}")
    lines.append(f"assembly {assembly}")
    lines.append(f"threads {threads}")
    lines.append(f"output {output}")
    lines.append("")
    lines.append("[materials]")
    lines.append("# id  E(Pa)  nu  thickness(m)")
    lines.append(f"1 {young:.6e} {poisson:.2f} {thickness:.4f}")
    lines.append("")
    lines.append("[nodes]")
    lines.append("# id  x(m)  y(m)")
    for j in range(ny + 1):
        y = j * dy
        for i in range(nx1):
            lines.append(f"{nid(i, j)} {i * dx:.6f} {y:.6f}")
    lines.append("")
    lines.append("[elements]")
    lines.append("# id  n1 n2 n3 n4  material_id（逆时针）")
    eid = 0
    for j in range(ny):
        for i in range(nx):
            eid += 1
            lines.append(
                f"{eid} {nid(i, j)} {nid(i + 1, j)} {nid(i + 1, j + 1)} {nid(i, j + 1)} 1"
            )
    lines.append("")
    lines.append("[boundary]")
    lines.append("# node_id ux uy：左端整边固支")
    for j in range(ny + 1):
        lines.append(f"{nid(0, j)} 0.0 0.0")
    lines.append("")
    lines.append("[loads]")
    lines.append("# node_id fx(N) fy(N)：右端均布剪力，按节点从属长度分配")
    # 右端每个节点承担 dy/2 的从属长度（角点再按一半计），总剪力 = traction*height*thickness
    node_len = [dy] * (ny + 1)
    node_len[0] = dy / 2.0
    node_len[-1] = dy / 2.0
    for j in range(ny + 1):
        fy = traction * thickness * node_len[j]
        lines.append(f"{nid(nx, j)} 0.0 {fy:.6f}")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description="生成 LightFEM 结构化 Q4 网格输入文件")
    ap.add_argument("--n", type=int, default=200, help="x、y 方向单元数（正方形网格），默认 200")
    ap.add_argument("--nx", type=int, default=None, help="x 方向单元数，覆盖 --n")
    ap.add_argument("--ny", type=int, default=None, help="y 方向单元数，覆盖 --n")
    ap.add_argument("--length", type=float, default=3.0, help="域长，默认 3.0 m")
    ap.add_argument("--height", type=float, default=1.0, help="域高，默认 1.0 m")
    ap.add_argument("--thickness", type=float, default=0.01, help="厚度，默认 0.01 m")
    ap.add_argument("--young", type=float, default=2.10e11, help="弹性模量，默认 210 GPa")
    ap.add_argument("--poisson", type=float, default=0.30, help="泊松比，默认 0.30")
    ap.add_argument("--traction", type=float, default=-1.0e5, help="右端均布剪力，默认 -1e5 Pa")
    ap.add_argument("--assembly", default="serial", choices=["serial", "atomic", "coloring", "thread_local"])
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--reordering", default="rcm", choices=["none", "rcm"])
    ap.add_argument("--output", default="result.vtk")
    ap.add_argument("-o", "--out", default="mesh.fem", help="输出的 .fem 路径")
    args = ap.parse_args()

    nx = args.nx if args.nx is not None else args.n
    ny = args.ny if args.ny is not None else args.n
    if nx < 1 or ny < 1:
        print("nx/ny 必须为正整数", file=sys.stderr)
        return 2

    text = build(nx, ny, args.length, args.height, args.thickness, args.young,
                 args.poisson, args.traction, args.assembly, args.threads,
                 args.reordering, args.output)
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print(f"已生成 {args.out}：{nx}x{ny} 网格，{nx * ny} 单元，"
          f"{(nx + 1) * (ny + 1)} 节点，{2 * (nx + 1) * (ny + 1)} 自由度")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
