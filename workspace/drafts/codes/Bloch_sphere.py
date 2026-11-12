# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 13:16:26 2026

@author: kunpu
"""

"""
绘制 Bloch 球 (Bloch Sphere)
依赖: numpy, matplotlib
运行: python bloch_sphere.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


# ---------------------------------------------------------------- #
# 1. 一个能在 3D 坐标轴里画“真箭头”的小工具类
# ---------------------------------------------------------------- #
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return float(np.min(zs))


# ---------------------------------------------------------------- #
# 2. 量子态 -> Bloch 矢量坐标
# ---------------------------------------------------------------- #
def bloch_coords(alpha, beta):
    """
    由态 |psi> = alpha|0> + beta|1> 计算 Bloch 矢量 (x, y, z)。
    会自动归一化。
    """
    amp = np.array([alpha, beta], dtype=complex)
    amp = amp / np.linalg.norm(amp)
    x = 2 * np.real(np.conj(amp[0]) * amp[1])
    y = 2 * np.imag(np.conj(amp[0]) * amp[1])
    z = np.abs(amp[0]) ** 2 - np.abs(amp[1]) ** 2
    return np.array([x, y, z])


def bloch_coords_angles(theta, phi):
    """由极角 theta、方位角 phi 得到单位球面上的 Bloch 矢量。"""
    return np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(theta)])


# ---------------------------------------------------------------- #
# 3. 画球体本身
# ---------------------------------------------------------------- #
def draw_bloch_sphere(ax=None, n=60, sphere_color="#4C9BE8", alpha=0.10,
                      fill=True):
    """画一个半透明的 Bloch 球（含经纬网、坐标轴、基矢标注）。"""
    if ax is None:
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = ax.figure

    # --- 球面（经纬网格）---
    u = np.linspace(0, 2 * np.pi, n)
    v = np.linspace(0, np.pi, n)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    if fill:
        ax.plot_surface(xs, ys, zs, rstride=2, cstride=2,
                        color=sphere_color, alpha=alpha,
                        linewidth=0.35, edgecolor="#9ec8f0", shade=False)
    else:
        grid_color = "#9ec8f0"
        for i in range(0, n, 2):
            ax.plot(np.cos(u[i]) * np.sin(v), np.sin(u[i]) * np.sin(v),
                    np.cos(v), color=grid_color, lw=0.35)
        for j in range(2, n, 2):
            ax.plot(np.cos(u) * np.sin(v[j]), np.sin(u) * np.sin(v[j]),
                    np.full_like(u, np.cos(v[j])), color=grid_color, lw=0.35)

    # --- 三条坐标轴（带箭头）---
    axis_lens = (1.8, 1.35, 1.35)
    for i, axis_len in enumerate(axis_lens):
        seg = np.zeros((2, 3))
        seg[:, i] = [-axis_len, axis_len]
        arrow = Arrow3D(seg[:, 0], seg[:, 1], seg[:, 2],
                        mutation_scale=14, lw=0.8, ls="--",
                        arrowstyle="-|>", color="0.45")
        ax.add_artist(arrow)

    ax.text(axis_lens[0] + 0.15, 0, 0, r"$d_x$", fontsize=14,
            ha="center", va="center")
    ax.text(0, axis_lens[1] + 0.15, 0, r"$d_y$", fontsize=14,
            ha="center", va="center")
    ax.text(0, 0, axis_lens[2] + 0.1, r"$d_z$", fontsize=14, ha="center")

    # --- 视觉设置 ---
    lim = max(axis_lens) + 0.25
    ax.set_box_aspect((1, 1, 1))          # 保证球不变形
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_axis_off()
    ax.view_init(elev=20, azim=35)
    return ax, fig


# ---------------------------------------------------------------- #
# 4. 添加一个 Bloch 矢量（箭头 + 投影虚线）
# ---------------------------------------------------------------- #
def add_bloch_vector(ax, vec, color="crimson", label=None, lw=2.2):
    x, y, z = vec
    arrow = Arrow3D([0, x], [0, y], [0, z], mutation_scale=18,
                    lw=lw, arrowstyle="-|>", color=color)
    ax.add_artist(arrow)
    ax.scatter([x], [y], [z], color=color, s=25)

    # 到三个坐标平面的投影虚线
    ax.plot([x, x], [y, y], [0, z], color=color, lw=1, ls=":")
    ax.plot([x, x], [0, y], [0, 0], color=color, lw=1, ls=":")
    ax.plot([0, x], [y, y], [0, 0], color=color, lw=1, ls=":")

    if label:
        ax.text(1.1 * x, 1.1 * y, 1.1 * z, label,
                color=color, fontsize=13, ha="center")


# ---------------------------------------------------------------- #
# 5. 示例
# ---------------------------------------------------------------- #
if __name__ == "__main__":
    ax, fig = draw_bloch_sphere()

    show_vectors = False
    samples = [
        (r"$|0\rangle$",  bloch_coords(1, 0),      "tab:red"),
        (r"$|+\rangle$",  bloch_coords(1, 1),      "tab:blue"),
        (r"$|+i\rangle$", bloch_coords(1, 1j),     "tab:green"),
        (r"$|-\rangle$",  bloch_coords(1, -1),     "tab:purple"),
    ]

    if show_vectors:
        for label, vec, color in samples:
            add_bloch_vector(ax, vec, color=color, label=label)

    plt.tight_layout()
    plt.savefig("bloch_sphere.png", dpi=200, bbox_inches="tight")
    plt.show()