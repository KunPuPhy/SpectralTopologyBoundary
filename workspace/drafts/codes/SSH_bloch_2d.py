# -*- coding: utf-8 -*-
"""
SSH 模型 d 矢量轨迹在二维 d_x-d_y 平面上的投影。

Bloch 球在俯视图中退化为灰色虚线的单位圆。
拓扑相: t1 = 0.5, t2 = 1.0
平庸相: t1 = 1.0, t2 = 0.5
运行: python SSH_bloch_2d.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch


def ssh_d(k, t1, t2):
    k = np.asarray(k, dtype=float)
    return np.column_stack([t1 + t2 * np.cos(k),
                            t2 * np.sin(k),
                            np.zeros_like(k)])


def draw_reference(ax, lim=3.1, axis_lens=(1.8, 1.35)):
    th = np.linspace(0.0, 2.0 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color="0.65", ls="--", lw=1.8, zorder=1)

    ax.add_patch(FancyArrowPatch(
        (-axis_lens[0], 0.0), (axis_lens[0], 0.0),
        arrowstyle="-|>", mutation_scale=18, lw=1.4, ls="-",
        color="black", zorder=2))
    ax.add_patch(FancyArrowPatch(
        (0.0, -axis_lens[1]), (0.0, axis_lens[1]),
        arrowstyle="-|>", mutation_scale=18, lw=1.4, ls="-",
        color="black", zorder=2))

    ax.text(axis_lens[0] + 0.12, 0.0, r"$d_x$", fontsize=14,
            ha="center", va="center")
    ax.text(0.0, axis_lens[1] + 0.12, r"$d_y$", fontsize=14,
            ha="center", va="center")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")


def plot_ssh_bloch_2d(t1, t2, outfile):
    k = np.linspace(-np.pi, np.pi, 801)
    d = ssh_d(k, t1, t2)
    d_hat = d[:, :2] / np.linalg.norm(d, axis=1, keepdims=True)

    fig, ax = plt.subplots(figsize=(7, 7))
    draw_reference(ax)

    pts = d_hat[:, None, :]
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    seg_len = np.linalg.norm(np.diff(d_hat, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg_len)])
    s = 0.5 * (s[:-1] + s[1:]) / s[-1]
    fade = np.clip((s - 0.75) / 0.25, 0.0, 1.0)
    alpha = 1.0 - fade
    colors = np.column_stack([np.ones_like(s), np.zeros_like(s),
                              np.zeros_like(s), alpha])
    traj = LineCollection(segs, colors=colors, linewidths=3.6, zorder=10)
    ax.add_collection(traj)

    start = d_hat[0]
    ax.add_patch(FancyArrowPatch(
        (0.0, 0.0), tuple(start), arrowstyle="-|>", mutation_scale=24,
        lw=3.6, color="red", zorder=12))
    ax.scatter([start[0]], [start[1]], color="red", s=36, zorder=13)

    cosang = d_hat @ start
    ang_max = np.arccos(np.clip(cosang.min(), -1.0, 1.0))
    candidates = np.where(cosang <= cosang.min() + 1e-6)[0]
    max_vecs = []
    for i in candidates:
        v = d_hat[i]
        if all(np.linalg.norm(v - u) > 1e-3 for u in max_vecs):
            max_vecs.append(v)

    for v in max_vecs:
        ax.add_patch(FancyArrowPatch(
            (0.0, 0.0), tuple(v), arrowstyle="-|>", mutation_scale=24,
            lw=3.6, color="tab:blue", zorder=14))

    print(f"{outfile}: max angle = {np.degrees(ang_max):.2f} deg")

    handles = [
        Line2D([0], [0], color="red", lw=3.6,
               label=r"$\hat{\mathbf{d}}(k)$"),
        Line2D([0], [0], color="tab:blue", lw=3.6,
               label=r"$\hat{\mathbf{d}}(k_{\max})$"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=11)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    plot_ssh_bloch_2d(0.5, 1.0, "ssh_bloch_2d_topo.png")
    plot_ssh_bloch_2d(1.0, 0.5, "ssh_bloch_2d_trivial.png")
    plt.show()
