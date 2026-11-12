# -*- coding: utf-8 -*-
"""
SSH 模型 d 矢量在 Bloch 球上的投影轨迹。

拓扑相: t1 = 0.5, t2 = 1.0
平庸相: t1 = 1.0, t2 = 0.5
运行: python SSH_bloch.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from Bloch_sphere import Arrow3D, draw_bloch_sphere


def ssh_d(k, t1, t2):
    k = np.asarray(k, dtype=float)
    return np.column_stack([t1 + t2 * np.cos(k),
                            t2 * np.sin(k),
                            np.zeros_like(k)])


def plot_ssh_bloch(t1, t2, outfile):
    k = np.linspace(-np.pi, np.pi, 801)
    d = ssh_d(k, t1, t2)
    d_hat = d / np.linalg.norm(d, axis=1, keepdims=True)

    ax, fig = draw_bloch_sphere()

    pts = d_hat[:, None, :]
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    seg_len = np.linalg.norm(np.diff(d_hat, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg_len)])
    s = 0.5 * (s[:-1] + s[1:]) / s[-1]
    fade = np.clip((s - 0.75) / 0.25, 0.0, 1.0)
    alpha = 1.0 - fade
    colors = np.column_stack([np.ones_like(s), np.zeros_like(s),
                              np.zeros_like(s), alpha])
    traj = Line3DCollection(segs, colors=colors, linewidths=2.2, zorder=10)
    ax.add_collection3d(traj)

    start = d_hat[0]
    ax.add_artist(Arrow3D([0.0, start[0]], [0.0, start[1]], [0.0, start[2]],
                          mutation_scale=18, lw=2.2, arrowstyle="-|>",
                          color="red", zorder=12))
    ax.scatter([start[0]], [start[1]], [start[2]], color="red",
               s=18, zorder=13)

    handles = [
        Line2D([0], [0], color="red", lw=2.2,
               label=r"$\hat{\mathbf{d}}(k)$"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=11)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    plot_ssh_bloch(0.5, 1.0, "ssh_bloch_topo.png")
    plot_ssh_bloch(1.0, 0.5, "ssh_bloch_trivial.png")
    plt.show()
