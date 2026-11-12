"""Plot I = 1 + cos(theta) and I = 1 - cos(theta) on the same figure.

theta ranges over [-2*pi, 2*pi]; the figure is saved to ../fig/.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

theta = np.linspace(-2 * np.pi, 2 * np.pi, 2000)
I_plus = 1 + np.cos(theta)
I_minus = 1 - np.cos(theta)

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(theta, I_plus, color="tab:blue", lw=2, label=r"$I_+ = 1+\cos\theta$")
ax.plot(theta, I_minus, color="tab:red", lw=2, ls="--", label=r"$I_- = 1-\cos\theta$")

ax.set_xlim(-2 * np.pi, 2 * np.pi)
ax.set_ylim(0, 2)
ax.set_xlabel(r"$\theta$", fontsize=14)
ax.set_ylabel(r"$I$", fontsize=14)
ax.set_xticks(
    [-2 * np.pi, -3 * np.pi / 2, -np.pi, -np.pi / 2, 0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]
)
ax.set_xticklabels(
    [r"$-2\pi$", r"$-3\pi/2$", r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"]
)
ax.set_yticks([0, 0.5, 1, 1.5, 2])
ax.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.7)
ax.grid(alpha=0.3)
ax.legend(fontsize=12, loc="upper right", framealpha=0.9)

fig.tight_layout()

fig_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fig")
os.makedirs(fig_dir, exist_ok=True)
out_png = os.path.join(fig_dir, "I_vs_theta.png")
out_pdf = os.path.join(fig_dir, "I_vs_theta.pdf")
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print(f"saved: {os.path.normpath(out_png)}")
print(f"saved: {os.path.normpath(out_pdf)}")
