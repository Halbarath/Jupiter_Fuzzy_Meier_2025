#!/usr/bin/env python3

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 13})
plt.rcParams.update({"lines.linewidth": 1.75})

fs_label = 16
fs_legend = 13
alpha = 0.9

cmap = matplotlib.cm.get_cmap("tab10")
colors = [cmap(i) for i in range(3)]
labels = ["Post evolution (hot)", "Post evolution (warm)", "Post evolution (cool)"]
titles = {
    "headon_intermediate": "Head-on impact (intermediate state)",
    "headon_final": "Head-on impact (final state)",
    "oblique_intermediate": "Oblique impact (intermediate state)",
    "oblique_final": "Oblique impact (final state)",
}

def plot_profiles(ax: plt.Axes, impact: str, x: str = "radius") -> None:
    """Plots the profiles of the heavy-element abundance for different impacts.

    Args:
        ax (plt.Axes): Instance of the matplotlib Axes class to plot on.
        impact (str): The impact scenario to plot, e.g., "headon_final".
        x (str, optional): Variable to use for the x-axis. Options are "radius"
            or "mass". Note that both values are normalised. Defaults to "radius".
    """
    assert x in ["radius", "mass"], "x must be either 'radius' or 'mass'"

    file = f"pre_evolution_{impact}.data"
    src = os.path.join("data", impact, file)
    data = np.loadtxt(src, skiprows=1)
    x_vals = data[:, 0] if x == "radius" else data[:, 1]
    ax.plot(
        x_vals,
        data[:, 2],
        color="black",
        ls="--",
        alpha=alpha,
        label="Post impact (this work)",
    )

    file = "liu_impact_headon.data" if "headon" in impact else "liu_impact_oblique.data"
    src = os.path.join("data", file)
    data = np.loadtxt(src, skiprows=1)
    x_vals = data[:, 0] if x == "radius" else data[:, 1]
    ax.plot(
        x_vals,
        data[:, 2],
        color="tab:grey",
        ls=":",
        alpha=alpha,
        label="Post impact (Liu et al. 2019)",
    )

    x_left = 0.45 if x == "radius" else 0.2
    x_right = 0.7 if x == "radius" else 0.6
    ax.axvspan(
        x_left,
        x_right,
        color="tab:green",
        alpha=0.1,
        hatch="//",
        label="Upper dilute-core boundary \nfrom interior models",
    )

    folder = os.path.join("data", impact)
    files = os.listdir(folder)
    for file in files:
        if "post_evolution" not in file:
            continue
        src = os.path.join(folder, file)
        i_col = int(file[-6])
        data = np.loadtxt(src, skiprows=1)
        x_vals = data[:, 0] if x == "radius" else data[:, 1]
        ax.plot(
            x_vals,
            data[:, 2],
            color=colors[i_col],
            label=labels[i_col],
            alpha=alpha,
        )
    ax.set_title(titles[impact], fontsize=fs_label)


for prefix in ["headon", "oblique"]:
    fig, axes = plt.subplots(2, 1, figsize=(5, 8), sharex=True, sharey=True)
    plot_profiles(ax=axes[0], impact=f"{prefix}_intermediate")
    plot_profiles(ax=axes[1], impact=f"{prefix}_final")
    axes[0].set_xlim(left=0, right=1)
    axes[1].set_xlabel("Normalised radius", fontsize=fs_label)
    axes[1].set_ylabel("Heavy-element fraction", fontsize=fs_label)
    axes[0].set_ylabel("Heavy-element fraction", fontsize=fs_label)
    axes[1].legend(
        ncol=2,
        fontsize=9,
        frameon=False,
        bbox_to_anchor=(1, -0.2),
    )
    fig.tight_layout()
    plt.savefig("composition_profiles_{}_radius.pdf".format(prefix),bbox_inches='tight')
