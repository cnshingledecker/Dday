from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


def add_species_cluster(ax, center_x: float, center_y: float, labels: list[tuple[str, str]]) -> None:
    offsets = [(-0.09, 0.08), (0.0, -0.02), (0.09, 0.08)]
    for (label, color), (dx, dy) in zip(labels, offsets):
        circle = Circle((center_x + dx, center_y + dy), 0.06, facecolor=color, edgecolor="black", lw=1.4)
        ax.add_patch(circle)
        ax.text(center_x + dx, center_y + dy, label, ha="center", va="center", fontsize=13, weight="bold")


def main() -> None:
    figure_dir = Path(__file__).resolve().parent
    paper_dir = figure_dir.parent

    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    bg = FancyBboxPatch(
        (0.02, 0.08),
        0.96,
        0.84,
        boxstyle="round,pad=0.02,rounding_size=0.04",
        facecolor="#f5f8fc",
        edgecolor="#9fb6cf",
        lw=1.8,
    )
    ax.add_patch(bg)

    left_box = FancyBboxPatch(
        (0.06, 0.22),
        0.22,
        0.56,
        boxstyle="round,pad=0.03,rounding_size=0.03",
        facecolor="#e8f0fb",
        edgecolor="#5b84c4",
        lw=1.6,
    )
    ax.add_patch(left_box)
    ax.text(0.17, 0.67, r"pure O$_2$ ice", ha="center", va="center", fontsize=18, weight="bold")
    ax.text(0.17, 0.46, r"O$_2$ + O$_2$ + O$_2$", ha="center", va="center", fontsize=16)
    ax.text(0.17, 0.30, r"10 K bulk mantle", ha="center", va="center", fontsize=12)

    ax.text(0.38, 0.69, r"$h\nu$", fontsize=22, weight="bold", color="#8c2d04")
    ax.text(0.38, 0.58, r"($< 11.3$ eV)", fontsize=11, color="#8c2d04")
    arrow1 = FancyArrowPatch((0.30, 0.50), (0.47, 0.50), arrowstyle="simple", mutation_scale=22, color="#8c2d04")
    ax.add_patch(arrow1)

    add_species_cluster(
        ax,
        0.58,
        0.50,
        [(r"O$_2^+$", "#ffd6a5"), (r"e$^-$", "#d9d9d9"), (r"O$^-$", "#cdeac0")],
    )
    ax.text(0.58, 0.27, "persistent ions and electrons", ha="center", va="center", fontsize=12)

    arrow2 = FancyArrowPatch((0.67, 0.50), (0.82, 0.50), arrowstyle="simple", mutation_scale=22, color="#2a6f97")
    ax.add_patch(arrow2)
    ax.text(0.745, 0.62, "ion-neutral / ion-ion", ha="center", va="center", fontsize=11, color="#2a6f97")
    ax.text(0.745, 0.38, "bulk chemistry", ha="center", va="center", fontsize=11, color="#2a6f97")

    right_box = FancyBboxPatch(
        (0.84, 0.22),
        0.10,
        0.56,
        boxstyle="round,pad=0.03,rounding_size=0.03",
        facecolor="#edf7ed",
        edgecolor="#4b8f5a",
        lw=1.6,
    )
    ax.add_patch(right_box)
    ax.text(0.89, 0.60, r"O$_3$", ha="center", va="center", fontsize=24, weight="bold", color="#216e39")
    ax.text(0.89, 0.38, "formation", ha="center", va="center", fontsize=12, color="#216e39")

    fig.tight_layout(pad=0)
    fig.savefig(paper_dir / "toc_ionice.png", dpi=200)
    fig.savefig(paper_dir / "toc_ionice.pdf")


if __name__ == "__main__":
    main()
