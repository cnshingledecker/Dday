from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


INITIAL_O2 = 5.7e22
FIG1_WIDTH_IN = 278.571013 / 72.0
FIG1_HEIGHT_IN = 216.737467 / 72.0
MODEL1_BLUE = "#098ec3"
IONOFF_RED = "#cf4520"


def read_bo3_curve(path: Path) -> tuple[list[float], list[float]]:
    fluence: list[float] = []
    pct: list[float] = []
    with path.open() as handle:
        reader = csv.reader(handle)
        next(reader)  # comment header
        for row in reader:
            if not row:
                continue
            fl = float(row[0])
            abundance = float(row[1])
            fluence.append(fl)
            pct.append(abundance / INITIAL_O2 * 100.0)
    return fluence, pct


def read_eval_points(path: Path) -> tuple[list[float], list[float]]:
    fluence: list[float] = []
    exp_value: list[float] = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            fluence.append(float(row["Fluence"]))
            exp_value.append(float(row["Exp_Value"]))
    return fluence, exp_value


def main() -> None:
    figure_dir = Path(__file__).resolve().parent
    paper_dir = figure_dir.parent
    repo_dir = paper_dir.parent
    paired_dir = repo_dir / "reviewer_response" / "comment_1.9_paired_kristen"

    on_x, on_y = read_bo3_curve(paired_dir / "bO3_ion_on.csv")
    off_x, off_y = read_bo3_curve(paired_dir / "bO3_ion_off.csv")
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "axes.linewidth": 1.2,
            "axes.labelsize": 12.5,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 10.5,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 1.2,
            "ytick.major.width": 1.2,
            "xtick.minor.width": 0.9,
            "ytick.minor.width": 0.9,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.minor.size": 3.5,
            "ytick.minor.size": 3.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(figsize=(FIG1_WIDTH_IN, FIG1_HEIGHT_IN))
    ax.plot(on_x, on_y, color=MODEL1_BLUE, lw=2.2, label=r"Model 1")
    ax.plot(
        off_x,
        off_y,
        color=IONOFF_RED,
        lw=2.2,
        ls=(0, (5.5, 3.0)),
        label=r"ion-off control",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(8.0e14, 8.0e17)
    ax.set_ylim(4.0e-1, 3.0e1)
    ax.set_xlabel(r"Fluence $\left(\frac{\mathit{particles}}{\mathit{cm}^{2}}\right)$")
    ax.set_ylabel(r"$[\mathrm{O}_{3}]/[\mathrm{O}_{2\,\mathit{initial}}]\times 100\%$")

    legend = ax.legend(
        loc="upper left",
        frameon=True,
        facecolor="white",
        edgecolor="#bfbfbf",
        framealpha=0.88,
        fancybox=True,
        handlelength=2.0,
    )
    legend.get_frame().set_linewidth(1.0)

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.17, top=0.97)
    fig.savefig(paper_dir / "f7.png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(paper_dir / "f7.pdf", bbox_inches="tight", pad_inches=0.02)


if __name__ == "__main__":
    main()
