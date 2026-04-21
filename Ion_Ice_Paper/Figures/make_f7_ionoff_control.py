from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


INITIAL_O2 = 5.7e22


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
    exp_x, exp_y = read_eval_points(paired_dir / "eval_points.csv")

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 11,
            "legend.fontsize": 10,
        }
    )

    fig, ax = plt.subplots(figsize=(6.3, 4.5))
    ax.plot(on_x, on_y, color="#1f77b4", lw=2.4, label="ion-on (Model 1)")
    ax.plot(off_x, off_y, color="#c23b22", lw=2.4, ls="--", label="ion-off control")
    ax.plot(
        exp_x,
        exp_y,
        linestyle="none",
        marker="o",
        ms=6.5,
        mfc="white",
        mec="black",
        mew=1.2,
        label="experiment",
    )

    ax.set_xscale("log")
    ax.set_xlabel(r"Fluence (photons cm$^{-2}$)")
    ax.set_ylabel(r"O$_3$ / initial O$_2$ (%)")
    ax.grid(True, which="both", ls=":", lw=0.7, alpha=0.55)
    ax.legend(loc="upper left", frameon=False)

    inset = inset_axes(ax, width="42%", height="38%", loc="lower left", borderpad=1.5)
    inset.plot(on_x, on_y, color="#1f77b4", lw=1.5)
    inset.plot(off_x, off_y, color="#c23b22", lw=1.8, ls="--")
    inset.plot(exp_x, exp_y, linestyle="none", marker="o", ms=4.0, mfc="white", mec="black", mew=0.9)
    inset.set_xscale("log")
    inset.set_xlim(7.0e13, 1.2e18)
    inset.set_ylim(0.0, 4.4)
    inset.set_title("zoom: ion-off regime", fontsize=9)
    inset.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
    inset.tick_params(labelsize=8)

    fig.tight_layout()
    fig.savefig(paper_dir / "f7.png", dpi=300)
    fig.savefig(paper_dir / "f7.pdf")


if __name__ == "__main__":
    main()
