"""
Overlay plot: bO3 abundance (% of initial O2) vs fluence for
 - ion-on   (Model 1, Table 5 deltas, all 26 photoprocess rows active)
 - ion-off  (delta=0 on the 12 charged-product rows; neutral rows unchanged)
 - experimental helper points (Gerakines data; carried inside singleRMSD via
   exportable_custom_functions.setup_experimental_data)

This script reads the per-branch bO3.csv outputs snapshotted into this folder
and regenerates the experimental points so the figure is self-contained.
"""
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, "/Users/cnshingledecker/Research/Dday-newnetwork-porting")
from exportable_custom_functions import setup_experimental_data

INITIAL_O2 = 5.7e22

exp = setup_experimental_data()

def load_run(path):
    df = pd.read_csv(path, header=1, names=["Fluence", "Abundance"])
    df["Pct"] = df["Abundance"] / INITIAL_O2 * 100.0
    return df

on = load_run("bO3_ion_on.csv")
off = load_run("bO3_ion_off.csv")

fig, ax = plt.subplots(figsize=(7.5, 5.5))
ax.plot(on["Fluence"], on["Pct"], "-", color="C0", lw=2,
        label=r"ion-on (Model 1, Table 5 $\delta$)")
ax.plot(off["Fluence"], off["Pct"], "--", color="C3", lw=2,
        label=r"ion-off ($\delta=0$ on charged-product rows)")
ax.plot(exp["expX"], exp["expY"], "ko", ms=7, mfc="white",
        mew=1.5, label="experimental (Gerakines helper pts)")

ax.set_xscale("log")
ax.set_xlabel("Fluence (photons cm$^{-2}$)")
ax.set_ylabel(r"bO$_3$ / initial O$_2$ (%)")
ax.set_title("Reviewer Comment 1.9 paired comparison")
ax.grid(True, which="both", ls=":", alpha=0.5)
ax.legend(loc="upper left", frameon=False)

# Annotate RMSD summary
on_sum = pd.read_csv("summary_ion_on.csv").set_index("Metric")["Value"]
off_sum = pd.read_csv("summary_ion_off.csv").set_index("Metric")["Value"]
txt = (f"Unweighted RMSD  (ion-on): {on_sum['Unweighted_RMSD']:.3f}\n"
       f"Unweighted RMSD (ion-off): {off_sum['Unweighted_RMSD']:.3f}\n"
       f"Weighted RMSD   (ion-on): {on_sum['Weighted_RMSD']:.3f}\n"
       f"Weighted RMSD  (ion-off): {off_sum['Weighted_RMSD']:.3f}")
ax.text(0.98, 0.05, txt, transform=ax.transAxes, ha="right", va="bottom",
        family="monospace", fontsize=9,
        bbox=dict(boxstyle="round", fc="white", alpha=0.85, ec="0.7"))

plt.tight_layout()
plt.savefig("fig_paired_overlay.png", dpi=150)
plt.savefig("fig_paired_overlay.pdf")
print("Wrote fig_paired_overlay.png, fig_paired_overlay.pdf")
