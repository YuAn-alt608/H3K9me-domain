# -*- coding: utf-8 -*-
"""Figure 1, redrawn from the source tables."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from data import r, style, FAM, TRACKS, KEY, ORANGE, BLUE, GREY, DGREY

style(plt)

rho = r("d7_within_rho.tsv").pivot(index="family", columns="track", values="rho").loc[FAM]
rho = rho[[KEY[t] for t in TRACKS]]
rho.columns = TRACKS

dd = r("d10_continuous.tsv")
FOCAL = ["tnrc18", "NFYB", "POLR2A", "INPUT", "ETS1", "CBX1", "CBX5", "CBX3"]
null = dd[~dd.track.isin(FOCAL)].groupby("track").rho.mean()
foc = dd[dd.track.isin(FOCAL)].groupby("track").rho.mean()

orc = r("d7_matched_or.tsv").pivot(index="family", columns="track", values="OR_matched").loc[FAM]
orc = orc[[KEY[t] for t in TRACKS]]
orc.columns = TRACKS

pk = r("d7_peak_rho.tsv").pivot(index="family", columns="track", values="rho_peak").loc[FAM]

fig = plt.figure(figsize=(6.0, 5.12))
gs = fig.add_gridspec(2, 2, left=0.105, right=0.985, top=0.915, bottom=0.135,
                      hspace=0.78, wspace=0.62, height_ratios=[1, 0.78])

# ---------------- A heatmap ----------------
axA = fig.add_subplot(gs[0, 0])
M = rho.values
im = axA.imshow(M, cmap="Blues", vmin=0, vmax=0.73, aspect="auto")
axA.set_xticks(range(len(TRACKS)))
axA.set_xticklabels([t if t != "INPUT" else "ChIP-seq\ninput" for t in TRACKS],
                    rotation=90, ha="center", va="top", fontsize=6.2)
axA.set_yticks(range(len(FAM)))
axA.set_yticklabels(FAM, fontsize=6.1)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        axA.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=5.0,
                 color="white" if v > 0.42 else "#333333")
axA.tick_params(length=1.5, pad=1.5)
for s in axA.spines.values():
    s.set_visible(False)
axA.set_xticks(np.arange(-.5, len(TRACKS), 1), minor=True)
axA.set_yticks(np.arange(-.5, len(FAM), 1), minor=True)
axA.grid(which="minor", color="white", linewidth=0.9)
axA.tick_params(which="minor", length=0)
axA.text(-0.30, 1.135, "A", transform=axA.transAxes, fontsize=9, fontweight="bold", va="top")
axA.text(-0.21, 1.135, "Within-family correlation with H3K9me3",
         transform=axA.transAxes, fontsize=7.1, fontweight="bold", va="top")
cax = fig.add_axes([0.456, 0.660, 0.0075, 0.200])
cb = fig.colorbar(im, cax=cax)
cb.set_ticks([0, 0.2, 0.4, 0.6])
cb.ax.tick_params(labelsize=5.4, length=1.4, pad=1.0)
cb.outline.set_visible(False)
cax.set_title("ρ", fontsize=6.0, pad=2.5)

# ---------------- B null distribution ----------------
axB = fig.add_subplot(gs[0, 1])
rng = np.random.default_rng(7)
axB.scatter(0.30 + rng.uniform(-0.055, 0.055, len(null)), null.values, s=9,
            facecolor="#CFCFCF", edgecolor="#B0B0B0", linewidth=0.3, zorder=2)
axB.boxplot([null.values], positions=[0.30], widths=0.13, showfliers=False,
            patch_artist=True, medianprops=dict(color="#3A3A3A", lw=1.0),
            boxprops=dict(facecolor="none", edgecolor="#7A7A7A", lw=0.6),
            whiskerprops=dict(color="#7A7A7A", lw=0.6),
            capprops=dict(color="#7A7A7A", lw=0.6))
axB.axhline(null.median(), color="#9A9A9A", lw=0.6, ls=(0, (3, 2.5)), zorder=1)

order = [("CBX3", BLUE, "CBX3 (HP1γ)"), ("CBX5", BLUE, "CBX5 (HP1α)"),
         ("CBX1", BLUE, "CBX1 (HP1β)"), ("INPUT", "#4D4D4D", "ChIP-seq input"),
         ("ETS1", "#4D4D4D", "ETS1"), ("POLR2A", "#4D4D4D", "POLR2A"),
         ("NFYB", "#4D4D4D", "NFYB"), ("TNRC18", ORANGE, "TNRC18")]
lab_y = {"CBX3": 0.512, "CBX5": 0.455, "CBX1": 0.417, "ETS1": 0.330,
         "INPUT": 0.285, "POLR2A": 0.215, "NFYB": 0.168, "TNRC18": 0.060}
for t, col, lab in order:
    y = foc[KEY[t]]
    yl = lab_y[t]
    axB.plot([0.86, 1.02], [yl, yl], color=col, lw=1.2, solid_capstyle="butt")
    axB.text(1.05, yl, lab, fontsize=6.0, va="center", ha="left", color=col,
             fontweight="bold" if t == "TNRC18" else "normal")
    if abs(yl - y) > 0.004:
        axB.annotate("", xy=(1.015, y), xytext=(1.015, yl),
                     arrowprops=dict(arrowstyle="-", color=col, lw=0.45, alpha=0.75))
axB.set_xlim(0.10, 1.98)
axB.set_ylim(0.030, 0.575)
axB.set_xticks([0.30])
axB.set_xticklabels(["36 random\nTF tracks"], fontsize=6.0)
axB.set_ylabel("mean within-family Spearman\ncorrelation with H3K9me3", fontsize=6.2, labelpad=2)
axB.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5])
axB.tick_params(axis="y", labelsize=6.0)
axB.text(-0.185, 1.135, "B", transform=axB.transAxes, fontsize=9, fontweight="bold", va="top")
axB.text(-0.085, 1.135, "Placement against 36 random proteins",
         transform=axB.transAxes, fontsize=7.1, fontweight="bold", va="top")

# ---------------- C prevalence-matched odds ratios ----------------
axC = fig.add_subplot(gs[1, 0])
x = np.arange(len(FAM))
for t in ["CBX5", "CBX1", "ETS1", "POLR2A", "NFYB", "INPUT"]:
    axC.plot(x, orc[t].values, color=GREY, lw=0.8, marker="o", ms=1.7,
             mfc=GREY, mec=GREY, zorder=1)
axC.plot(x, orc["TNRC18"].values, color=ORANGE, lw=1.5, marker="o", ms=3.3,
         mfc=ORANGE, mec=ORANGE, zorder=4, label="TNRC18")
axC.plot(x, orc["CBX3"].values, color=BLUE, lw=1.5, marker="o", ms=3.3,
         mfc=BLUE, mec=BLUE, zorder=3, label="CBX3 (HP1γ)")
axC.axhline(1, color="#9A9A9A", lw=0.6, ls=(0, (3, 2.5)), zorder=0)
axC.set_yscale("log")
axC.set_ylim(0.55, 900)
axC.set_yticks([1, 10, 100])
axC.set_yticklabels(["$10^0$", "$10^1$", "$10^2$"], fontsize=6.0)
axC.set_xticks(x)
axC.set_xticklabels(FAM, rotation=52, ha="right", fontsize=5.9)
axC.set_ylabel("odds ratio at matched prevalence", fontsize=6.5, labelpad=2)
h, l = axC.get_legend_handles_labels()
h.append(Line2D([0], [0], color=GREY, lw=0.8))
l.append("six further control tracks")
axC.legend(h, l, loc="upper center", bbox_to_anchor=(0.60, 1.01), fontsize=5.8,
           handlelength=1.5, borderpad=0.1, labelspacing=0.3, ncol=1)
axC.text(-0.215, 1.165, "C", transform=axC.transAxes, fontsize=9, fontweight="bold", va="top")
axC.text(-0.135, 1.165, "Odds ratios at matched positive rate",
         transform=axC.transAxes, fontsize=7.1, fontweight="bold", va="top")

# ---------------- D 300 bp window maximum ----------------
axD = fig.add_subplot(gs[1, 1])
w = 0.38
axD.bar(x - w / 2, pk["tnrc18"].values, width=w, color=ORANGE, label="TNRC18", lw=0)
axD.bar(x + w / 2, pk["CBX3"].values, width=w, color=BLUE, label="CBX3 (HP1γ)", lw=0)
axD.axhline(0, color="#8C8C8C", lw=0.6)
axD.set_xticks(x)
axD.set_xticklabels(FAM, rotation=52, ha="right", fontsize=5.9)
axD.set_ylabel("within-family Spearman\ncorrelation, 300 bp maximum", fontsize=6.2, labelpad=2)
axD.set_ylim(-0.11, 0.86)
axD.set_yticks([0.0, 0.2, 0.4, 0.6])
axD.legend(loc="upper left", bbox_to_anchor=(0.02, 1.00), ncol=1, fontsize=5.8,
           handlelength=1.1, borderpad=0.1, labelspacing=0.3)
axD.text(-0.215, 1.165, "D", transform=axD.transAxes, fontsize=9, fontweight="bold", va="top")
axD.text(-0.135, 1.165, "Window-maximum readout",
         transform=axD.transAxes, fontsize=7.1, fontweight="bold", va="top")

fig.savefig("fig1.png", dpi=450, facecolor="white")
print("fig1 written")
