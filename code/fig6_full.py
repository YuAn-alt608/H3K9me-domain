"""Figure 6, rebuilt to match the v11 text.

A: paired 5' vs 3' medians for the mark and for DNase, with the raw share and
   the composition-controlled share.
B: observed/expected as the composition control is refined (120,000-window pool).
C: the paired 5'-3' difference on two assemblies, six values.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = pd.read_csv(os.path.join(ROOT, "data", "tables", "l1pa5_pool_expansion.tsv"), sep="\t")
T = T[T.cell == "H1 DNase (ENCODE)"]

INK, GREY, BLUE, ORANGE = "#1a1a1a", "#8a8a8a", "#2f6f9f", "#c8813a"

f, axes = plt.subplots(1, 3, figsize=(170 / 25.4, 58 / 25.4),
                       gridspec_kw=dict(width_ratios=[0.95, 1.30, 0.92], wspace=0.46))
f.subplots_adjust(left=0.062, right=0.978, top=0.80, bottom=0.21)

# ---------------- A ----------------
ax = axes[0]
ax.set_title("Paired 5\u2032\u20133\u2032 signal", fontsize=8.6, pad=6, loc="left", x=0.11)
vals = [(0.364, 0.000, "H3K9me3", ORANGE, "57.7%"), (0.0707, 0.0115, "DNase", BLUE, "98.6%")]
for i, (five, three, lab, c, share) in enumerate(vals):
    ax.bar(i - 0.18, five, width=0.34, color=c)
    ax.bar(i + 0.18, three, width=0.34, color="white", edgecolor=c, lw=1.0)
    fmt = "%.4f" if i == 1 else "%.3f"
    ax.text(i - 0.18, five + 0.012, fmt % five, ha="center", fontsize=7.0)
    ax.text(i + 0.18, three + 0.012, fmt % three, ha="center", fontsize=7.0)
ax.set_xticks([0, 1], ["H3K9me3\n5\u2032>3\u2032 in %s" % vals[0][4],
                      "DNase\n5\u2032>3\u2032 in %s" % vals[1][4]])
ax.set_ylabel("normalised signal", fontsize=7.6)
ax.set_xlim(-0.55, 1.55)
ax.set_ylim(0, 0.50)
ax.tick_params(labelsize=6.0, length=2, pad=1.6)
ax.text(0.98, 0.97, "after GC control,\nDNase 5\u2032>3\u2032 = 24.0%",
        transform=ax.transAxes, fontsize=6.4, color=BLUE, va="top", ha="right")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---------------- B ----------------
ax = axes[1]
ax.set_title("Composition-matched expectation", fontsize=8.6, pad=6, loc="left", x=0.07)
ORDER = [("k5xg3 (pool quantiles)", "GC 3\nbins"), ("k5xg8 (pool quantiles)", "GC 8\nbins"),
         ("k5 x GC 0.02 bins", "GC 0.02"), ("k5 x GC 0.01 bins", "GC 0.01"),
         ("1:1 balanced (caliper 0.10 SD)", "1:1\nmatched")]
for i, (key, lab) in enumerate(ORDER):
    for part, filled, col, dy in [("5'", True, INK, 0.15), ("3' (negative control)", False, GREY, 0.0)]:
        row = T[(T.design == key) & (T.part == part)]
        if not len(row):
            continue
        r = row.iloc[0]
        x = i + dy
        if np.isfinite(r.ci_lo):
            ax.plot([x, x], [r.ci_lo, r.ci_hi], color=col, lw=1.0)
        ax.plot([x], [r.ratio_medians], marker="o", ms=3.8,
                mfc=col if filled else "white", mec=col, mew=0.9, zorder=3)
ax.axhline(1.0, color=GREY, lw=0.7, ls=(0, (3, 2)))
ax.set_xticks(range(len(ORDER)), [l for _, l in ORDER])
ax.set_xlabel("composition control", fontsize=7.6, labelpad=2)
ax.set_ylabel("observed / expected", fontsize=7.6)
ax.set_xlim(-0.5, len(ORDER) - 0.35)
ax.set_ylim(0, 1.55)
ax.tick_params(labelsize=6.0, length=2, pad=1.6)
handles = [Line2D([], [], marker="o", ls="", ms=3.8, mfc=INK, mec=INK, label="5\u2032 UTR"),
           Line2D([], [], marker="o", ls="", ms=3.8, mfc="white", mec=GREY, label="3\u2032 end (control)")]
ax.legend(handles=handles, fontsize=6.6, frameon=False, loc="upper right",
          handletextpad=0.4, borderpad=0.1, labelspacing=0.3)

for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---------------- C ----------------
ax = axes[2]
ax.set_title("Reference assembly", fontsize=8.6, pad=6, loc="left", x=0.12)
vals = [(1.94, 1.59), (0.83, 0.48), (0.51, 0.02)]
for i, (g, h) in enumerate(vals):
    ax.bar(i - 0.18, g, width=0.34, color=BLUE)
    ax.bar(i + 0.18, h, width=0.34, color="white", edgecolor=BLUE, lw=1.0, hatch="////")
    ax.text(i - 0.18, g + 0.05, "%.2f" % g, ha="center", fontsize=6.8)
    ax.text(i + 0.18, h + 0.05, "%.2f" % h, ha="center", fontsize=6.8)
ax.set_xticks(range(3), ["H1", "HepG2", "K562"])
ax.set_ylabel("paired 5\u2032 \u2212 3\u2032 DNase", fontsize=7.2)
ax.set_ylim(-0.07, 2.35)
ax.set_xlim(-0.55, 2.55)
ax.tick_params(labelsize=6.0, length=2, pad=1.6)
handles = [Line2D([], [], marker="s", ls="", ms=4, mfc=BLUE, mec=BLUE, label="GRCh38"),
           Line2D([], [], marker="s", ls="", ms=4, mfc="white", mec=BLUE, label="hg19")]
ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right",
          handletextpad=0.4, borderpad=0.1, labelspacing=0.3)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

for ax_, lab in zip(axes, ["A", "B", "C"]):
    ax_.text(0.0, 1.005, lab, transform=ax_.transAxes, fontsize=10.5,
             fontweight="bold", va="bottom", ha="left")

out = os.path.join(ROOT, "figures", "fig6_composition.png")
f.savefig(out, dpi=600)
print("wrote", out)
