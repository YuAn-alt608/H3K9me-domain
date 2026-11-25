# -*- coding: utf-8 -*-
"""Shared loaders and style for the redrawn manuscript figures."""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "data", "tables")
TAB232 = TAB
BM = TAB

FAM = ["L1PA3", "L1PA4", "L1PA5", "L1PA6", "LTR12C",
       "MER11A", "MER11B", "MER11C", "MER41A", "MER41B"]
TRACKS = ["CBX3", "CBX5", "CBX1", "INPUT", "ETS1", "POLR2A", "NFYB", "TNRC18"]
TRACK_LAB = {"CBX3": "CBX3 (HP1γ)", "CBX5": "CBX5 (HP1α)", "CBX1": "CBX1 (HP1β)",
             "INPUT": "ChIP input", "ETS1": "ETS1", "POLR2A": "POLR2A",
             "NFYB": "NFYB", "TNRC18": "TNRC18"}
# internal table names are lower case for tnrc18
KEY = {"CBX3": "CBX3", "CBX5": "CBX5", "CBX1": "CBX1", "INPUT": "INPUT",
       "ETS1": "ETS1", "POLR2A": "POLR2A", "NFYB": "NFYB", "TNRC18": "tnrc18"}

# palette
ORANGE = "#C95B0C"
BLUE = "#17669E"
GREEN = "#1F8A76"
GREY = "#B9B9B9"
DGREY = "#5A5A5A"
PALE = "#EDEDED"
SHADE = "#F2F2F2"


def style(plt):
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7.5,
        "xtick.labelsize": 6.6,
        "ytick.labelsize": 6.6,
        "legend.fontsize": 6.6,
        "axes.linewidth": 0.6,
        "axes.edgecolor": "#8C8C8C",
        "axes.labelcolor": "#222222",
        "text.color": "#222222",
        "xtick.color": "#8C8C8C",
        "ytick.color": "#8C8C8C",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.2,
        "ytick.major.size": 2.2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 450,
        "figure.dpi": 130,
        "legend.frameon": False,
        "lines.linewidth": 1.1,
    })


def panel(ax, letter, title, x=-0.13, y=1.10):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")
    ax.text(x + 0.055, y, title, transform=ax.transAxes, fontsize=7.6,
            fontweight="bold", va="top", ha="left")


def r(label, **kw):
    return pd.read_csv(os.path.join(TAB, label), sep="\t", **kw)
