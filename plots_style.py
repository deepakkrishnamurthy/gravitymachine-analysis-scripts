# plot_style.py

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# --- Fluorescence palette (LUT-matching) ---
FL_PALETTE = {
    "FITC":     "#00FFFF",   # cyan
    "TRITC":    "#FF00FF",   # magenta
    "Jacalin":  "#00FF00",   # green
    "FarRed":   "#FFD700",   # gold/yellow
    "Hoechst":  "#1E90FF",   # bright blue
    "Membrane": "#FF8C00",   # orange
}

# Optional: slightly deeper tones for plots
FL_PALETTE_DARK = {
    "FITC":    "#00AACC",
    "TRITC":   "#CC00CC",
    "Jacalin": "#00AA00",
    "FarRed":  "#CCAA00",
}

# FL colormaps
FL_cmaps = {
    "FITC":    LinearSegmentedColormap.from_list("FITC_cmap", ["black", FL_PALETTE["FITC"]]),
    "TRITC":   LinearSegmentedColormap.from_list("TRITC_cmap", ["black", FL_PALETTE["TRITC"]]),
    "Jacalin": LinearSegmentedColormap.from_list("Jacalin_cmap", ["black", FL_PALETTE["Jacalin"]]),
    "FarRed":  LinearSegmentedColormap.from_list("FarRed_cmap", ["black", FL_PALETTE["FarRed"]]),
    "Hoechst": LinearSegmentedColormap.from_list("Hoechst_cmap", ["black", FL_PALETTE["Hoechst"]]),
    "Membrane":LinearSegmentedColormap.from_list("Membrane_cmap", ["black", FL_PALETTE["Membrane"]]),
}
REPLICATE_PALETTE = {
    1: "#0072B2",  # blue
    2: "#D55E00",  # vermillion
    3: "#009E73",  # green
}
# --- Global publication-quality style ---
FL_STYLE = {
    "figure.figsize": (6, 4),
    "figure.dpi": 150,
    "savefig.dpi": 300,

    # ---- fonts (THIS is the key part) ----
    "font.family": "sans-serif",
    "font.sans-serif": [
        "Arial",
        "Helvetica",
        "Helvetica Neue",
        "Nimbus Sans",
        "DejaVu Sans",  # final fallback (always present)
    ],

    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,

    # axes
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.0,
    "axes.titleweight": "bold",
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "axes.titlecolor": "black",

    # grid
    "axes.grid": True,
    "grid.linestyle": "--",
    "grid.color": "0.85",
    "grid.linewidth": 0.7,

    # lines
    "lines.linewidth": 2.5,
    "lines.markersize": 8,

    # ticks
    "xtick.major.size": 6,
    "xtick.major.width": 1.2,
    "ytick.major.size": 6,
    "ytick.major.width": 1.2,

    # legend
    "legend.frameon": False,

    # Illustrator / PDF safety
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

FL_STYLE_DARK = {
    "figure.figsize": (6, 4),
    "figure.dpi": 150,
    "savefig.dpi": 300,

    # background
    "figure.facecolor": "#000000",
    "axes.facecolor": "#000000",

    # fonts
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,

    # axes + spine colors
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.0,
    "axes.titleweight": "bold",
    "axes.edgecolor": "white",      # ← white spines
    "axes.labelcolor": "white",     # ← white labels
    "axes.titlecolor": "white",     # ← white title

    # ticks
    "xtick.color": "white",         # ← white tick labels
    "ytick.color": "white",
    "xtick.major.size": 6,
    "xtick.major.width": 1.2,
    "ytick.major.size": 6,
    "ytick.major.width": 1.2,

    # grid (subtle, non-distracting)
    "axes.grid": True,
    "grid.linestyle": "--",
    "grid.color": "0.3",            # ← dark gray grid line
    "grid.linewidth": 0.7,

    # lines
    "lines.linewidth": 2.5,
    "lines.markersize": 8,

    # legend
    "legend.frameon": False,
    "legend.edgecolor": "white",
    "text.color": "white",

    # Keep text as TEXT (not paths) in Illustrator
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

FL_STYLE_DARK_LARGE = {
    "figure.figsize": (8, 6),
    "figure.dpi": 150,
    "savefig.dpi": 300,

    # background
    "figure.facecolor": "#000000",
    "axes.facecolor": "#000000",

    # fonts: BIG & CLEAR
    "font.family": "sans-serif",
    "font.size": 18,          # base text
    "axes.labelsize": 20,     # axis labels
    "axes.titlesize": 22,     # title
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,

    # axes + spines
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 2.5,
    "axes.titleweight": "bold",
    "axes.edgecolor": "white",
    "axes.labelcolor": "white",
    "axes.titlecolor": "white",

    # ticks: bigger + thicker
    "xtick.color": "white",
    "ytick.color": "white",
    "xtick.major.size": 10,
    "xtick.major.width": 2.5,
    "ytick.major.size": 10,
    "ytick.major.width": 2.5,

    # grid
    "axes.grid": True,
    "grid.linestyle": "--",
    "grid.color": "0.4",
    "grid.linewidth": 1.2,

    # lines
    "lines.linewidth": 4.0,
    "lines.markersize": 12,

    # legend
    "legend.frameon": False,
    "legend.edgecolor": "white",
    "legend.fontsize": 24,
    "text.color": "white",

    # Illustrator-friendly
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def use_fl_style():
    """Apply the fluorescence plotting style globally."""
    plt.style.use("default")
    mpl.rcParams.update(FL_STYLE)

def use_fl_style_dark():
    """Apply the dark fluorescence plotting style globally."""
    plt.style.use("default")
    mpl.rcParams.update(FL_STYLE_DARK)

def use_fl_style_dark_large():
    """Apply the large dark fluorescence plotting style globally."""
    plt.style.use("default")
    mpl.rcParams.update(FL_STYLE_DARK_LARGE)