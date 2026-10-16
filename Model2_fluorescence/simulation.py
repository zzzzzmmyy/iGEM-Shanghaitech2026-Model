"""Generate the eight final fluorescence figures as PNG and SVG."""
from pathlib import Path
import argparse
import os

def main():
    import os
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--input", type=Path, default=Path("output"), help="Results from model.py")
    cli.add_argument("--out", type=Path, default=Path("figures"))
    args = cli.parse_args()
    INPUT = str(args.input.resolve())
    args.out.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(args.out.resolve() / ".mplconfig"))
    import os
    import json
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"

    HERE = os.path.dirname(os.path.abspath(__file__))
    OUT = str(args.out.resolve())
    os.makedirs(OUT, exist_ok=True)

    C_GFP = "#2a9d3f"       # green: GFP-positive
    C_MCH = "#d6336c"       # red: mCherry-positive
    C_CTRL = "#6c757d"      # gray: controls (PBS/blank)
    C_CROSS_G = "#74c69d"
    C_CROSS_M = "#f4a6c0"

    df = pd.read_csv(os.path.join(INPUT, "df920_tidy.csv"))


    def classify(label):
        if label == "PBS":
            return "PBS"
        if label == "空白":
            return "Blank"
        return label


    df["group"] = df["label"].apply(classify)

    order = ["Blank", "PBS",
             "mCherry(OD600=0.5)", "mCherry(OD600=1)", "mCherry(OD600=5)", "mCherry(OD600=10)",
             "GFP(OD600=0.5)", "GFP(OD600=1)", "GFP(OD600=5)", "GFP(OD600=10)"]


    def group_color(g):
        if g in ("Blank", "PBS"):
            return C_CTRL
        if g.startswith("mCherry"):
            return C_MCH
        if g.startswith("GFP"):
            return C_GFP
        return "black"


    # ============ Figure 1: channel specificity by sample group ============
    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 10, "axes.titlesize": 12, "axes.linewidth": 0.8,
        "figure.facecolor": "#000000", "axes.facecolor": "#000000",
        "savefig.facecolor": "#000000", "axes.edgecolor": "#9CA3AF",
        "text.color": "#F3F4F6", "axes.labelcolor": "#E5E7EB",
        "xtick.color": "#D1D5DB", "ytick.color": "#D1D5DB",
        "svg.fonttype": "none",
    }):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))
        colors = {"Control": "#9CA3AF", "GFP": "#56B4E9", "mCherry": "#CC79A7"}
        short_labels = [g if g in ("Blank", "PBS") else
                        g.split("(")[0] + "\nOD " + g.split("=")[1].rstrip(")") for g in order]
        for ax, col, title, yscale in [
            (axes[0], "GFP_raw", "(A) GFP", "log"),
            (axes[1], "mCherry_raw", "(B) mCherry", "linear"),
        ]:
            for i, g in enumerate(order):
                sub = df[df["group"] == g][col].dropna().values
                if len(sub) == 0:
                    continue
                color = colors["Control" if g in ("Blank", "PBS") else g.split("(")[0]]
                jitter = np.random.default_rng(0).normal(0, 0.06, len(sub))
                ax.scatter(np.full(len(sub), i) + jitter, sub, color=color, s=42, alpha=0.85,
                           zorder=3, edgecolor="#000000", linewidth=0.4)
                mean = sub.mean()
                sd = sub.std(ddof=1) if len(sub) > 1 else 0
                ax.errorbar([i], [mean], yerr=[sd], fmt="_", color="#F3F4F6", markersize=18,
                            markeredgewidth=1.6, elinewidth=1.2, capsize=3, zorder=4)
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels(short_labels, rotation=45, ha="right", fontsize=8)
            ax.set_ylabel("Fluorescence (a.u.)" + (" · log scale" if yscale == "log" else ""))
            ax.set_title(title, loc="left", fontweight="bold", pad=12, y=1.04)
            if yscale == "log":
                ax.set_yscale("log")
                ax.tick_params(axis="y", which="minor", length=0)
            else:
                ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
            ax.set_axisbelow(True)
            ax.grid(axis="y", color="#FFFFFF", alpha=0.12, linewidth=0.65)
            ax.tick_params(axis="x", length=0, pad=8)
            ax.tick_params(axis="y", labelsize=9, width=0.7, length=3)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
        legend_el = [
            Line2D([0], [0], marker="o", linestyle="none", color=colors["Control"], markersize=6, label="Controls"),
            Line2D([0], [0], marker="o", linestyle="none", color=colors["mCherry"], markersize=6, label="mCherry-positive"),
            Line2D([0], [0], marker="o", linestyle="none", color=colors["GFP"], markersize=6, label="GFP-positive"),
            Line2D([0], [0], marker="_", linestyle="none", color="#F3F4F6", markersize=12, markeredgewidth=1.6, label="Mean ± SD"),
        ]
        fig.suptitle("Channel specificity", fontsize=15, fontweight="bold", y=0.985)
        fig.legend(handles=legend_el, loc="upper center", bbox_to_anchor=(0.5, 0.915),
                   ncol=4, fontsize=10, frameon=False, columnspacing=2)
        fig.subplots_adjust(left=0.065, right=0.985, bottom=0.24, top=0.73, wspace=0.22)
        for extension in ("png", "svg"):
            fig.savefig(os.path.join(OUT, f"figure1_channel_specificity.{extension}"),
                        dpi=300, bbox_inches="tight")
        plt.close(fig)
    print("saved figure1_channel_specificity (PNG + SVG)")

    # ============ Figure 2: dose-response + crosstalk (4 panels) ============
    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 10, "axes.titlesize": 12, "axes.linewidth": 0.8,
        "figure.facecolor": "#000000", "axes.facecolor": "#000000",
        "savefig.facecolor": "#000000", "axes.edgecolor": "#9CA3AF",
        "text.color": "#F3F4F6", "axes.labelcolor": "#E5E7EB",
        "xtick.color": "#D1D5DB", "ytick.color": "#D1D5DB",
        "svg.fonttype": "none",
    }):
        fig, axes = plt.subplots(2, 2, figsize=(12, 9))

        gfp_pos = df[df["label"].astype(str).str.startswith("GFP")].copy()
        mch_pos = df[df["label"].astype(str).str.startswith("mCherry")].copy()
        ctrl = df[df["group"].isin(["PBS", "Blank"])].copy()

        ax = axes[0][0]
        ax.scatter(gfp_pos["OD600"], gfp_pos["GFP_raw"], color="#56B4E9", s=42, zorder=3, edgecolor="#000000", linewidth=0.4)
        ax.scatter(ctrl["OD600"], ctrl["GFP_raw"], color="#9CA3AF", s=32, zorder=2, alpha=0.7, marker="x", label="Controls")
        ax.set_xlabel("OD600")
        ax.set_ylabel("GFP fluorescence (a.u.)")
        ax.set_title("(A) GFP response", loc="left", fontsize=12, fontweight="bold", pad=12, y=1.07)
        ax.set_yscale("log")
        ax.set_ylabel("GFP fluorescence (a.u.) · log scale")
        ax.tick_params(axis="y", which="minor", length=0)
        ax.set_axisbelow(True)
        ax.grid(color="#FFFFFF", alpha=0.12, linewidth=0.65)
        ax.tick_params(labelsize=9, width=0.7, length=3)
        ax.legend(fontsize=9, frameon=False)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        ax = axes[0][1]
        ax.scatter(mch_pos["OD600"], mch_pos["mCherry_raw"], color="#CC79A7", s=42, zorder=3, edgecolor="#000000", linewidth=0.4)
        ax.scatter(ctrl["OD600"], ctrl["mCherry_raw"], color="#9CA3AF", s=32, zorder=2, alpha=0.7, marker="x", label="Controls")
        ax.set_xlabel("OD600")
        ax.set_ylabel("mCherry fluorescence (a.u.)")
        ax.set_title("(B) mCherry response", loc="left", fontsize=12, fontweight="bold", pad=12, y=1.07)
        ax.set_axisbelow(True)
        ax.grid(color="#FFFFFF", alpha=0.12, linewidth=0.65)
        ax.tick_params(labelsize=9, width=0.7, length=3)
        ax.legend(fontsize=9, frameon=False)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        ax = axes[1][0]
        ax.scatter(mch_pos["OD600"], mch_pos["GFP_raw"], color="#CC79A7", s=42, zorder=3,
                   edgecolor="#CC79A7", linewidth=1.2, label="mCherry-positive")
        blank_gfp_mean = ctrl["GFP_raw"].mean()
        ax.axhline(blank_gfp_mean, color="#9CA3AF", ls="--", lw=1.3, label="Control mean")
        ax.set_xlabel("OD600")
        ax.set_ylabel("GFP fluorescence (a.u.)")
        ax.set_title("(C) mCherry → GFP", loc="left", fontsize=12, fontweight="bold", pad=12, y=1.07)
        ax.set_axisbelow(True)
        ax.grid(color="#FFFFFF", alpha=0.12, linewidth=0.65)
        ax.tick_params(labelsize=9, width=0.7, length=3)
        ax.legend(fontsize=9, frameon=False)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        ax = axes[1][1]
        ax.scatter(gfp_pos["OD600"], gfp_pos["mCherry_raw"], color="#56B4E9", s=42, zorder=3,
                   edgecolor="#56B4E9", linewidth=1.2, label="GFP-positive")
        blank_mch_mean = ctrl["mCherry_raw"].mean()
        ax.axhline(blank_mch_mean, color="#9CA3AF", ls="--", lw=1.3, label="Control mean")
        ax.set_xlabel("OD600")
        ax.set_ylabel("mCherry fluorescence (a.u.)")
        ax.set_title("(D) GFP → mCherry", loc="left", fontsize=12, fontweight="bold", pad=12, y=1.07)
        ax.set_axisbelow(True)
        ax.grid(color="#FFFFFF", alpha=0.12, linewidth=0.65)
        ax.tick_params(labelsize=9, width=0.7, length=3)
        ax.legend(fontsize=9, frameon=False)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        for ax in (axes[0][1], axes[1][0], axes[1][1]):
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        fig.suptitle("Dose response and crosstalk", fontsize=15, fontweight="bold", y=0.985)
        fig.subplots_adjust(left=0.08, right=0.985, bottom=0.08, top=0.85, wspace=0.27, hspace=0.48)
        for extension in ("png", "svg"):
            fig.savefig(os.path.join(OUT, f"figure2_doseresponse_crosstalk.{extension}"),
                        dpi=300, bbox_inches="tight")
        plt.close(fig)
    print("saved figure2_doseresponse_crosstalk (PNG + SVG)")



    import os
    import pandas as pd
    import numpy as np
    from scipy import stats
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"

    HERE = os.path.dirname(os.path.abspath(__file__))
    OUT = str(args.out.resolve())
    os.makedirs(OUT, exist_ok=True)

    C_GFP = "#2a9d3f"
    C_MCH = "#d6336c"
    C_CTRL = "#6c757d"
    C_EXTRA = "#f0a000"

    df = pd.read_csv(os.path.join(INPUT, "df921_tidy.csv"))


    def is_numeric_label(x):
        try:
            float(x)
            return True
        except (ValueError, TypeError):
            return False


    df["is_ctrl"] = df["label"].isin(["PBS", "空白"])
    df["is_frac"] = df["label"].apply(is_numeric_label)
    main_grid = {round(0.05 * i, 2) for i in range(21)}
    df.loc[df["is_frac"], "fraction"] = df.loc[df["is_frac"], "label"].astype(float)
    df["is_extra"] = df["is_frac"] & (~df["fraction"].isin(main_grid))

    frac_df = df[df["is_frac"]].copy()

    # ============ Fig 1: dose-response (GFP / mCherry channel vs GFP fraction) ============
    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 10, "axes.titlesize": 12, "axes.linewidth": 0.8,
        "figure.facecolor": "#000000", "axes.facecolor": "#000000",
        "savefig.facecolor": "#000000", "axes.edgecolor": "#9CA3AF",
        "text.color": "#F3F4F6", "axes.labelcolor": "#E5E7EB",
        "xtick.color": "#D1D5DB", "ytick.color": "#D1D5DB",
        "svg.fonttype": "none",
    }):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
        main_color, test_color = "#56B4E9", "#CC79A7"
        for ax, col, title in [
            (axes[0], "GFP_raw", "(A) GFP"),
            (axes[1], "mCherry_raw", "(B) mCherry"),
        ]:
            main = frac_df[~frac_df["is_extra"]]
            extra = frac_df[frac_df["is_extra"]]
            g = main.groupby("fraction")[col].agg(["mean", "std"]).reset_index().sort_values("fraction")
            ax.plot(g["fraction"], g["mean"], "-", color=main_color, lw=1.6, alpha=0.8, zorder=2)
            ax.fill_between(g["fraction"], g["mean"] - g["std"], g["mean"] + g["std"],
                            color=main_color, alpha=0.16, zorder=1)
            ax.scatter(main["fraction"], main[col], color=main_color, s=38, alpha=0.85,
                       zorder=3, edgecolor="#000000", linewidth=0.4, label="Calibration")
            ax.scatter(extra["fraction"], extra[col], color=test_color, s=70, alpha=0.95,
                       zorder=4, marker="^", edgecolor="#000000", linewidth=0.5, label="Held-out test")
            ax.set_xlabel("GFP fraction")
            ax.set_ylabel("Fluorescence (a.u.)")
            ax.set_title(title, loc="left", fontweight="bold", pad=12)
            ax.set_axisbelow(True)
            ax.grid(color="#FFFFFF", alpha=0.12, linewidth=0.65)
            ax.tick_params(labelsize=9, width=0.7, length=3)
            ax.set_xticks(np.linspace(0, 1, 6))
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

        handles, labels = axes[0].get_legend_handles_labels()
        fig.suptitle("Fluorescence response", fontsize=15, fontweight="bold", y=0.985)
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.915),
                   ncol=2, frameon=False, fontsize=10, columnspacing=2.5)
        fig.subplots_adjust(left=0.07, right=0.985, bottom=0.16, top=0.76, wspace=0.25)
        for extension in ("png", "svg"):
            fig.savefig(os.path.join(OUT, f"figure1_dose_response.{extension}"), dpi=300, bbox_inches="tight")
        plt.close(fig)
    print("saved figure1_dose_response (PNG + SVG)")

    # ============ Fig 2: blank baseline check + OD confound check ============
    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 10, "axes.titlesize": 12, "axes.linewidth": 0.8,
        "figure.facecolor": "#000000", "axes.facecolor": "#000000",
        "savefig.facecolor": "#000000", "axes.edgecolor": "#9CA3AF",
        "text.color": "#F3F4F6", "axes.labelcolor": "#E5E7EB",
        "xtick.color": "#D1D5DB", "ytick.color": "#D1D5DB",
        "svg.fonttype": "none",
    }):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
        main_color, test_color = "#56B4E9", "#CC79A7"
        for ax in axes:
            ax.set_axisbelow(True)
            ax.grid(axis="y", color="#FFFFFF", alpha=0.12, linewidth=0.65)
            ax.tick_params(labelsize=9, width=0.7, length=3)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

        ax = axes[0]
        blank_v = df[df["label"] == "空白"]["GFP_raw"].values
        pbs_v = df[df["label"] == "PBS"]["GFP_raw"].values
        f0_v = frac_df[frac_df["fraction"] == 0.0]["GFP_raw"].values
        groups = [("Blank", blank_v, main_color),
                  ("PBS", pbs_v, main_color),
                  ("mCherry only", f0_v, test_color)]
        for i, (name, vals, color) in enumerate(groups):
            jitter = np.random.default_rng(1).normal(0, 0.05, len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals, color=color, s=55,
                       zorder=3, edgecolor="#000000", linewidth=0.4)
            ax.errorbar([i], [vals.mean()], yerr=[vals.std(ddof=1) if len(vals) > 1 else 0],
                        fmt="_", color="#F3F4F6", markersize=20, markeredgewidth=1.6,
                        elinewidth=1.2, capsize=4, zorder=4)
        t1, p1 = stats.ttest_ind(blank_v, f0_v, equal_var=False)
        t2, p2 = stats.ttest_ind(pbs_v, f0_v, equal_var=False)
        all_values = np.concatenate([blank_v, pbs_v, f0_v])
        span = np.ptp(all_values)
        for left, right, level, pvalue in [(1, 2, 0.14, p2), (0, 2, 0.34, p1)]:
            height = all_values.max() + level * span
            cap = 0.035 * span
            ax.plot([left, left, right, right], [height-cap, height, height, height-cap],
                    color="#D1D5DB", linewidth=0.8)
            ax.text((left+right)/2, height+0.025*span, f"p = {pvalue:.4f}",
                    ha="center", va="bottom", fontsize=9)
        ax.set_ylim(all_values.min()-0.12*span, all_values.max()+0.57*span)
        ax.set_xlim(-0.45, 2.45)
        ax.set_xticks([0, 1, 2], labels=[g[0] for g in groups])
        ax.set_ylabel("GFP fluorescence (a.u.)")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        ax.set_title("(A) GFP background", loc="left", fontweight="bold", pad=12, y=1.06)

        ax = axes[1]
        non_ctrl = frac_df.copy()
        main2 = non_ctrl[~non_ctrl["is_extra"]]
        extra2 = non_ctrl[non_ctrl["is_extra"]]
        ax.scatter(main2["fraction"], main2["OD600"], color=main_color, s=38, alpha=0.85,
                   zorder=3, edgecolor="#000000", linewidth=0.4, label="Calibration")
        ax.scatter(extra2["fraction"], extra2["OD600"], color=test_color, s=70, marker="^",
                   zorder=4, edgecolor="#000000", linewidth=0.5, label="Held-out test")
        blank_od = df[df["is_ctrl"]]["OD600"]
        ax.axhspan(blank_od.min(), blank_od.max(), color="#9CA3AF", alpha=0.18,
                   label="Blank/PBS range")
        r, p = stats.pearsonr(non_ctrl["fraction"], non_ctrl["OD600"])
        ax.text(0.04, 0.95, f"r = {r:.3f}, p = {p:.3f}", transform=ax.transAxes,
                ha="left", va="top", fontsize=9)
        ax.set_xlabel("GFP fraction")
        ax.set_ylabel("OD600")
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.set_title("(B) Cell density", loc="left", fontweight="bold", pad=12, y=1.06)
        ax.legend(fontsize=9, frameon=False, loc="lower right")
        fig.suptitle("Background and cell density", fontsize=15, fontweight="bold", y=0.985)
        fig.subplots_adjust(left=0.08, right=0.985, bottom=0.16, top=0.78, wspace=0.28)
        for extension in ("png", "svg"):
            fig.savefig(os.path.join(OUT, f"figure2_baseline_od_confound.{extension}"),
                        dpi=300, bbox_inches="tight")
        plt.close(fig)
    print("saved figure2_baseline_od_confound (PNG + SVG)")
    print(f"blank mean={blank_v.mean():.0f} sd={blank_v.std(ddof=1):.0f}")
    print(f"pbs mean={pbs_v.mean():.0f} sd={pbs_v.std(ddof=1):.0f}")
    print(f"f0 mean={f0_v.mean():.0f} sd={f0_v.std(ddof=1):.0f}")
    print(f"t-test blank vs f0: t={t1:.3f} p={p1:.5f}")
    print(f"t-test pbs vs f0: t={t2:.3f} p={p2:.5f}")
    print(f"OD vs fraction pearson r={r:.4f} p={p:.4f}")

    # ============ Fig 4: ratio-transform preview GFP/(GFP+mCherry) vs fraction ============
    work = frac_df.copy()
    work["ratio_raw"] = work["GFP_raw"] / (work["GFP_raw"] + work["mCherry_raw"])
    work["denom"] = work["GFP_raw"] + work["mCherry_raw"]

    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 10, "axes.titlesize": 12, "axes.linewidth": 0.8,
        "figure.facecolor": "#000000", "axes.facecolor": "#000000",
        "savefig.facecolor": "#000000", "axes.edgecolor": "#9CA3AF",
        "text.color": "#F3F4F6", "axes.labelcolor": "#E5E7EB",
        "xtick.color": "#D1D5DB", "ytick.color": "#D1D5DB",
        "svg.fonttype": "none",
    }):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
        main_color, test_color = "#56B4E9", "#CC79A7"
        for ax in axes:
            ax.set_axisbelow(True)
            ax.grid(color="#FFFFFF", alpha=0.12, linewidth=0.65)
            ax.tick_params(labelsize=9, width=0.7, length=3)
            ax.set_xlabel("GFP fraction")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

        ax = axes[0]
        valid = work[work["fraction"] > 0]
        main3 = valid[~valid["is_extra"]]
        extra3 = valid[valid["is_extra"]]
        ax.scatter(main3["fraction"], main3["ratio_raw"], color=main_color, s=38, alpha=0.85,
                   zorder=3, edgecolor="#000000", linewidth=0.4, label="Calibration (f > 0)")
        ax.scatter(extra3["fraction"], extra3["ratio_raw"], color=test_color, s=70, marker="^",
                   zorder=4, edgecolor="#000000", linewidth=0.5, label="Held-out test")
        ax.plot([0, 1], [0, 1], color="#9CA3AF", lw=1, ls=(0, (4, 4)),
                label="y = x", zorder=2)
        r_valid, p_valid = stats.pearsonr(valid["fraction"], valid["ratio_raw"])
        rho_valid, prho_valid = stats.spearmanr(valid["fraction"], valid["ratio_raw"])
        ax.text(0.04, 0.65, f"Pearson r = {r_valid:.3f}\nSpearman ρ = {rho_valid:.3f}",
                transform=ax.transAxes, ha="left", va="top", fontsize=9, linespacing=1.6)
        ax.set_ylabel("GFP / (GFP + mCherry)")
        ax.set_title("(A) Signal ratio", loc="left", fontweight="bold", pad=12, y=1.06)
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.legend(fontsize=9, frameon=False, loc="lower right")

        ax = axes[1]
        allw = work.sort_values("fraction")
        # Preserve the low-fraction x range; scale y to the data visible in this window.
        for is_zero, color, label in [(False, main_color, "f > 0"),
                                      (True, test_color, "mCherry only")]:
            sub = allw[(allw["fraction"] == 0) == is_zero]
            ax.scatter(sub["fraction"], sub["denom"], color=color, s=45, alpha=0.85,
                       zorder=3, edgecolor="#000000", linewidth=0.4, label=label)
        ax.set_xlim(-0.03, 0.35)
        visible = allw[allw["fraction"] <= 0.35]["denom"]
        margin = 0.07 * (visible.max() - visible.min())
        ax.set_ylim(visible.min() - margin, visible.max() + margin)
        ax.set_xticks([0, 0.1, 0.2, 0.3])
        ax.set_ylabel("GFP + mCherry (a.u.)")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        ax.set_title("(B) Denominator", loc="left", fontweight="bold", pad=12, y=1.06)
        ax.legend(fontsize=9, frameon=False, loc="upper left")
        fig.suptitle("Fluorescence ratio", fontsize=15, fontweight="bold", y=0.985)
        fig.subplots_adjust(left=0.08, right=0.985, bottom=0.16, top=0.78, wspace=0.28)
        for extension in ("png", "svg"):
            fig.savefig(os.path.join(OUT, f"figure3_ratio_preview.{extension}"),
                        dpi=300, bbox_inches="tight")
        plt.close(fig)
    print("saved figure3_ratio_preview (PNG + SVG), ratio r (f>0) =", r_valid, "rho=", rho_valid)
    f0_denom = work[work["fraction"] == 0]["denom"].values
    print("f=0 denom values:", f0_denom)


    import os
    import json
    import pandas as pd
    import numpy as np
    from sklearn.linear_model import LinearRegression
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.metrics import r2_score
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"

    HERE = os.path.dirname(os.path.abspath(__file__))
    OUT = str(args.out.resolve())

    C_TRAIN = "#2a78d6"
    C_TEST = "#eb6834"

    results_df = pd.read_pickle(os.path.join(INPUT, "results_df.pkl"))
    train = pd.read_pickle(os.path.join(INPUT, "train_921.pkl"))
    test = pd.read_pickle(os.path.join(INPUT, "test_921.pkl"))


    PAPER_STYLE = {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 12,
        "axes.linewidth": 0.8,
        "figure.facecolor": "#000000",
        "axes.facecolor": "#000000",
        "savefig.facecolor": "#000000",
        "axes.edgecolor": "#9CA3AF",
        "text.color": "#F3F4F6",
        "axes.labelcolor": "#E5E7EB",
        "xtick.color": "#D1D5DB",
        "ytick.color": "#D1D5DB",
        "svg.fonttype": "none",
    }
    MODEL_COLORS = ("#56B4E9", "#CC79A7")

    # ============ model comparison ============
    name_short = {
        "M1_GFP_only_linear": "M1\nGFP\nlinear",
        "M2_mCherry_only_linear": "M2\nmCherry\nlinear",
        "M3_two_channel_linear": "M3\nTwo-channel\nlinear",
        "M4_ratio_linear": "M4\nRatio\nlinear",
        "M5_ratio_quadratic": "M5\nRatio\nquadratic",
        "M6_GFP_quadratic": "M6\nGFP\nquadratic",
        "M7_random_forest_two_channel": "M7\nRandom forest\n(two-channel)",
    }
    results_df["short"] = results_df["model"].map(name_short)

    # Apply the shared paper style locally to each figure.
    with plt.rc_context(PAPER_STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
        # Sky blue and reddish purple from the Okabe-Ito scientific palette.
        colors = MODEL_COLORS
        x = np.arange(len(results_df))
        w = 0.34

        for ax in axes:
            ax.set_axisbelow(True)
            ax.grid(axis="y", color="#FFFFFF", alpha=0.12, linewidth=0.65)
            ax.set_xticks(x)
            ax.set_xticklabels(results_df["short"], fontsize=7.5)
            ax.tick_params(axis="x", length=0, pad=9)
            ax.tick_params(axis="y", labelsize=9, width=0.7, length=3)
            ax.set_xlim(-0.65, len(results_df) - 0.35)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

        ax = axes[0]
        r2_cv_clipped = results_df["R2_CV"].clip(lower=-0.3)
        r2_test_clipped = results_df["R2_test"].clip(lower=-0.3)
        ax.bar(x - w/2, r2_cv_clipped, width=w, color=colors[0], label="Cross-validation (LOGO)")
        ax.bar(x + w/2, r2_test_clipped, width=w, color=colors[1], label="External test")
        ax.axhline(0, color="#9CA3AF", lw=0.8)
        ax.axhline(1, color="#9CA3AF", lw=0.8, ls=(0, (3, 3)))
        ax.set_ylabel("R²")
        ax.set_ylim(-0.4, 1.08)
        ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_title("(A) R²", loc="left", fontweight="bold", pad=12)

        ax = axes[1]
        ax.bar(x - w/2, results_df["RMSE_CV"], width=w, color=colors[0])
        ax.bar(x + w/2, results_df["RMSE_test"], width=w, color=colors[1])
        # A log axis reveals small errors while retaining every original RMSE value.
        ax.set_yscale("log")
        ax.set_ylim(0.01, 15)
        ax.set_yticks([0.01, 0.1, 1, 10], labels=["0.01", "0.1", "1", "10"])
        ax.tick_params(axis="y", which="minor", length=0)
        ax.set_ylabel("RMSE (log scale)")
        ax.set_title("(B) RMSE", loc="left", fontweight="bold", pad=12)

        handles, labels = axes[0].get_legend_handles_labels()
        fig.suptitle("Model comparison", fontsize=15, fontweight="bold", y=0.985)
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.915),
                   ncol=2, frameon=False, fontsize=10, handlelength=1.4, columnspacing=2.5)
        fig.subplots_adjust(left=0.07, right=0.985, bottom=0.16, top=0.76, wspace=0.22)
        for extension in ("png", "svg"):
            filename = f"figure2-1_model_comparison.{extension}"
            fig.savefig(os.path.join(OUT, filename), dpi=300, bbox_inches="tight")
        plt.close(fig)
    print("saved figure2-1_model_comparison (PNG + SVG)")
    print(results_df[["model", "R2_CV", "RMSE_CV", "R2_test", "RMSE_test"]].to_string(index=False))

    # ============ Recompute best-model (M3) OOF + test predictions for diagnostics ============
    Xtrain = train[["GFPc", "mCherryc"]].values
    ytrain = train["fraction"].values
    groups = train["fraction"].values
    Xtest = test[["GFPc", "mCherryc"]].values
    ytest = test["fraction"].values

    logo = LeaveOneGroupOut()
    oof_pred_m3 = np.zeros(len(ytrain))
    for tr_idx, va_idx in logo.split(Xtrain, ytrain, groups):
        m = LinearRegression()
        m.fit(Xtrain[tr_idx], ytrain[tr_idx])
        oof_pred_m3[va_idx] = m.predict(Xtrain[va_idx])

    m3_full = LinearRegression().fit(Xtrain, ytrain)
    pred_test_m3 = m3_full.predict(Xtest)
    print(f"\nM3 two-channel linear regression coefficients: GFPc={m3_full.coef_[0]:.6e}, "
          f"mCherryc={m3_full.coef_[1]:.6e}, intercept={m3_full.intercept_:.4f}")

    # ============ best-model diagnostics ============
    resid_train = oof_pred_m3 - ytrain
    resid_test = pred_test_m3 - ytest
    band = np.std(resid_train)

    with plt.rc_context(PAPER_STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
        colors = MODEL_COLORS

        for ax in axes:
            ax.set_axisbelow(True)
            ax.grid(color="#FFFFFF", alpha=0.12, linewidth=0.65)
            ax.tick_params(labelsize=9, width=0.7, length=3)
            ax.set_xlabel("True GFP fraction")
            ax.set_xlim(-0.08, 1.08)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)

        ax = axes[0]
        ideal, = ax.plot([0, 1], [0, 1], color="#9CA3AF", lw=1,
                        ls=(0, (4, 4)), zorder=1, label="Ideal prediction")
        cv_points = ax.scatter(ytrain, oof_pred_m3, color=colors[0], s=38, alpha=0.85,
                               zorder=3, edgecolor="#000000", linewidth=0.4,
                               label="Cross-validation (LOGO)")
        test_points = ax.scatter(ytest, pred_test_m3, color=colors[1], s=70, marker="^",
                                 zorder=4, edgecolor="#000000", linewidth=0.5,
                                 label="External test")
        ax.set_ylabel("Predicted GFP fraction")
        ax.set_ylim(-0.08, 1.08)
        ax.set_title("(A) Predicted vs. true", loc="left", fontweight="bold", pad=12)
        ax.legend(handles=[ideal], frameon=False, fontsize=9, loc="upper left")

        ax = axes[1]
        ax.axhline(0, color="#9CA3AF", lw=1, ls=(0, (4, 4)), zorder=2)
        ax.scatter(ytrain, resid_train, color=colors[0], s=38, alpha=0.85,
                   zorder=3, edgecolor="#000000", linewidth=0.4)
        ax.scatter(ytest, resid_test, color=colors[1], s=70, marker="^",
                   zorder=4, edgecolor="#000000", linewidth=0.5)
        ax.axhspan(-band, band, color=colors[0], alpha=0.14, zorder=0, label="CV ±1 SD")
        ax.set_ylabel("Residual (predicted − true)")
        ax.set_title("(B) Residuals", loc="left", fontweight="bold", pad=12)
        ax.legend(frameon=False, fontsize=9, loc="upper right")

        fig.suptitle("Best-model diagnostics (M3)", fontsize=15, fontweight="bold", y=0.985)
        fig.legend(handles=[cv_points, test_points], loc="upper center",
                   bbox_to_anchor=(0.5, 0.915), ncol=2, frameon=False,
                   fontsize=10, handlelength=1.4, columnspacing=2.5)
        fig.subplots_adjust(left=0.07, right=0.985, bottom=0.16, top=0.76, wspace=0.25)
        for extension in ("png", "svg"):
            filename = f"figure2-2_best_model_diagnostics.{extension}"
            fig.savefig(os.path.join(OUT, filename), dpi=300, bbox_inches="tight")
        plt.close(fig)
    print("saved figure2-2_best_model_diagnostics (PNG + SVG)")
    print(f"train residual: mean={resid_train.mean():.4f} sd={resid_train.std(ddof=1):.4f}")
    print(f"test residual: mean={resid_test.mean():.4f} sd={resid_test.std(ddof=1):.4f}")

    # ============ calibration curve ============
    resid_sd = resid_train.std(ddof=1)
    xs = np.linspace(0, 1, 100)
    rmse_all_test = np.sqrt(np.mean((pred_test_m3 - ytest) ** 2))
    mae_all_test = np.mean(np.abs(pred_test_m3 - ytest))

    with plt.rc_context(PAPER_STYLE):
        fig, ax = plt.subplots(figsize=(7.2, 6.4))
        colors = MODEL_COLORS
        ideal, = ax.plot([0, 1], [0, 1], color="#9CA3AF", lw=1,
                        ls=(0, (4, 4)), label="Ideal calibration", zorder=1)
        cv_points = ax.scatter(ytrain, oof_pred_m3, color=colors[0], s=38, alpha=0.85,
                               zorder=3, edgecolor="#000000", linewidth=0.4,
                               label="Cross-validation (LOGO)")
        test_points = ax.scatter(ytest, pred_test_m3, color=colors[1], s=70, marker="^",
                                 zorder=4, edgecolor="#000000", linewidth=0.5,
                                 label="External test")
        residual_band = ax.fill_between(xs, xs - 1.96 * resid_sd, xs + 1.96 * resid_sd,
                                        color=colors[0], alpha=0.14,
                                        label="CV ±1.96 SD", zorder=0)
        ax.set_xlabel("True GFP fraction")
        ax.set_ylabel("Predicted GFP fraction")
        ax.set_axisbelow(True)
        ax.grid(color="#FFFFFF", alpha=0.12, linewidth=0.65)
        ax.tick_params(labelsize=9, width=0.7, length=3)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.1, 1.1)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.legend(handles=[ideal, residual_band], fontsize=9, frameon=False, loc="lower right")
        fig.suptitle("Calibration curve (M3)", fontsize=15, fontweight="bold", y=0.985)
        fig.legend(handles=[cv_points, test_points], loc="upper center",
                   bbox_to_anchor=(0.5, 0.915), ncol=2, frameon=False,
                   fontsize=10, handlelength=1.4, columnspacing=2)
        fig.subplots_adjust(left=0.12, right=0.975, bottom=0.12, top=0.79)
        for extension in ("png", "svg"):
            filename = f"figure2-3_calibration_curve.{extension}"
            fig.savefig(os.path.join(OUT, filename), dpi=300, bbox_inches="tight")
        plt.close(fig)
    print("saved figure2-3_calibration_curve (PNG + SVG)")

    with open(os.path.join(INPUT, "model_blank_baseline.json")) as f:
        baseline = json.load(f)

    summary = {
        "blank_baseline_GFP": baseline["blank_gfp"],
        "blank_baseline_mCherry": baseline["blank_mch"],
        "M3_coef_GFPc": float(m3_full.coef_[0]),
        "M3_coef_mCherryc": float(m3_full.coef_[1]),
        "M3_intercept": float(m3_full.intercept_),
        "M3_train_cv_R2": float(r2_score(ytrain, oof_pred_m3)),
        "M3_test_R2": float(r2_score(ytest, pred_test_m3)),
        "M3_test_RMSE": float(rmse_all_test),
        "M3_test_MAE": float(mae_all_test),
    }
    with open(os.path.join(INPUT, "final_model_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
