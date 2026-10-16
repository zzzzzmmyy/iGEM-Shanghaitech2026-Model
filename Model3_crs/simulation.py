"""Reproduce four C:R:S Wiki figures.

Run: python simulation.py --out output
Check existing outputs: python simulation.py --out output --check
"""
import argparse
import os
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import model as crs_model
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from model import DEFAULT_PARAMS, simulate_batch
BATCH_OUTPUT_STEM = 'figure3-1_trajectories'
BATCH_BACKGROUND = '#000000'
BATCH_FOREGROUND = '#E8ECF1'
BATCH_MUTED = '#ABB4C0'
BATCH_RATIOS = (('Equal inoculum', '1:1:1', (0.01 / 3, 0.01 / 3, 0.01 / 3)), ('High producer', '5:1:1', (0.01 * 5 / 7, 0.01 / 7, 0.01 / 7)), ('Low producer', '0.1:1:1', (0.01 * 0.1 / 2.1, 0.01 / 2.1, 0.01 / 2.1)), ('Very low producer', '0.02:1:1', (0.01 * 0.02 / 2.02, 0.01 / 2.02, 0.01 / 2.02)))
BATCH_STRAINS = (('C (colicin producer)', '#F3AC69', '-'), ('R (colicin resistant)', '#76B7F2', '--'), ('S (colicin sensitive)', '#79D3AF', '-.'))
BATCH_STYLE = {'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 11, 'text.color': BATCH_FOREGROUND, 'axes.facecolor': BATCH_BACKGROUND, 'figure.facecolor': BATCH_BACKGROUND, 'savefig.facecolor': BATCH_BACKGROUND, 'axes.labelcolor': BATCH_FOREGROUND, 'axes.edgecolor': '#67717E', 'axes.linewidth': 0.8, 'axes.titleweight': 'normal', 'xtick.color': BATCH_MUTED, 'ytick.color': BATCH_MUTED, 'xtick.labelsize': 10, 'ytick.labelsize': 10, 'grid.color': '#56606E', 'grid.linewidth': 0.6, 'svg.fonttype': 'none'}

def batch_create_figure():
    """Plot the original four scenarios without changing model calculations."""
    with plt.rc_context(BATCH_STYLE):
        (fig, axes) = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.095, right=0.98, bottom=0.09, top=0.79, hspace=0.34, wspace=0.13)
        for (panel, ax, (label, ratio, (c0, r0, s0))) in zip('ABCD', axes.flat, BATCH_RATIOS):
            sol = simulate_batch((c0, s0, r0), t_max=30, params=DEFAULT_PARAMS)
            (C, S, R, T) = sol.y
            total = C + S + R
            for (values, (legend, color, linestyle)) in zip((C, R, S), BATCH_STRAINS):
                ax.plot(sol.t, values / total, color=color, linestyle=linestyle, label=legend, linewidth=2.2)
            ax.set_title(f'{label}  |  C:R:S = {ratio}', loc='left', fontsize=11, pad=11)
            ax.text(-0.065, 1.045, panel, transform=ax.transAxes, fontsize=12, fontweight='bold', ha='right', va='bottom')
            ax.set_xlim(0, 30)
            ax.set_ylim(-0.02, 1.02)
            ax.xaxis.set_major_locator(MultipleLocator(5))
            ax.yaxis.set_major_locator(MultipleLocator(0.2))
            ax.set_axisbelow(True)
            ax.grid(axis='y', alpha=0.28)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(length=3, width=0.8, pad=5)
        for ax in axes[1]:
            ax.set_xlabel('Time (model units)', labelpad=9)
        for ax in axes[:, 0]:
            ax.set_ylabel('Relative abundance', labelpad=10)
        fig.suptitle('Single-batch dynamics in well-mixed culture', y=0.965, fontsize=17)
        fig.text(0.5, 0.915, 'Effect of initial C:R:S composition', ha='center', color=BATCH_MUTED, fontsize=11)
        (handles, labels) = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.884), ncol=3, frameon=False, handlelength=3.2, columnspacing=2.4, fontsize=10)
    return fig

def batch_main(output_dir=None):
    """Save CSV data, a 300 dpi PNG and an editable vector SVG."""
    output_dir = Path(output_dir) if output_dir else ROOT / 'output'
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / f'{BATCH_OUTPUT_STEM}.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(('initial_ratio_C_R_S', 'time_model_units', 'C_fraction', 'R_fraction', 'S_fraction'))
        for (_, ratio, (c0, r0, s0)) in BATCH_RATIOS:
            sol = simulate_batch((c0, s0, r0), t_max=30, params=DEFAULT_PARAMS)
            (C, S, R, T) = sol.y
            total = C + S + R
            writer.writerows(zip([ratio] * len(sol.t), sol.t, C / total, R / total, S / total))
    fig = batch_create_figure()
    try:
        with plt.rc_context(BATCH_STYLE):
            for extension in ('png', 'svg'):
                output = output_dir / f'{BATCH_OUTPUT_STEM}.{extension}'
                fig.savefig(output, dpi=300, facecolor=BATCH_BACKGROUND, transparent=False)
                print(f'Saved {output}')
    finally:
        plt.close(fig)
import numpy as np
from matplotlib.colors import ListedColormap
from model import DEFAULT_PARAMS, simulate_serial_dilution
from model import SQRT3_2, bary_to_xy, simplex_grid
TERNARY_N_CYCLES = 60
TERNARY_GRID_DIVISIONS = 28
TERNARY_OUTPUT_STEM = 'figure3-2_deterministic_ternary'
TERNARY_BACKGROUND = '#000000'
TERNARY_FOREGROUND = '#E8ECF1'
TERNARY_MUTED = '#ABB4C0'
TERNARY_COLORS = {'C': '#F3AC69', 'R': '#76B7F2', 'S': '#79D3AF', 'none': '#C4A1ED'}
TERNARY_STYLE = {'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 11, 'text.color': TERNARY_FOREGROUND, 'axes.facecolor': TERNARY_BACKGROUND, 'figure.facecolor': TERNARY_BACKGROUND, 'savefig.facecolor': TERNARY_BACKGROUND, 'axes.labelcolor': TERNARY_FOREGROUND, 'axes.edgecolor': '#67717E', 'xtick.color': TERNARY_MUTED, 'ytick.color': TERNARY_MUTED, 'svg.fonttype': 'none'}

def ternary_compute_results():
    """Retain the original simulation and persistence/classification logic."""
    (xs, ys, persistence, first_lost) = ([], [], [], [])
    grid = simplex_grid(TERNARY_GRID_DIVISIONS)
    for (index, (rC, rR, rS)) in enumerate(grid, start=1):
        if rC + rR + rS == 0:
            continue
        (_, ext) = simulate_serial_dilution((rC, rR, rS), n_cycles=TERNARY_N_CYCLES, params=DEFAULT_PARAMS)
        cycles = {k: v if v is not None else TERNARY_N_CYCLES + 1 for (k, v) in ext.items()}
        first = min(cycles, key=cycles.get)
        pers = cycles[first] - 1
        (x, y) = bary_to_xy(rC, rR, rS)
        xs.append(x)
        ys.append(y)
        persistence.append(min(pers, TERNARY_N_CYCLES))
        first_lost.append(first if cycles[first] <= TERNARY_N_CYCLES else 'none')
        if index % 50 == 0 or index == len(grid):
            print(f'Simulated {index}/{len(grid)} initial compositions', flush=True)
    return tuple((np.array(values) for values in (xs, ys, persistence, first_lost)))

def ternary_draw_frame(ax):
    """Local English frame; shared ternary utilities and other figures stay intact."""
    ax.plot([0, 1, 0.5, 0], [0, 0, SQRT3_2, 0], color='#ABB4C0', linewidth=0.85, zorder=3)
    ax.text(0.5, SQRT3_2 + 0.045, 'C (colicin producer)', ha='center', va='bottom')
    ax.text(0, -0.055, 'R (resistant)', ha='center', va='top')
    ax.text(1, -0.055, 'S (sensitive)', ha='center', va='top')
    ax.set_xlim(-0.14, 1.14)
    ax.set_ylim(-0.105, SQRT3_2 + 0.105)
    ax.set_aspect('equal')
    ax.axis('off')

def ternary_create_figure(results):
    (xs, ys, persistence, first_lost) = results
    with plt.rc_context(TERNARY_STYLE):
        fig = plt.figure(figsize=(12.8, 7.4))
        axes = [fig.add_axes((0.035, 0.22, 0.445, 0.575)), fig.add_axes((0.52, 0.22, 0.445, 0.575))]
        for ax in axes:
            ternary_draw_frame(ax)
        cmap = ListedColormap(plt.get_cmap('viridis')(np.linspace(0.18, 1, 256)))
        dots = axes[0].scatter(xs, ys, c=persistence, cmap=cmap, s=32, vmin=0, vmax=TERNARY_N_CYCLES, linewidths=0, zorder=2)
        colorbar_ax = fig.add_axes((0.105, 0.145, 0.305, 0.024))
        colorbar = fig.colorbar(dots, cax=colorbar_ax, orientation='horizontal', ticks=np.arange(0, TERNARY_N_CYCLES + 1, 10))
        colorbar.set_label('Three-strain coexistence (passages)', labelpad=8, fontsize=10)
        colorbar.outline.set_visible(False)
        colorbar.ax.tick_params(labelsize=9, length=3)
        labels = {'C': 'C below threshold first', 'R': 'R below threshold first', 'S': 'S below threshold first', 'none': f'No loss within {TERNARY_N_CYCLES} passages'}
        for strain in ('C', 'R', 'S', 'none'):
            mask = first_lost == strain
            if mask.any():
                axes[1].scatter(xs[mask], ys[mask], color=TERNARY_COLORS[strain], s=32, label=labels[strain], linewidths=0, zorder=2)
        (handles, legend_labels) = axes[1].get_legend_handles_labels()
        fig.legend(handles, legend_labels, loc='upper center', bbox_to_anchor=(0.7425, 0.18), ncol=3 if len(handles) <= 3 else 2, frameon=False, columnspacing=1.4, handletextpad=0.4, fontsize=10)
        fig.suptitle('Initial composition shapes coexistence', y=0.965, fontsize=17)
        fig.text(0.5, 0.91, 'Deterministic serial dilution in well-mixed culture | 60 passages', ha='center', color=TERNARY_MUTED, fontsize=11)
        fig.text(0.065, 0.827, 'A', fontsize=13, fontweight='bold')
        fig.text(0.095, 0.827, 'Coexistence duration', fontsize=12)
        fig.text(0.55, 0.827, 'B', fontsize=13, fontweight='bold')
        fig.text(0.58, 0.827, 'First strain below threshold', fontsize=12)
    return fig

def ternary_main(output_dir=None):
    output_dir = Path(output_dir) if output_dir else ROOT / 'output'
    output_dir.mkdir(parents=True, exist_ok=True)
    results = ternary_compute_results()
    with (output_dir / f'{TERNARY_OUTPUT_STEM}.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(('r_C', 'r_R', 'r_S', 'x', 'y', 'coexistence_passages', 'first_lost'))
        for (ratio, x, y, persistence, first) in zip(simplex_grid(TERNARY_GRID_DIVISIONS), *results):
            writer.writerow((*ratio, x, y, int(persistence), first))
    fig = ternary_create_figure(results)
    try:
        with plt.rc_context(TERNARY_STYLE):
            for extension in ('png', 'svg'):
                output = output_dir / f'{TERNARY_OUTPUT_STEM}.{extension}'
                fig.savefig(output, dpi=300, facecolor=TERNARY_BACKGROUND, transparent=False)
                print(f'Saved {output}', flush=True)
    finally:
        plt.close(fig)
    persistence = results[2]
    print('median persistence:', np.median(persistence), 'max:', persistence.max())
    print('fraction with persistence >= 30 cycles:', (persistence >= 30).mean())
import matplotlib.patheffects as path_effects
from matplotlib.ticker import MaxNLocator
POTENCY_N_CYCLES = 40
POTENCY_RC_VALUES = np.linspace(0.01, 0.5, 22)
POTENCY_POTENCY_VALUES = np.linspace(0.5, 10, 18)
POTENCY_OUTPUT_STEM = 'figure3-3_biodiversity_law'
POTENCY_BACKGROUND = '#000000'
POTENCY_FOREGROUND = '#E8ECF1'
POTENCY_MUTED = '#ABB4C0'
POTENCY_STYLE = {'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 11, 'text.color': POTENCY_FOREGROUND, 'axes.facecolor': POTENCY_BACKGROUND, 'figure.facecolor': POTENCY_BACKGROUND, 'savefig.facecolor': POTENCY_BACKGROUND, 'axes.labelcolor': POTENCY_FOREGROUND, 'axes.edgecolor': '#67717E', 'axes.linewidth': 0.8, 'xtick.color': POTENCY_MUTED, 'ytick.color': POTENCY_MUTED, 'xtick.labelsize': 10, 'ytick.labelsize': 10, 'mathtext.fontset': 'dejavusans', 'svg.fonttype': 'none'}

def potency_compute_results():
    """Run the unchanged original parameter scan and persistence calculation."""
    Z = np.zeros((len(POTENCY_POTENCY_VALUES), len(POTENCY_RC_VALUES)))
    for (i, pot) in enumerate(POTENCY_POTENCY_VALUES):
        params = dict(DEFAULT_PARAMS)
        params['kappa'] = pot
        params['kill'] = pot
        for (j, rC) in enumerate(POTENCY_RC_VALUES):
            rR = rS = (1 - rC) / 2
            (_, ext) = simulate_serial_dilution((rC, rR, rS), n_cycles=POTENCY_N_CYCLES, params=params)
            cycles = {k: v if v is not None else POTENCY_N_CYCLES + 1 for (k, v) in ext.items()}
            pers = min(cycles.values()) - 1
            Z[i, j] = min(pers, POTENCY_N_CYCLES)
        print(f'Simulated {(i + 1) * len(POTENCY_RC_VALUES)}/{Z.size} parameter combinations', flush=True)
    return Z

def potency_create_figure(Z):
    with plt.rc_context(POTENCY_STYLE):
        (fig, ax) = plt.subplots(figsize=(9.5, 7.3))
        fig.subplots_adjust(left=0.115, right=0.795, bottom=0.175, top=0.805)
        cmap = ListedColormap(plt.get_cmap('viridis')(np.linspace(0.18, 1, 256)))
        heatmap = ax.pcolormesh(POTENCY_RC_VALUES, POTENCY_POTENCY_VALUES, Z, shading='auto', cmap=cmap, rasterized=False)
        colorbar_ax = fig.add_axes((0.83, 0.175, 0.024, 0.63))
        colorbar = fig.colorbar(heatmap, cax=colorbar_ax)
        colorbar.set_label('Three-strain coexistence (passages)', labelpad=13, fontsize=11)
        colorbar.locator = MaxNLocator(nbins=6, integer=True)
        colorbar.update_ticks()
        colorbar.outline.set_visible(False)
        colorbar.ax.tick_params(length=3, labelsize=10)
        contours = ax.contour(POTENCY_RC_VALUES, POTENCY_POTENCY_VALUES, Z, levels=[5, 15], colors=POTENCY_FOREGROUND, linewidths=1.15)
        contour_labels = ax.clabel(contours, fmt=lambda value: f'{value:g} passages', fontsize=9, inline=True, inline_spacing=7)
        for label in contour_labels:
            label.set_path_effects([path_effects.withStroke(linewidth=2.5, foreground='#18212B')])
        ax.set_xlabel('Initial producer fraction, $r_C$', labelpad=11)
        ax.set_ylabel('Joint toxin parameter, $p$ (a.u.)', labelpad=12)
        ax.set_xticks(np.arange(0, 0.51, 0.1))
        ax.set_yticks(np.arange(2, 11, 2))
        ax.tick_params(length=3, width=0.8, pad=6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.suptitle('Initial composition and toxin activity', y=0.96, fontsize=17)
        fig.text(0.5, 0.905, 'Well-mixed serial dilution | 40 passages', ha='center', fontsize=11, color=POTENCY_MUTED)
        fig.text(0.455, 0.055, '$r_R = r_S = (1-r_C)/2$    |    $\\kappa = k_{\\mathrm{kill}} = p$', ha='center', fontsize=10, color=POTENCY_MUTED)
    return fig

def potency_main(output_dir=None):
    output_dir = Path(output_dir) if output_dir else ROOT / 'output'
    output_dir.mkdir(parents=True, exist_ok=True)
    Z = potency_compute_results()
    with (output_dir / f'{POTENCY_OUTPUT_STEM}.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(('r_C', 'r_R', 'r_S', 'potency_p', 'kappa', 'kill', 'coexistence_passages'))
        for (i, pot) in enumerate(POTENCY_POTENCY_VALUES):
            for (j, rC) in enumerate(POTENCY_RC_VALUES):
                writer.writerow((rC, (1 - rC) / 2, (1 - rC) / 2, pot, pot, pot, int(Z[i, j])))
    fig = potency_create_figure(Z)
    try:
        with plt.rc_context(POTENCY_STYLE):
            for extension in ('png', 'svg'):
                output = output_dir / f'{POTENCY_OUTPUT_STEM}.{extension}'
                fig.savefig(output, dpi=300, facecolor=POTENCY_BACKGROUND, transparent=False)
                print(f'Saved {output}', flush=True)
    finally:
        plt.close(fig)
    print(f'Coexistence range: {Z.min():g}-{Z.max():g} passages across {Z.size} conditions')
import hashlib
import importlib.util
import json
import platform
import time

def robust_scenarios(base):
    """One factor at a time. R is fixed as the growth-rate reference."""
    cases = [dict(name='baseline', label='Baseline', family='baseline', params=dict(base), duration=12.0)]
    for key in ('mu_C', 'mu_S'):
        for factor in (0.9, 1.1):
            params = dict(base)
            params[key] *= factor
            cases.append(dict(name=f'{key}_{factor:.1f}', label=f'{key} x {factor:.1f}', family=key, params=params, duration=12.0))
    for duration in (9.0, 15.0):
        cases.append(dict(name=f'duration_{duration:g}', label=f'Duration = {duration:g}', family='duration', params=dict(base), duration=duration))
    return cases

def robust_persistence(fractions, threshold):
    """Rows are post-passage C,R,S fractions; first crossing is retained."""
    import numpy as np
    crossings = np.flatnonzero(np.any(fractions < threshold, axis=1))
    if not len(crossings):
        return (len(fractions), 'none')
    index = int(crossings[0])
    first = int(np.flatnonzero(fractions[index] < threshold)[0])
    return (index, ('C', 'R', 'S')[first])

def robust_compute(model, c_values, n_cycles):
    import numpy as np
    rows = []
    for case in robust_scenarios(model.DEFAULT_PARAMS):
        for c in c_values:
            (trajectory, _) = model.simulate_serial_dilution((float(c), float((1 - c) / 2), float((1 - c) / 2)), n_cycles=n_cycles, t_cycle=case['duration'], params=case['params'])
            fractions = trajectory[1:]
            if fractions.shape != (n_cycles, 3) or not np.isfinite(fractions).all() or np.any(fractions < 0) or (not np.allclose(fractions.sum(axis=1), 1)):
                raise RuntimeError(f"Invalid trajectory: {case['name']}, C={c}")
            thresholds = (1e-07, 1e-06, 1e-05) if case['name'] == 'baseline' else (1e-06,)
            for threshold in thresholds:
                (count, lost) = robust_persistence(fractions, threshold)
                name = case['name'] if threshold == 1e-06 else f'threshold_{threshold:.0e}'
                family = case['family'] if threshold == 1e-06 else 'threshold'
                rows.append(dict(scenario=name, family=family, r_C=float(c), r_R=float((1 - c) / 2), r_S=float((1 - c) / 2), mu_C=case['params']['mu_C'], mu_R=case['params']['mu_R'], mu_S=case['params']['mu_S'], t_cycle=case['duration'], dilution=100.0, threshold=threshold, horizon_passages=n_cycles, coexistence_passages=count, coexistence_time=count * case['duration'], first_lost=lost, right_censored=count == n_cycles))
        print(f"Completed {case['name']}: {len(c_values)} compositions", flush=True)
    return rows

def robust_summarize(rows):
    import numpy as np
    records = []
    for name in dict.fromkeys((row['scenario'] for row in rows)):
        group = [r for r in rows if r['scenario'] == name]
        baseline = next((r for r in group if np.isclose(r['r_C'], 1 / 3, atol=1e-12, rtol=0)))
        best = max((r['coexistence_passages'] for r in group))
        winners = [r['r_C'] for r in group if r['coexistence_passages'] == best]
        band = [r for r in group if 0.05 - 1e-12 <= r['r_C'] <= 0.2 + 1e-12]
        records.append(dict(scenario=name, equal_ratio_passages=baseline['coexistence_passages'], max_sampled_passages=best, maximizing_C_fractions=winners, band_points=len(band), band_points_better_than_equal=sum((r['coexistence_passages'] > baseline['coexistence_passages'] for r in band)), band_best_passages=max((r['coexistence_passages'] for r in band)), band_worst_passages=min((r['coexistence_passages'] for r in band))))
    return records

def robust_plot(rows, output):
    os.environ.setdefault('MPLCONFIGDIR', str(output / '.mplconfig'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    style = {'figure.facecolor': 'black', 'axes.facecolor': 'black', 'savefig.facecolor': 'black', 'text.color': '#E8ECF1', 'axes.labelcolor': '#E8ECF1', 'axes.edgecolor': '#718096', 'xtick.color': '#BCC5D0', 'ytick.color': '#BCC5D0', 'font.size': 10, 'font.family': 'DejaVu Sans', 'svg.fonttype': 'none'}
    labels = {'baseline': 'Baseline', 'mu_C_0.9': 'Producer growth x 0.9', 'mu_C_1.1': 'Producer growth x 1.1', 'mu_S_0.9': 'Sensitive growth x 0.9', 'mu_S_1.1': 'Sensitive growth x 1.1', 'duration_9': '9 time units / passage', 'duration_15': '15 time units / passage', 'threshold_1e-07': 'Threshold = 1e-7', 'threshold_1e-05': 'Threshold = 1e-5'}
    with plt.rc_context(style):
        (fig, axes) = plt.subplots(2, 2, figsize=(12, 8.5), sharex=True, sharey=True)
        for (ax, family, title) in zip(axes.flat, ('mu_C', 'mu_S', 'duration', 'threshold'), ('A  Producer growth', 'B  Sensitive growth', 'C  Passage duration', 'D  Detection threshold')):
            names = ['baseline'] + list(dict.fromkeys((r['scenario'] for r in rows if r['family'] == family)))
            ax.axvspan(0.05, 0.2, color='#A6B5C8', alpha=0.12, zorder=0)
            ax.axvline(1 / 3, color='#A6B5C8', ls=':', lw=1)
            for (name, color) in zip(names, ('#E8ECF1', '#5BC8F2', '#F4B66A')):
                group = sorted((r for r in rows if r['scenario'] == name), key=lambda r: r['r_C'])
                ax.plot([r['r_C'] for r in group], [r['coexistence_passages'] for r in group], color=color, marker='.', markersize=3, lw=1.4, label=labels[name])
            ax.set_title(title, loc='left', pad=10)
            ax.legend(frameon=False, fontsize=9)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.spines[['top', 'right']].set_visible(False)
            ax.set_xlim(0, 0.51)
            ax.set_ylim(-0.5, max((r['coexistence_passages'] for r in rows)) + 2)
        for ax in axes[1]:
            ax.set_xlabel('Initial producer fraction, C')
        for ax in axes[:, 0]:
            ax.set_ylabel('Three-strain coexistence (passages)')
        fig.suptitle('Sensitivity of initial-composition guidance', fontsize=16, y=0.98)
        fig.text(0.5, 0.02, 'R = S initially | Shading: proposed 5-20% C band | Dotted line: C:R:S = 1:1:1', ha='center', fontsize=9, color='#BCC5D0')
        fig.tight_layout(rect=(0, 0.05, 1, 0.95))
        for extension in ('png', 'svg'):
            fig.savefig(output / f'figure3-4_robustness.{extension}', dpi=220)
        plt.close(fig)

def robust_write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

def robust_main(output_dir=None):
    output = Path(output_dir).resolve() if output_dir else ROOT / 'output'
    output.mkdir(parents=True, exist_ok=True)
    import numpy as np
    import scipy
    import matplotlib
    model_path = ROOT / 'model.py'
    model = crs_model
    start = time.perf_counter()
    c_values = np.unique(np.r_[np.linspace(0.01, 0.5, 31), [0.05, 0.1, 0.15, 0.2, 1 / 3]])
    rows = robust_compute(model, c_values, 60)
    summary = robust_summarize(rows)
    robust_write_csv(output / 'robustness_scan.csv', rows)
    with (output / 'summary.json').open('w', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=2)
    robust_plot(rows, output)
    metadata = dict(model_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(), python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__, c_values=c_values.tolist(), simulated_trajectories=7 * len(c_values), rows=len(rows), elapsed_seconds=round(time.perf_counter() - start, 2), design='OAT scenarios, not a probability distribution or global sensitivity analysis', params=model.DEFAULT_PARAMS, n_cycles=60, scenarios=[dict(name=c['name'], params=c['params'], duration=c['duration']) for c in robust_scenarios(model.DEFAULT_PARAMS)], thresholds=[1e-07, 1e-06, 1e-05], solver='Original solve_ivp defaults; max_step=0.1')
    (output / 'RUN_METADATA.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    lines = ['# 计算结果记录（自动生成）', '', '以下是情景分析，不是实验验证或全局最优搜索。', '', '| 情景 | 等比例共存轮数 | 扫描最大轮数 | 5–20% 中优于等比例的点数/该区间采样点数 | 达到最大值的 C 比例 |', '| --- | ---: | ---: | ---: | --- |']
    for s in summary:
        winners = ', '.join((f'{100 * c:.2f}%' for c in s['maximizing_C_fractions']))
        lines.append(f"| {s['scenario']} | {s['equal_ratio_passages']} | {s['max_sampled_passages']} | {s['band_points_better_than_equal']}/{s['band_points']} | {winners} |")
    lines += ['', '最大值只针对所采样的点；离散轮数可能产生多个并列点。不同传代时长的轮数并非相同培养时间。', 'CSV 同时记录完成共存轮次对应的培养时间；不包含稀释和操作用时，也不表示精确阈值跨越时刻。']
    (output / 'RESULTS_zh.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(f'Complete: {output}', flush=True)
from scipy.integrate import solve_ivp

def check_tighter_run(model, c, params, duration):
    """Same rhs, protocol, clipping and detection; tighter solver tolerances only."""
    state = np.array([c, (1 - c) / 2, (1 - c) / 2, 0.0]) * 0.01
    fractions = []
    for _ in range(60):
        sol = solve_ivp(model.rhs, (0, duration), state, args=(params,), max_step=0.1, rtol=1e-08, atol=1e-12)
        if not sol.success:
            raise RuntimeError(sol.message)
        state = np.r_[np.maximum(sol.y[:3, -1], 0) / 100, 0.0]
        total = state[:3].sum()
        if total <= 0:
            raise RuntimeError('Nonpositive total population')
        fractions.append(state[[0, 2, 1]] / total)
    return robust_persistence(np.asarray(fractions), 1e-06)

def check_main(output_dir=None):
    model = crs_model
    output = Path(output_dir) if output_dir else ROOT / 'output'
    with (output / 'robustness_scan.csv').open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    metadata = json.loads((output / 'RUN_METADATA.json').read_text(encoding='utf-8'))
    assert len(rows) == 9 * len(metadata['c_values'])
    assert len({(r['scenario'], r['r_C']) for r in rows}) == len(rows)
    for row in rows:
        assert np.isclose(sum((float(row[k]) for k in ('r_C', 'r_R', 'r_S'))), 1)
        assert np.isclose(float(row['r_R']), float(row['r_S']))
        assert 0 <= int(row['coexistence_passages']) <= 60
        assert np.isclose(float(row['coexistence_time']), int(row['coexistence_passages']) * float(row['t_cycle']))
    for c in metadata['c_values']:
        group = [r for r in rows if float(r['r_C']) == c and r['scenario'] in ('baseline', 'threshold_1e-07', 'threshold_1e-05')]
        counts = [int(r['coexistence_passages']) for r in sorted(group, key=lambda r: float(r['threshold']))]
        assert len(counts) == 3 and counts[0] >= counts[1] >= counts[2]
    assert robust_persistence(np.array([[0.3, 0.3, 0.4], [1e-08, 0.5, 0.5]]), 1e-06) == (1, 'C')
    assert robust_persistence(np.array([[0.5, 0.5, 1e-08]]), 1e-06) == (0, 'S')
    assert robust_persistence(np.array([[1e-08, 1e-08, 1.0]]), 1e-06) == (0, 'C')
    assert robust_persistence(np.array([[0.3, 0.3, 0.4]]), 1e-06) == (1, 'none')
    for threshold in (1e-07, 1e-06, 1e-05):
        (trajectory, loss) = model.simulate_serial_dilution((0.1, 0.45, 0.45), n_cycles=60, extinct_thresh=threshold)
        expected = min((v if v is not None else 61 for v in loss.values())) - 1
        (actual, _) = robust_persistence(trajectory[1:], threshold)
        assert actual == expected
        name = 'baseline' if threshold == 1e-06 else f'threshold_{threshold:.0e}'
        exported = next((r for r in rows if r['scenario'] == name and float(r['r_C']) == 0.1))
        assert actual == int(exported['coexistence_passages'])
    numerical_checks = []
    cases = {case['name']: case for case in robust_scenarios(model.DEFAULT_PARAMS)}
    selected = [('baseline', c) for c in (0.01, 0.1, 1 / 3)]
    for name in cases:
        if name != 'baseline':
            best = max((r for r in rows if r['scenario'] == name), key=lambda r: int(r['coexistence_passages']))
            selected.append((name, float(best['r_C'])))
    for (name, c) in selected:
        case = cases[name]
        (strict, strict_lost) = check_tighter_run(model, c, case['params'], case['duration'])
        original = next((r for r in rows if r['scenario'] == name and float(r['r_C']) == c))
        numerical_checks.append(dict(scenario=name, r_C=c, original_passages=int(original['coexistence_passages']), tighter_passages=strict, original_first_lost=original['first_lost'], tighter_first_lost=strict_lost, matched=strict == int(original['coexistence_passages']) and strict_lost == original['first_lost']))
    report = dict(metric_and_export_checks='PASS', original_threshold_detection='PASS', tighter_tolerances=dict(rtol=1e-08, atol=1e-12, checks=numerical_checks), scope='Three baseline points and one sampled maximizer in each of six changed dynamical scenarios; not full-grid convergence verification')
    (output / 'CHECK_RESULT.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))
    if not all((c['matched'] for c in numerical_checks)):
        print('REVIEW: tolerance-sensitive outcomes; do not interpret these as robust predictions.')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'output')
    parser.add_argument('--figure', choices=('all', '3-1', '3-2', '3-3', '3-4'), default='all')
    parser.add_argument('--check', action='store_true', help='Verify existing outputs and run nine tighter-tolerance spot checks.')
    args = parser.parse_args()
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.check:
        check_main(output)
        return
    for (number, runner) in (('3-1', batch_main), ('3-2', ternary_main), ('3-3', potency_main), ('3-4', robust_main)):
        if args.figure in ('all', number):
            print(f'Running Figure {number}', flush=True)
            runner(output)
    print(f'Done: {output}', flush=True)
if __name__ == '__main__':
    main()
