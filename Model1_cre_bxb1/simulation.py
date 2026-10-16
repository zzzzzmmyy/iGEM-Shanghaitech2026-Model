"""Reproduce the 11 numerical Wiki figures as PNG and SVG.

Default analysis uses the reference parameters. --input changes only Figures
1-11 to 1-15; the six basis plots retain their reference settings.
Figures 1-1, 1-2, 1-4, and 1-5 are maintained separately and are not generated.
"""
from pathlib import Path
import argparse
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from design_model import Model, Settings, run_analysis

ROOT = Path(__file__).resolve().parent


def plot_basis_figures(destination):
    """Regenerate the six wiki figures from the standard-library prediction model."""
    from dataclasses import replace
    from pathlib import Path

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    from design_model import Model, cre_outcomes, population_at, _site_at

    ROOT = Path(__file__).resolve().parent
    DESTINATION = ROOT / 'figures'
    BACKGROUND, TEXT, MUTED = '#000000', '#F2F4F7', '#B7BEC8'
    LIME, BLUE, ORANGE = '#B7D66D', '#6BBFFF', '#FFB86B'
    SITE_COLORS = ['#D8EBA8', LIME, '#91C957', '#63B87E', '#39A99F']


    def style_axes(ax):
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='#252A30', linewidth=.6)
        ax.tick_params(length=3, pad=6)


    def save(fig, name):
        destination = DESTINATION
        destination.mkdir(exist_ok=True)
        fig.savefig(destination / f'{name}.png', dpi=300, bbox_inches='tight', facecolor=BACKGROUND)
        fig.savefig(destination / f'{name}.svg', bbox_inches='tight', facecolor=BACKGROUND)
        plt.close(fig)


    def single(model, left_arm, right_arm, time=24):
        return _site_at(replace(model, left_arm_bp=left_arm), right_arm, time)['type_a_fraction']


    def main(destination=None):
        nonlocal DESTINATION
        DESTINATION = Path(destination) if destination is not None else ROOT / 'figures'
        model = Model()
        sites = cre_outcomes(model)
        times = [24 * i / 299 for i in range(300)]
        sweep = [300 + (6000 - 300) * i / 119 for i in range(120)]
        plt.rcParams.update({
            'font.family': 'DejaVu Sans', 'font.size': 10,
            'figure.facecolor': BACKGROUND, 'axes.facecolor': BACKGROUND,
            'savefig.facecolor': BACKGROUND, 'text.color': TEXT,
            'axes.labelcolor': TEXT, 'axes.edgecolor': MUTED,
            'xtick.color': MUTED, 'ytick.color': MUTED,
            'axes.linewidth': .8, 'axes.spines.top': False,
            'axes.spines.right': False, 'svg.fonttype': 'none',
            'xtick.direction': 'out', 'ytick.direction': 'out',
        })

        fig, ax = plt.subplots(figsize=(7.6, 4.8), layout='constrained')
        bars = ax.bar([s['site'] for s in sites],
                      [s['selection_probability'] * 100 for s in sites],
                      width=.58, color=LIME, zorder=3)
        ax.bar_label(bars, fmt='%.1f%%', padding=7, color=TEXT, fontsize=11)
        ax.set_xticks([s['site'] for s in sites],
                      [f'lox71_{s["site"]}\n{s["synapsis_distance_bp"]:,.0f} bp' for s in sites])
        ax.set_ylim(0, 40)
        ax.set_yticks(range(0, 41, 10))
        ax.set_xlim(.4, model.n_sites + .6)
        ax.set_xlabel('Cre site and pre-recombination sequence distance', labelpad=13)
        ax.set_ylabel('Selection probability (%)', labelpad=12)
        ax.set_title('Distance-weighted Cre site selection', loc='center', fontsize=12, pad=20)
        style_axes(ax)
        ax.tick_params(axis='x', length=0, pad=10)
        save(fig, 'figure1-3_Layer1_Cre_selection')

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout='constrained')
        fig.suptitle('Arm-length-dependent differentiation', fontsize=13)
        for arm, color in zip([1000, 2000, 4000], [LIME, BLUE, ORANGE]):
            axes[0].plot(times, [None if t == 0 else single(model, arm, 773, t) * 100 for t in times],
                         color=color, lw=2, label=f'Left arm: {arm:,} bp')
        axes[0].set_xlim(0, 24)
        axes[0].set_xlabel('Time (h)', labelpad=10)
        axes[0].set_title('A  Differentiation time courses', loc='center', fontsize=11, pad=15)
        axes[1].plot([h / 1000 for h in sweep], [single(model, h, 773) * 100 for h in sweep],
                     color=LIME, lw=2, label='24 h')
        axes[1].plot([h / 1000 for h in sweep],
                     [_site_at(replace(model, left_arm_bp=h), 773, 24)['limiting_type_a_fraction'] * 100
                      for h in sweep], color=BLUE, lw=2, ls='--', label='Limit as time approaches infinity')
        axes[1].set_xlim(.3, 6)
        axes[1].set_xlabel('Left arm length (kb)', labelpad=10)
        axes[1].set_title('B  Finite-time and limiting fractions', loc='center', fontsize=11, pad=15)
        for ax in axes:
            ax.set_ylim(0, 100)
            ax.set_ylabel('Venus among differentiated cells (%)', labelpad=10)
            style_axes(ax)
            ax.axhline(50, color=MUTED, ls=':', lw=.8, alpha=.6)
            ax.legend(frameon=False, fontsize=9, loc='upper right')
        save(fig, 'figure1-6_ODE_validation')

        ratios = [1, 2.5, 4, 7, 10]
        for shorter, name in [(True, 'figure1-7_1to10'), (False, 'figure1-8_10to1')]:
            arms = [773 / r if shorter else 773 * r for r in ratios]
            finite = [single(model, arm, 773) * 100 for arm in arms]
            infinite = [_site_at(replace(model, left_arm_bp=arm), 773, 24)
                        ['limiting_type_a_fraction'] * 100 for arm in arms]
            fig, ax = plt.subplots(figsize=(7.6, 4.8), layout='constrained')
            ax.plot(range(5), finite, 'o-', color=LIME, lw=2, ms=6, label='24 h')
            ax.plot(range(5), infinite, 's--', color=BLUE, lw=2, ms=5,
                    label='Limit as time approaches infinity')
            for j, value in enumerate(finite):
                ax.annotate(f'{value:.1f}%', (j, value), xytext=(0, 10),
                            textcoords='offset points', ha='center')
            ax.set_xticks(range(5), [f'1:{r:g}' if shorter else f'{r:g}:1' for r in ratios])
            ax.set_xlim(-.25, 4.25)
            ax.set_ylim(0, 100)
            ax.set_xlabel('Left arm : right arm', labelpad=12)
            ax.set_ylabel('Venus among differentiated cells (%)', labelpad=12)
            end = '1:10' if shorter else '10:1'
            ax.set_title(f'Arm ratio and Venus fraction (1:1 to {end})', loc='center', fontsize=12, pad=20)
            style_axes(ax)
            ax.axhline(50, color=MUTED, ls=':', lw=.8, alpha=.6)
            ax.legend(frameon=False, fontsize=9, loc='upper left')
            save(fig, name)

        curve = [population_at(model, t) for t in times]
        fig, ax = plt.subplots(figsize=(11, 5), layout='constrained')
        for site, color in zip(sites, SITE_COLORS):
            values = [_site_at(model, site['right_arm_bp'], t)['type_a_fraction'] for t in times]
            ax.plot(times, [None if value is None else value * 100 for value in values],
                    color=color, lw=1.4,
                    label=f'lox71_{site["site"]}: {site["right_arm_bp"]:,.0f} bp, P={site["selection_probability"]*100:.1f}%')
        for key, color, line, label in [
            ('weighted_site_fraction', TEXT, '--', 'Expected mean across picked wells'),
        ]:
            ax.plot(times, [None if r[key] is None else r[key] * 100 for r in curve],
                    color=color, lw=2.5, ls=line, zorder=5, label=label)
        ax.set_xlim(0, 24)
        ax.set_ylim(0, 100)
        ax.set_xlabel('Time (h)', labelpad=12)
        ax.set_ylabel('Venus among differentiated cells (%)', labelpad=12)
        ax.set_title('Within-well differentiation by Cre outcome', loc='center', fontsize=12, pad=20)
        style_axes(ax)
        ax.axhline(50, color=MUTED, ls=':', lw=.8, alpha=.6)
        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=8.5, frameon=False)
        save(fig, 'figure1-9_Layer3_timecourse')

        n_values, spacers = list(range(1, 11)), list(range(200, 2001, 200))
        finite, infinite = [], []
        for spacer in spacers:
            predictions = [population_at(replace(model, n_sites=n, spacer_bp=spacer), 24) for n in n_values]
            finite.append([r['weighted_site_fraction'] * 100 for r in predictions])
            infinite.append([r['limiting_type_a_fraction'] * 100 for r in predictions])
        cmap = LinearSegmentedColormap.from_list(
            'wiki_green', ['#101C23', '#17433E', '#337B62', '#78AE68', '#D8EBA8'])
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), layout='constrained')
        fig.suptitle('Expected mean across randomly picked wells', fontsize=13)
        for ax, values, title in zip(axes, [finite, infinite],
                                     ['A  Mean well fraction at 24 h', 'B  Limiting mean well fraction']):
            im = ax.pcolormesh([.5 + i for i in range(11)], list(range(100, 2200, 200)),
                               values, cmap=cmap, vmin=0, vmax=100, shading='flat')
            lo, hi = min(map(min, values)), max(map(max, values))
            levels = [v for v in [10, 20, 30, 50, 70, 90] if lo < v < hi]
            if levels:
                contours = ax.contour(n_values, spacers, values, levels=levels,
                                      colors=TEXT, linewidths=.8, alpha=.9)
                ax.clabel(contours, fmt='%d%%', fontsize=8)
            ax.set_xticks(n_values)
            ax.set_yticks(spacers)
            ax.set_xlabel('Number of candidate lox71 sites', labelpad=10)
            ax.set_ylabel('NS spacer length (bp)', labelpad=10)
            ax.set_title(title, loc='center', fontsize=11, pad=15)
            ax.tick_params(length=3, pad=6)
        bar = fig.colorbar(im, ax=axes, shrink=.85, pad=.03, ticks=range(0, 101, 20))
        bar.set_label('Expected mean A fraction across wells (%)', labelpad=12)
        bar.ax.tick_params(colors=MUTED)
        bar.outline.set_edgecolor(MUTED)
        bar.outline.set_linewidth(.6)
        save(fig, 'figure1-10_Layer3_heatmap')
        print(f'Saved six basis PNG/SVG pairs to {DESTINATION.resolve()}')
    main(destination)

def plot_library_figures(source, destination=None):
    """Draw independent-well library figures from exported analysis tables."""
    import argparse
    import csv
    import json
    import math
    from pathlib import Path

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    from matplotlib.ticker import PercentFormatter

    ROOT=Path(__file__).resolve().parents[1]
    BG,FG,MUTED='#000000','#F2F4F7','#B7BEC8'
    GREEN,BLUE,ORANGE,PINK='#B7D66D','#6BBFFF','#FFB86B','#DF9BC1'
    SITE_COLORS=['#D8EBA8',GREEN,'#91C957','#63B87E','#39A99F']
    CMAP=LinearSegmentedColormap.from_list('wiki_green',['#14242B','#2B6658','#6DAB74','#D8EBA8'])


    def read_table(path):
        with Path(path).open(encoding='utf-8-sig',newline='') as stream:
            rows=list(csv.DictReader(stream))
        for row in rows:
            for key,value in row.items():
                if value=='':
                    row[key]=None
                elif value in ('True','False'):
                    row[key]=value=='True'
                else:
                    try:
                        row[key]=float(value)
                    except ValueError:
                        pass
        return rows


    def style():
        plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                             'figure.facecolor':BG,'axes.facecolor':BG,'savefig.facecolor':BG,
                             'text.color':FG,'axes.labelcolor':FG,'axes.edgecolor':MUTED,
                             'xtick.color':MUTED,'ytick.color':MUTED,'axes.linewidth':.8,
                             'axes.spines.top':False,'axes.spines.right':False,
                             'svg.fonttype':'none','xtick.direction':'out','ytick.direction':'out'})


    def grid(ax):
        ax.set_axisbelow(True)
        ax.grid(axis='y',color='#252A30',lw=.6)
        ax.tick_params(length=3,pad=5)


    def save(fig,destination,name):
        destination.mkdir(parents=True,exist_ok=True)
        fig.savefig(destination/f'{name}.png',dpi=300,bbox_inches='tight',facecolor=BG)
        fig.savefig(destination/f'{name}.svg',bbox_inches='tight',facecolor=BG)
        plt.close(fig)


    def generate_figures(source,destination=None):
        source=Path(source)
        destination=Path(destination) if destination is not None else source/'figures'
        summary=json.loads((source/'summary.json').read_text(encoding='utf-8'))
        settings=summary['settings']
        metrics=summary['metrics']
        sites=read_table(source/'site_outcomes.csv')
        wells=read_table(source/'example_wells.csv')
        curve=read_table(source/'coverage_by_picks.csv')
        targets=read_table(source/'target_picking.csv')
        designs=read_table(source/'design_scan.csv')
        n, time=settings['colonies'],settings['time_h']
        confidence=settings['confidence']
        style()

        # Expected counts, exact marginal count intervals, and within-well fractions.
        fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
        fig.suptitle(f'Expected representation after picking {n} colonies',fontsize=14)
        x=list(range(len(sites)))
        expected=[r['expected_wells'] for r in sites]
        # An equal-tailed discrete interval can exclude the mean near probability extremes.
        axes[0].bar(x,expected,color=GREEN,width=.6,zorder=3)
        for j,row in enumerate(sites):
            axes[0].vlines(j,row['count_interval_low'],row['count_interval_high'],color=FG,lw=1)
            axes[0].hlines([row['count_interval_low'],row['count_interval_high']],j-.08,j+.08,color=FG,lw=1)
            axes[0].annotate(f'{row["expected_wells"]:.1f}',(j,max(row['count_interval_high'],row['expected_wells'])),
                             xytext=(0,6),textcoords='offset points',ha='center',fontsize=9)
        axes[0].set_ylabel('Number of wells')
        axes[0].set_title(f'A  Expected counts and {confidence:.0%} marginal intervals',fontsize=11,pad=12)
        upper=max([r['count_interval_high'] for r in sites]+expected)
        axes[0].set_ylim(0,max(1,upper)*1.2)
        axes[1].bar(x,[r['type_a_fraction']*100 for r in sites],color=BLUE,width=.6)
        for j,row in enumerate(sites):
            axes[1].annotate(f'{row["type_a_fraction"]:.1%}',(j,row['type_a_fraction']*100),
                             xytext=(0,6),textcoords='offset points',ha='center',fontsize=9)
        axes[1].set_ylim(0,100)
        axes[1].set_ylabel('A fraction among differentiated cells (%)')
        axes[1].set_title(f'B  Expected composition within a well at {time:g} h',fontsize=11,pad=12)
        for ax in axes:
            ax.set_xticks(x,[f'Site {int(r["site"])}' for r in sites])
            ax.set_xlabel('Clonal Cre structure',labelpad=10)
            grid(ax)
        save(fig,destination,'figure1-11_expected_well_counts')

        # The plate is a seeded draw, not a predicted spatial pattern.
        fig,(ax,ax2)=plt.subplots(1,2,figsize=(13,5.5),layout='constrained',
                                gridspec_kw={'width_ratios':[1.75,1]})
        fig.suptitle(f'One simulated {n}-well library (seed={settings["seed"]})',fontsize=14)
        first_plate=[row for row in wells if row['plate']==1]
        matrix=[[math.nan for _ in range(12)] for _ in range(8)]
        for row in first_plate:
            position=int(row['colony'])-1
            matrix[position//12][position%12]=row['type_a_fraction']
        color=CMAP.copy()
        color.set_bad('#11151A')
        im=ax.imshow(matrix,cmap=color,vmin=0,vmax=1,aspect='equal')
        for row in first_plate:
            position=int(row['colony'])-1
            fraction=row['type_a_fraction']
            ax.text(position%12,position//12,f'S{int(row["site"])}\n{fraction:.1%}',ha='center',va='center',
                    fontsize=7.5,color=BG if fraction>.6 else FG)
        ax.set_xticks(range(12),range(1,13))
        ax.set_yticks(range(8),list('ABCDEFGH'))
        ax.set_title('A  First 96-well plate; row order is arbitrary',fontsize=11,pad=12)
        bar=fig.colorbar(im,ax=ax,pad=.02,shrink=.75,format=PercentFormatter(1))
        bar.set_label('Within-well A fraction',fontsize=9)
        x=list(range(len(sites)))
        ax2.bar([v-.18 for v in x],[r['expected_wells'] for r in sites],width=.35,color=GREEN,label='Expected')
        ax2.bar([v+.18 for v in x],[r['example_wells'] for r in sites],width=.35,color=BLUE,label='This simulation')
        ax2.set_xticks(x,[f'S{int(r["site"])}' for r in sites])
        ax2.set_xlabel('Clonal structure')
        ax2.set_ylabel(f'Counts across all {n} wells')
        ax2.set_title('B  Random counts vary around expectations',fontsize=11,pad=12)
        ax2.legend(frameon=False,fontsize=9)
        grid(ax2)
        save(fig,destination,'figure1-12_plate_example')

        # Compare coverage with diversity; a one-type library is not a broad library.
        fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
        fig.suptitle('How many colonies are needed to cover the library?',fontsize=14)
        picks=[int(row['picks']) for row in curve]
        a=[row['probability_all_structures'] for row in curve]
        b=[row['probability_all_ratio_types'] for row in curve]
        same=all(abs(u-v)<1e-10 for u,v in zip(a,b))
        axes[0].plot(picks,a,color=GREEN,lw=2,label='All structures / usable ratio types' if same else 'All structures')
        if not same:
            axes[0].plot(picks,b,color=BLUE,lw=2,ls='--',label='All usable ratio types')
        axes[0].plot(picks,[r['probability_all_library_targets'] for r in curve],
                     color=ORANGE,lw=2,ls='--',label='All requested target intervals')
        axes[0].axhline(confidence,color=MUTED,ls=':',lw=1)
        axes[0].set_ylim(0,1.04)
        axes[0].yaxis.set_major_formatter(PercentFormatter(1))
        axes[0].set_ylabel('Probability of covering the required set')
        axes[0].set_title('A  Joint coverage probability',fontsize=11,pad=12)
        axes[0].legend(frameon=False,fontsize=8,loc='lower right')
        required=metrics['required_all_ratio_types']['picks']
        note=f'{confidence:.0%} coverage of usable ratio types: {required} picks' if required is not None else 'Coverage threshold not reached within search cap'
        axes[0].text(.98,.43,note,transform=axes[0].transAxes,ha='right',fontsize=8.5,color=MUTED)
        axes[1].plot(picks,[r['expected_distinct_structures'] for r in curve],color=GREEN,lw=2,label='Distinct structures')
        axes[1].plot(picks,[r['expected_distinct_ratio_types'] for r in curve],color=BLUE,lw=2,ls='--',label='Distinct usable ratio bins')
        axes[1].set_ylabel('Expected number of represented types')
        axes[1].set_ylim(0,max(1,metrics['number_structures'])*1.1)
        axes[1].set_title('B  Expected diversity',fontsize=11,pad=12)
        axes[1].legend(frameon=False,fontsize=9,loc='lower right')
        for ax in axes:
            ax.axvline(n,color=MUTED,ls=':',lw=1)
            ax.set_xlabel('Number of colonies picked',labelpad=10)
            ax.set_xlim(0,max(picks))
            grid(ax)
        save(fig,destination,'figure1-13_coverage_probability')

        # Target bands are stated explicitly; an unreachable target is not extrapolated.
        try:
            from design_model import hit_probability
        except ImportError:
            from design_model import hit_probability
        fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
        fig.suptitle(f'Finding a target-ratio well (tolerance = {settings["target_tolerance"]*100:g} percentage points)',fontsize=13)
        palette=[GREEN,BLUE,ORANGE,PINK,'#A998DC']
        for j,row in enumerate(targets):
            axes[0].plot(picks,[hit_probability(row['single_pick_probability'],v) for v in picks],
                         color=palette[j%len(palette)],lw=2,
                         label=f'Target {row["target"]:.0%}; one-pick chance {row["single_pick_probability"]:.1%}')
        axes[0].set_ylim(0,1.04)
        axes[0].yaxis.set_major_formatter(PercentFormatter(1))
        axes[0].set_xlabel('Number of colonies picked',labelpad=10)
        axes[0].set_ylabel('Probability of at least one matching well')
        axes[0].set_title('A  At least one hit',fontsize=11,pad=12)
        axes[0].axhline(confidence,color=MUTED,ls=':',lw=1)
        axes[0].axvline(n,color=MUTED,ls=':',lw=1)
        axes[0].legend(frameon=False,fontsize=8,loc='lower right')
        finite=[row['required_picks'] for row in targets if row['required_picks'] is not None]
        scale=max(finite,default=1)
        for j,row in enumerate(targets):
            value=row['required_picks']
            axes[1].bar(j,value or 0,color=palette[j%len(palette)],width=.6)
            label=str(int(value)) if value is not None else row['requirement_status'].replace('_','\n')
            axes[1].text(j,(value or 0)+scale*.04,label,ha='center',va='bottom',fontsize=9,color=FG)
        axes[1].set_xticks(range(len(targets)),[f'{row["target"]:.0%}' for row in targets])
        axes[1].set_ylim(0,scale*1.4)
        axes[1].set_xlabel('Target A fraction within a well',labelpad=10)
        axes[1].set_ylabel('Minimum colonies to pick')
        axes[1].set_title(f'B  Picking requirement at {confidence:.0%} success',fontsize=11,pad=12)
        for ax in axes:
            grid(ax)
        save(fig,destination,'figure1-14_target_picking')

        # Library objectives are explicit, rather than optimizing only a plate mean.
        counts=sorted({int(row['n_sites']) for row in designs})
        spacers=sorted({row['spacer_bp'] for row in designs})
        lookup={(int(row['n_sites']),row['spacer_bp']):row for row in designs}
        fig,axes=plt.subplots(1,3,figsize=(15,5.3),layout='constrained')
        goals=', '.join(f'{v:.0%}' for v in settings['library_targets'])
        fig.suptitle(f'Library design: {n} picks, {time:g} h; requested ratios {goals}',fontsize=13)
        keys=['expected_distinct_ratio_types','required_picks_library','single_pick_target_probability']
        titles=['A  Expected usable ratio types','B  Required picks (lower is better)',f'C  One-pick chance near {settings["target"]:.0%}']
        labels=['Expected type count',f'Picks at {confidence:.0%} success','Matching probability (%)']
        for ax,key,title,cbar_label in zip(axes,keys,titles,labels):
            values=[[lookup[(n_sites,spacer)][key] for n_sites in counts] for spacer in spacers]
            if key=='single_pick_target_probability':
                values=[[value*100 for value in row] for row in values]
            values=[[math.nan if value is None else value for value in row] for row in values]
            finite=[v for row in values for v in row if math.isfinite(v)]
            color=CMAP.copy()
            color.set_bad('#12171D')
            im=ax.imshow(values,origin='lower',aspect='auto',cmap=color,
                          vmin=0 if key!='required_picks_library' else 1,
                          vmax=max(finite,default=1) or 1)
            for i,row in enumerate(values):
                for j,value in enumerate(row):
                    if not math.isfinite(value):
                        ax.text(j,i,'×',ha='center',va='center',color=MUTED,fontsize=9)
            best=summary['best_design']
            if best and int(best['n_sites']) in counts and best['spacer_bp'] in spacers:
                ax.plot(counts.index(int(best['n_sites'])),spacers.index(best['spacer_bp']),
                         marker='o',markersize=8,markerfacecolor='none',markeredgecolor=FG,markeredgewidth=1.2)
            ax.set_xticks(range(len(counts)),counts)
            positions=list(range(len(spacers)))
            if len(positions)>12:
                positions=positions[::max(1,len(positions)//10)]
            ax.set_yticks(positions,[f'{spacers[i]:g}' for i in positions])
            ax.set_xlabel('Candidate lox71 sites',labelpad=10)
            ax.set_ylabel('NS spacer length (bp)',labelpad=10)
            ax.set_title(title,fontsize=10.5,pad=12)
            bar=fig.colorbar(im,ax=ax,pad=.025,shrink=.85)
            bar.set_label(cbar_label,fontsize=8.5)
            bar.ax.tick_params(labelsize=8)
        fig.supxlabel('Circles: first-ranked design for the selected objective.  ×: unreachable or above the pick-count search cap.',
                       fontsize=9,color=MUTED)
        save(fig,destination,'figure1-15_library_design_space')
        print(f'Saved five independent-well figure pairs to {destination.resolve()}')
    generate_figures(source, destination)


FIGURE_STEMS = [
    'figure1-3_Layer1_Cre_selection',

    'figure1-6_ODE_validation', 'figure1-7_1to10', 'figure1-8_10to1',
    'figure1-9_Layer3_timecourse', 'figure1-10_Layer3_heatmap',
    'figure1-11_expected_well_counts', 'figure1-12_plate_example',
    'figure1-13_coverage_probability', 'figure1-14_target_picking',
    'figure1-15_library_design_space',
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('figures'),help='PNG/SVG output directory')
    parser.add_argument('--input',type=Path,help='Existing design_model.py output for Figures 1-11 to 1-15')
    args = parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    plot_basis_figures(args.out)
    if args.input:
        plot_library_figures(args.input,args.out)
    else:
        with tempfile.TemporaryDirectory(dir=args.out.resolve(),prefix='cre_bxb1_figures_') as temporary:
            run_analysis(Model(),Settings(),Path(temporary))
            plot_library_figures(Path(temporary),args.out)
    for stem in FIGURE_STEMS:
        for extension in ('png','svg'):
            if not (args.out/f'{stem}.{extension}').is_file():
                raise RuntimeError(f'Missing figure: {stem}.{extension}')
    print(f'Saved all 11 numerical figures (22 PNG/SVG files) to {args.out.resolve()}')


if __name__ == '__main__':
    main()
