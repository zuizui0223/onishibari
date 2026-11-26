#!/usr/bin/env python3
"""Render the v51 comparison from aggregate results, without individual records."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


def forest(ax, rows, null, title):
    ax.axvline(null, color='#9aa4af', ls='--', lw=1.2, zorder=0)
    for y, (label, point, bounds, color) in enumerate(reversed(rows)):
        lo, hi = bounds
        ax.errorbar(point, y, xerr=[[point-lo], [hi-point]], fmt='o',
                    color=color, ms=7.5, capsize=4, lw=2.1)
        ax.text(1.04, y, f'{point:.2f}  [{lo:.2f}, {hi:.2f}]',
                transform=ax.get_yaxis_transform(), va='center', fontsize=10)
    ax.set_yticks(range(len(rows)), [r[0] for r in reversed(rows)])
    ax.set_ylim(-.55, len(rows)-.45)
    ax.set_title(title, loc='left', pad=15, fontsize=12, fontweight='bold')
    ax.grid(axis='x', color='#e9edf1', lw=.8)
    ax.set_axisbelow(True)
    for side in ['top','right','left']:
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#b8c1ca')
    ax.tick_params(axis='y', length=0, pad=10)


def main(mercurialis, pulsatilla, output):
    m = json.loads(mercurialis.read_text())['within_flowering_mothers']
    p = json.loads(pulsatilla.read_text())['typed_arm_contrast_jecol_60']
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10,
                         'svg.fonttype':'none', 'svg.hashsalt':'onishibari-v51'})
    fig, axes = plt.subplots(2, 1, figsize=(12, 9.4),
                             gridspec_kw={'height_ratios':[4,3]})
    fig.subplots_adjust(left=.32, right=.75, top=.83, bottom=.27, hspace=.98)
    fig.suptitle('Male allocation, selfing, and the choice of endpoint',
                 x=.035, ha='left', y=.975, fontsize=17, fontweight='bold')
    fig.text(.035, .935, 'Secondary analyses of two published plant systems | aggregate estimates and 95% intervals',
             color='#526171', fontsize=11)

    rows = [('Pooled association\n164 mothers', m['OR_per_doubling'],
             m['ci95_within_array_mother_bootstrap'], '#006e79')]
    for i, record in enumerate(m['by_population_varying_male_slope_shared_adjusters'].values(), 1):
        rows.append((f'Array {i}  ({record["mothers"]} mothers)', record['OR_per_doubling'],
                     record['ci95_mother_HC3'], '#6c8994'))
    forest(axes[0], rows, 1, 'A. Wind pollination: Mercurialis annua')
    axes[0].set_xscale('log')
    axes[0].set_xlim(.93, 1.9)
    axes[0].set_xticks([1, 1.25, 1.5, 1.75])
    axes[0].xaxis.set_major_formatter(ScalarFormatter())
    axes[0].minorticks_off()
    axes[0].set_xlabel('Selfing odds ratio per doubling of male flower count', labelpad=9)
    c = m['standardized_100_200_male_flowers']
    lo, hi = c['ci95_difference_within_array_mother_bootstrap']
    axes[0].text(-.66, -.51,
        f"Model-standardized selfing: {100*c['probabilities'][0]:.1f}% at 100 male flowers; "
        f"{100*c['probabilities'][1]:.1f}% at 200\n"
        f"Difference: +{100*c['difference']:.1f} percentage points [{100*lo:.1f}, {100*hi:.1f}]",
        transform=axes[0].transAxes, fontsize=10.5, color='#00636f', linespacing=1.5)

    labels = [('Mature seed harvest','Mature_seed_number','#697884'),
              ('Estimated outcrossed harvest','outcrossed_harvest','#007f80'),
              ('Weighted maternal output\nSelfed seeds weighted 0.05','discounted_harvest','#7960a0')]
    rows = []
    for label, key, color in labels:
        r = p['all_minus_intact'][key]
        rows.append((label, r['estimate'], r['flower_percentile_95'], color))
    forest(axes[1], rows, 0, 'B. Insect pollination: Pulsatilla alpina')
    axes[1].set_xlim(-40, 80)
    axes[1].set_xticks([-40, 0, 40, 80])
    axes[1].set_xlabel('All stamens removed minus intact (seed units per retained flower)', labelpad=9)
    axes[1].text(-.66, -.70,
        'Retained, seed-typed families: 6 with all stamens removed versus 43 intact.\n'
        'The comparison is conditional on retention and typing; it is not an intention-to-treat effect.',
        transform=axes[1].transAxes, fontsize=10.5, color='#526171', linespacing=1.5)

    fig.text(.035, .058,
        'A: female flower count, above-ground biomass and array adjusted. Pooled interval: 1,999 whole-mother resamples within arrays;\n'
        'array-specific intervals: mother-level HC3, with shared female-count and biomass coefficients. No flower-number intervention was fitted.\n'
        'B: 10,000 whole-flower resamples within arms; estimated selfing fractions and delta = 0.95 held fixed. Assignment and finite-seed\n'
        'uncertainty are not propagated. Neither panel measures offspring survival or individual floral construction biomass.',
        fontsize=9, color='#526171', linespacing=1.65)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix('.svg'), metadata={'Date':None}, facecolor='white')
    fig.savefig(output.with_suffix('.png'), dpi=160, facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mercurialis', type=Path,
                        default=Path('data/2026-11-22_mercurialis_selfing_refit_v51.json'))
    parser.add_argument('--pulsatilla', type=Path,
                        default=Path('data/2026-11-22_pulsatilla_selfing_refit_v51.json'))
    parser.add_argument('--output', type=Path,
                        default=Path('figures/2026-11-22_selfing_pathway_v51'))
    args = parser.parse_args()
    main(args.mercurialis, args.pulsatilla, args.output)
