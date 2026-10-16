# MOSAIC Cre–Bxb1 model

Predict the A/B composition of individual colony-derived wells, estimate colony-picking requirements, and compare circuit designs. This release contains three files:

- `design_model.py`: the prediction model and design command-line interface.
- `simulation.py`: Wiki figures, saved as PNG and SVG.
- `README.md`: installation, commands, inputs, and outputs.

## Run

Use Python 3.10 or newer. Open a terminal in this directory.

Predictions require only the Python standard library:

```bash
python design_model.py --out output
```

To reproduce Figure 1-3 and Figures 1-6 through 1-15:

```bash
python -m pip install matplotlib==3.10.9
python simulation.py --out figures
```

This generates 11 plots, each saved as PNG and SVG. 

To explore more ratio types with 96 colony picks:

```bash
python design_model.py --n-sites 10 --spacer-bp 600 --colonies 96 --time 24 --objective diversity --out exploration
python simulation.py --input exploration --out exploration/figures
```

`--input` changes Figures 1-11 through 1-15 only. To draw only the five analysis figures for a custom run, add `--plots` to the prediction command.

To search for a specified ratio:

```bash
python design_model.py --target 0.60 --target-tolerance 0.05 --objective target --out target_search
```

The current circuit prediction and the recommended circuit from the design scan are separate results. An unavailable target returns `unreachable`; increasing colony picks cannot create a ratio absent from that circuit.

## Inputs

Time starts at Bxb1 induction after single-colony isolation. Lengths are in bp; all fractions use 0–1, so `0.05` is five percentage points and `0.0001` is 0.01%.

| Option | Meaning | Default |
| --- | --- | --- |
| `--n-sites` | Number of candidate lox71 sites | 5 |
| `--spacer-bp` | Length of one NS sequence, excluding lox | 773 |
| `--cre-base-bp` | Cre distance baseline; alias: `--base-right-bp` | 773 |
| `--bxb1-base-bp` | Invariant attP-to-fixed-lox66 interval | 773 |
| `--lox-site-bp` | Length contribution of one lox site | 34 |
| `--left-arm-bp` | Fixed Bxb1 left arm | 773 |
| `--gamma` | Cre inverse-distance exponent | 1 |
| `--k-a0`, `--k-b0` | Raw-bp cleavage-rate coefficients, h⁻¹ | 0.38, 0.14 |
| `--beta` | Arm-length exponent | 1.3 |
| `--r-a`, `--r-b` | Transcription-machinery binding rates, h⁻¹ | 0.56, 0.01 |
| `--alpha` | Blocked/unblocked rate ratio | 0.01 |
| `--colonies` | Colony picks; alias: `--wells` | 96 |
| `--time`, `--horizon` | Prediction time and time-course end, h | 24, 72 |
| `--target` | Single target A fraction | 0.50 |
| `--library-targets` | Space-separated targets to cover together | 0.40 0.50 0.60 0.70 |
| `--target-tolerance` | Absolute half-width shared by all target intervals | 0.05 |
| `--confidence` | Required picking success probability and marginal interval level | 0.95 |
| `--bin-width` | Fixed ratio-bin width, starting at zero | 0.05 |
| `--min-yield` | Minimum differentiated fraction for a usable structure | 0.0001 |
| `--objective` | `library`, `diversity`, or `target` | `library` |
| `--scan-sites` | Space-separated site counts to scan | 1–10 |
| `--scan-spacers` | Space-separated NS lengths to scan | 200–2,000, step 200 |
| `--max-picks` | Search cap for minimum colony counts | 10,000 |
| `--seed` | Random seed for the example wells only | 2026 |
| `--out` | Prediction output directory | `output` |
| `--plots` | Also draw five custom independent-well figures | Off |

Use `python design_model.py --help` for the complete command reference. Exact coverage supports at most 16 structures and 16 requested target events. The default scan evaluates 100 circuit designs.

## Outputs

| File | Contents |
| --- | --- |
| `site_outcomes.csv` | Each structure's probability, distances, arm lengths, five states, A/B fractions, differentiated yield, expected counts, and marginal count intervals |
| `ratio_types.csv` | Usable ratio bins, their structures, total probabilities, expected counts, and occurrence probabilities |
| `example_wells.csv` | One seeded colony-picking example: plate/well positions, structures, expected compositions, and target matches |
| `coverage_by_picks.csv` | Coverage and target-hit probabilities, and expected type counts, across picking budgets |
| `target_picking.csv` | Matching structures, single-pick probability, hit probability, and required picks for each target |
| `design_scan.csv` | Candidate circuits, ratio ranges, diversity, coverage, target hits, required picks, and objective-specific ranks |
| `per_site_timecourse.csv` | Five-state trajectories, A fractions, differentiated yields, and limiting fractions for each structure |
| `summary.json` | Inputs, metrics, minimum picks, recommended designs, auxiliary pooled results, and assumptions |

`type_a_fraction` describes A among differentiated A/B cells; `differentiated_yield` describes differentiated cells among the normalized whole well. At time zero the conditional A/B fraction is undefined and exported as an empty CSV field. Unusable structures still consume picks; remaining probabilities are not renormalized.

The `library` objective ranks by joint coverage of requested targets, then expected usable diversity and shorter total NS length. `diversity` prioritizes expected usable type count, then ratio span and NS length. `target` prioritizes the chance of at least one target hit, then single-pick probability, target error, and NS length.

## Figure names

| Wiki figure | Output stem; both `.png` and `.svg` are saved |
| --- | --- |
| 1-3 | `figure1-3_Layer1_Cre_selection` |
| 1-6 | `figure1-6_ODE_validation` |
| 1-7 | `figure1-7_1to10` |
| 1-8 | `figure1-8_10to1` |
| 1-9 | `figure1-9_Layer3_timecourse` |
| 1-10 | `figure1-10_Layer3_heatmap` |
| 1-11 | `figure1-11_expected_well_counts` |
| 1-12 | `figure1-12_plate_example` |
| 1-13 | `figure1-13_coverage_probability` |
| 1-14 | `figure1-14_target_picking` |
| 1-15 | `figure1-15_library_design_space` |

References: [Shoura et al., 2012](https://doi.org/10.1093/nar/gks430); [Weber et al., 2016](https://doi.org/10.1186/s12918-016-0290-3); [Zhang and Lutz, 2002](https://doi.org/10.1093/nar/gnf089); [An et al., 2026](https://doi.org/10.1038/s41586-026-10259-3). The project-specific geometry and distance weighting remain modeling assumptions requiring calibration.
