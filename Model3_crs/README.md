# C:R:S application model

A theoretical application of MOSAIC's programmable cell-composition library: compare initial producer (C), resistant (R) and sensitive (S) fractions by their subsequent trajectories and transient coexistence in well-mixed culture. 

It contains three files:

- `model.py`: original ODE and serial-passage solver, plus simplex utilities.
- `simulation.py`: reproduce all four figures and simulation CSVs; optional numerical checks.
- `README.md`: installation, parameters, interpretation and reuse terms.

## Run

Verified with Python 3.10.20. Open a terminal in this folder:

```bash
python -m pip install numpy==2.2.5 scipy==1.15.3 matplotlib==3.10.9
python simulation.py --out output
python simulation.py --out output --check
```

For one figure:

```bash
python simulation.py --figure 3-4 --out output
```

## Analyses and outputs

| Figure | Settings | Output stem |
| --- | --- | --- |
| 3-1 | C:R:S=1:1:1, 5:1:1, 0.1:1:1, 0.02:1:1; duration 30; 400 samples | `figure3-1_trajectories` |
| 3-2 | Simplex resolution 28; 435 compositions; 60 passages | `figure3-2_deterministic_ternary` |
| 3-3 | 22 C fractions from 0.01 to 0.50; 18 p values from 0.5 to 10; 40 passages | `figure3-3_biodiversity_law` |
| 3-4 | 36 C fractions; initial R=S; seven dynamical scenarios, 60 passages | `figure3-4_robustness` |
