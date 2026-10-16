"""Well-mixed C:R:S model used by Figures 3-1 to 3-3.

State order: C, S, R, T. Public initial ratios use C, R, S order.
Shared logistic resource limitation; globally acting toxin kills S only.
Growth rates are relative to mu_R=1. A physical hour interpretation requires
calibration of the reference growth rate. See README.md for assumptions and parameter provenance.
Literature motivates mechanisms, not quantitative parameter validation.
"""

import numpy as np
from scipy.integrate import solve_ivp

# 1. Illustrative parameters: mu_S > mu_R > mu_C; not fitted to this project.
MU_S = 56.2 / 43.6   # ~1.29
MU_R = 1.0
MU_C = 39.3 / 43.6   # ~0.90

K = 1.0               # shared carrying capacity (logistic, resource competition)
KAPPA = 6.0           # toxin production rate per C cell (well-mixed / global toxin)
DELTA = 6.0           # toxin degradation rate
KILL = 8.0            # kill-rate constant of toxin acting on S (mass action)

DEFAULT_PARAMS = dict(mu_S=MU_S, mu_R=MU_R, mu_C=MU_C, K=K,
                       kappa=KAPPA, delta=DELTA, kill=KILL)


def rhs(t, y, p):
    C, S, R, T = y
    N = C + S + R
    dC = p['mu_C'] * C * max(0.0, 1 - N / p['K'])
    dS = p['mu_S'] * S * max(0.0, 1 - N / p['K']) - p['kill'] * T * S
    dR = p['mu_R'] * R * max(0.0, 1 - N / p['K'])
    dT = p['kappa'] * C - p['delta'] * T
    return [dC, dS, dR, dT]


def simulate_batch(y0, t_max=24.0, params=DEFAULT_PARAMS, n_eval=400):
    sol = solve_ivp(rhs, [0, t_max], list(y0) + [0.0], args=(params,),
                     dense_output=True, max_step=0.05,
                     t_eval=np.linspace(0, t_max, n_eval))
    return sol


def simulate_serial_dilution(r0, n_cycles=15, t_cycle=12.0, dilution=100.0,
                              params=DEFAULT_PARAMS, extinct_thresh=1e-6):
    """
    Repeated batch-growth + dilution cycles = the standard well-mixed analogue of a
    long-term coexistence assay (serial passage / 'mixed treatment' as in Kerr et al. 2002).
    r0: (rC, rR, rS) initial fractions (sum=1) at cycle 0, seeded at low total density.
    Returns (trajectory, extinction_cycle). Subsequent rows of trajectory are
    C:R:S fractions; the first row retains the original inoculum densities.
    Extinction is a numerical relative-abundance threshold, not measured loss.
    Detection is checked after each passage; missing-at-inoculation strains are
    first recorded at passage 1. Ties retain C, R, S dictionary order.
    Threshold crossing does not remove a strain from the ODE state.
    """
    C, R, S = r0[0] * 0.01, r0[1] * 0.01, r0[2] * 0.01  # start at 1% of K total
    T = 0.0
    frac_traj = [(C, R, S)]
    extinction_cycle = {'C': None, 'R': None, 'S': None}
    for cyc in range(1, n_cycles + 1):
        sol = solve_ivp(rhs, [0, t_cycle], [C, S, R, T], args=(params,),
                         max_step=0.1)
        C, S, R, T = sol.y[:, -1]
        C, S, R = max(C, 0), max(S, 0), max(R, 0)
        # dilute into fresh medium (toxin washed out, resets to 0)
        C, S, R, T = C / dilution, S / dilution, R / dilution, 0.0
        total = C + S + R
        if total > 0:
            frac_traj.append((C / total, R / total, S / total))
        else:
            frac_traj.append((0, 0, 0))
        for name, val in zip(['C', 'R', 'S'], [C, R, S]):
            if extinction_cycle[name] is None and total > 0 and val / total < extinct_thresh:
                extinction_cycle[name] = cyc
    return np.array(frac_traj), extinction_cycle


if __name__ == "__main__":
    print("module loaded OK")


SQRT3_2 = np.sqrt(3) / 2

def bary_to_xy(rC, rR, rS):
    """Barycentric (C top, R bottom-left, S bottom-right) -> Cartesian for plotting."""
    x = 0.5 * (2 * rS + rC) / (rC + rR + rS)
    y = SQRT3_2 * rC / (rC + rR + rS)
    return (x, y)

def simplex_grid(n):
    """Generate a triangular grid of (rC, rR, rS) with rC+rR+rS=1, resolution n."""
    pts = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            pts.append((i / n, j / n, k / n))
    return np.array(pts)
