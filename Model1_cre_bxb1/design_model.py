"""Predict colony-derived Cre-Bxb1 ratio libraries using Python 3.10+.

Cre precedes colony isolation. Bxb1 starts separately in each clonal well.
Prediction requires only the standard library; --plots also requires simulation.py
and matplotlib. Fractions use a 0-1 scale and time uses hours.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path


@dataclass(frozen=True)
class Model:
    n_sites: int = 5
    spacer_bp: float = 773.0
    base_right_bp: float = 773.0
    bxb1_base_bp: float = 773.0
    lox_site_bp: float = 34.0
    left_arm_bp: float = 773.0
    gamma: float = 1.0
    k_a0: float = 0.38
    k_b0: float = 0.14
    beta: float = 1.3
    r_a: float = 0.56
    r_b: float = 0.01
    alpha: float = 0.01

    def __post_init__(self):
        if isinstance(self.n_sites, bool) or not isinstance(self.n_sites, int) or self.n_sites < 1:
            raise ValueError('n_sites must be a positive integer')
        for name in ('spacer_bp', 'base_right_bp', 'left_arm_bp',
                     'k_a0', 'k_b0', 'beta', 'alpha'):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be finite and positive')
        for name in ('gamma', 'r_a', 'r_b', 'lox_site_bp', 'bxb1_base_bp'):
            value = getattr(self, name)
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'{name} must be finite and nonnegative')
        if self.alpha > 1:
            raise ValueError('alpha must not exceed 1')


def _nonnegative(value, name):
    if not math.isfinite(value) or value < 0:
        raise ValueError(f'{name} must be finite and nonnegative')


def _decay_integral(rate, t):
    """Integral of exp(-rate*s) from zero to t, evaluated stably."""
    return -math.expm1(-rate * t) / rate


def _binding_state(binding_rate, exit_rate, founder_exit, t):
    """Return the binding-state occupancy and its time integral."""
    gap = founder_exit - exit_rate
    if abs(gap) <= 1e-10 * max(founder_exit, exit_rate):
        rate = (founder_exit + exit_rate) / 2
        x = rate * t
        occupancy = binding_rate * t * math.exp(-x)
        if x < 0.1:
            # Stable series for the integral of s*exp(-rate*s).
            term, total = 1.0, 0.5
            for m in range(1, 30):
                term *= -x / m
                total += term / (m + 2)
            integral = binding_rate * t * t * total
        else:
            integral = binding_rate * (-math.expm1(-x) - x * math.exp(-x)) / rate**2
        return occupancy, integral
    # Choose the slower exponential to avoid exponential overflow.
    slower = min(founder_exit, exit_rate)
    occupancy = (binding_rate * math.exp(-slower * t)
                 * -math.expm1(-abs(gap) * t) / abs(gap))
    integral = binding_rate * (_decay_integral(exit_rate, t)
                              - _decay_integral(founder_exit, t)) / gap
    return occupancy, integral


def cre_outcomes(model: Model):
    """Map fixed lox66 / candidate lox71_i into lengths and Cre weights.

    bxb1_base_bp is the invariant interval from attP to the left fixed
    lox66, with boundaries chosen consistently for the post-Cre arm.
    Each added unit includes NS and one lox-length contribution. Under
    this bookkeeping convention the retained recombined site is counted
    once, not in both the baseline and the repeat contribution.
    base_right_bp retains the original effective pre-Cre distance offset.
    Cre synapsis distances include i NS units and i-1 intervening lox71
    sites. The Bxb1 baseline is independent of this pre-Cre offset.
    """
    unit = model.spacer_bp + model.lox_site_bp
    distances = [model.base_right_bp + i * model.spacer_bp + (i - 1) * model.lox_site_bp
                 for i in range(1, model.n_sites + 1)]
    weights = [(distance / distances[0]) ** (-model.gamma)
               for distance in distances]
    total = math.fsum(weights)
    return [dict(site=i, left_arm_bp=model.left_arm_bp,
                 right_arm_bp=model.bxb1_base_bp + i * unit,
                 repeat_unit_bp=unit,
                 synapsis_distance_bp=distance,
                 selection_probability=weight / total)
            for i, (distance, weight) in enumerate(zip(distances, weights), 1)]


def _site_at(model, h_b, t):
    """Solve all five states exactly for a single post-Cre structure."""
    k_a = model.k_a0 * model.left_arm_bp ** (-model.beta)
    k_b = model.k_b0 * h_b ** (-model.beta)
    v_a, v_b = model.alpha * k_a, model.alpha * k_b
    q = k_a + k_b + model.r_a + model.r_b
    m_a, int_m_a = _binding_state(model.r_a, v_a + k_b, q, t)
    m_b, int_m_b = _binding_state(model.r_b, v_b + k_a, q, t)
    int_d = _decay_integral(q, t)
    l_a = k_a * int_d + v_a * int_m_a + k_a * int_m_b
    l_b = k_b * int_d + v_b * int_m_b + k_b * int_m_a
    differentiated = l_a + l_b
    limiting = (k_a + model.r_a * v_a / (v_a + k_b)
                + model.r_b * k_a / (v_b + k_a)) / q
    return dict(D=math.exp(-q * t), L_A=l_a, L_B=l_b, M_A=m_a, M_B=m_b,
                k_a=k_a, k_b=k_b, differentiated_yield=differentiated,
                type_a_fraction=l_a / differentiated if differentiated > 0 else None,
                limiting_type_a_fraction=limiting)


def population_at(model: Model, t: float = 24):
    """Return pooled composition, yield, and the weighted site-level fraction.

    Fractions use a 0-1 scale. Composition is None at t=0 because there
    are no differentiated cells. Yield is relative to the initial population.
    """
    _nonnegative(t, 't')
    outcomes = cre_outcomes(model)
    states = [_site_at(model, row['right_arm_bp'], t) for row in outcomes]
    pooled = {key: math.fsum(row['selection_probability'] * state[key]
                            for row, state in zip(outcomes, states))
              for key in ('D', 'L_A', 'L_B', 'M_A', 'M_B')}
    differentiated = pooled['L_A'] + pooled['L_B']
    fraction = pooled['L_A'] / differentiated if differentiated > 0 else None
    weighted = (math.fsum(row['selection_probability'] * state['type_a_fraction']
                         for row, state in zip(outcomes, states)) if t > 0 else None)
    return dict(time_h=t, type_a_fraction=fraction,
                type_b_fraction=1 - fraction if fraction is not None else None,
                differentiated_yield=differentiated,
                type_a_yield=pooled['L_A'], type_b_yield=pooled['L_B'],
                weighted_site_fraction=weighted,
                limiting_type_a_fraction=math.fsum(
                    row['selection_probability'] * state['limiting_type_a_fraction']
                    for row, state in zip(outcomes, states)), **pooled)


def sampling_report(model: Model, times, min_yield=0.0001, tolerance=0.01,
                    window_h=6, total_events=100_000, min_per_type=1,
                    window_steps=60):
    """Screen times using yield, expected events, and sampled composition drift.

    Future drift is evaluated on a uniform grid over the window. This is a
    numerical screening criterion, not a proof of asymptotic steady state.
    """
    for value, name in [(min_yield, 'min_yield'), (tolerance, 'tolerance'),
                        (min_per_type, 'min_per_type')]:
        _nonnegative(value, name)
    if min_yield > 1 or tolerance > 1:
        raise ValueError('min_yield and tolerance must not exceed 1')
    if not math.isfinite(window_h) or window_h <= 0:
        raise ValueError('window_h must be finite and positive')
    if not isinstance(window_steps, int) or window_steps < 1:
        raise ValueError('window_steps must be a positive integer')
    if not isinstance(total_events, int) or total_events < 1:
        raise ValueError('total_events must be a positive integer')
    reports = []
    for t in sorted(set(times)):
        endpoint = population_at(model, t)
        fraction = endpoint['type_a_fraction']
        drift = (max(abs(population_at(model, t + window_h * j / window_steps)
                         ['type_a_fraction'] - fraction)
                     for j in range(1, window_steps + 1))
                 if fraction is not None else None)
        events_a = total_events * endpoint['type_a_yield']
        events_b = total_events * endpoint['type_b_yield']
        yield_ok = endpoint['differentiated_yield'] >= min_yield
        events_ok = min(events_a, events_b) >= min_per_type
        drift_ok = drift is not None and drift <= tolerance
        reports.append(dict(time_h=t, type_a_fraction=fraction,
                            differentiated_yield=endpoint['differentiated_yield'],
                            expected_type_a_events=events_a, expected_type_b_events=events_b,
                            max_sampled_fraction_change=drift, yield_pass=yield_ok,
                            events_pass=events_ok, drift_pass=drift_ok,
                            eligible=yield_ok and events_ok and drift_ok))
    earliest = next((row['time_h'] for row in reports if row['eligible']), None)
    return reports, earliest


def design_scan(model: Model, t=24, target=0.50, min_yield=0.0001,
                site_counts=range(1, 11), spacer_lengths=range(200, 2001, 200)):
    """Return feasible designs sorted by pooled target error and spacer burden."""
    _nonnegative(t, 't')
    if not math.isfinite(target) or not 0 <= target <= 1:
        raise ValueError('target must be between 0 and 1')
    if not math.isfinite(min_yield) or not 0 <= min_yield <= 1:
        raise ValueError('min_yield must be between 0 and 1')
    counts, lengths = list(site_counts), list(spacer_lengths)
    designs = []
    for n in counts:
        for spacer in lengths:
            result = population_at(replace(model, n_sites=n, spacer_bp=spacer), t)
            if result['type_a_fraction'] is None or result['differentiated_yield'] < min_yield:
                continue
            designs.append(dict(n_sites=n, spacer_bp=spacer, time_h=t,
                                type_a_fraction=result['type_a_fraction'],
                                differentiated_yield=result['differentiated_yield'],
                                target_error=abs(result['type_a_fraction'] - target),
                                total_spacer_bp=n * spacer))
    designs.sort(key=lambda row: (row['target_error'], row['total_spacer_bp'],
                                  row['n_sites'], row['spacer_bp']))
    return [dict(rank=i, **row) for i, row in enumerate(designs, 1)]


def _write_csv(path, rows, columns):
    with path.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)



# Independent-well sampling, exact coverage, and design ranking.
import argparse
import csv
import json
import math
import random
from dataclasses import asdict, dataclass, replace
from pathlib import Path

MAX_EXACT_EVENTS = 16


@dataclass(frozen=True)
class Settings:
    colonies: int = 96
    time_h: float = 24.0
    horizon_h: float = 72.0
    target: float = 0.50
    target_tolerance: float = 0.05
    confidence: float = 0.95
    bin_width: float = 0.05
    min_yield: float = 0.0001
    library_targets: tuple = (0.40, 0.50, 0.60, 0.70)
    max_picks: int = 10000
    seed: int = 2026
    objective: str = 'library'

    def __post_init__(self):
        object.__setattr__(self, 'library_targets', tuple(sorted(set(self.library_targets))))
        for name in ('colonies', 'max_picks'):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f'{name} must be a positive integer')
        for name in ('target', 'target_tolerance', 'min_yield'):
            value = getattr(self, name)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f'{name} must be between zero and one')
        if not math.isfinite(self.confidence) or not 0 < self.confidence < 1:
            raise ValueError('confidence must be strictly between zero and one')
        if not math.isfinite(self.bin_width) or not 0 < self.bin_width <= 1:
            raise ValueError('bin_width must be greater than zero and at most one')
        if not math.isfinite(self.time_h) or self.time_h <= 0:
            raise ValueError('time_h must be positive; composition is undefined at zero')
        if not math.isfinite(self.horizon_h) or self.horizon_h < self.time_h:
            raise ValueError('horizon_h must be finite and at least time_h')
        if not self.library_targets or len(self.library_targets) > MAX_EXACT_EVENTS:
            raise ValueError(f'Provide 1-{MAX_EXACT_EVENTS} library targets')
        if any(not math.isfinite(t) or not 0 <= t <= 1 for t in self.library_targets):
            raise ValueError('Library targets must be between zero and one')
        if self.objective not in ('library', 'diversity', 'target'):
            raise ValueError('objective must be library, diversity, or target')


def hit_probability(q, picks):
    """Probability of at least one success in independent colony picks."""
    if picks < 0 or int(picks) != picks or not math.isfinite(q) or not 0 <= q <= 1:
        raise ValueError('Require integer picks >= 0 and a probability between zero and one')
    if picks == 0 or q == 0:
        return 0.0
    if q == 1:
        return 1.0
    return -math.expm1(picks * math.log1p(-q))


class Coverage:
    """Exact inclusion-exclusion for possibly overlapping success events.

    Each event is a set of structure indices. An event is covered when at
    least one picked colony belongs to that set. Unusable structures retain
    their original probability and are never renormalized away.
    """
    def __init__(self, probabilities, events):
        self.probabilities = tuple(probabilities)
        if not probabilities or any(not math.isfinite(p) or p < 0 for p in probabilities):
            raise ValueError('Require nonnegative finite structure probabilities')
        if not math.isclose(math.fsum(probabilities), 1, abs_tol=1e-10):
            raise ValueError('Structure probabilities must sum to one')
        if len(events) > MAX_EXACT_EVENTS:
            raise ValueError(f'Exact coverage supports at most {MAX_EXACT_EVENTS} events')
        self.events = tuple(frozenset(event) for event in events)
        if any(any(i < 0 or i >= len(probabilities) for i in event) for event in self.events):
            raise ValueError('Event structure index is out of range')
        self.possible = bool(self.events) and all(
            math.fsum(probabilities[i] for i in event) > 0 for event in self.events)
        self.terms = []
        # Merge equal union masses to reduce cancellation and repeated work.
        if self.possible:
            union_masks = [0] * (1 << len(self.events))
            coefficients = {}
            masks = [sum(1 << i for i in event) for event in self.events]
            for subset in range(1, 1 << len(self.events)):
                low = subset & -subset
                bit = low.bit_length() - 1
                union = union_masks[subset ^ low] | masks[bit]
                union_masks[subset] = union
                sign = -1 if subset.bit_count() % 2 else 1
                coefficients[union] = coefficients.get(union, 0) + sign
            for union, coefficient in coefficients.items():
                mass = min(1.0, math.fsum(probabilities[i] for i in range(len(probabilities))
                                         if union & (1 << i)))
                if coefficient:
                    self.terms.append((coefficient, mass))

    def probability(self, picks):
        if isinstance(picks, bool) or not isinstance(picks, int) or picks < 0:
            raise ValueError('picks must be a nonnegative integer')
        if picks == 0 or not self.possible:
            return 0.0
        result = math.fsum([1.0] + [
            coefficient * (0.0 if mass == 1 else math.exp(picks * math.log1p(-mass)))
            for coefficient, mass in self.terms])
        return max(0.0, min(1.0, result))

    def required_picks(self, confidence=.95, limit=10000):
        """Return the smallest tested integer or an explicit non-success status."""
        if not 0 < confidence < 1 or not isinstance(limit, int) or limit < 1:
            raise ValueError('Require confidence in (0,1) and a positive integer limit')
        if not self.possible:
            return dict(picks=None, status='unreachable')
        upper = 1
        while upper < limit and self.probability(upper) < confidence:
            upper = min(limit, upper * 2)
        if self.probability(upper) < confidence:
            return dict(picks=None, status='above_max_picks')
        lower = 0
        while upper - lower > 1:
            middle = (upper + lower) // 2
            if self.probability(middle) >= confidence:
                upper = middle
            else:
                lower = middle
        return dict(picks=upper, status='reached')


def binomial_interval(n, p, confidence=.95):
    """Equal-tailed marginal predictive interval for one category's well count."""
    if p <= 0:
        return 0, 0
    if p >= 1:
        return n, n
    low_tail, high_tail = (1 - confidence) / 2, (1 + confidence) / 2
    accumulated, lower, upper = 0.0, None, n
    log_p, log_q, log_n = math.log(p), math.log1p(-p), math.lgamma(n + 1)
    for count in range(n + 1):
        probability = math.exp(log_n - math.lgamma(count + 1) - math.lgamma(n - count + 1)
                               + count * log_p + (n - count) * log_q)
        accumulated += probability
        if lower is None and accumulated >= low_tail:
            lower = count
        if accumulated >= high_tail:
            upper = count
            break
    return lower if lower is not None else n, upper


def _bin(fraction, width):
    number = math.ceil(1 / width)
    index = min(number - 1, max(0, int(math.floor(fraction / width + 1e-12))))
    lower, upper = index * width, min(1.0, (index + 1) * width)
    label = f'{lower:.3f}-{upper:.3f}'
    return index, lower, upper, label


def structure_predictions(model, settings, intervals=True):
    """One row per clonal structure; each row predicts a class of wells."""
    if model.n_sites > MAX_EXACT_EVENTS:
        raise ValueError(f'Exact all-structure coverage supports up to {MAX_EXACT_EVENTS} sites')
    rows = []
    for row in cre_outcomes(model):
        state = _site_at(model, row['right_arm_bp'], settings.time_h)
        p, fraction = row['selection_probability'], state['type_a_fraction']
        index, lower, upper, bin_name = _bin(fraction, settings.bin_width)
        usable = state['differentiated_yield'] >= settings.min_yield
        low, high = binomial_interval(settings.colonies, p, settings.confidence) if intervals else (None, None)
        rows.append(dict(**row, **state, type_b_fraction=1-fraction,
                         expected_wells=settings.colonies*p,
                         count_sd=math.sqrt(settings.colonies*p*(1-p)),
                         count_interval_low=low, count_interval_high=high,
                         probability_present=hit_probability(p, settings.colonies),
                         usable=usable, ratio_bin=index, bin_lower=lower, bin_upper=upper,
                         bin_label=bin_name, target_error=abs(fraction-settings.target),
                         target_match=usable and abs(fraction-settings.target) <= settings.target_tolerance+1e-12))
    return rows


def library_metrics(rows, settings, include_curves=True):
    p = [row['selection_probability'] for row in rows]
    usable = [i for i, row in enumerate(rows) if row['usable']]
    bin_indices = sorted({rows[i]['ratio_bin'] for i in usable})
    bin_events = [{i for i in usable if rows[i]['ratio_bin'] == b} for b in bin_indices]
    structure_cover = Coverage(p, [{i} for i in range(len(rows))])
    ratio_cover = Coverage(p, bin_events)
    target_values = sorted(set(settings.library_targets + (settings.target,)))
    target_events = {t: {i for i in usable
                         if abs(rows[i]['type_a_fraction']-t) <= settings.target_tolerance+1e-12}
                     for t in target_values}
    library_cover = Coverage(p, [target_events[t] for t in settings.library_targets])
    target_rows = []
    for target in target_values:
        event = target_events[target]
        q = math.fsum(p[i] for i in event)
        required = Coverage(p, [event]).required_picks(settings.confidence, settings.max_picks)
        target_rows.append(dict(target=target, tolerance=settings.target_tolerance,
                                matching_sites=';'.join(str(rows[i]['site']) for i in sorted(event)),
                                single_pick_probability=q,
                                expected_matching_wells=settings.colonies*q,
                                probability_at_least_one=hit_probability(q, settings.colonies),
                                required_picks=required['picks'], requirement_status=required['status']))
    type_rows = []
    for b, event in zip(bin_indices, bin_events):
        q = math.fsum(p[i] for i in event)
        representative = rows[min(event)]
        low, high = binomial_interval(settings.colonies, q, settings.confidence) if include_curves else (None,None)
        type_rows.append(dict(ratio_bin=b, bin_lower=representative['bin_lower'],
                              bin_upper=representative['bin_upper'], bin_label=representative['bin_label'],
                              sites=';'.join(str(rows[i]['site']) for i in sorted(event)),
                              probability=q, expected_wells=settings.colonies*q,
                              probability_present=hit_probability(q, settings.colonies),
                              count_interval_low=low, count_interval_high=high))
    single_target = next(row for row in target_rows if row['target'] == settings.target)
    fractions = [rows[i]['type_a_fraction'] for i in usable]
    summary = dict(number_structures=len(rows), number_usable_structures=len(usable),
                   number_ratio_types=len(type_rows),
                   usable_probability=math.fsum(p[i] for i in usable),
                   ratio_min=min(fractions) if fractions else None,
                   ratio_max=max(fractions) if fractions else None,
                   ratio_span=max(fractions)-min(fractions) if fractions else 0,
                   expected_distinct_structures=math.fsum(hit_probability(q,settings.colonies) for q in p),
                   expected_distinct_ratio_types=math.fsum(r['probability_present'] for r in type_rows),
                   probability_all_structures=structure_cover.probability(settings.colonies),
                   probability_all_ratio_types=ratio_cover.probability(settings.colonies),
                   probability_all_library_targets=library_cover.probability(settings.colonies),
                   probability_target_hit=single_target['probability_at_least_one'],
                   single_pick_target_probability=single_target['single_pick_probability'],
                   best_usable_target_error=min((rows[i]['target_error'] for i in usable), default=None),
                   expected_mean_well_fraction=math.fsum(p[i]*rows[i]['type_a_fraction'] for i in range(len(rows))),
                   mean_well_differentiated_yield=math.fsum(p[i]*rows[i]['differentiated_yield'] for i in range(len(rows))),
                   required_all_structures=structure_cover.required_picks(settings.confidence, settings.max_picks),
                   required_all_ratio_types=ratio_cover.required_picks(settings.confidence, settings.max_picks),
                   required_all_library_targets=library_cover.required_picks(settings.confidence, settings.max_picks),
                   required_target=Coverage(p,[target_events[settings.target]]).required_picks(
                       settings.confidence,settings.max_picks))
    curve = []
    if include_curves:
        end = max(192, settings.colonies*2)
        # Limit table size while retaining all small pick counts and the plate size.
        grid = set(range(min(end,512)+1))
        if end > 512:
            grid.update(round(end*j/256) for j in range(257))
        grid.add(settings.colonies)
        for value in [summary['required_all_structures'], summary['required_all_ratio_types'],
                      summary['required_all_library_targets'], summary['required_target']]:
            if value['picks'] is not None and value['picks'] <= end:
                grid.add(value['picks'])
        for picks in sorted(grid):
            curve.append(dict(picks=picks, probability_all_structures=structure_cover.probability(picks),
                               probability_all_ratio_types=ratio_cover.probability(picks),
                               probability_all_library_targets=library_cover.probability(picks),
                               probability_target_hit=hit_probability(single_target['single_pick_probability'],picks),
                               expected_distinct_structures=math.fsum(hit_probability(q,picks) for q in p),
                               expected_distinct_ratio_types=math.fsum(hit_probability(r['probability'],picks)
                                                                      for r in type_rows)))
    return summary, type_rows, target_rows, curve


def simulate_wells(rows, settings):
    """One reproducible plate realization; it does not model measurement noise."""
    generator = random.Random(settings.seed)
    indices = generator.choices(range(len(rows)), weights=[r['selection_probability'] for r in rows],
                                k=settings.colonies)
    result = []
    for index, structure_index in enumerate(indices):
        site = rows[structure_index]
        position = index % 96
        result.append(dict(colony=index+1, plate=index//96+1,
                           well=f'{chr(65+position//12)}{position%12+1:02d}',
                           site=site['site'], right_arm_bp=site['right_arm_bp'],
                           type_a_fraction=site['type_a_fraction'], type_b_fraction=site['type_b_fraction'],
                           differentiated_yield=site['differentiated_yield'], usable=site['usable'],
                           ratio_bin=site['ratio_bin'], target_match=site['target_match']))
    return result


def scan_designs(model, settings, counts=range(1,11), spacers=range(200,2001,200)):
    """Compare designs by library coverage, expected diversity, or a single target."""
    result = []
    counts, spacers = list(counts), list(spacers)
    for n in counts:
        for spacer in spacers:
            candidate = replace(model,n_sites=n,spacer_bp=spacer)
            sites = structure_predictions(candidate,settings,intervals=False)
            metrics,_,_,_ = library_metrics(sites,settings,include_curves=False)
            result.append(dict(n_sites=n,spacer_bp=spacer,
                               total_spacer_bp=n*spacer,
                               number_usable_structures=metrics['number_usable_structures'],
                               number_ratio_types=metrics['number_ratio_types'],
                               usable_probability=metrics['usable_probability'],
                               ratio_min=metrics['ratio_min'],ratio_max=metrics['ratio_max'],
                               ratio_span=metrics['ratio_span'],
                               expected_distinct_ratio_types=metrics['expected_distinct_ratio_types'],
                               probability_all_ratio_types=metrics['probability_all_ratio_types'],
                               probability_all_library_targets=metrics['probability_all_library_targets'],
                               probability_target_hit=metrics['probability_target_hit'],
                               single_pick_target_probability=metrics['single_pick_target_probability'],
                               best_usable_target_error=metrics['best_usable_target_error'],
                               required_picks_all_types=metrics['required_all_ratio_types']['picks'],
                               all_types_requirement_status=metrics['required_all_ratio_types']['status'],
                               required_picks_target=metrics['required_target']['picks'],
                               target_requirement_status=metrics['required_target']['status'],
                               required_picks_library=metrics['required_all_library_targets']['picks'],
                               library_requirement_status=metrics['required_all_library_targets']['status']))
    def ordering(row):
        valid = row['number_usable_structures'] > 0
        error = row['best_usable_target_error'] if valid else math.inf
        if settings.objective=='diversity':
            return (not valid,-row['expected_distinct_ratio_types'],-row['ratio_span'],row['total_spacer_bp'])
        if settings.objective=='target':
            return (not valid,-row['probability_target_hit'],-row['single_pick_target_probability'],error,row['total_spacer_bp'])
        return (not valid,-row['probability_all_library_targets'],
                -row['expected_distinct_ratio_types'],row['total_spacer_bp'])
    result.sort(key=ordering)
    return [dict(rank=i,**row) for i,row in enumerate(result,1)]


def write_csv(path, rows, columns=None):
    if columns is None:
        if not rows:
            raise ValueError('Specify columns for an empty table')
        columns=list(rows[0])
    with Path(path).open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def run_analysis(model,settings,out,counts=range(1,11),spacers=range(200,2001,200)):
    out=Path(out)
    counts,spacers=list(counts),list(spacers)
    if not counts or not spacers:
        raise ValueError('The design scan requires nonempty site and spacer lists')
    rows=structure_predictions(model,settings)
    metrics,types,targets,coverage=library_metrics(rows,settings)
    wells=simulate_wells(rows,settings)
    for row in rows:
        row['example_wells']=sum(w['site']==row['site'] for w in wells)
    for row in types:
        row['example_wells']=sum(w['usable'] and w['ratio_bin']==row['ratio_bin'] for w in wells)
    designs=scan_designs(model,settings,counts,spacers)
    def objective_possible(row):
        if settings.objective == 'library':
            return row['library_requirement_status'] != 'unreachable'
        if settings.objective == 'target':
            return row['single_pick_target_probability'] > 0
        return row['number_ratio_types'] > 0
    feasible=[row for row in designs if objective_possible(row)]
    target_feasible=[row for row in designs if row['single_pick_target_probability'] > 0]
    diverse_feasible=[row for row in designs if row['number_ratio_types'] > 0]
    times=sorted({settings.horizon_h*i/144 for i in range(145)}|{settings.time_h})
    timecourse=[]
    for site in rows:
        for time in times:
            state=_site_at(model,site['right_arm_bp'],time)
            timecourse.append(dict(site=site['site'],time_h=time,
                                   selection_probability=site['selection_probability'],**state))
    summary=dict(model=asdict(model),settings=asdict(settings),metrics=metrics,
                  targets=targets,scan_settings=dict(site_counts=counts,spacer_lengths_bp=spacers),
                  best_design=feasible[0] if feasible else None,
                  best_by_single_target=min(target_feasible,key=lambda r:(-r['probability_target_hit'],
                                                                  -r['single_pick_target_probability'],
                                                                  r['best_usable_target_error'],r['total_spacer_bp'])) if target_feasible else None,
                  best_by_expected_diversity=max(diverse_feasible,key=lambda r:(r['expected_distinct_ratio_types'],
                                                                       r['ratio_span'],-r['total_spacer_bp'])) if diverse_feasible else None,
                  auxiliary_pooled_population=population_at(model,settings.time_h),
                  assumptions=[
                      'Cre rearrangement is complete and stable before colony picking.',
                      'One picked colony carries one modeled structure; Bxb1 is induced after transfer to its own well.',
                      'Colony picks are independent and unbiased; no differential plating or growth weights are included.',
                      'Fixed lox66 and candidate lox71_i: h_B(i)=L_B,0+i*(L_NS+L_lox71), where L_B,0 is the attP-to-fixed-lox66 interval. Cre distances include intervening lox71 sites: d_i=L_base,Cre+i*L_NS+(i-1)*L_lox71. The two baselines are independent.',
                      'One genotype has one deterministic expected within-well fraction at a given time; biological and detection noise are not modeled.',
                      'Ratio types are fixed bins starting at zero, not experimentally validated resolution classes.',
                      'Only structures meeting min_yield contribute to usable bins and target hits; their probabilities are not renormalized.',
                      'The example plate is seeded sampling, not experimental data or guaranteed future well positions.',
                      'All-type coverage uses exact inclusion-exclusion up to 16 events; floating-point roundoff is clipped to [0,1].'])
    out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'site_outcomes.csv',rows)
    write_csv(out/'ratio_types.csv',types,
               list(types[0]) if types else ['ratio_bin','bin_lower','bin_upper','bin_label','sites','probability',
                                             'expected_wells','probability_present','count_interval_low','count_interval_high','example_wells'])
    write_csv(out/'example_wells.csv',wells)
    write_csv(out/'coverage_by_picks.csv',coverage)
    write_csv(out/'target_picking.csv',targets)
    write_csv(out/'design_scan.csv',designs)
    write_csv(out/'per_site_timecourse.csv',timecourse)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False),encoding='utf-8')
    return summary


def main():
    parser=argparse.ArgumentParser(description='Predict independent-well ratios, colony coverage, and candidate designs.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--out',type=Path,default=Path('output'),help='Directory for seven CSV tables and summary.json')
    parser.add_argument('--colonies','--wells',dest='colonies',type=int,default=96,help='Number of picked colonies; one per well')
    parser.add_argument('--time',type=float,default=24,help='Hours after Bxb1 induction within each well')
    parser.add_argument('--horizon',type=float,default=72,help='End of exported per-structure time courses (h)')
    parser.add_argument('--target',type=float,default=.50,help='Single-well target type A fraction')
    parser.add_argument('--target-tolerance',type=float,default=.05,help='Absolute half-width shared by single-target and library-target intervals')
    parser.add_argument('--confidence',type=float,default=.95,help='Desired picking success probability and marginal count interval level')
    parser.add_argument('--bin-width',type=float,default=.05,help='Width of fixed usable-ratio bins on [0,1]')
    parser.add_argument('--min-yield',type=float,default=.0001,help='Minimum differentiated fraction for a usable well type')
    parser.add_argument('--library-targets',type=float,nargs='+',default=[.4,.5,.6,.7],help='Desired well-ratio targets; intervals may overlap')
    parser.add_argument('--max-picks',type=int,default=10000,help='Search cap for minimum colony counts')
    parser.add_argument('--seed',type=int,default=2026,help='Random seed used only for the example plate')
    parser.add_argument('--objective',choices=['library','diversity','target'],default='library',help='Design ranking objective')
    parser.add_argument('--scan-sites',type=int,nargs='+',default=list(range(1,11)),help='Candidate site counts; at most 16')
    parser.add_argument('--scan-spacers',type=float,nargs='+',default=list(range(200,2001,200)),help='Candidate spacer lengths (bp)')
    parser.add_argument('--plots',action='store_true',help='Generate five independent-well PNG/SVG pairs; requires matplotlib and simulation.py')
    descriptions=dict(n_sites='Number of candidate lox71 sites',spacer_bp='NS spacer length (bp)',
                       base_right_bp='Original Cre distance offset (bp); not the Bxb1 baseline',left_arm_bp='Fixed left-arm length (bp)',
                       lox_site_bp='Length of each lox71 / retained lox contribution (bp)',
                       bxb1_base_bp='Invariant attP-to-fixed-lox66 interval in h_B (bp)',
                       gamma='Cre inverse-distance exponent',k_a0='A-arm raw-bp rate coefficient (h^-1)',
                       k_b0='B-arm raw-bp rate coefficient (h^-1)',beta='Arm-length rate exponent',
                       r_a='A-side binding rate (h^-1)',r_b='B-side binding rate (h^-1)',alpha='Blocked/unblocked rate factor')
    for name,field in Model.__dataclass_fields__.items():
        flags=['--'+name.replace('_','-')]
        if name=='base_right_bp':
            flags.append('--cre-base-bp')
        parser.add_argument(*flags,type=int if name=='n_sites' else float,
                            default=field.default,help=descriptions[name])
    args=parser.parse_args()
    try:
        model=Model(**{name:getattr(args,name) for name in Model.__dataclass_fields__})
        settings=Settings(colonies=args.colonies,time_h=args.time,horizon_h=args.horizon,
                           target=args.target,target_tolerance=args.target_tolerance,
                           confidence=args.confidence,bin_width=args.bin_width,min_yield=args.min_yield,
                           library_targets=tuple(args.library_targets),max_picks=args.max_picks,
                           seed=args.seed,objective=args.objective)
        summary=run_analysis(model,settings,args.out,args.scan_sites,args.scan_spacers)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(dict(settings=summary['settings'],metrics=summary['metrics'],
                          best_design=summary['best_design']),indent=2,allow_nan=False))
    print(f'Saved seven CSV files and summary.json to {args.out.resolve()}')
    if args.plots:
        try:
            from simulation import plot_library_figures as generate_figures
        except ModuleNotFoundError as error:
            parser.error(f'Prediction tables were saved; plotting requires matplotlib. '
                         f'Run python -m pip install matplotlib==3.10.9. Missing module: {error.name}')
        generate_figures(args.out)


if __name__=='__main__':
    main()
