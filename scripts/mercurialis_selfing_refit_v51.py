#!/usr/bin/env python3
"""Reproduce the assigned-seed summary and estimate the within-array selfing association.

Usage: python scripts/mercurialis_selfing_refit_v51.py --output output.json [--archive private.zip]
Dependencies: numpy, pandas, scipy. The source ZIP is read in memory. All outputs
are aggregate. The association is conditional on both flower sexes being recorded,
mature seed assignment and known above-ground biomass. No survival model is used.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import zipfile

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm

VERSION = '10.5281/zenodo.15644955'
URL = 'https://zenodo.org/api/records/15644955/files/Geito_NewPhyto_v2.zip/content'
SHA = 'a18b7d296110be1e9817d6e1a3994365a4e7ae0aab50435e7e6883e1f3de4042'
MD5 = 'ea198a76d60cdf412638c77e6e0e3d80'
FILES = ['3G.csv', '3G_ID01.csv', 'Genotype_3P.csv', 'Selfing.R',
         'Fitness.R', 'Trade-off.R', 'README.md.txt']
DRAWS = 1999
SEED = 70356


def source(local_path):
    if local_path is None:
        req = urllib.request.Request(URL, headers={'User-Agent': 'Published-data-reproduction/0.1'})
        with urllib.request.urlopen(req, timeout=90) as response:
            blob = response.read()
    else:
        blob = local_path.read_bytes()
    if hashlib.sha256(blob).hexdigest() != SHA or hashlib.md5(blob).hexdigest() != MD5:
        raise ValueError('Source archive hash mismatch')
    with zipfile.ZipFile(io.BytesIO(blob)) as archive:
        files = {name: archive.read('Geito_NewPhyto_v2/' + name) for name in FILES}
    hashes = {name: hashlib.sha256(data).hexdigest() for name, data in files.items()}
    frames = {name: pd.read_csv(io.BytesIO(data)) for name, data in files.items() if name.endswith('.csv')}
    return frames, hashes


def genotype_summary(genotypes):
    offspring = genotypes.loc[genotypes.Parent_Offspring.eq('Offspring')]
    alleles = offspring.iloc[:, 3:]
    if alleles.shape[1] != 18:
        raise ValueError('Expected paired allele columns at nine loci')
    # Each locus requires two recorded allele entries, including homozygotes.
    typed = alleles.notna().to_numpy().reshape(len(offspring), 9, 2).all(axis=2).sum(axis=1)
    return {
        'genotype_parent_rows': int(genotypes.Parent_Offspring.eq('Parent').sum()),
        'genotype_progeny_rows': len(offspring),
        'genotype_progeny_at_least_five_complete_loci': int((typed >= 5).sum()),
        'genotype_progeny_more_than_five_complete_loci': int((typed > 5).sum()),
    }


def design(data, separate_slopes=False, male_count=None):
    pop = pd.get_dummies(data.Population).to_numpy(dtype=float)
    male = np.log2(data.Male_flower_total.to_numpy() if male_count is None
                   else np.full(len(data), male_count))
    predictor = pop * male[:, None] if separate_slopes else male
    return np.column_stack([pop, predictor, np.log2(data.Female_flower_total),
                            np.log2(data.Biomass_whole)])


def fit(X, k, n):
    beta = np.zeros(X.shape[1])
    for _ in range(50):
        p = expit(X @ beta)
        w = n * p * (1 - p)
        info = X.T @ (w[:, None] * X)
        delta = np.linalg.solve(info, X.T @ (k - n * p))
        beta += delta
        if np.max(np.abs(delta)) < 1e-9:
            return beta
    raise ArithmeticError('Grouped binomial fit failed to converge')


def hc3(X, k, n, beta):
    p = expit(X @ beta)
    w = n * p * (1 - p)
    bread = np.linalg.inv(X.T @ (w[:, None] * X))
    h = w * np.einsum('ij,jk,ik->i', X, bread, X)
    scores = X * ((k - n * p) / (1 - h))[:, None]
    se = np.sqrt(np.diag(bread @ (scores.T @ scores) @ bread))
    return np.exp(beta[:, None] + se[:, None] * np.array([-1., 1.]) * norm.ppf(.975))


def estimate(data):
    k = data.Selfed_seed_number.to_numpy()
    n = data.Genotyped_number.to_numpy()
    X = design(data)
    group_X = design(data, separate_slopes=True)
    beta = fit(X, k, n)
    group_beta = fit(group_X, k, n)
    pop_groups = data.groupby('Population', sort=True).indices
    offset = len(pop_groups)
    interval = hc3(X, k, n, beta)
    group_interval = hc3(group_X, k, n, group_beta)
    scenarios = [design(data, male_count=m) for m in [100, 200]]
    probabilities = [float(expit(Z @ beta).mean()) for Z in scenarios]

    # Sample whole maternal rows within each of the three fixed observed arrays.
    # This quantifies sampling variation among measured mothers; it does not
    # represent new garden realizations or correction for paternity call errors.
    rng = np.random.default_rng(SEED)
    boot = np.empty((DRAWS, 2))
    for i in range(DRAWS):
        ix = np.concatenate([rng.choice(rows, size=len(rows), replace=True)
                             for rows in pop_groups.values()])
        b = fit(X[ix], k[ix], n[ix])
        boot[i, 0] = np.exp(b[offset])
        prob = [expit(Z[ix] @ b).mean() for Z in scenarios]
        boot[i, 1] = prob[1] - prob[0]

    per_pop = {}
    for i, (name, rows) in enumerate(pop_groups.items()):
        sample = data.iloc[rows]
        per_pop[name] = {
            'mothers': len(sample), 'assigned': int(sample.Genotyped_number.sum()),
            'selfed': int(sample.Selfed_seed_number.sum()),
            'OR_per_doubling': float(np.exp(group_beta[offset + i])),
            'ci95_mother_HC3': group_interval[offset + i].tolist(),
            'male_flower_count_q25_q75': sample.Male_flower_total.quantile([.25, .75]).tolist(),
        }
    return {
        'model': 'Grouped binomial logit of selfed among assigned mature seeds, log2 male flowers + log2 female flowers + log2 end aboveground biomass + population intercepts',
        'mothers': len(data), 'assigned_seeds': int(n.sum()), 'selfed_seeds': int(k.sum()),
        'OR_per_doubling': float(np.exp(beta[offset])),
        'ci95_mother_HC3': interval[offset].tolist(),
        'ci95_within_array_mother_bootstrap': np.quantile(boot[:, 0], [.025, .975]).tolist(),
        'standardized_100_200_male_flowers': {
            'probabilities': probabilities, 'difference': probabilities[1] - probabilities[0],
            'ci95_difference_within_array_mother_bootstrap': np.quantile(boot[:, 1], [.025, .975]).tolist(),
            'standardization': 'each measured mother equal weight; observed female flowers, biomass and array retained',
        },
        'by_population_varying_male_slope_shared_adjusters': per_pop,
        'bootstrap_draws': DRAWS, 'bootstrap_seed': SEED,
    }


def analyze(frames):
    field = frames['3G.csv']
    if field.duplicated(['Population', 'ID']).any():
        raise ValueError('Nonunique field unit')
    have_seeds = field.Genotyped_number.gt(0)
    if not np.all(field.Selfing_rate.notna() == have_seeds):
        raise ValueError('Family denominator and selfing-rate availability disagree')
    assayed = field.loc[have_seeds]
    rates = assayed.Selfed_seed_number / assayed.Genotyped_number
    if not np.allclose(rates, assayed.Selfing_rate, rtol=0, atol=1e-8):
        raise ValueError('Reported selfing fractions disagree with counts')
    known_size = field.loc[field.Biomass_whole.notna()]
    eligible = assayed.loc[assayed.Biomass_whole.notna()]
    dual = eligible.loc[eligible.Male_flower_total.gt(0) & eligible.Female_flower_total.gt(0)]
    scen = frames['3G_ID01.csv']
    if scen.duplicated(['Population', 'ID', 'Inbreeding']).any():
        raise ValueError('Duplicate scenario record')
    if set(scen.Inbreeding) != {'ID=0','ID=1'}:
        raise ValueError('Unexpected fitness scenarios')
    scenario_rows = scen.groupby(['Population', 'ID']).Inbreeding.nunique()
    if len(scenario_rows) != len(field) or not scenario_rows.eq(2).all():
        raise ValueError('Expected paired scenarios on the same plants')
    return {
        'reproduction': {
            'parents_all': len(field), 'parents_by_population': field.groupby('Population').size().tolist(),
            'whole_biomass_observed': len(known_size),
            'average_female_and_male_flowers_180_size_known': known_size[['Female_flower_total','Male_flower_total']].mean().tolist(),
            'mothers_with_assigned_mature_seeds': len(assayed),
            'genotyped_number_sum_equals_paper_assigned': int(assayed.Genotyped_number.sum()),
            'all_assigned_selfed_calls': int(assayed.Selfed_seed_number.sum()),
            'source_model_mothers_with_size_and_progeny': len(eligible),
            'source_model_assigned_and_selfed': [int(eligible.Genotyped_number.sum()),int(eligible.Selfed_seed_number.sum())],
            'source_model_mean_of_maternal_fractions': float((eligible.Selfed_seed_number/eligible.Genotyped_number).mean()),
            'source_model_mothers_with_zero_male_count': int(eligible.Male_flower_total.eq(0).sum()),
            'source_model_mothers_with_zero_female_count': int(eligible.Female_flower_total.eq(0).sum()),
            'scenario_rows': len(scen), 'unique_paired_plants': len(scenario_rows),
            **genotype_summary(frames['Genotype_3P.csv']),
        },
        'within_flowering_mothers': estimate(dual.reset_index(drop=True)),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    frames, hashes = source(args.archive)
    result = {'concept_doi': '10.5281/zenodo.14542260', 'version_doi': VERSION,
              'source_url': URL, 'archive_sha256': SHA, 'archive_zenodo_md5': MD5,
              'member_sha256': hashes, **analyze(frames)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
