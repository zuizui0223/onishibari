"""Reanalyze male allocation in Aljiboury & Friedman (2022); output aggregates only.

python scripts/ambrosia_allocation_refit_v50.py --input-dir /private/data --output results.json
If --input-dir is omitted, work in memory with the public version archive.
Python dependencies: numpy, scipy, pandas. No individual rows are written.
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
from scipy import stats

VERSION_URL = 'https://datadryad.org/api/v2/versions/178493/download'
FIELD = 'AljibouryFriedman_field_data.csv'
PROGENY = 'AljibouryFriedman_progeny_data.csv'
DIGESTS = {
    FIELD: '9c703da516bed53363d1678d2fb95d228d11956790fd03d6777e5b3ef4631aa2',
    PROGENY: 'a9612ecbb6a6f57acdd4b96462b030f8295fc18170c2d3fd2d5cb783bdc5bbdd',
    'README_file.txt': '9530359c9617f22468ac331c34f5728c9e1f568fcd4c169d58696bac5ba6c0ff',
}


def inputs(folder):
    if folder:
        blobs = {name: (folder / name).read_bytes() for name in DIGESTS}
    else:
        req = urllib.request.Request(VERSION_URL, headers={
            'User-Agent': 'Published-data-reproduction/0.1', 'Accept': 'application/zip'})
        with urllib.request.urlopen(req, timeout=90) as response:
            with zipfile.ZipFile(io.BytesIO(response.read())) as archive:
                blobs = {name: archive.read(name) for name in DIGESTS}
    for name, blob in blobs.items():
        if hashlib.sha256(blob).hexdigest() != DIGESTS[name]:
            raise ValueError(f'Unexpected SHA-256 for {name}')
    return pd.read_csv(io.BytesIO(blobs[FIELD])), pd.read_csv(io.BytesIO(blobs[PROGENY]))


def contrast(y, category):
    lower = y[category.eq(1)].dropna()
    upper = y[category.eq(4)].dropna()
    return {
        'donors': int(y.count()), 'positive': int(y.gt(0).sum()), 'zeros': int(y.eq(0).sum()),
        'sum': int(y.sum()),
        'category_means': y.groupby(category).mean().tolist(),
        'category_n': y.groupby(category).count().tolist(),
        'high_minus_low': float(upper.mean() - lower.mean()),
        'high_over_low': float(upper.mean() / lower.mean()),
    }


def mass_design(data, within_cap):
    factors = ['block', 'edge']
    if within_cap:
        factors.append('cat')
    dummies = pd.get_dummies(data[factors].astype(str), drop_first=True)
    return np.column_stack([np.ones(len(data)), np.log(data.mf_wt), dummies])


def elasticity(X, y, covariance=False):
    # Poisson pseudo-likelihood. Zeros contribute directly to the count mean.
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    for _ in range(60):
        mu = np.exp(X @ beta)
        hessian = X.T @ (mu[:, None] * X)
        step = np.linalg.solve(hessian, X.T @ (y - mu))
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    else:
        raise ArithmeticError('Poisson mean fit did not converge')
    if not covariance:
        return float(beta[1])
    mu = np.exp(X @ beta)
    bread = np.linalg.inv(X.T @ (mu[:, None] * X))
    leverage = mu * np.einsum('ij,jk,ik->i', X, bread, X)
    scores = X * ((y - mu) / (1 - leverage))[:, None]
    variance = bread @ (scores.T @ scores) @ bread
    se = np.sqrt(variance[1, 1])
    return {
        'beta_log_mass': float(beta[1]), 'se_donor_hc3': float(se),
        'ci95_donor_hc3': (beta[1] + np.array([-1, 1]) * stats.norm.ppf(.975) * se).tolist(),
        'doubling_ratio': float(2**beta[1]),
    }


def mother_bootstrap(field, progeny, donors, designs):
    # Keep the observed donor plants fixed and resample whole maternal families.
    # Respect maternal selection by block and allocation category (all are interior).
    counts = pd.crosstab(progeny.Mother_ID, progeny.Candidate_father_ID)
    counts = counts.reindex(columns=donors.index, fill_value=0)
    if not np.array_equal(counts.sum(axis=0).to_numpy(), donors.sire_total.to_numpy()):
        raise ValueError('Progeny rows do not reconstruct sampled donor counts')
    mother = field.loc[counts.index]
    family = counts.to_numpy()
    rng = np.random.default_rng(14544)
    repetitions = 1999
    by_stratum = mother.groupby(['block', 'cat']).indices.values()
    weight = np.zeros((repetitions, len(mother)), dtype=int)
    for rows in by_stratum:
        size = len(rows)
        weight[:, rows] = rng.multinomial(size, np.full(size, 1 / size), size=repetitions)
    totals = weight @ family
    return {
        name: np.quantile([elasticity(X, y) for y in totals], [.025, .975]).tolist()
        for name, X in designs.items()
    }


def analyze(field, progeny):
    field = field.set_index('id', verify_integrity=True)
    if not progeny.Offspring_ID.is_unique:
        raise ValueError('Duplicate progeny record')
    for column in ['Mother_ID', 'Candidate_father_ID']:
        if not progeny[column].isin(field.index).all():
            raise ValueError(f'Unlinked IDs in {column}')
    archive_y = progeny.Candidate_father_ID.value_counts().reindex(field.index, fill_value=0)
    maternal = progeny.Mother_ID.value_counts().reindex(field.index, fill_value=0)
    focal = field.seeds_geno.notna()
    if not np.array_equal(field.loc[focal,'seeds_geno'], maternal[focal]) or maternal[~focal].any():
        raise ValueError('Maternal tally mismatch')
    primary = field.dropna(subset=['sire_total'])
    known = field.sire_total.notna()
    if not np.array_equal(field.sire_total[known], archive_y[known]):
        raise ValueError('Observed sire tallies disagree')
    documented = field.sire_total.fillna(archive_y.where(archive_y.gt(0)))
    definitions = {
        'field_observed': field.sire_total,
        'observed_plus_documented_assignments': documented,
        'appearances_as_sire_in_progeny_file': archive_y,
    }
    fits = {}
    for name, y in definitions.items():
        data = field.loc[y.notna()]
        fits[name] = {
            'pooled': elasticity(mass_design(data, False), y.dropna().to_numpy(), True),
            'within_cap': elasticity(mass_design(data, True), y.dropna().to_numpy(), True),
        }
    designs = {name: mass_design(primary, within) for name, within in [('pooled',False),('within_cap',True)]}
    boot = mother_bootstrap(field, progeny, primary, designs)
    for name, interval in boot.items():
        fits['field_observed'][name]['ci95_fixed_donors_maternal_family_resampling'] = interval
    missing = ~known
    return {
        'status': 'original_public_data_secondary_refit_v50',
        'article_doi': '10.1111/evo.14544',
        'data_doi': '10.5061/dryad.4tmpg4fcp',
        'source_version_zip': VERSION_URL, 'sha256': DIGESTS,
        'interpretation': {
            'allocation_contrasts': 'Descriptive treatment-category means among archived outcomes; not a complete-data intention-to-treat estimator.',
            'mass_models': 'Observational Poisson pseudo-likelihood log-mean models, not the published negative-binomial mixed models or published power fit.',
            'mass_endpoint': 'Remaining male inflorescence dry mass at senescence; not cumulative seasonal construction or maintenance cost.',
            'donor_hc3': 'Working donor-level independence, not a full covariance model for the shared mating network.',
            'maternal_resampling': 'Conditional maternal-family sampling variation with the measured donor set fixed; not total biological uncertainty.',
            'missing_outcomes': 'Three donors remain unclassified after recovery of two positive assignments for one missing donor. Archive appearances are not biological zero success.',
        },
        'audit': {
            'field_rows': len(field), 'assigned_category_n': field.groupby('cat').size().tolist(),
            'block_x_category_counts': pd.crosstab(field.block, field.cat).to_numpy().tolist(),
            'source_populations': int(field['pop'].nunique()),
            'unique_progeny': int(progeny.Offspring_ID.nunique()),
            'distinct_archive_fathers': int(progeny.Candidate_father_ID.nunique()),
            'focal_mothers': int(focal.sum()), 'mothers_per_block': field[focal].groupby('block').size().tolist(),
            'maternal_assigned_count_range': [int(maternal[focal].min()), int(maternal[focal].max())],
            'maternal_assignment_total': int(maternal.sum()), 'maternal_sums_match': True,
            'observed_field_donor_totals_match_progeny': True,
            'field_missing_donor_totals': int(missing.sum()),
            'field_missing_with_positive_documentation': int((missing & archive_y.gt(0)).sum()),
            'documented_assignments_in_missing': int(archive_y[missing].sum()),
            'mass_positive_and_complete': bool(field.mf_wt.gt(0).all()),
            'mean_mass_all_assigned_by_cat': field.groupby('cat').mf_wt.mean().tolist(),
        },
        'allocation_contrasts': {name: contrast(y, field.cat) for name, y in definitions.items()},
        'observational_mass_models': fits,
        'maternal_bootstrap': {'draws':1999,'seed':14544,'strata':'block by maternal allocation category',
            'population':'resample measured maternal families conditional on fixed field-observed donor set'},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = analyze(*inputs(args.input_dir))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, allow_nan=False, indent=2) + '\n', encoding='utf-8')
