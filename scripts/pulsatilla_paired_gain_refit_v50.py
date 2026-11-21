#!/usr/bin/env python3
"""Reproduce and compare seasonal organ-count gain curves in Pulsatilla alpina.

Run with Python 3, NumPy, pandas, SciPy: python scripts/pulsatilla_paired_gain_refit_v50.py
Archives are downloaded from versioned public records if absent. Output is
aggregate JSON only. The resampling unit is the focal plant. Intervals condition on this population,
its observed mating network and deposited reproductive estimates; they do not
propagate paternity-assignment or inbreeding-depression estimation uncertainty.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import tempfile
from zipfile import ZipFile

import numpy as np
import pandas as pd

from scipy.optimize import minimize_scalar
from scipy.stats import t


def fit_power(allocation, success):
    x = np.asarray(allocation, dtype=float)
    y = np.asarray(success, dtype=float)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    assert np.all(x >= 0) and np.all(y >= 0)
    positive = x > 0
    # The fitted prediction at x=0 is zero for b>0. Its residual (even if
    # nonzero) is constant in the optimization. Use every row in residual SE.
    log_x = np.log(x[positive] / np.median(x[positive]))
    obs = y[positive]
    min_b = 1e-8 if not positive.all() else -12.0
    def loss(b):
        h = np.exp(b * log_x)
        residual = obs - h * (np.dot(h, obs) / np.dot(h, h))
        return np.dot(residual, residual)
    fit = minimize_scalar(loss, bounds=(min_b,12.), method='bounded', options={'xatol':1e-9})
    assert fit.success
    b = fit.x
    h = np.exp(b * log_x)
    a = np.dot(h, obs) / np.dot(h, h)
    derivative = np.column_stack([h, a * h * log_x])
    rss = fit.fun + np.dot(y[~positive], y[~positive])
    se = np.sqrt(np.linalg.inv(derivative.T @ derivative)[1,1] * rss / (len(x)-2))
    return float(b), float(se), float(2*t.sf(abs(b-1)/se,len(x)-2))


ARCHIVES = {
    'amnat': ('15549445', 'AmNat_R3_Zenodo.zip',
              'c31152fb988d5602113c91628806691519e3700351df3b79a403822ca942a463'),
    'cub': ('15183584', 'Fitness_landscape_Zenodo.zip',
            '2cb1c78efe86128255dabed72d01799cfe71e10208386f55319f708a8535c141'),
}


def get_archive(label, data_dir):
    record, filename, checksum = ARCHIVES[label]
    data_dir.mkdir(parents=True, exist_ok=True)
    path = data_dir / filename
    if not path.exists():
        url = f'https://zenodo.org/records/{record}/files/{filename}?download=1'
        payload = urllib.request.urlopen(url, timeout=180).read()
        assert hashlib.sha256(payload).hexdigest() == checksum
        path.write_bytes(payload)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum
    return ZipFile(path)


def table(archive, filename, hashes):
    data = archive.read(filename)
    hashes[filename] = hashlib.sha256(data).hexdigest()
    return pd.read_csv(io.BytesIO(data), keep_default_na=False, na_values=['NA',''])


def curves_for_outcomes(organ_counts, outcomes, replications, generator):
    # All outcomes in a dataset use the same index vector on each draw.
    names = list(outcomes)
    arrays = [(organ_counts[sex].to_numpy(), y.to_numpy()) for sex,y in outcomes.values()]
    n = len(arrays[0][0])
    actual = {name:fit_power(x, y) for name,(x,y) in zip(names,arrays)}
    boot = np.empty((replications,len(arrays)))
    for i in range(replications):
        index = generator.integers(n,size=n)
        for j,(x,y) in enumerate(arrays):
            boot[i,j] = fit_power(x[index],y[index])[0]
    boot_table = pd.DataFrame(boot, columns=names)
    estimates = {name:x[0] for name,x in actual.items()}
    point_and_ci = {}
    for condition in sorted(set(name.rsplit('_',1)[1] for name in names)):
        name = 'female_minus_male_'+condition
        f,m = 'female_'+condition, 'male_'+condition
        boot_table[name] = boot_table[f]-boot_table[m]
        estimates[name] = estimates[f]-estimates[m]
    if 'female_mature' in names and 'female_discounted' in names:
        for component in ['female','male','female_minus_male']:
            name = component+'_discounted_minus_mature'
            f,m = component+'_discounted',component+'_mature'
            boot_table[name] = boot_table[f]-boot_table[m]
            estimates[name] = estimates[f]-estimates[m]
    for name,estimate in estimates.items():
        ci = list(boot_table[name].quantile([.025,.975]))
        point_and_ci[name] = {'b_or_contrast':estimate,'paired_percentile_95':ci}
        if name in actual:
            point_and_ci[name]['nls_se'] = actual[name][1]
            x = arrays[names.index(name)][0]
            lower = 0.0 if np.any(x==0) else -12.0
            point_and_ci[name]['draws_near_lower_bound'] = int((boot_table[name]<=lower+1e-5).sum())
            point_and_ci[name]['draws_near_upper_bound'] = int((boot_table[name]>=12-1e-5).sum())
    return {'n':n, 'replicates':replications, 'results':point_and_ci}


def run(data_dir, reps, seed):
    files = {}
    with get_archive('amnat',data_dir) as z:
        d = table(z,'AmNat_R3_Zenodo/data_Ind_fitness_SexFunction.csv',files)
        rm = 'AmNat_R3_Zenodo/code_Fig.1 Gain curve_revision.R'
        files[rm] = hashlib.sha256(z.read(rm)).hexdigest()
    with get_archive('cub',data_dir) as z:
        c = table(z,'Fitness_landscape_Zenodo/Data_RS.csv',files).set_index('Ind_ID')
    f = d[d.Sex_function=='female function'].set_index('Ind_ID').sort_index()
    m = d[d.Sex_function=='male function'].set_index('Ind_ID').sort_index()
    assert f.index.is_unique and f.index.equals(m.index)
    predicate = (f.Mani=='N') & f.Biomass.notna() & f.Flowering_date.notna()
    assert len(f)==135 and predicate.sum()==89
    variables = ['Mani','Biomass','Flowering_date','Ind_Pistil_number','Ind_Stamen_number_func','Flower_number']
    assert f[variables].equals(m[variables])
    raw_n = len(f)
    excluded = {'manipulated':int((f.Mani!='N').sum()),
                'unmanipulated_incomplete':int(((f.Mani=='N')&~predicate).sum())}
    cub_pop_overlap = len(c.index.intersection(f.index))
    f,m = f.loc[predicate],m.loc[predicate]
    pistils,stamens,fitness = 'Ind_Pistil_number','Ind_Stamen_number_func','RS_AD_Ind_ID95_r'
    data89 = {'female':f[pistils], 'male':m[stamens]}
    assert f[fitness].notna().all() and m[fitness].notna().all()
    assert f.loc[f[pistils]==0,fitness].eq(0).all()
    common = f.index.intersection(c.index)
    assert len(common)==59 and (c.loc[common,'Stamen_removal']=='None').all()
    assert c.loc[common,'Pistil_number'].eq(f.loc[common,pistils]).all()
    assert c.loc[common,'Stamen_number_func'].eq(m.loc[common,stamens]).all()
    c59 = c.loc[common]
    assert (f.loc[common,fitness]-c59.F_RS_ID95).abs().max()<=0.5
    assert (m.loc[common,fitness]-c59.M_RS_ID95).abs().max()<=0.5
    rng = np.random.default_rng(seed)
    out = {'status':'original_public_data_secondary_refit_v50',
           'article_dois':['10.1086/737646','10.1016/j.cub.2025.04.022'],
           'archive_provenance':{name:{'record':record,'filename':file,'sha256':sha}
                                 for name,(record,file,sha) in ARCHIVES.items()},
           'interpretation':{
               'resampling_unit':'Same focal-plant indices for all sex functions and outcome definitions in each draw.',
               'uncertainty':'Conditional on the observed mating network, estimated paternity, rounded primary outcomes and fixed inbreeding depression delta=0.95.',
               'sampling_scope':'One population in 2022; companion paper is overlapping, not independent replication.',
               'investment':'Pistil and calibrated stamen counts; not individually measured floral dry mass.',
               'success':'Archived seasonal reproductive-success estimates; not observed survival or lifetime fitness.',
               'covariates':'Univariate organ-count curves; no causal interpretation or adjustment for plant size and flowering time.'},
           'archive_members_sha256':files,'population_n':raw_n,'exclusions':excluded,
           'seed':seed,'n89_allocation':{sex:{'positive':int((v>0).sum()),'min_positive':float(v[v>0].min()),'max':float(v.max())} for sex,v in data89.items()},
           'n89_zero_success':dict(female=int(f[fitness].eq(0).sum()),male=int(m[fitness].eq(0).sum())),
           'cub_n':len(c),'cub_population_overlap':cub_pop_overlap,'cub_analysis_overlap':len(common),
           'n59_positive_pistils':int((c59.Pistil_number>0).sum()),
           'n89_flower_numbers':{int(k):int(v) for k,v in f.Flower_number.value_counts().items()},
           'n89_selfing_missing_with_nonzero_female_fitness':int((f.Selfing_rate.isna()&(f[fitness]>0)).sum())}
    out['primary_n89'] = curves_for_outcomes(data89, {'female_discounted':('female',f[fitness]),
                                            'male_discounted':('male',m[fitness])}, reps, rng)
    out['outcome_sensitivity_shared_n59'] = curves_for_outcomes(
        {'female':c59.Pistil_number,'male':c59.Stamen_number_func},
        {'female_discounted':('female',c59.F_RS_ID95),'male_discounted':('male',c59.M_RS_ID95),
         'female_mature':('female',c59.F_RS_ID0),'male_mature':('male',c59.M_RS_ID0)}, reps,rng)
    # Whole companion trial includes removals. Point estimates verify source coding only.
    out['companion_n87_replication'] = {f'{sex}_{d}':fit_power(x,y)[0] for sex,x,base in [
        ('female',c.Pistil_number,'F_RS_ID'),('male',c.Stamen_number_func,'M_RS_ID')]
        for d,y in [(0,c[base+'0']),(95,c[base+'95'])]}
    foc = []
    for i in range(len(f)):
        ix = np.arange(len(f))!=i
        bf = fit_power(f[pistils][ix],f[fitness][ix])[0]
        bm = fit_power(m[stamens][ix],m[fitness][ix])[0]
        foc.append(bf-bm)
    out['n89_leave_one_out_difference_range'] = [min(foc),max(foc)]
    return out


if __name__=='__main__':
    args=argparse.ArgumentParser()
    args.add_argument('--data-dir','--cache-dir',type=Path,default=Path(tempfile.gettempdir())/'onishibari_gain_curves_v50')
    args.add_argument('--replicates',type=int,default=10000)
    args.add_argument('--seed',type=int,default=20261121)
    args.add_argument('--output',type=Path,default=Path('pulsatilla_aggregate.json'))
    opts=args.parse_args()
    if opts.replicates < 100:
        raise ValueError('Use at least 100 paired resamples; final report uses 10000.')
    opts.output.parent.mkdir(parents=True, exist_ok=True)
    opts.output.write_text(json.dumps(run(opts.data_dir,opts.replicates,opts.seed),indent=2,allow_nan=False)+'\n')
