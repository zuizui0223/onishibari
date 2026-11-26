#!/usr/bin/env python3
"""Analyze two questions in published Pulsatilla alpina single-flower families.

1. Within intact, typed bisexual flowers, how do organ counts and opening date
   relate to directly observed selfing among the typed seeds?
2. How do the retained bisexual flower groups compare after stamen removal?

Download three exact Dryad files from official Zenodo preservation records and
the Current Biology archive. The output contains only aggregates. Requires
Python 3.12, numpy 2.3.5, scipy 1.17.0 and pandas 2.2.3.

python scripts/pulsatilla_selfing_refit_v51.py --cache-dir PATH --output result.json
"""
import argparse
from hashlib import sha256
import io
import json
from pathlib import Path
import tempfile
import time
from urllib.request import urlopen
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy.special import expit


# The CSV hashes match Dryad's published SHA-256 checksums.
FILES = {
    'census': ('10121671','population.csv','f03cfd1f7f7d4c829dc21691298fa7cbb0930f3719cf5cadeb4b8fc6bd85cbc2'),
    'selfing': ('12748590','Data_selfing.csv','6c6d15b649e157930242f73075a6c77fd93192aa2aba527f320f7a8176f1272e'),
    'phenotypes': ('12748590','Data_phenotypic_selection.csv','e49fef29cefa6e1ed419500a0e312a49524f04876cca5eb21dda14657fe11a98'),
    'fitness': ('15183584','Fitness_landscape_Zenodo.zip','2cb1c78efe86128255dabed72d01799cfe71e10208386f55319f708a8535c141'),
}
CB_MEMBER = 'Fitness_landscape_Zenodo/Data_RS.csv'
CB_SHA256 = '49b6c5147aadf0d3c9de5bacda59b1cd8f25df0d49f521cea5d577668d6a27ba'
ARMS = ['intact','half','all']
ARM_MAP = {'N':'intact','None':'intact','Intact':'intact','Half':'half','SR50':'half','All':'all','SR100':'all'}
TERMS = ['50_more_stamens','100_more_pistils','7_days_later']
SEED = 20261122


def download(cache, label):
    record,name,digest = FILES[label]
    path = cache / name
    if path.exists():
        if sha256(path.read_bytes()).hexdigest()!=digest:
            raise ValueError(f'Checksum mismatch in cached file {name}')
        return path
    url = f'https://zenodo.org/records/{record}/files/{name}?download=1'
    for attempt in range(3):
        try:
            with urlopen(url, timeout=90) as response:
                content = response.read()
            if sha256(content).hexdigest()!=digest:
                raise ValueError(f'Downloaded checksum mismatch for {name}')
            with tempfile.NamedTemporaryFile(dir=cache,delete=False) as temp:
                temp.write(content)
                tmp = Path(temp.name)
            tmp.replace(path)
            return path
        except OSError:
            if attempt==2:
                raise
            time.sleep(1+attempt)


def read_frame(content):
    # 'None' is the name of the intact arm in Current Biology.
    df = pd.read_csv(content,keep_default_na=False,na_values=['NA',''],dtype={'Ind_ID':str})
    df['arm'] = df.Stamen_removal.map(ARM_MAP)
    if df.arm.isna().any():
        raise ValueError('Unrecognized treatment name')
    return df


def interval(point, draws):
    return {'estimate':float(point),'flower_percentile_95':np.quantile(draws,[.025,.975]).tolist()}


def logit(X, typed, selfed):
    # Grouped binomial score. Each bootstrap draw resamples whole flowers.
    b = np.zeros(X.shape[1])
    for _ in range(80):
        p = expit(X@b)
        information = X.T @ ((typed*p*(1-p))[:,None]*X)
        step = np.linalg.solve(information,X.T@(selfed-typed*p))
        b += step
        if np.max(np.abs(step))<1e-9:
            return b
    raise RuntimeError('Logistic regression did not converge')


def intact_models(t, reps, rng):
    g = t[t.arm=='intact'].sort_values('Ind_ID')
    X = np.column_stack([np.ones(len(g)),(g.Stamen_number_func-180)/50,
                         (g.Pistil_number-220)/100,(g.Flowering_date-160)/7])
    typed,selfed = (g[col].to_numpy() for col in ['Total_typed_seed','Total_selfed_seed'])
    def estimates(ind):
        return np.exp(logit(X[ind],typed[ind],selfed[ind])[1:])
    b = estimates(np.arange(len(g)))
    boot = np.stack([estimates(rng.integers(len(g),size=len(g))) for _ in range(reps)])
    return {
        'flowers':len(g),'typed_seeds':int(typed.sum()),'selfed_typed_seeds':int(selfed.sum()),
        'organ_date_support':{k:[float(g[k].min()),float(g[k].max())] for k in ['Stamen_number_func','Pistil_number','Flowering_date']},
        'selfing_odds_ratio':{v:interval(b[i],boot[:,i]) for i,v in enumerate(TERMS)},
    }


def by_arm_means(g, cols, reps, rng):
    draws={}
    means={}
    for arm in ARMS:
        rows = g[g.arm==arm].sort_values('Ind_ID')
        values = rows[cols].to_numpy()
        means[arm] = values.mean(axis=0)
        ix = rng.integers(len(rows), size=(reps,len(rows)))
        draws[arm] = values[ix].mean(axis=1)
    return {
        'n':{arm:int((g.arm==arm).sum()) for arm in ARMS},
        'group_means':{arm:dict(zip(cols,means[arm].tolist())) for arm in ARMS},
        'all_minus_intact':{col:interval(means['all'][i]-means['intact'][i],draws['all'][:,i]-draws['intact'][:,i]) for i,col in enumerate(cols)},
    }


def audit_and_analyze(cache, reps):
    cache.mkdir(parents=True,exist_ok=True)
    census,typed_file,phen = (read_frame(download(cache,k)) for k in ['census','selfing','phenotypes'])
    with ZipFile(download(cache,'fitness')) as z:
        data = z.read(CB_MEMBER)
    assert sha256(data).hexdigest()==CB_SHA256
    cb = read_frame(io.BytesIO(data))
    p = census[census.Flower_number==1]
    assert len(p)==104 and len(cb)==87 and p.Ind_ID.is_unique and cb.Ind_ID.is_unique
    matched = p.merge(cb,on='Ind_ID',suffixes=('_p','_c'),validate='one_to_one')
    assert (matched.arm_p==matched.arm_c).all() and (matched.Floral_gender_p==matched.Floral_gender_c).all()
    # Count the published arm totals, including two outcome rows with identifiers
    # unmatched in the census. Do not assign an inferred identity to either.
    def counts(g,sex):
        return {a:int(((g.arm==a)&(g.Floral_gender==sex)).sum()) for a in ARMS}
    cohort={sex:{'census_104':counts(p,sex),'retained_87':counts(cb,sex)} for sex in ['H','M']}
    for sub in cohort.values():
        sub['census_minus_retained']={a:sub['census_104'][a]-sub['retained_87'][a] for a in ARMS}

    t = phen[phen.Inbreeding_depression==0].copy()
    assert len(t)==60 and len(typed_file)==60 and t.Ind_ID.is_unique and typed_file.Ind_ID.is_unique
    # One label differs between the two 2024 CSV files. Verify the complete
    # one-to-one match on observed organ, treatment, date and count fields.
    keys = ['arm','Pistil_number','Stamen_number_func','Flowering_date','Total_typed_seed','Total_selfed_seed','Petal_length','Stalk_height']
    aligned = typed_file.merge(t,on=keys,validate='one_to_one',suffixes=('_s','_p'))
    assert len(aligned)==60
    assert np.max(abs(aligned.Selfing_rate-aligned.Total_selfed_seed/aligned.Total_typed_seed)) < 1e-8
    assert (t.Total_typed_seed>0).all() and (t.Total_typed_seed>=t.Total_selfed_seed).all()
    t['typed_self_fraction'] = t.Total_selfed_seed/t.Total_typed_seed
    t['outcrossed_harvest'] = t.Mature_seed_number*(1-t.typed_self_fraction)
    t['discounted_harvest'] = t.Mature_seed_number*(1-.95*t.typed_self_fraction)
    bisexual = cb[cb.Floral_gender=='H'].copy()
    assert len(bisexual)==66 and (bisexual.Pistil_number>0).all() and (cb.loc[cb.Floral_gender=='M','Pistil_number']==0).all()
    overlap = t.merge(bisexual,on='Ind_ID',suffixes=('_j','_c'),validate='one_to_one')
    assert len(overlap)==60
    assert (overlap.arm_j==overlap.arm_c).all()
    for col in ['Pistil_number','Stamen_number_func']:
        assert np.array_equal(overlap[col+'_j'],overlap[col+'_c'])
    assert np.array_equal(overlap.F_RS_ID0,overlap.Mature_seed_number)
    # Exact decomposition of the archived endpoints, not a paternity observation.
    bisexual['implied_selfed']=(bisexual.F_RS_ID0-bisexual.F_RS_ID95)/.95
    bisexual['implied_outcrossed']=bisexual.F_RS_ID0-bisexual.implied_selfed
    assert np.allclose(bisexual.implied_outcrossed+.05*bisexual.implied_selfed,bisexual.F_RS_ID95)
    abs_disagreement = abs(overlap.discounted_harvest-overlap.F_RS_ID95)
    implied_fraction = (overlap.F_RS_ID0-overlap.F_RS_ID95)/(.95*overlap.F_RS_ID0)
    authors_t = t[t.Total_typed_seed>=6]
    rng = np.random.default_rng(SEED)
    return {
        'design':{'population_seasons':1,'census_flowers':len(census),'census_flowering_plants':census.Ind_ID.nunique(),
                  'single_flowered_census':len(p),'single_flowered_outcomes':len(cb),'exact_census_outcome_id_matches':len(matched),
                  'unmatched_retained_by_arm':counts(cb[~cb.Ind_ID.isin(p.Ind_ID)],'H'),'by_flower_gender':cohort,
                  'phenotype_file_families_linked_to_cb_by_id_and_traits':len(overlap),
                  'typed_vs_phenotype_exact_id_matches':int((aligned.Ind_ID_s==aligned.Ind_ID_p).sum()),
                  'typed_vs_phenotype_unique_observation_matches':len(aligned)},
        'seed':SEED,'bootstrap_draws':reps,
        'typed_seed_counts':{a:{'flowers':int((t.arm==a).sum()),'selfed':int(t.loc[t.arm==a,'Total_selfed_seed'].sum()),
                                'typed':int(t.loc[t.arm==a,'Total_typed_seed'].sum())} for a in ARMS},
        'authors_six_seed_threshold':{'flowers':len(authors_t),'typed':int(authors_t.Total_typed_seed.sum()),'selfed':int(authors_t.Total_selfed_seed.sum())},
        'intact_adjusted_associations':intact_models(t,reps,rng),
        'typed_arm_contrast_jecol_60':by_arm_means(t,['Mature_seed_number','outcrossed_harvest','discounted_harvest'],reps,rng),
        'decomposition_audit':{'cb_n':len(bisexual),'cb_implied_outcrossed_below_zero':int((bisexual.implied_outcrossed< -1e-8).sum()),
                               'cb_implied_outcrossed_min':float(bisexual.implied_outcrossed.min()),
                               'cb_implied_selfed_below_zero':int((bisexual.implied_selfed< -1e-8).sum()),
                               'cb_all_removed_positive_implied_selfed':int((bisexual.loc[bisexual.arm=='all','implied_selfed']>1e-8).sum()),
                               'cb_vs_typed_rate_absolute_difference_median':float(abs(implied_fraction-overlap.typed_self_fraction).median()),
                               'cb_vs_typed_rate_absolute_difference_max':float(abs(implied_fraction-overlap.typed_self_fraction).max()),
                               'cb_vs_typed_weighted_absolute_difference_median':float(abs_disagreement.median()),
                               'cb_vs_typed_weighted_absolute_difference_max':float(abs_disagreement.max())},
        'source_files':{k:{'zenodo_record':r,'file':f,'sha256':h} for k,(r,f,h) in FILES.items()},
        'cb_member':{'file':CB_MEMBER,'sha256':CB_SHA256},
    }


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replicates',type=int,default=10000)
    args=parser.parse_args()
    if args.replicates<1:
        parser.error('--replicates must be positive')
    result=audit_and_analyze(args.cache_dir,args.replicates)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
