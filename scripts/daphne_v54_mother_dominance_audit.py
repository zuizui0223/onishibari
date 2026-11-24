#!/usr/bin/env python3
"""v0.54 supplementary audit: the mother-unit distribution behind pair counts.
No raw mother/father identifiers or source rows printed or persisted.
"""
import json, numpy as np
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data
from daphne_jezoensis_2018_2020_flower_to_paternity_choice_v30 import build
f=data()
for year in (2018,2020):
    q=build(year,f)
    a=q["counts"].astype(int)
    m=a.sum(axis=1)
    pairs=m*(m-1)//2
    same=np.sum(a*(a-1)//2,axis=1)
    ratio=np.divide(same,pairs,out=np.full(len(same),np.nan),where=pairs>0)
    eligible=pairs>0
    weighted=same.sum()/pairs.sum()
    equal=np.nanmean(ratio)
    leaveout=[float((same.sum()-same[i])/(pairs.sum()-pairs[i])) for i in range(len(m))
              if pairs.sum()>pairs[i]]
    top=np.argsort(pairs)[::-1]
    top_same=np.argsort(same)[::-1]
    # Original input CSV columns and duplicate pattern: report only column headings.
    name=(f"NF{year}_2_primer16_res1.2.3_var1_150000_220702.csv")
    raw=f[name]
    c=[k for k in raw.columns if "offspring" in k.lower() or
       "seed" in k.lower() or "sample" in k.lower()]
    result={
      "year":year,"n_mothers":len(m),"offspring":int(m.sum()),
      "source_original_assigned_file_columns":list(map(str,raw.columns)),
      "source_offspring_id_columns":c,
      "total_pairs":int(pairs.sum()),"total_same_father_pairs":int(same.sum()),
      "pair_weighted_fraction":float(weighted),
      "mother_equal_fraction":float(equal),
      "mother_equal_median":float(np.nanmedian(ratio)),
      "mothers_with_two_or_more_offspring":int(eligible.sum()),
      "mothers_with_five_or_more_offspring":int((m>=5).sum()),
      "mother_seed_count_quantiles":[float(v) for v in np.quantile(m,[0,.25,.5,.75,.9,1])],
      "largest_mother_family_size":int(m.max()),
      "top1_mother_fraction_of_all_within_mother_pairs":float(pairs[top[0]]/pairs.sum()),
      "top3_mothers_fraction_of_all_within_mother_pairs":float(pairs[top[:3]].sum()/pairs.sum()),
      "top1_mother_fraction_of_all_same_father_pairs":float(same[top_same[0]]/same.sum()),
      "top3_mothers_fraction_of_all_same_father_pairs":float(same[top_same[:3]].sum()/same.sum()),
      "leave_one_mother_out_pooled_ratio_range":[float(min(leaveout)),float(max(leaveout))],
      "remove_largest_sampled_mother_pooled_ratio":float((same.sum()-same[np.argmax(m)])/(pairs.sum()-pairs[np.argmax(m)])),
      "share_of_mothers_with_all_offspring_same_father_if_n_ge_2":float((ratio[eligible]==1).mean()),
    }
    print("DAPHNE_MATERNAL_CLUSTER_AUDIT_V54",json.dumps(result,ensure_ascii=False,allow_nan=False))
