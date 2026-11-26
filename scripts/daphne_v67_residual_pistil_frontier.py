#!/usr/bin/env python3
"""v67: explicit conditional break-even frontier for retained pistils.

Use previously source-verified v55 natural-only paternity assignments (female
mothers, strict H pollen donors, assignment probability >=.95) and v64 natural
fruit observations (strict H). The ratio fruits / assigned sample offspring
is NOT a lifetime fitness ratio. theta = expected maternal gene-copy value per
recorded fruit / effective total external paternal gene-copy value per sampled
assigned seed. The latter includes coverage/representativeness, survival,
other recipient sex classes, and reproductive value. No empirical theta.
For a hypothetical small-effect pistil-loss variant with complete loss of
female function, unchanged mate environment, common fractional male gain g,
and additional future net benefit h in units of current paternal contribution:
    retain-minus-loss / P = theta * R - g - h
For fixed floral budget, c=fraction of per-flower resource cost saved,
eta=fraction of savings allocated to extra flowers at reduced flower cost,
and assumed male fitness elasticity beta:
    extra flowers factor = 1 + eta*c/(1-c)
    g = factor**beta - 1.
This is a conditional accounting model, not an ESS/invasion analysis,
not a causal estimate of beta, c, h, theta, or selection.
"""
import json,hashlib,math
from pathlib import Path

P55=Path("data/2026-11-24_daphne_treatment_specific_sire_v55.json")
P64=Path("data/2026-11-24_daphne_fruit_activation_v64.json")
OUT=Path("analysis_results_daphne_residual_pistil_frontier_v67")
ZIP_SHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
THETA=(0.025,0.05,0.1,0.25,0.5,1.,2.)
ETA=(0.25,0.5,1.)
H=(0.,0.02,0.05,0.1)
C=(0.,0.01,0.02,0.03,0.05,0.1,0.2)
def frontier(r,theta,beta,eta,h=0.):
    target=theta*r-h
    if target<=0:
        return {"required_additional_current_paternal_gain":max(0.,target),
                "break_even_cost_fraction":0.,
                "future_gain_alone_weakly_suffices":True}
    dr=math.expm1(math.log1p(target)/beta)
    return {"required_additional_current_paternal_gain":target,
            "break_even_cost_fraction":dr/(eta+dr),
            "future_gain_alone_weakly_suffices":False}
def evaluate(r,theta,beta,eta,c,h):
    g=math.expm1(beta*math.log1p(eta*c/(1-c)))
    return {"paternal_gain_fraction":g,"retain_minus_loss_over_current_paternity":theta*r-g-h}
def main():
    a=json.loads(P55.read_text());b=json.loads(P64.read_text())
    assert a["zip_sha256"]==b["original_zip_sha256"]==ZIP_SHA
    results=[]
    for y in (2018,2020):
        p=next(x for x in a["years"] if x["year"]==y)
        q=next(x for x in p["by_treatment"] if x["treatment"]=="Natural")
        m=next(x for x in b["sex_year"] if x["year"]==y and x["morph_exact_label"]=="H")
        n,fruits,seeds=int(m["n_original_genets_valid"]),int(m["fruits_total"]),int(q["offspring"])
        assert (n,fruits,seeds)=={2018:(105,53,110),2020:(71,28,131)}[y]
        assert p["strict_H_father_p095_by_treatment"]["Natural"]==seeds
        r=fruits/seeds; beta_obs=float(q["model_beta_10m"])
        scenarios=[]
        for theta in THETA:
            for beta,kind in [(0.25,"hypothetical"),(0.5,"hypothetical"),(1.,"hypothetical"),(2.,"hypothetical"),(beta_obs,"observational_elasticity_sensitivity_only")]:
                for eta in ETA:
                    for h in H:
                        z=frontier(r,theta,beta,eta,h)
                        c=z["break_even_cost_fraction"]
                        if not z["future_gain_alone_weakly_suffices"]:
                            assert abs(evaluate(r,theta,beta,eta,c,h)["retain_minus_loss_over_current_paternity"])<1e-10
                        scenarios.append({"theta":theta,"beta":beta,"beta_type":kind,"eta":eta,"future_gain_h":h,**z})
        results.append({"year":y,"strict_H_flowering_genets_with_valid_fruits":n,
          "flowers":int(m["flowers_total"]),"natural_H_fruits":fruits,
          "sampled_natural_female_mother_seed_assignments_to_H_p095":seeds,
          "source_natural_only_sampled_mother_genets_with_H_assignments":int(q["mothers"]),
          "raw_fruit_to_assigned_seed_count_ratio_not_fitness":r,
          "observational_paternity_choice_flower_elasticity_lambda10":beta_obs,
          "frontier_grid":scenarios,
          "example_theta0_1_beta1_eta1_future0":{str(c):evaluate(r,.1,1.,1.,c,0.) for c in C}})
    # Selfing contribution: outcross descendants inherit one copy of a
    # parental nuclear allele in this gene-copy accounting; selfed inherit
    # two parental gametic copies in expectation under equal Mendelian
    # transmission. Counts are relative to allele frequency, not genotypic
    # fixation. P is external siring, excluding own selfed seed double count.
    selfing=[]
    for s in (0.,.5,.6,.8,1.):
        for d in (0.,.5,.9,.95,1.):
            selfing.append({"selfed_fraction_at_mature_seed":s,"assumed_relative_postmature_inbreeding_depression":d,
              "relative_nuclear_gene_copy_value_per_mature_seed":(1-s)+2*s*(1-d)})
    return {"status":"v67_conditional_accounting_frontier_not_estimated_selection",
       "source_study":"10.1002/ajb2.16239","original_zip_sha256":ZIP_SHA,
       "source_aggregate_sha256":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (P55,P64)},
       "by_year":results,"selfing_sensitivity":selfing,
       "limits":["same originating Daphne jezoensis population; not new independent field data",
         "fruit count is not surviving offspring or a calibrated maternal contribution",
         "sampled female-mother paternity is not total external siring; same-cohort source ID audit separately required",
         "theta includes unknown conversion, sampling expansion and reproductive values",
         "source father assignments depend on spatial prior; p095 not sire truth",
         "observational flower-size/crop elasticity is not causal resource response",
         "future gain h and c and eta unmeasured; no frequency-dependent fitness or demographic ESS",
         "no inference for focal Daphne pseudomezereum"]}
if __name__=="__main__":
    OUT.mkdir(exist_ok=True)
    x=main()
    (OUT/"daphne_residual_pistil_frontier_v67.json").write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    summary=[]
    for z in x["by_year"]:
        summary.append({k:v for k,v in z.items() if k!="frontier_grid"})
    print("DAPHNE_V67_FRONTIER",json.dumps(summary,allow_nan=False),flush=True)
