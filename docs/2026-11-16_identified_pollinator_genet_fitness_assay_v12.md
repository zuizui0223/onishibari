# Preregisterable field assay: unique pollinator arrivals versus repeated probes, pollen quality and maternal/paternal fitness

Date: 2026-11-16 JST; design protocol v0.12; **not new observation data, a pre-data analysis plan**. Target focal system: functional flower-morph variation in *Daphne pseudomezereum* (Onishibari). Its morphological sex system, retained pistil functionality, pollen self-compatibility, actual floral construction biomass, genetic siring and reliable season-long pollinator identities remain **unverified in this repository**. Never label all apparently male flowers non-female-functional until controlled crosses.

## Novelty boundary
- Makino, Ohashi & Sakai (2007), *Functional Ecology* DOI 10.1111/j.1365-2435.2006.01211.x already followed **individually tagged Bombus diversus** at 27 *Cirsium purpuratum* plants and split unique visiting foragers from repeated returns. Large floral displays increased regular and occasional visitor densities; surrounding flower density increased occasional visitors but suppressed regular revisit rates. Distinguishing visitors vs returns is not a new discovery.
- Mitchell et al. 2004 DOI 10.1111/j.1365-2435.2004.00812.x manipulated *Mimulus ringens* real flower number and recorded within-plant consecutive probes; Karron et al. 2009 DOI 10.1093/aob/mcp030 linked visit sequence to actual fruit-level selfing.
- Joffard et al. 2026 DOI 10.1002/ajb2.70276, original Dataverse 10.57745/CF1VNW tracks bee-level successive stigma pollen deposition (25 tagged pollinators) **and separately** offspring father diversity (21 mothers); these are two different experiments. The new contribution would be **same-genet and same-season simultaneously** observing distinct pollinator identity, effective compatible pollen, viable maternal offspring and male fathering, then testing context-dependent limits, reproductive capacity and future survival.
- Barbot et al. 2022 DOI 10.1098/rspb.2022.1987, original Dryad 10.5061/dryad.0rxwdbs4d provides plant-level visitor number, total visits, seeds, genetic mate counts and male MEMM siring, but lacks verified unique bee IDs. v0.12 new source fit: visitor vs visits Spearman rho .983 females/.835 males; adjusted male siring slope positive for each separately but not jointly, LOOCV favors neither jointly. It cannot replace marked visitation.

## Prespecified mechanistic hypotheses
H1 **Detection vs repeated probing**: At fixed natural flower number and observation effort, increasing visual advertisement may increase (a) newly arriving visitor identities per 2h and/or (b) same-individual returns and consecutive within-genet probes. H1 is two outcomes; direction of each depends on surrounding floral density (per Makino 2007).

H2 **Pollen replenishment and sequence**: For a known individually identified pollinator/return bout, the probability and amount of **conspecific genetically compatible** pollen deposited per stigma declines with the number of previous flowers it probed since an external donor; independent newly pollen-loaded bouts replenish that opportunity. Do NOT infer compatible deposition merely from total visible grains: verify compatible germination/tubes or donor assignment on separate flower cohorts.

H3 **Sex-function specific payoffs**: A natural true pollen-producing flower may raise genetically assigned paternal offspring even if dummy visually added flowers have the same apparent display, whereas additional apparent display may increase arrivals without increasing female viable seeds in already pollen-sufficient contexts. Measure true pollen grains and actual flower biomass per donor—not only number of open flowers or visual area.

H4 **Temporal developmental gate or insurance**: Early seed/fruit resource allocation can depress development of later fertilized flowers; conversely subsequent flowering may rescue fruit set when early pollen receipt was unsuccessful. *Brassica* 2026 and *Sagittaria trifolia* 2026 demonstrate opposite original-source temporal responses, not a shared mechanistic treatment effect.

H5 **Future compensation**: A high current-year seed or paternal payoff can carry future survival/flowering costs, which must be estimated on independently identified genets over multiple years. *Silene latifolia* 2004–2007 observational association was sensitive to final missing survival and endpoint scan; it does not prove floral dry-mass costs.

## Conservation-first staged implementation
Phase 0 **verify morphological and genetic identities / plant safety**: record an individual **genet ID** distinct from stems/ramets; phenotype both presumptive sex morphs; minimal viability controls (bagging, hand self, hand outcross using genetically distinct donors, and sham controls), appropriate permits and contamination controls. Positive hand outcross and stigma pollen-tube controls distinguish failed pollen delivery from nonfunctional gynoecium. *Self-incompatibility must never be asserted from flower shape alone.*

Phase 1 **observational non-destructive tagging**: matched observation sessions of the same available plants for nominal 2 h per session (effort recorded as actual active minutes, not imputed if interrupted); plant, full date, start, observer and camera IDs; true open flower count and branch; counted visits; pollinator individual identities only when identifiable with ethical tags or unambiguous continuous video tracking; explicit UNKNOWN for unmarked/unresolved insects; per-visit same-genet flower sequence and between-plant donor movement if observable; background floral density and pollinator activity. Record "no visit with monitored effort" as zero, missing effort as missing.

Phase 2 **visual-only manipulation when sufficient safe replication**: random assign sham versus dummy visual display within comparable genet/day/morph blocks, making sure fake displays do not supply viable pollen/nectar/ovules. To avoid a site-wide confounded before–after, counterbalance order when repeated within a plant and include period/date controls and observer/camera-detection controls. Outcomes unique identified visitor arrival and within-genet repeated probing must be modeled separately.

Phase 3 **reproductive functions / experimental cohorts only if Phase 0 works**: individually tagged buds randomized across (i) natural open access, (ii) controlled hand outcross, (iii) hand-self and sham when permitted; pollen count/deposition/pollen tube assay on separate destructively sampled flowers, embryo status, mature viable seeds, genetic sire assignments. A flower removed for stigma assays cannot later set seed; never treat these as single-flower longitudinal transitions. Fruit-load thinning requires abundant common-garden material and permits, **not a default on rare wild genets**. Removal affects both resource and hormonal signals and does not isolate one mechanism.

Phase 4 **future fitness and sex**: repeat exactly identified genotype monitoring next reproductive seasons; viability, survival, recruitment, maternal seed, paternal seed genotyping. No claiming complete lifetime fitness without adult survival / descendant establishment.

Existing ethically frozen ceiling for rare wild plants: at most **50% of confirmed accessible buds and at most 10 confirmed buds per independently identified genet** across combined experimental plus destructive allocations. Do not relax to obtain sample sizes; exclude underpowered or sensitive genets, and use common garden or alternate abundant species as a functional demonstration.

## Minimal anonymized structured tables (no coordinates or original individual codes in public artifact)

**genet_season**: study_site_id_private, genet_id_private, calendar_year, documented_morph, reproductive_morph_function_verified, genet_size_pretreatment, aboveground_biomass_if_measured, independent_genet_confirmed, subsequent_survival_and_display_measured.

**flower**: genet_id_private, flower_id_private, branch_id_private, true_or_dummy_flower, unique_new_flower_indicator, date_opened, duration_open_hours, dry_flower_mass_mg_if_measured, pollen_grains_per_true_flower, functional_pistil_test_group, treatment_randomization_block, pollen_assay_cohort_vs_seed_fate_cohort, anther_and_stigma_function, fruit_and_seed_outcome_if_non_destructive.

**observation**: site_id_private, genet_id_private, start_and_end_JST, exact_observer_effort_minutes, camera_identity, observer_id_private, observer_presence, natural_open_flower_number, dummy_display_number, neighboring_bloom_density, measured_sun_and_weather, individual_identification_completeness, zero_visits_with_positive_effort.

**visitor_event**: observation_id_private, visit_timestamp_JST, real_identified_insect_id_private_or_UNKNOWN, insect_taxon, bee_tag_confirmed, fresh_arrival_to_focal_genet, reentry_after_leaving_focal_genet, probed_flower_id_private, within_current_bout_probe_order, previous_external_plant_donor_confirmed_or_UNKNOWN, physical_reproductive_contact_recorded, visitor_identity_detection_confidence.

**pollen_stigma_assay**: flower_id_private, visitor_id_private_or_UNKNOWN, deposited_conspecific_grains, pollen_donor_genet_if_genotyped, compatible_tube_number_if_assayed, effective_pollination_indicator_only_if_verified, sampling_date, source_cohort, methods_QC.

**seed_fate**: flower_id_private, preexisting_ovules_if_counted, fertilized_ovules_if_assayed, aborted_embryos_if_assayed, fruit_present, number_mature_seeds, number_viable_seeds, maternal_genet, paternal_genet_if_genotyped, followup_established_offspring_known.

## Analysis and identifiability specification

**Unit of independent inference**: genet × independent session/flower cohort, clustered within genet/site/date (not each flower, grain or visit). If individual visitor identity is unreliable, **do not estimate new-bee versus same-bee return effects**; report only total visits and visitor-identity coverage. For time-to-event models, define time since independently verified prior external donor, but not from camera gaps.

Primary advertisement response:
- count of **identified unique pollen-bearing visitor individuals arriving** per exactly observed 120 min, with log-effort offset, zero-inflation checked, standardized surrounding floral density and biologically defined genotype block;
- number of **same-identified-individual plant returns** (not number of flowers probed by one continuous visit);
- number of consecutive within-genet flower probes per identified arrival.

Primary functional reproductive endpoint:
- maternal genetically viable seeds **per genet across all eligible flowers**, as well as fruit/ovule probabilities separately;
- paternal genetically assigned offspring **per focal genet**, with parentage detection probabilities and mother's genotype availability;
- first compatible pollen on a flower with positive exposure vs seed formation and eventually future survival.

Two covariate sets are needed: an intention-to-treat effect of dummy display with **pretreatment-only adjustment**, and a *descriptive* conditioning on visitor composition to explore potential mediation. Post-treatment visit number must not be adjusted in the first causal effect model, as this blocks part of the intervention. In the second model, the quantity of visits and number of visitor identities need enough independent variation, otherwise coefficients are unidentifiable (v0.12 source observed 2022 collinearity rho .983 in females/.835 in males).

Falsification examples:
- dummy display increases both distinct arrivals and seed output, contradicts a simple "advertisement never produces fitness" rule;
- visual display increases only returns with no extra compatible pollen, supports a return-only bottleneck under the monitored configuration;
- independent arrivals increase but genetic fathers do not, then **paternal donor diversity, compatibility or postzygotic capacity** remains a candidate gate;
- natural flower production increases paternal seed siring but dummy display does not, after a controlled genet and flower-mass comparison, supports reproductive substrate beyond pure visual signal;
- larger true floral resource expenditure increases same-year success but reduces future survival under *randomized resource allocation*, supports a true trade-off rather than observational life-history covariance.

## Stopping rules and objective limits
- No site-level inference when insufficient independent identified genets; do not inflate n by counting flowers or repeated bouts as independent plants.
- Report ID-identification rate and explicitly estimate bounds under UNKNOWN insects instead of assigning them arbitrarily to fresh or returning status.
- No "proven compatibility" without pollen-tube/hand-cross controls.
- No perfect pollinator tag identity from species classification alone (YOLO or generic bumblebee class is not an individual tag).
- No universal fitness elasticity in mg of floral mass without paired true dry biomass and genetically verified descendant measures.
- No assertion that repeated visits produce gene-flow across genets without independent pollen donor genotype evidence.

**Desired central publishable claim**: ecological conditions change which stage of floral investment (advertisement, independent pollen replenishment, gamete supply, compatible fertilization, development and temporal insurance) limits **genet-wide female and male offspring**, with cross-clade replication. Existing global 1241-species flower/pollen result and 2007 tagged-bee work establish prior theory, so the novelty must be a new matched causal design or independent source-backed comparative estimates, not merely restating "multiple stages".
