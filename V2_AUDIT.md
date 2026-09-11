# Forensic audit and v2 correction plan

Audit date: 2026-09-12

Scope: bioRxiv v1 (doi: 10.1101/2025.11.14.688421), the repository at commit
`f316ab2`, the serialized model/result objects, and authoritative GDC/GEO
records. The original scripts and results remain in place as a record; they are
not valid evidence for the manuscript's predictive or biomarker claims.

## Bottom line

The published analysis does not establish a prognostic classifier, external
validation, or five validated prognostic biomarkers. The main result is an
exploratory demonstration that ComBat can force two highly confounded data
collections to overlap in PCA space. That result is not, by itself, evidence
that biological signal was preserved or that a model generalizes.

The useful scientific goal can be salvaged by treating the five genes as
hypotheses, using time-to-event outcomes with censoring, keeping every learned
preprocessing and selection step inside resampling, and reserving entire
outcome-labelled cohorts for external validation.

## Confirmed errors

### Fatal to the original conclusions

1. **GSE71729 is not RNA-seq.** It is expression profiling by Agilent array
   (GPL20769). The manuscript and repository repeatedly call it Illumina HiSeq
   RNA-seq and use that false premise to claim a platform-matched analysis.
2. **The 357 samples are not 357 patients eligible for prognosis.** The GEO
   record contains 145 primary PDACs, 61 metastases, 17 cell lines, 46 normal
   pancreas samples, and 88 normal samples from distant sites. The old pipeline
   used all 357 as a single validation cohort without specimen filtering.
3. **There was no external validation.** No GSE71729 outcome was joined and no
   discrimination, calibration, or survival statistic was calculated. A
   histogram of unlabeled predictions cannot demonstrate correctness,
   confidence, uncertainty, transport, or generalization.
4. **Vital status was substituted for survival.** `Alive` and `Dead` at last
   contact discard follow-up time and censoring. A patient alive after two days
   and one alive after years receive the same label; early-censored patients
   can be treated as favourable outcomes. This is not a valid prognostic
   endpoint.
5. **The external cohort influenced the correction applied to training and
   test samples.** ComBat was fit once to the combined TCGA and GSE71729 matrix.
   This is a transductive operation, not a deployable training-only transform,
   and invalidates the claim of true held-out validation.

### Leakage and validation defects

- TCGA VST, the 14,137-gene intersection, ComBat, and the top-500 variance
  filter were computed before Random Forest OOB evaluation. An OOB patient
  therefore influenced preprocessing and feature selection used to predict
  that patient. OOB is not a complete pipeline-level validation here.
- The 500 features were selected from the corrected full TCGA matrix and then
  reused for the raw model. This makes the raw-versus-corrected comparison
  asymmetric. The reported 0.5 percentage-point difference has no uncertainty
  estimate and is not evidence of improvement.
- XGBoost was fit on all 178 patients with 500 predictors, depth 6, learning
  rate 0.3 and 100 rounds, with no tuning and no internal or external
  performance estimate. Its importance ranking is not independent validation.
- Batch, cohort, platform, specimen type, and disease setting are confounded.
  `mod = NULL` does not prove biology is preserved. A near-zero batch
  silhouette can be achieved by removing real cohort/specimen biology.
- Method choice used only batch mixing. There was no biology-retention metric,
  negative control, outcome-aware validation, or uncertainty interval.
- Harmony and fastMNN are primarily integration methods for cell-level/single-
  cell settings. Their use here on bulk, cross-platform cohort matrices is not
  justified, and Harmony coordinates are not a gene-expression transform for
  a prognostic assay.
- Tree models do not require conventional feature scaling, so lack of scaling
  is not the central defect. The actual problem is incomparable platform
  measurement plus a correction learned from the held-out cohort.

### Reporting and reproducibility errors

- `Training Accuracy` is calculated as `1 - OOB error`; it is OOB accuracy,
  not resubstitution training accuracy.
- `Alive: 52.9%, Dead: 74.2%` are corrected-model class-specific OOB recalls,
  not class balance. Class balance is 85/178 versus 93/178.
- The raw model actually has 47.1% Alive recall and 78.5% Dead recall, but the
  manuscript prints the corrected model's 52.9%/74.2% values for both rows.
- The Shiny overview hard-codes `92.6%` model accuracy. The saved corrected RF
  object has 64.0% OOB accuracy. No repository result supports 92.6%.
- Calling a broader unlabeled prediction distribution “appropriate
  uncertainty” or “better generalization” is unsupported. It can equally
  indicate miscalibration or sample-type confounding.
- Gini/Gain rankings from models trained and interpreted on the same cohort do
  not establish prognostic association, novelty, robustness, or independent
  validation. The five genes must be described as exploratory candidates.
- Processed input data, a package lockfile, a session manifest, data checksums,
  and an executable end-to-end results target are absent. One script contains
  an absolute Windows path and documentation names do not match several files.

## What v2 changes

1. Filters to one primary tumour sample per patient and requires positive
   follow-up time plus an explicit event indicator.
2. Uses `Surv(time, event)` and penalized Cox regression. Performance is
   Harrell's concordance index with a patient-bootstrap confidence interval.
3. Performs variance filtering, outcome-guided screening, centering/scaling,
   penalty tuning, and model fitting inside the training partition of each
   outer split.
4. Provides leave-one-cohort-out validation. For cross-platform use, the
   default transform is within-sample percentile rank, which can be computed
   for a new sample without borrowing the held-out cohort distribution.
5. Keeps cohorts separate for association discovery and meta-analyzes
   cohort-specific Cox estimates instead of ComBat-pooling outcomes.
6. Emits machine-readable sample audits, out-of-fold predictions, performance
   estimates, selected-gene stability, and candidate-specific Cox results.

## Corrected fast-run results

These results are a retrospective reanalysis, not a confirmatory study. TCGA
uses the cBioPortal PanCancer Atlas RSEM snapshot rather than a pinned raw-count
GDC reprocessing.

| Analysis | Patients / events | C-index (95% bootstrap CI) | Interpretation |
|---|---:|---:|---|
| TCGA repeated nested 5-fold CV | 176 / 92 | 0.630 (0.567-0.694) | Modest internal discrimination |
| Train TCGA, hold out GSE85916 | 79 / 57 | 0.527 (0.435-0.607) | No demonstrated transport |
| Train GSE85916, hold out TCGA | 176 / 92 | 0.576 (0.508-0.647) | Weak reverse-direction transport |
| Original fixed five-gene Cox, train TCGA, test GSE85916 | 79 / 57 | 0.569 (0.476-0.657) | Confidence interval includes chance |

All five original candidates have adverse univariate associations in TCGA, but
the effect sizes attenuate in GSE85916. In a two-cohort random-effects analysis,
only LAMC2 is nominally significant (HR 1.58 per cohort SD, 95% CI 1.04-2.42,
unadjusted p=0.034); none survives BH correction across the five, and
heterogeneity is high (I-squared 75%-93%). These are hypotheses, not validated
biomarkers.

The TCGA nested-CV feature-selection stability also does not reproduce the
original panel: DKK1 appears in 12 of 25 outer fits, while LAMC2, ITGB6, GPRC5A
and MAL2 appear in none. The most frequent genes are FAM83A (25/25), UCA1
(24/25), LY6D (19/25), COL17A1 (18/25), and HMGA2 (16/25). This stability list
is exploratory and itself requires external confirmation.

## Interpretation gates

- **Current status:** corrected code, provisional nested-CV TCGA reanalysis,
  and a retrospective cross-platform test in GSE85916. The external result is
  negative and must not be relabelled as successful validation.
- **External-validation go gate:** at least one independent primary-PDAC cohort
  with expression, patient-level overall-survival time, event/censoring status,
  and an auditable sample-ID join.
- **Biomarker go gate:** directionally consistent cohort-specific effects,
  acceptable heterogeneity, multiplicity control, and pre-specified external
  performance with confidence intervals.
- **Clinical-use gate:** locked assay and model, independent clinical cohort,
  calibration and clinical-utility analysis. V2 does not meet this gate.

## Authoritative records

- GSE71729: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71729>
- GDC TCGA-PAAD: <https://portal.gdc.cancer.gov/projects/TCGA-PAAD>
- Moffitt et al. source study: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4912058/>
- GSE85916 external survival cohort:
  <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85916>
- cBioPortal Datahub TCGA PanCancer Atlas snapshot used for the fast v2 run:
  <https://github.com/cBioPortal/datahub/tree/master/public/paad_tcga_pan_can_atlas_2018>
