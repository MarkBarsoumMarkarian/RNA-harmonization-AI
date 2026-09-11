# Leakage-Safe Cross-Cohort Transcriptomic Survival Analysis in PDAC

This repository contains a corrected reanalysis of the project originally
reported as *Batch-Harmonized Machine Learning Framework for Cross-Cohort RNA
Biomarker Discovery in Pancreatic Adenocarcinoma*.

> **Version 2 correction notice (2026-09-12):** The original analysis does not
> establish a prognostic classifier, external validation, or five validated
> biomarkers. GSE71729 is an Agilent microarray containing mixed primary,
> metastatic, normal, and cell-line specimens; the original outcome also
> reduced survival to alive/dead status and ignored follow-up time and
> censoring. The v1 scripts and artifacts are preserved for provenance but must
> not be used as evidence of predictive performance.

Read the complete [forensic audit and interpretation](V2_AUDIT.md) before using
the results.

## Corrected scientific question

Can a transcriptomic survival model developed in one pancreatic ductal
adenocarcinoma cohort retain discrimination in a completely held-out cohort
measured on another platform, without letting the held-out samples influence
normalization, feature selection, tuning, or model fitting?

The current answer is **not convincingly**. The corrected analysis finds modest
internal TCGA discrimination but no demonstrated cross-platform transport.

## Headline results

| Analysis | Patients / events | C-index (95% patient-bootstrap CI) | Interpretation |
|---|---:|---:|---|
| TCGA repeated nested 5-fold CV | 176 / 92 | 0.630 (0.567-0.694) | Modest internal discrimination |
| Train TCGA, hold out GSE85916 | 79 / 57 | 0.527 (0.435-0.607) | No demonstrated transport |
| Train GSE85916, hold out TCGA | 176 / 92 | 0.576 (0.508-0.647) | Weak reverse-direction transport |
| Original fixed five-gene Cox, train TCGA, test GSE85916 | 79 / 57 | 0.569 (0.476-0.657) | Confidence interval includes chance |

The original five genes—LAMC2, DKK1, ITGB6, GPRC5A, and MAL2—remain exploratory
hypotheses. All have adverse univariate associations in TCGA, but this is the
original discovery cohort. In the two-cohort random-effects analysis, none
survives BH correction across the five and heterogeneity is high. LAMC2 is only
nominally significant: HR 1.58 per cohort SD, 95% CI 1.04-2.42, unadjusted
p=0.034, BH FDR=0.170.

## Cohorts and eligibility

| Dataset | Platform | Prognostic use in v2 | Patients / events |
|---|---|---|---:|
| TCGA-PAAD PanCancer Atlas | RNA-seq; cBioPortal RSEM snapshot | Development and internal validation | 176 / 92 |
| GSE85916 | GPL13667 Affymetrix Human Genome U219 array | Retrospective whole-cohort holdout | 79 / 57 |
| GSE71729 | GPL20769 Agilent microarray | Provenance/sample-composition audit only | Not eligible |

GSE71729 contains 145 primary tumors, 61 metastases, 17 cell lines, 46 normal
pancreas samples, and 88 distant-site normal samples. It is excluded from
prognostic modelling because the repository has no auditable patient-level
survival-time/event join for its primary tumors.

The fast TCGA reanalysis uses the cBioPortal PanCancer Atlas RSEM snapshot. A
manuscript-grade final analysis should pin and reprocess a GDC STAR-count
release and preserve raw-file provenance.

## Leakage controls in v2

- One primary tumor per patient.
- Positive overall-survival time and explicit death/censoring indicator.
- `Surv(time, event)` rather than alive/dead classification.
- Repeated nested patient-level cross-validation for model development.
- MAD filtering, imputation, centering/scaling, outcome-guided screening,
  hyperparameter tuning, and penalized Cox fitting learned only from training
  partitions.
- Entire cohorts held out from model development.
- No joint ComBat, cohort-wise test scaling, or fitting to the held-out cohort
  distribution.
- Sample-wise percentile-rank transformation for the cross-platform analysis.
- Cohort-specific Cox estimates and random-effects meta-analysis for candidate
  biomarkers.
- Machine-readable sample flows, predictions, selected features, input URLs,
  byte sizes, MD5 checksums, package metadata, and confidence intervals.

## Reproduce the fast v2 analysis

Run from the repository root with R. Base R plus the recommended `survival`
package are required.

```bash
Rscript scripts/v2/01_audit_gse71729.R
Rscript scripts/v2/02_prepare_tcga.R
Rscript scripts/v2/03_nested_cox_tcga.R
Rscript scripts/v2/06_prepare_gse85916.R
Rscript scripts/v2/04_loco_cox.R
Rscript scripts/v2/05_meta_analyze_candidates.R
Rscript scripts/v2/07_validate_fixed_panel.R
```

Downloads are cached under `data/v2/cache/` and are excluded from version
control. See [scripts/v2/README.md](scripts/v2/README.md) for the cohort schema,
outputs, and extension instructions.

## Main v2 outputs

- `results/v2/gse71729_sample_audit.csv`
- `results/v2/tcga_sample_flow.csv`
- `results/v2/gse85916_sample_flow.csv`
- `results/v2/tcga_nested_cv_performance.csv`
- `results/v2/tcga_oof_predictions.csv`
- `results/v2/tcga_selected_gene_stability.csv`
- `results/v2/loco_performance.csv`
- `results/v2/loco_predictions.csv`
- `results/v2/fixed_five_external_performance.csv`
- `results/v2/candidate_cohort_cox.csv`
- `results/v2/candidate_random_effects_meta.csv`
- `results/v2/tcga_input_manifest.csv`
- `results/v2/gse85916_input_manifest.csv`
- `results/v2/run_manifest.txt`

## What the original analysis got wrong

The detailed evidence is in [V2_AUDIT.md](V2_AUDIT.md). In brief:

- GSE71729 was incorrectly described as platform-matched RNA-seq.
- All 357 heterogeneous GSE71729 samples were combined without specimen
  filtering.
- Unlabelled predictions were described as external validation.
- Alive/dead status replaced time-to-event analysis.
- ComBat was fitted jointly to training and purported external samples.
- Full-data preprocessing and top-variance selection preceded Random Forest
  OOB scoring.
- “Training Accuracy” was actually OOB accuracy.
- “Class Balance” values were class-specific recalls, and the raw-model recalls
  were reported incorrectly.
- The Shiny app's hard-coded 92.6% accuracy is unsupported; the serialized
  corrected Random Forest has approximately 64.0% OOB accuracy.

## Interpretation and next phase

Version 2 is a retrospective correction, not a clinical model. GSE85916 was
selected during the repair and should be treated as exploratory external
testing rather than a prospectively locked confirmation cohort.

The next phase is a multi-cohort transportability study that evaluates locked
PDAC transcriptomic signatures across eligible survival cohorts and compares
clinical-only, expression-only, and clinical-plus-expression models. It should
report discrimination, calibration, time-dependent performance, biological
signal preservation, and uncertainty under whole-cohort holdout validation.

## Data sources

- [GSE71729](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71729)
- [GSE85916](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85916)
- [GDC TCGA-PAAD](https://portal.gdc.cancer.gov/projects/TCGA-PAAD)
- [cBioPortal PanCancer Atlas snapshot](https://github.com/cBioPortal/datahub/tree/master/public/paad_tcga_pan_can_atlas_2018)
- [Original bioRxiv preprint](https://doi.org/10.1101/2025.11.14.688421)

## Repository history

The original scripts, figures, serialized models, and results remain in their
existing directories as an auditable record. They are v1 artifacts and are
superseded by `scripts/v2/`, `results/v2/`, and this correction notice. The
original repository state is preserved at tag `v1-original-preprint`.

## License

MIT License. See [LICENSE](LICENSE).
