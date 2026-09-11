# Fast corrected v2

The v2 scripts are intentionally separate from the invalid v1 artifacts.

## Fast reproducible run

```r
source("scripts/v2/01_audit_gse71729.R")
source("scripts/v2/02_prepare_tcga.R")
source("scripts/v2/03_nested_cox_tcga.R")
source("scripts/v2/06_prepare_gse85916.R")
source("scripts/v2/04_loco_cox.R")
source("scripts/v2/05_meta_analyze_candidates.R")
source("scripts/v2/07_validate_fixed_panel.R")
```

Only base R and the recommended `survival` package are required. Downloads are
cached under `data/v2/cache/` and are not committed.

The fast TCGA analysis uses the cBioPortal PanCancer Atlas RSEM snapshot so the
corrected statistical workflow can run quickly. It is a provisional v2 result,
not an exact reprocessing of GDC STAR counts. A final manuscript analysis should
pin a GDC release, retain raw-count provenance and checksums, and perform
training-only normalization inside resampling where applicable.

## Outputs

- `results/v2/gse71729_sample_audit.csv`
- `results/v2/gse71729_manifest.csv`
- `results/v2/tcga_sample_flow.csv`
- `results/v2/tcga_nested_cv_performance.csv`
- `results/v2/tcga_oof_predictions.csv`
- `results/v2/tcga_selected_gene_stability.csv`
- `results/v2/original_five_gene_cox.csv`
- `results/v2/run_manifest.txt`
- `results/v2/gse85916_sample_flow.csv`
- `results/v2/tcga_input_manifest.csv`
- `results/v2/gse85916_input_manifest.csv`
- `results/v2/loco_performance.csv`
- `results/v2/loco_predictions.csv`
- `results/v2/candidate_random_effects_meta.csv`
- `results/v2/fixed_five_external_performance.csv`

## External validation / LOCO input

`scripts/v2/06_prepare_gse85916.R` creates `data/v2/cohorts.rds` from TCGA and
the GSE85916 primary-tumour microarray cohort. The generic LOCO input is a
named list. Each cohort contains:

```r
list(
  expression = expression_matrix, # genes x samples; unique gene symbols
  clinical = data.frame(
    sample_id = colnames(expression_matrix),
    patient_id = ...,
    time = ...,                  # positive, one declared unit
    event = ...,                 # 1 death/event, 0 censored
    specimen_class = "primary_tumor",
    stringsAsFactors = FALSE
  ),
  platform = "declared platform and preprocessing"
)
```

The script refuses unlabeled cohorts, non-primary specimens, duplicated
patients, non-positive times, or fewer than ten events. It never ComBat-fits the
held-out cohort with the training cohorts. The cross-platform default is a
per-sample percentile-rank transform followed by training-only gene scaling.
