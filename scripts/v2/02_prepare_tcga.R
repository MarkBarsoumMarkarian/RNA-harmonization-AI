options(stringsAsFactors = FALSE, timeout = max(600, getOption("timeout")))

dir.create("data/v2/cache", recursive = TRUE, showWarnings = FALSE)
dir.create("results/v2", recursive = TRUE, showWarnings = FALSE)

base_url <- paste0(
  "https://media.githubusercontent.com/media/cBioPortal/datahub/master/",
  "public/paad_tcga_pan_can_atlas_2018/"
)
files <- c(
  "data_mrna_seq_v2_rsem.txt",
  "data_clinical_patient.txt",
  "data_clinical_sample.txt"
)

for (filename in files) {
  destination <- file.path("data/v2/cache", filename)
  if (!file.exists(destination)) {
    message("Downloading ", filename, "...")
    download.file(paste0(base_url, filename), destination, mode = "wb")
  }
}

input_paths <- file.path("data/v2/cache", files)
input_info <- file.info(input_paths)
write.csv(
  data.frame(
    file = files,
    source_url = paste0(base_url, files),
    bytes = input_info$size,
    md5 = unname(tools::md5sum(input_paths)),
    stringsAsFactors = FALSE
  ),
  "results/v2/tcga_input_manifest.csv",
  row.names = FALSE
)

read_cbio <- function(path) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE,
    quote = "",
    na.strings = c("", "NA", "NaN")
  )
}

patient <- read_cbio("data/v2/cache/data_clinical_patient.txt")
sample <- read_cbio("data/v2/cache/data_clinical_sample.txt")
expression_file <- "data/v2/cache/data_mrna_seq_v2_rsem.txt"
expression <- read.delim(
  expression_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  na.strings = c("", "NA", "NaN")
)

required_patient <- c("PATIENT_ID", "OS_STATUS", "OS_MONTHS")
required_sample <- c("PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE")
if (!all(required_patient %in% names(patient))) stop("Missing patient fields")
if (!all(required_sample %in% names(sample))) stop("Missing sample fields")
if (!all(c("Hugo_Symbol", "Entrez_Gene_Id") %in% names(expression))) {
  stop("Unexpected cBioPortal expression schema")
}

flow <- data.frame(step = character(), n = integer())
add_flow <- function(step, n) {
  flow <<- rbind(flow, data.frame(step = step, n = as.integer(n)))
}

add_flow("PanCancer Atlas clinical patients", nrow(patient))
sample <- sample[sample$SAMPLE_TYPE == "Primary", , drop = FALSE]
add_flow("Primary tumour sample records", nrow(sample))

clinical <- merge(sample, patient, by = "PATIENT_ID", all.x = TRUE, sort = FALSE)
clinical$time <- suppressWarnings(as.numeric(clinical$OS_MONTHS))
clinical$event <- ifelse(
  grepl("^1:", clinical$OS_STATUS), 1L,
  ifelse(grepl("^0:", clinical$OS_STATUS), 0L, NA_integer_)
)
clinical <- clinical[is.finite(clinical$time) & clinical$time > 0 & !is.na(clinical$event), ]
add_flow("Positive OS time and known event status", nrow(clinical))

clinical <- clinical[order(clinical$PATIENT_ID, clinical$SAMPLE_ID), ]
clinical <- clinical[!duplicated(clinical$PATIENT_ID), ]
add_flow("One primary tumour per patient", nrow(clinical))

sample_columns <- intersect(clinical$SAMPLE_ID, names(expression))
clinical <- clinical[match(sample_columns, clinical$SAMPLE_ID), ]
add_flow("Matched expression and survival", nrow(clinical))

gene_symbol <- trimws(as.character(expression$Hugo_Symbol))
keep_gene <- !is.na(gene_symbol) & nzchar(gene_symbol)
values <- as.matrix(expression[keep_gene, sample_columns, drop = FALSE])
storage.mode(values) <- "double"
gene_symbol <- gene_symbol[keep_gene]

# The cBioPortal file contains non-negative RSEM abundance estimates.
# A log2(x + 1) transform is fixed before any resampling. Every data-adaptive
# filter, scaler, selector and model fit occurs later inside training folds.
values <- log2(values + 1)

if (anyDuplicated(gene_symbol)) {
  index <- split(seq_along(gene_symbol), gene_symbol)
  values <- do.call(rbind, lapply(index, function(i) {
    if (length(i) == 1L) values[i, ] else colMeans(values[i, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(values) <- names(index)
} else {
  rownames(values) <- gene_symbol
}

clinical_out <- data.frame(
  sample_id = clinical$SAMPLE_ID,
  patient_id = clinical$PATIENT_ID,
  cohort = "TCGA_PAAD_PanCancer_Atlas_2018",
  platform = "Illumina RNA-seq; cBioPortal RSEM abundance snapshot",
  specimen_class = "primary_tumor",
  time = clinical$time,
  event = clinical$event,
  age = suppressWarnings(as.numeric(clinical$AGE)),
  stage = clinical$AJCC_PATHOLOGIC_TUMOR_STAGE,
  stringsAsFactors = FALSE
)
rownames(clinical_out) <- clinical_out$sample_id

if (!identical(colnames(values), clinical_out$sample_id)) stop("Expression/clinical order failure")
if (anyDuplicated(clinical_out$patient_id)) stop("Patient duplication after filtering")
if (sum(clinical_out$event) < 10L) stop("Too few events")

tcga <- list(
  expression = values,
  clinical = clinical_out,
  platform = unique(clinical_out$platform),
  provenance = list(
    source = "cBioPortal Datahub paad_tcga_pan_can_atlas_2018",
    purpose = "fast provisional v2 reanalysis; not exact GDC STAR-count reproduction"
  )
)

saveRDS(tcga, "data/v2/tcga_cbioportal.rds")
write.csv(clinical_out, "results/v2/tcga_clinical_used.csv", row.names = FALSE)
write.csv(flow, "results/v2/tcga_sample_flow.csv", row.names = FALSE)

message(
  "TCGA prepared: ", ncol(values), " patients, ", sum(clinical_out$event),
  " events, ", nrow(values), " genes."
)
print(flow, row.names = FALSE)
