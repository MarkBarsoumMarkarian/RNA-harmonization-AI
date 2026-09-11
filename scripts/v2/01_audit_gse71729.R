options(stringsAsFactors = FALSE)

dir.create("data/v2/cache", recursive = TRUE, showWarnings = FALSE)
dir.create("results/v2", recursive = TRUE, showWarnings = FALSE)

url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/",
  "GSE71729/matrix/GSE71729_series_matrix.txt.gz"
)
cache_file <- "data/v2/cache/GSE71729_series_matrix.txt.gz"

if (!file.exists(cache_file)) {
  message("Downloading authoritative GSE71729 Series Matrix metadata...")
  download.file(url, cache_file, mode = "wb", quiet = FALSE)
}

con <- gzfile(cache_file, open = "rt")
on.exit(close(con), add = TRUE)
meta <- character()
repeat {
  block <- readLines(con, n = 100L, warn = FALSE)
  if (!length(block)) break
  meta <- c(meta, block[startsWith(block, "!")])
  if (any(startsWith(block, "!series_matrix_table_begin"))) break
}
close(con)

sample_values <- function(prefix) {
  hit <- meta[startsWith(meta, prefix)]
  if (!length(hit)) stop("Missing Series Matrix field: ", prefix)
  values <- strsplit(hit[[1L]], "\t", fixed = TRUE)[[1L]][-1L]
  sub('"$', "", sub('^"', "", values))
}

titles <- sample_values("!Sample_title")
gsm <- sample_values("!Sample_geo_accession")
platform <- unique(sample_values("!Sample_platform_id"))

if (length(titles) != length(gsm)) stop("GSE71729 title/accession mismatch")
if (length(gsm) != 357L) stop("Expected 357 GSE71729 samples; found ", length(gsm))
if (!identical(platform, "GPL20769")) {
  stop("Expected GSE71729 platform GPL20769; found ", paste(platform, collapse = ", "))
}

specimen_class <- ifelse(
  grepl("-Primary-Pancreas$", titles), "primary_tumor",
  ifelse(
    grepl("-Met-", titles), "metastasis",
    ifelse(
      grepl("-CellLine$", titles), "cell_line",
      ifelse(
        grepl("-Normal-Pancreas$", titles), "normal_pancreas",
        ifelse(grepl("-Normal-", titles), "normal_distant", "unresolved")
      )
    )
  )
)

manifest <- data.frame(
  sample_id = gsm,
  title = titles,
  dataset = "GSE71729",
  platform_id = "GPL20769",
  assay = "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F",
  specimen_class = specimen_class,
  eligible_primary_tumor = specimen_class == "primary_tumor",
  survival_in_geo_series_matrix = FALSE,
  stringsAsFactors = FALSE
)

if (anyDuplicated(manifest$sample_id)) stop("Duplicated GSE71729 GSM identifiers")
if (any(manifest$specimen_class == "unresolved")) stop("Unresolved specimen titles")

audit <- as.data.frame(table(manifest$specimen_class), stringsAsFactors = FALSE)
names(audit) <- c("specimen_class", "n")
audit$eligible_for_primary_tumor_survival <- audit$specimen_class == "primary_tumor"
audit$has_survival_in_geo_series_matrix <- FALSE
audit <- audit[order(match(
  audit$specimen_class,
  c("primary_tumor", "metastasis", "cell_line", "normal_pancreas", "normal_distant")
)), ]

expected <- c(
  primary_tumor = 145L,
  metastasis = 61L,
  cell_line = 17L,
  normal_pancreas = 46L,
  normal_distant = 88L
)
observed <- setNames(audit$n, audit$specimen_class)[names(expected)]
if (!identical(as.integer(observed), as.integer(expected))) {
  stop("Specimen counts do not match the authoritative GEO design")
}

write.csv(manifest, "results/v2/gse71729_manifest.csv", row.names = FALSE)
write.csv(audit, "results/v2/gse71729_sample_audit.csv", row.names = FALSE)

message("GSE71729 audit complete:")
print(audit, row.names = FALSE)
message(
  "Gate: expression is available for 145 primary tumours, but the GEO Series ",
  "Matrix does not contain patient-level survival time/status. Do not call ",
  "unlabelled predictions external validation."
)
