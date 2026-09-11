options(stringsAsFactors = FALSE, timeout = max(1200, getOption("timeout")))

dir.create("data/v2/cache", recursive = TRUE, showWarnings = FALSE)
dir.create("results/v2", recursive = TRUE, showWarnings = FALSE)

series_file <- "data/v2/cache/GSE85916_series_matrix.txt.gz"
platform_file <- "data/v2/cache/GPL13667.txt"
series_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85916/matrix/",
  "GSE85916_series_matrix.txt.gz"
)
platform_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?",
  "targ=self&acc=GPL13667&form=text&view=full"
)

if (!file.exists(series_file)) {
  message("Downloading GSE85916 Series Matrix...")
  download.file(series_url, series_file, mode = "wb")
}
if (!file.exists(platform_file)) {
  message("Downloading GPL13667 annotation (large file)...")
  download.file(platform_url, platform_file, mode = "wb")
}

input_paths <- c(series_file, platform_file)
input_info <- file.info(input_paths)
write.csv(
  data.frame(
    file = basename(input_paths),
    source_url = c(series_url, platform_url),
    bytes = input_info$size,
    md5 = unname(tools::md5sum(input_paths)),
    stringsAsFactors = FALSE
  ),
  "results/v2/gse85916_input_manifest.csv",
  row.names = FALSE
)

strip_fields <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
  sub('"$', "", sub('^"', "", fields))
}

# Read the Series Matrix once to locate the table and extract the authoritative
# sample-level survival fields. The expression table is then reread directly.
series_lines <- readLines(gzfile(series_file), warn = FALSE)
table_begin <- match("!series_matrix_table_begin", series_lines)
table_end <- match("!series_matrix_table_end", series_lines)
if (is.na(table_begin) || is.na(table_end) || table_end <= table_begin + 1L) {
  stop("Malformed GSE85916 Series Matrix")
}
if (!any(grepl('^!Series_platform_id\\t"GPL13667"$', series_lines))) {
  stop("GSE85916 is not recorded on GPL13667")
}
if (!any(grepl('^!Series_type\\t"Expression profiling by array"$', series_lines))) {
  stop("GSE85916 is not recorded as an expression microarray")
}

sample_ids <- strip_fields(series_lines[grep("^!Sample_geo_accession\\t", series_lines)[1L]])
sample_titles <- strip_fields(series_lines[grep("^!Sample_title\\t", series_lines)[1L]])
sources <- strip_fields(series_lines[grep("^!Sample_source_name_ch1\\t", series_lines)[1L]])
characteristic_lines <- grep("^!Sample_characteristics_ch1\\t", series_lines, value = TRUE)
os_line <- characteristic_lines[grepl('"os\\.year:', characteristic_lines)][1L]
death_line <- characteristic_lines[grepl('"death:', characteristic_lines)][1L]
if (is.na(os_line) || is.na(death_line)) stop("Missing OS time or death fields")

os_years <- suppressWarnings(as.numeric(sub("^os\\.year:\\s*", "", strip_fields(os_line))))
event <- suppressWarnings(as.integer(sub("^death:\\s*", "", strip_fields(death_line))))
if (length(unique(c(length(sample_ids), length(sample_titles), length(sources),
                    length(os_years), length(event)))) != 1L) {
  stop("GSE85916 metadata fields have different lengths")
}
if (any(tolower(sources) != "pancreatic tumor")) {
  stop("Unexpected non-tumour sample in GSE85916")
}

n_expression_rows <- table_end - table_begin - 2L
rm(series_lines)
gc(verbose = FALSE)
expression_table <- read.delim(
  gzfile(series_file),
  header = TRUE,
  sep = "\t",
  skip = table_begin,
  nrows = n_expression_rows,
  check.names = FALSE,
  quote = '"',
  comment.char = "",
  na.strings = c("", "NA", "null")
)
if (!identical(names(expression_table)[-1L], sample_ids)) {
  stop("Series Matrix sample order does not match metadata")
}

# Stream the large platform annotation and retain only the ID and Gene Symbol
# fields. Ambiguous multi-gene probes are excluded, and repeated probes for one
# gene are collapsed by a fixed mean that does not use survival outcomes.
annotation_connection <- file(platform_file, open = "rt")
inside_table <- FALSE
id_index <- NA_integer_
symbol_index <- NA_integer_
probe_id <- character()
gene_symbol <- character()
finished <- FALSE
while (!finished) {
  chunk <- readLines(annotation_connection, n = 2500L, warn = FALSE)
  if (!length(chunk)) break
  for (line in chunk) {
    if (!inside_table) {
      if (identical(line, "!platform_table_begin")) inside_table <- TRUE
      next
    }
    if (is.na(id_index)) {
      header <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      id_index <- match("ID", header)
      symbol_index <- match("Gene Symbol", header)
      if (is.na(id_index) || is.na(symbol_index)) stop("GPL13667 annotation fields missing")
      next
    }
    if (identical(line, "!platform_table_end")) {
      finished <- TRUE
      break
    }
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(fields) >= max(id_index, symbol_index)) {
      probe_id <- c(probe_id, fields[id_index])
      gene_symbol <- c(gene_symbol, fields[symbol_index])
    }
  }
}
close(annotation_connection)

symbol_map <- setNames(trimws(gene_symbol), probe_id)
mapped_symbol <- unname(symbol_map[as.character(expression_table[[1L]])])
keep_probe <- !is.na(mapped_symbol) & nzchar(mapped_symbol) &
  mapped_symbol != "---" & !grepl("///", mapped_symbol, fixed = TRUE)
values <- as.matrix(expression_table[keep_probe, -1L, drop = FALSE])
storage.mode(values) <- "double"
mapped_symbol <- mapped_symbol[keep_probe]
if (any(!is.finite(values))) stop("Non-finite GSE85916 expression values")

probe_groups <- split(seq_along(mapped_symbol), mapped_symbol)
values <- do.call(rbind, lapply(probe_groups, function(i) {
  if (length(i) == 1L) values[i, ] else colMeans(values[i, , drop = FALSE])
}))
rownames(values) <- names(probe_groups)
colnames(values) <- sample_ids

clinical_all <- data.frame(
  sample_id = sample_ids,
  patient_id = sample_ids,
  cohort = "GSE85916",
  platform = "GPL13667 Affymetrix Human Genome U219 Array",
  specimen_class = "primary_tumor",
  time = os_years * 12,
  event = event,
  title = sample_titles,
  source = sources,
  stringsAsFactors = FALSE
)
keep_patient <- is.finite(clinical_all$time) & clinical_all$time > 0 &
  clinical_all$event %in% c(0L, 1L)
clinical <- clinical_all[keep_patient, , drop = FALSE]
values <- values[, clinical$sample_id, drop = FALSE]

if (!identical(colnames(values), clinical$sample_id)) stop("Expression/clinical order failure")
if (anyDuplicated(clinical$patient_id)) stop("Duplicated patient in GSE85916")
if (sum(clinical$event) < 10L) stop("Too few GSE85916 events")

gse85916 <- list(
  expression = values,
  clinical = clinical,
  platform = unique(clinical$platform),
  provenance = list(
    series = series_url,
    platform = platform_url,
    expression = "submitter-normalized Series Matrix",
    probe_collapse = "unambiguous Gene Symbol only; arithmetic mean across probes",
    os_unit = "months (GEO os.year multiplied by 12)"
  )
)

tcga_file <- "data/v2/tcga_cbioportal.rds"
if (!file.exists(tcga_file)) stop("Run scripts/v2/02_prepare_tcga.R first")
tcga <- readRDS(tcga_file)
cohorts <- list(TCGA_PAAD = tcga, GSE85916 = gse85916)

saveRDS(gse85916, "data/v2/gse85916.rds")
saveRDS(cohorts, "data/v2/cohorts.rds")
write.csv(clinical, "results/v2/gse85916_clinical_used.csv", row.names = FALSE)
write.csv(
  data.frame(
    step = c("Series Matrix samples", "Primary tumours", "Positive OS time and known event", "Unique patients used"),
    n = c(length(sample_ids), sum(tolower(sources) == "pancreatic tumor"),
          sum(keep_patient), nrow(clinical))
  ),
  "results/v2/gse85916_sample_flow.csv",
  row.names = FALSE
)

message(
  "GSE85916 prepared: ", ncol(values), " patients, ", sum(clinical$event),
  " events, ", nrow(values), " uniquely mapped genes."
)
