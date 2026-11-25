#!/usr/bin/env Rscript
# Byg en visuel BAM-pool fra en brugerudvalgt liste.
if (sys.nframe() == 0L) {
    script <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(script) != 1L) stop("Run with Rscript.")
    script_dir <- dirname(normalizePath(sub("^--file=", "", script), mustWork = TRUE))
} else {
    script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = TRUE))
}
source(file.path(script_dir, "clinical_core.R"))
source(file.path(script_dir, "reference_core.R"))
run_pool <- function(args) {
opt <- reference_options(args,
    c("--list", "--output-dir", "--condition", "--target-reads", "--seed", "--threads"),
    c("--list", "--output-dir", "--condition"))
requested_condition <- opt[["--condition"]]
condition <- reference_condition(requested_condition)
number <- function(x, default, label) {
    x <- if (is.null(x)) default else x
    if (!grepl("^[0-9]+$", x) || is.na(suppressWarnings(as.integer(x))) || as.integer(x) < 1L) stop_reference("Invalid ", label)
    as.integer(x)
}
target <- number(opt[["--target-reads"]], "20000000", "--target-reads")
seed <- number(opt[["--seed"]], "100", "--seed")
threads <- number(opt[["--threads"]], "2", "--threads")
if (!file.exists(opt[["--list"]])) stop_reference("BAM list not found.")
lines <- trimws(readLines(opt[["--list"]], warn = FALSE))
lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
if (!length(lines)) stop_reference("BAM list is empty.")
bams <- vapply(lines, reference_bam, character(1))
if (anyDuplicated(bams)) stop_reference("BAM list contains duplicate paths.")
out <- opt[["--output-dir"]]
if (file.exists(out)) stop_reference("Output directory already exists: ", out)
if (!dir.create(out, recursive = FALSE)) stop_reference("Cannot create output directory: ", out)
out <- normalizePath(out, mustWork = TRUE)
writeLines("INCOMPLETE", file.path(out, "STATUS.txt"))
for (folder in c("Logs", "subsamples")) if (!dir.create(file.path(out, folder))) stop_reference("Cannot create ", folder)
logs <- file.path(out, "Logs")
writeLines(condition, file.path(out, "condition.txt"))
reference_run("samtools", "--version", file.path(out, "samtools_version.txt"), file.path(logs, "samtools_version.err"))
pool_files <- character(length(bams))
details <- vector("list", length(bams))
# Kontrollér hver BAM og subsample ved behov før merge.
dictionary <- NULL
for (i in seq_along(bams)) {
    prefix <- file.path(logs, sprintf("input_%04d", i))
    sq <- reference_header(bams[i], prefix)
    if (is.null(dictionary)) {
        dictionary <- sq
        writeLines(dictionary, file.path(out, "sequence_dictionary.tsv"))
    } else if (!identical(dictionary, sq)) stop_reference("Sequence dictionary/order differs for: ", bams[i])
    idx <- paste0(prefix, ".idxstats")
    reference_run("samtools", c("idxstats", bams[i]), idx, paste0(prefix, ".idxstats.err"))
    stats <- read.delim(idx, header = FALSE, colClasses = "character")
    if (ncol(stats) != 4L) stop_reference("Invalid idxstats for: ", bams[i])
    mapped <- suppressWarnings(as.numeric(stats[[3L]]))
    if (anyNA(mapped) || any(mapped < 0)) stop_reference("Invalid idxstats for: ", bams[i])
    count <- sum(mapped)
    if (count < 1) stop_reference("No mapped alignments in: ", bams[i])
    fraction <- min(1, target / count)
    actual_fraction <- 1
    pool_files[i] <- bams[i]
    if (fraction < 1) {
        digits <- floor(fraction * 1e8)
        if (digits < 1) stop_reference("Sampling fraction too small for 8-digit precision: ", bams[i])
        actual_fraction <- digits / 1e8
        pool_files[i] <- file.path(out, "subsamples", sprintf("input_%04d.bam", i))
        reference_run("samtools", c("view", "-@", as.character(threads), "-b", "-s",
            paste0(seed, ".", sprintf("%08d", digits)), "-o", pool_files[i], bams[i]),
            paste0(prefix, ".subsample.out"), paste0(prefix, ".subsample.err"))
        reference_run("samtools", c("quickcheck", pool_files[i]), paste0(prefix, ".subsample_check.out"), paste0(prefix, ".subsample_check.err"))
    }
    info <- file.info(bams[i])
    details[[i]] <- data.frame(input_bam = bams[i], condition = condition, mapped_alignments_idxstats = count,
                               target = target, sampled_fraction = actual_fraction, seed = seed,
                               size_bytes = info$size, modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE))
}
# Gem provenance, merge, indeksér og markér først derefter output som fuldført.
write.table(do.call(rbind, details), file.path(out, "inputs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(pool_files, file.path(out, "merge_inputs.txt"))
merged <- file.path(out, paste0("reference_", requested_condition, ".bam"))
reference_run("samtools", c("merge", "-@", as.character(threads), "-b", file.path(out, "merge_inputs.txt"), merged),
              file.path(logs, "merge.out"), file.path(logs, "merge.err"))
reference_run("samtools", c("quickcheck", merged), file.path(logs, "merged_check.out"), file.path(logs, "merged_check.err"))
reference_run("samtools", c("index", "-@", as.character(threads), merged), file.path(logs, "index.out"), file.path(logs, "index.err"))
writeLines(c("COMPLETE - visual reference only; not individual controls.",
             "Inputs and treatment membership must be independently verified.",
             "quickcheck is not a full integrity scan; idxstats target is approximate."),
           file.path(out, "STATUS.txt"))
cat("Pool ready: ", merged, "\n", sep = "")
}
if (sys.nframe() == 0L) run_pool(commandArgs(trailingOnly = TRUE))
