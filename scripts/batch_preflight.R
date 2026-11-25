#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript batch_preflight.R.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
source(file.path(script_dir, "clinical_core.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop_input("Usage: Rscript batch_preflight.R OUTPUT_PARENT REFERENCE_DIR SAI_MODULE < manifest")
output <- args[1L]
refs <- args[2L]
module <- args[3L]
con <- file("stdin", open = "r")
rows <- readLines(con, warn = FALSE)
close(con)
errors <- character()
check <- function(label, expr) {
    tryCatch(force(expr), error = function(e) {
        errors <<- c(errors, paste(label, conditionMessage(e)))
    })
}
check("Output:", {
    if (!dir.exists(output) || file.access(output, 2L) != 0L)
        stop("Output parent missing or not writable: ", output)
})
check("Reference BAM:", {
    if (!dir.exists(refs)) stop("Directory missing: ", refs)
    controls <- list.files(refs, pattern = "[.]bam$", full.names = TRUE)
    if (!length(controls)) stop("No reference BAM files: ", refs)
    for (bam in controls) {
        if (file.access(bam, 4L) != 0L) stop("Reference BAM not readable: ", bam)
        check_bam_index(bam)
    }
})
check("SAI module:", {
    if (!file.exists(module) || dir.exists(module) || file.info(module)$size == 0L ||
        file.access(module, 4L) != 0L) stop("Module missing, empty or unreadable: ", module)
})
if (!length(rows)) errors <- c(errors, "No patient records provided.")
for (row in rows) {
    fields <- strsplit(row, "|", fixed = TRUE)[[1L]]
    if (length(fields) != 6L || any(!nzchar(fields))) {
        errors <- c(errors, "Invalid manifest record: expected ID|variant|gene|transcript|BAM directory|single or paired.")
        next
    }
    id <- fields[1L]
    bamdir <- fields[5L]
    mode <- fields[6L]
    if (!mode %in% c("single", "paired")) {
        errors <- c(errors, paste("Invalid sample mode for", id))
        next
    }
    ids <- if (mode == "single") id else paste0(id, c("-B", "-A"))
    for (sample in ids) check(sample, {
        parse_sample_ids(sample)
        if (!dir.exists(bamdir)) stop("BAM directory missing: ", bamdir)
        bam <- find_sample_bam(bamdir, sample)
        if (file.access(bam, 4L) != 0L) stop("BAM not readable: ", bam)
        check_bam_index(bam)
        cat("OK:", sample, bam, "\n")
    })
}
if (length(errors)) {
    cat("\nPreflight errors:\n", paste(errors, collapse = "\n"), "\n", sep = "")
    quit(save = "no", status = 1L)
}
cat("\nPreflight passed; BAM files and indexes exist. Contents and sample identity were not verified.\n")
