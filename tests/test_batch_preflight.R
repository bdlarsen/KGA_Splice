#!/usr/bin/env Rscript
args <- grep("^--file=", commandArgs(), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)))
script <- file.path(base, "scripts", "batch_preflight.R")
tmp <- tempfile("batch_preflight_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
for (d in c("out", "refs", "bam")) dir.create(file.path(tmp, d))
for (name in c("G150-123-B-GENE_x.bam", "G150-123-A-GENE_x.bam", "G150-456_x.bam")) {
    bam <- file.path(tmp, "bam", name)
    file.create(bam, paste0(bam, ".bai"))
}
file.create(file.path(tmp, "refs", "pool.bam"), file.path(tmp, "refs", "pool.bam.bai"))
module <- file.path(tmp, "module.py")
writeLines("# test", module)
manifest <- file.path(tmp, "manifest.txt")
run_test <- function(rows) {
    writeLines(rows, manifest)
    output <- tempfile(tmpdir = tmp)
    status <- system2(file.path(R.home("bin"), "Rscript"),
        args = shQuote(c(script, file.path(tmp, "out"), file.path(tmp, "refs"), module)),
        stdin = manifest, stdout = output, stderr = output)
    list(status = status, output = paste(readLines(output, warn = FALSE), collapse = "\n"))
}
valid <- c(paste("123", "17-100-A-G", "GENE", "NM_1.1", file.path(tmp, "bam"), "paired", sep = "|"),
           paste("456", "17-100-A-G", "GENE", "NM_1.1", file.path(tmp, "bam"), "single", sep = "|"))
result <- run_test(valid)
stopifnot(result$status == 0L, grepl("Preflight passed", result$output, fixed = TRUE))
unlink(file.path(tmp, "bam", "G150-456_x.bam.bai"))
result <- run_test(valid)
stopifnot(result$status != 0L, grepl("Missing BAM index", result$output, fixed = TRUE),
          !grepl("Preflight passed", result$output, fixed = TRUE))
result <- run_test("invalid")
stopifnot(result$status != 0L, grepl("Invalid manifest record", result$output, fixed = TRUE))
cat("Batch preflight tests passed\n")
