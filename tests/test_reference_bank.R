#!/usr/bin/env Rscript
args <- grep("^--file=", commandArgs(), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)))
source(file.path(base, "scripts", "reference_bank.R"))
expect_error <- function(expr, expected) {
    err <- tryCatch({ force(expr); NULL }, error = identity)
    if (is.null(err) || !grepl(expected, conditionMessage(err), fixed = TRUE)) stop("Expected error: ", expected)
}
tmp <- tempfile("bank_test_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
fake_tool <- file.path(tmp, "regtools")
writeLines("mock version", fake_tool)
bank_regtools <- function() fake_tool
reference_run <- function(command, args, stdout, stderr) {
    writeLines(character(), stderr)
    if (command == "samtools" && args[1L] == "view") {
        writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:1000"), stdout)
    } else if (command == "regtools") {
        if (basename(tail(args, 1L)) == "two.bam") {
            writeLines("chr1\t140\t210\tJ1\t3\t+\t140\t210\t0\t2\t10,10\t0,60", stdout)
        } else {
            writeLines("chr1\t140\t310\tJ2\t4\t+\t140\t310\t0\t2\t10,10\t0,160", stdout)
        }
    } else writeLines(character(), stdout)
}
paths <- file.path(tmp, c("one.bam", "two.bam", "three.bam"))
for (p in paths) { file.create(p); file.create(paste0(p, ".bai")) }
manifest <- file.path(tmp, "manifest.tsv")
dat <- data.frame(sample_id = c("s1", "s2"), person_id = c("p1", "p2"), condition = c("minus", "minus"),
                  bam = paths[1:2], qc = "pass")
write.table(dat, manifest, sep = "\t", quote = FALSE, row.names = FALSE)
bank <- file.path(tmp, "bank")
init_bank(list("--bank" = bank, "--reference-id" = "GRCh38_test"))
add_samples(list("--bank" = bank, "--manifest" = manifest))
legacy_path <- file.path(bank, "samples", "s1", "metadata.rds")
legacy <- readRDS(legacy_path)
legacy$condition <- "minus"
saveRDS(legacy, legacy_path)
dat$condition <- "nt"
write.table(dat, manifest, sep = "\t", quote = FALSE, row.names = FALSE)
add_samples(list("--bank" = bank, "--manifest" = manifest))
stopifnot(nrow(bank_samples(bank)) == 2L, all(bank_samples(bank)$condition == "nt"))
extra <- data.frame(sample_id = "s3", person_id = "p3", condition = "plus", bam = paths[3], qc = "pass")
write.table(extra, manifest, sep = "\t", quote = FALSE, row.names = FALSE)
add_samples(list("--bank" = bank, "--manifest" = manifest))
stopifnot(nrow(bank_samples(bank)) == 3L)
bad <- extra
bad$sample_id <- "s4"
bad$bam <- paths[1]
write.table(bad, manifest, sep = "\t", quote = FALSE, row.names = FALSE)
expect_error(add_samples(list("--bank" = bank, "--manifest" = manifest)), "Person/condition already exists")
controls <- file.path(tmp, "controls.txt")
writeLines(c("s1", "s2"), controls)
patient_bed <- file.path(tmp, "patient.bed")
writeLines("chr1\t140\t210\tJ1\t10\t+\t140\t210\t0\t2\t10,10\t0,60", patient_bed)
options <- list("--bank" = bank, "--patient" = "p0", "--condition" = "minus",
                "--junctions" = patient_bed, "--controls" = controls, "--region" = "chr1:100-350",
                "--output-dir" = file.path(tmp, "compare"))
compare_samples(options)
stopifnot(identical(reference_condition(c("nt", "puro", "minus", "plus")), c("nt", "puro", "nt", "puro")))
results <- read.delim(file.path(tmp, "compare", "comparison.tsv"))
stopifnot(nrow(results) == 4L,
          results$control_reads[results$control_sample == "s2" & results$intron_end_0 == 200] == 3L,
          results$patient_status[results$control_sample == "s1" & results$intron_end_0 == 300] == "not_observed_coverage_unassessed",
          is.na(results$control_reads[results$control_sample == "s1" & results$intron_end_0 == 200]))
empty_bed <- file.path(tmp, "empty.bed")
file.create(empty_bed)
empty_options <- options
empty_options[["--junctions"]] <- empty_bed
empty_options[["--region"]] <- "chr1:700-800"
empty_options[["--output-dir"]] <- file.path(tmp, "empty")
compare_samples(empty_options)
stopifnot(nrow(read.delim(file.path(tmp, "empty", "comparison.tsv"))) == 0L)
options[["--patient"]] <- "p1"
options[["--output-dir"]] <- file.path(tmp, "self")
expect_error(compare_samples(options), "own control")
options[["--patient"]] <- "p0"
options[["--condition"]] <- "plus"
expect_error(compare_samples(options), "condition does not match")
options[["--condition"]] <- "minus"
options[["--output-dir"]] <- file.path(tmp, "corrupt")
inventory <- readRDS(file.path(bank, "samples", "s1", "inventory.rds"))
writeLines("damaged", file.path(bank, "samples", "s1", inventory$filename[1L]))
expect_error(compare_samples(options), "Missing or damaged chromosome file")
writeLines("changed", fake_tool)
expect_error(add_samples(list("--bank" = bank, "--manifest" = manifest)), "regtools differs")
cat("Reference bank mocked tests passed\n")
