#!/usr/bin/env Rscript
args <- grep("^--file=", commandArgs(), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)))
source(file.path(base, "scripts", "reference_pool.R"))
expect_error <- function(expr, expected) {
    err <- tryCatch({ force(expr); NULL }, error = identity)
    if (is.null(err) || !grepl(expected, conditionMessage(err), fixed = TRUE)) stop("Expected error: ", expected)
}
tmp <- tempfile("pool_test_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
reference_run <- function(command, args, stdout, stderr) {
    writeLines(character(), stderr)
    if (args[1L] == "--version") writeLines("samtools mock", stdout)
    else if (args[1L] == "view" && args[2L] == "-H")
        writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:1000"), stdout)
    else if (args[1L] == "idxstats")
        writeLines(c("chr1\t1000\t100\t0", "*\t0\t0\t0"), stdout)
    else if (args[1L] == "view") {
        stopifnot(args[which(args == "-s") + 1L] == "100.50000000")
        writeLines("mock bam", args[which(args == "-o") + 1L])
        writeLines(character(), stdout)
    } else if (args[1L] == "merge") {
        writeLines("mock pool", tail(args, 1L))
        writeLines(character(), stdout)
    } else if (args[1L] == "index") {
        writeLines("mock index", paste0(tail(args, 1L), ".bai"))
        writeLines(character(), stdout)
    } else writeLines(character(), stdout)
}
paths <- file.path(tmp, c("one.bam", "two.bam"))
for (p in paths) { file.create(p); file.create(paste0(p, ".bai")) }
list <- file.path(tmp, "list.txt")
writeLines(paths, list)
out <- file.path(tmp, "pool")
parameters <- c("--list", list, "--output-dir", out, "--condition", "minus", "--target-reads", "50")
run_pool(parameters)
info <- read.delim(file.path(out, "inputs.tsv"))
stopifnot(nrow(info) == 2L, all(info$sampled_fraction == 0.5),
          file.exists(file.path(out, "reference_minus.bam.bai")),
          startsWith(readLines(file.path(out, "STATUS.txt"))[1L], "COMPLETE"))
expect_error(run_pool(parameters), "already exists")
nt_out <- file.path(tmp, "pool_nt")
run_pool(c("--list", list, "--output-dir", nt_out, "--condition", "nt", "--target-reads", "50"))
stopifnot(file.exists(file.path(nt_out, "reference_nt.bam.bai")),
          identical(readLines(file.path(nt_out, "condition.txt")), "nt"))
previous_run <- reference_run
reference_run <- function(command, args, stdout, stderr) {
    if (args[1L] == "merge") stop_reference("mock merge failure")
    previous_run(command, args, stdout, stderr)
}
failed <- file.path(tmp, "failed")
expect_error(run_pool(c("--list", list, "--output-dir", failed, "--condition", "minus", "--target-reads", "50")), "mock merge failure")
stopifnot(identical(readLines(file.path(failed, "STATUS.txt")), "INCOMPLETE"))
writeLines(c(paths[1], paths[1]), list)
expect_error(run_pool(c("--list", list, "--output-dir", file.path(tmp, "duplicate"), "--condition", "plus")), "duplicate paths")
cat("Reference pool mocked tests passed\n")
