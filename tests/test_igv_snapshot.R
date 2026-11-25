#!/usr/bin/env Rscript
args <- grep("^--file=", commandArgs(), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)))
source(file.path(base, "scripts", "clinical_core.R"))
source(file.path(base, "scripts", "igv_snapshot_core.R"))
tmp <- tempfile("igv_snapshot_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
for (d in c("Logs", "Report", "IGV")) dir.create(file.path(tmp, d))
fasta <- file.path(tmp, "genome.fa")
writeLines(">chr17\nA", fasta)
session <- file.path(tmp, "IGV", "paired_IGV.xml")
writeLines('<Session genome="hg38" locus="chr17:1-100"/>', session)
launcher <- file.path(tmp, if (.Platform$OS.type == "windows") "igv.exe" else "igv.sh")
if (.Platform$OS.type == "windows") {
    stopifnot(file.copy(file.path(R.home("bin"), "Rscript.exe"), launcher))
} else {
    writeLines("#!/bin/sh", launcher)
    Sys.chmod(launcher, "0755")
}
finder <- function(cmd) "fake-xvfb"
output <- file.path(tmp, "Report", "igv_snapshot.svg")
mock <- function(command, args, stdout, stderr, timeout) {
    stopifnot(command == "fake-xvfb", timeout == 180L)
    batch <- readLines(file.path(tmp, "Logs", "igv_snapshot.batch"))
    stopifnot(identical(batch[1L], paste("genome", fasta)), "goto chr17:1-100" %in% batch,
              "snapshot igv_snapshot.svg" %in% batch)
    local <- xml2::read_xml(file.path(tmp, "Logs", "igv_snapshot_session.xml"))
    stopifnot(xml2::xml_attr(xml2::xml_root(local), "genome") == normalizePath(fasta))
    writeLines(character(), stdout)
    writeLines(character(), stderr)
    writeLines('<svg xmlns="http://www.w3.org/2000/svg"><text>test</text></svg>', output)
    0L
}
stopifnot(render_igv_snapshot(NULL, session, fasta, tmp, "chr17:1-100")$status == "not_enabled")
result <- render_igv_snapshot(launcher, session, fasta, tmp, "chr17:1-100", finder, mock)
stopifnot(result$status == "available", result$detail == "chr17:1-100")
result <- render_igv_snapshot(launcher, session, fasta, tmp, "chr17:1-100", finder,
                              function(command, args, stdout, stderr, timeout) {
                                  writeLines(character(), stdout)
                                  writeLines("SEVERE [test] Missing track", stderr)
                                  0L
                              })
stopifnot(result$status == "failed", grepl("reported an error", result$detail, fixed = TRUE))
unlink(output)
result <- render_igv_snapshot(launcher, session, fasta, tmp, "chr17:1-100", finder,
                              function(...) 1L)
stopifnot(result$status == "failed", grepl("failed or timed out", result$detail, fixed = TRUE))
result <- render_igv_snapshot(launcher, session, fasta, tmp, "chr17:1-100", finder,
                              function(command, args, stdout, stderr, timeout) {
                                  writeLines(character(), stdout)
                                  writeLines(character(), stderr)
                                  0L
                              })
stopifnot(result$status == "failed", grepl("did not create", result$detail, fixed = TRUE))
result <- render_igv_snapshot(launcher, session, fasta, tmp, "chr17:1-100", function(...) "", mock)
stopifnot(result$status == "failed", grepl("xvfb-run", result$detail, fixed = TRUE))
cat("IGV snapshot tests passed\n")
