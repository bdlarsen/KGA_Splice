#!/usr/bin/env Rscript
if (!requireNamespace("xml2", quietly = TRUE)) stop("xml2 is required")
args <- grep("^--file=", commandArgs(), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)))
source(file.path(base, "scripts", "clinical_core.R"))
tmp <- tempfile("igv_session_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
input <- file.path(tmp, "input.xml")
output <- file.path(tmp, "output.xml")
writeLines(c(
    '<Session genome="hg38" locus="chr1:1-100"><Resources>',
    '<Resource path="/drift/RNA/sample&amp;test.bam"/>',
    '<Resource path="/drift/RNA/out/variant.bed"/>',
    '<Resource path="/drift-other/control.bam"/>',
    '</Resources><Panel>',
    '<Track id="/drift/RNA/sample&amp;test.bam" name="/drift stays in name" visible="false"/>',
    '<Track id="/drift/RNA/sample&amp;test.bam_coverage"/>',
    '<Track id="/drift/RNA/sample&amp;test.bam_junctions"/>',
    '<Track id="/drift-other/control.bam"/>',
    '</Panel></Session>'
), input)
before <- readLines(input)
warning_text <- NULL
withCallingHandlers(write_vswarehouse_igv_session(input, output), warning = function(w) {
    warning_text <<- conditionMessage(w)
    invokeRestart("muffleWarning")
})
xml <- xml2::read_xml(output)
paths <- xml2::xml_attr(xml2::xml_find_all(xml, ".//Resource"), "path")
tracks <- xml2::xml_find_all(xml, ".//Track")
ids <- xml2::xml_attr(tracks, "id")
stopifnot(identical(readLines(input), before),
          identical(paths, c("/home/ghuser/Workspace/drift/RNA/sample&test.bam",
                             "/home/ghuser/Workspace/drift/RNA/out/variant.bed",
                             "/drift-other/control.bam")),
          identical(ids, c(paste0(paths[1L], c("", "_coverage", "_junctions")), paths[3L])),
          xml2::xml_attr(tracks[1L], "name") == "/drift stays in name",
          xml2::xml_attr(tracks[1L], "visible") == "false",
          xml2::xml_attr(xml2::xml_root(xml), "genome") == "hg38",
          grepl("/drift-other/control.bam", warning_text, fixed = TRUE))
cli_output <- file.path(tmp, "cli.xml")
script <- file.path(base, "scripts", "make_vswarehouse_session.R")
rscript <- file.path(R.home("bin"), "Rscript")
result <- system2(rscript, c(shQuote(script), shQuote(input), shQuote(cli_output)),
                  stdout = TRUE, stderr = TRUE)
stopifnot(is.null(attr(result, "status")), file.exists(cli_output),
          identical(xml2::xml_attr(xml2::xml_find_all(xml2::read_xml(cli_output), ".//Resource"), "path"), paths))
before_cli <- readLines(cli_output)
result <- suppressWarnings(system2(rscript, c(shQuote(script), shQuote(input), shQuote(cli_output)),
                                   stdout = TRUE, stderr = TRUE))
stopifnot(!is.null(attr(result, "status")), attr(result, "status") != 0L,
          any(grepl("Output already exists", result, fixed = TRUE)),
          identical(readLines(cli_output), before_cli))
writeLines("<Session><Resources/></Session>", input)
error <- tryCatch(write_vswarehouse_igv_session(input, output), error = identity)
stopifnot(inherits(error, "error"), grepl("no resource paths", conditionMessage(error), fixed = TRUE))
cat("IGV session tests passed\n")
