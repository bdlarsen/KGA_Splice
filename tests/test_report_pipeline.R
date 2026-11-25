#!/usr/bin/env Rscript
if (!requireNamespace("xml2", quietly = TRUE)) stop("xml2 is required for pipeline test")
args <- grep("^--file=", commandArgs(), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)))
tmp <- tempfile("report_pipeline_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
for (d in c("bam", "vep", "refs", "out", "vcf_input", "ab", "prefixed", "single", "igv_failed", "indel", "failed", "failed_plot", "sai", "sai_missing", "sai_failed")) dir.create(file.path(tmp, d))
for (id in c("1", "2")) {
    file.create(file.path(tmp, "bam", paste0(id, ".bam")))
    file.create(file.path(tmp, "bam", paste0(id, ".bam.bai")))
}
for (id in c("040966-2498-B", "040966-2498-A")) {
    file.create(file.path(tmp, "bam", paste0(id, "_illumina.bam")))
    file.create(file.path(tmp, "bam", paste0(id, "_illumina.bam.bai")))
}
for (id in c("G150-single", "G150-pair-B-PTEN", "G150-pair-A-PTEN")) {
    file.create(file.path(tmp, "bam", paste0(id, "_illumina.bam")))
    file.create(file.path(tmp, "bam", paste0(id, "_illumina.bam.bai")))
}
file.create(file.path(tmp, "refs", "pool.bam"))
file.create(file.path(tmp, "refs", "pool.bam.bai"))
fa <- file.path(tmp, "genome.fa")
writeLines(">chr17\nA", fa)
writeLines("chr17\t1000\t7\t1\t2", paste0(fa, ".fai"))
file.create(file.path(tmp, "genome.dict"))
module <- file.path(tmp, "sai10k_predictions.py")
writeLines("# simulated module", module)
gtf <- file.path(tmp, "gene.gtf")
writeLines(c('chr17\ttest\texon\t101\t150\t.\t+\t.\tgene_id "GENE"; transcript_id "NM_1.1";',
             'chr17\ttest\texon\t201\t250\t.\t+\t.\tgene_id "GENE"; transcript_id "NM_1.1";',
             'chr17\ttest\tCDS\t111\t150\t.\t+\t0\tgene_id "GENE"; transcript_id "NM_1.1";',
             'chr17\ttest\tCDS\t201\t240\t.\t+\t0\tgene_id "GENE"; transcript_id "NM_1.1";'), gtf)

run_one <- function(var, gt, output, fail_coverage = FALSE, fail_sashimi = FALSE, enable_sai = FALSE,
                    missing_sai = FALSE, fail_sai = FALSE, sample_ids = "1,2", fail_igv = FALSE,
                    reserved = FALSE) {
    env <- new.env(parent = globalenv())
    sai_calls <- 0L
    env$source <- function(file, ...) base::source(file, local = env)
    env$commandArgs <- function(trailingOnly = FALSE) {
        if (!trailingOnly) return(paste0("--file=", file.path(base, "scripts", "ClinV1.R")))
        c("-variant", var, "-gene", "GENE", "-tx", "NM_1.1",
          if (reserved) "-project-dir" else "-p", output,
          "-bam", file.path(tmp, "bam"), "-s", sample_ids, "-g", fa, "-gtf", gtf,
          "-vep", file.path(tmp, "vep"), "-gt", gt, "-refs", file.path(tmp, "refs"),
          if (enable_sai) c("-sai-module", module))
    }
    mock_run <- function(command, args, stdout, stderr) {
        writeLines(character(), stderr)
        if (command == "samtools" && args[1L] == "faidx") {
            writeLines(c(">chr17", if (var == "17-150-A-G" || grepl("\\.vcf$", var)) "A" else "AC"), stdout)
        } else if (command == "conda") {
            tool <- args[5L]
            if (tool == "python") {
                if (!("--check" %in% args)) {
                    if (fail_sai) stop("simulated SAI predictor failure")
                    data <- data.frame(transcript = "NM_1.1", status = "predicted", effect = "partial_exon_deletion",
                        description = "Donor shift", consequence = "exon deletion", size_bp = "1",
                        frameshift = "True", affects_coding = "True")
                    write.table(data, args[which(args == "--output") + 1L], sep = "\t", quote = FALSE, row.names = FALSE)
                    writeLines("{}", args[which(args == "--json") + 1L])
                }
                writeLines(character(), stdout)
                return(invisible(NULL))
            }
            if (tool == "sashimipy") {
                if ("--help" %in% args) {
                    writeLines("sashimipy help", stdout)
                    return(invisible(NULL))
                }
                stopifnot(args[which(args == "-t") + 1L] == "20", "--density-by-strand" %in% args,
                          args[which(args == "--rs") + 1L] == "XS:A:+")
                density <- read.delim(args[which(args == "--density") + 1L], header = FALSE)
                stopifnot(nrow(density) == length(strsplit(sample_ids, ",", fixed = TRUE)[[1L]]) + 1L,
                          all(density[[5L]] == "frf"))
                if (fail_sashimi) stop("Simulated sashimipy failure")
                writeLines("%PDF-1.4", args[which(args == "-o") + 1L])
                writeLines(character(), stdout)
                return(invisible(NULL))
            }
            out <- args[which(args == if (tool == "vep") "-o" else "-O") + 1L]
            vcf <- args[which(args == if (tool == "vep") "-i" else "-I") + 1L]
            input <- readLines(vcf)
            input <- input[!startsWith(input, "#")]
            if (tool == "vep") {
                header <- '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: SYMBOL|Consequence|Feature">'
                info <- "CSQ=GENE|splice_donor_variant|NM_1.1"
            } else {
                header <- "##fileformat=VCFv4.2"
                symbol <- if (args[which(args == "-D") + 1L] == "4999") "RefSeqTx-NM_1.1" else "GENE"
                if (symbol == "RefSeqTx-NM_1.1") sai_calls <<- sai_calls + 1L
                info <- paste0("SpliceAI=", strsplit(input, "\t", fixed = TRUE)[[1L]][5L], "|", symbol, "|0.1|0.0|0.0|0.0|1|0|0|0")
                if (symbol == "RefSeqTx-NM_1.1" && missing_sai) info <- "."
            }
            fields <- strsplit(input, "\t", fixed = TRUE)[[1L]]
            fields[8L] <- info
            writeLines(c(header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", paste(fields, collapse = "\t")), out)
            writeLines(character(), stdout)
        } else if (command == "samtools" && args[1L] == "view" && args[2L] == "-H") {
            writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr17\tLN:1000"), stdout)
        } else if (command == "samtools" && args[1L] == "idxstats") {
            n <- if (basename(args[2L]) == "1.bam") 1000000 else 2000000
            writeLines(c(paste("chr17", 1000, sprintf("%.0f", n), 0, sep = "\t"), "*\t0\t0\t0"), stdout)
        } else if (command == "samtools" && args[1L] == "view" && args[2L] == "-c") {
            writeLines("2", stdout)
        } else if (command == "samtools" && args[1L] == "depth") {
            positions <- if (var == "17-150-A-G" || grepl("\\.vcf$", var)) 150L else 150:151
            writeLines(paste("chr17", positions, 20, sep = "\t"), stdout)
        } else if (command == "sambamba" && args[1L] == "view") {
            writeLines("10", stdout)
        } else if (command == "bash") {
            if (fail_coverage) stop("Simulated coverage failure")
            writeLines("chr17\t100\t110\t2", args[4L])
            writeLines("chr17\t100\t110\t-2", args[5L])
            writeLines(character(), stdout)
        } else if (command == "regtools") {
            writeLines("chr17\t140\t210\tJ1\t10\t+\t140\t210\t0\t2\t10,10\t0,60", stdout)
        } else if (command == "gatk" && args[1L] == "ASEReadCounter") {
            out <- args[which(args == "-O") + 1L]
            writeLines(c("contig,position,refAllele,altAllele,refCount,altCount", "chr17,150,A,G,12,8"), out)
            writeLines(character(), stdout)
        } else writeLines(character(), stdout)
    }
    exprs <- parse(file.path(base, "scripts", "ClinV1.R"))
    for (expr in exprs) {
        if (grepl('^opt <- parse_options', paste(deparse(expr), collapse = ""))) {
            env$check_required_commands <- function(...) invisible(NULL)
            env$run_checked <- mock_run
            if (fail_igv) env$render_igv_snapshot <- function(...) list(status = "failed", detail = "simulated IGV failure")
        }
        eval(expr, env)
    }
    if (enable_sai) stopifnot(sai_calls == 1L)
    if (reserved) return(normalizePath(output))
    ids <- strsplit(sample_ids, ",", fixed = TRUE)[[1L]]
    file.path(output, if (length(ids) == 1L) paste0("GENE_", ids[1L], "_treatment_unknown") else
        paste0("GENE_", ids[1L], "_nt_", ids[2L], "_puro"))
}
reserved_project <- file.path(tmp, "1_GENE_unique")
dir.create(reserved_project)
dir.create(file.path(reserved_project, "tmp"))
reserved_result <- run_one("17-150-A-G", "het", reserved_project, reserved = TRUE)
stopifnot(file.exists(file.path(reserved_result, "STATUS.txt")),
          file.exists(file.path(reserved_result, "IGV", "paired_IGV.xml")),
          dir.exists(file.path(reserved_result, "tmp")),
          !dir.exists(file.path(reserved_result, "GENE_1_nt_2_puro")))
ab <- run_one("17-150-A-G", "het", file.path(tmp, "ab"), sample_ids = "040966-2498-B,040966-2498-A")
ab_samples <- read.delim(file.path(ab, "samples.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(ab_samples$label[1:2], c("040966-2498-B_nt", "040966-2498-A_puro")),
          identical(basename(ab_samples$input_bam[1:2]),
                    c("040966-2498-B_illumina.bam", "040966-2498-A_illumina.bam")),
          identical(ab_samples$condition[1:2], c("nt", "puro")))
prefixed <- run_one("17-150-A-G", "het", file.path(tmp, "prefixed"), sample_ids = "pair-B,pair-A")
prefixed_samples <- read.delim(file.path(prefixed, "samples.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(basename(prefixed_samples$input_bam[1:2]),
                    c("G150-pair-B-PTEN_illumina.bam", "G150-pair-A-PTEN_illumina.bam")),
          identical(prefixed_samples$condition[1:2], c("nt", "puro")))
single <- run_one("17-150-A-G", "het", file.path(tmp, "single"), sample_ids = "single")
single_samples <- read.delim(file.path(single, "samples.tsv"), stringsAsFactors = FALSE)
single_report <- paste(readLines(file.path(single, "Report", "index.html"), warn = FALSE), collapse = "")
single_status <- paste(readLines(file.path(single, "STATUS.txt")), collapse = "")
single_xml <- xml2::read_xml(file.path(single, "IGV", "paired_IGV.xml"))
single_names <- xml2::xml_attr(xml2::xml_find_all(single_xml, ".//Track"), "name")
stopifnot(identical(single_samples$condition, c("treatment_unknown", "control_unspecified")),
          basename(single_samples$input_bam[1]) == "G150-single_illumina.bam",
          grepl("Enkelt patientprøve med ukendt behandling", single_report, fixed = TRUE),
          grepl("no nt/puro comparison", single_status, fixed = TRUE),
          "single (treatment unknown) BAM" %in% single_names,
          !any(grepl("[+-]puro", single_names)))
project <- run_one("17-150-A-G", "het", file.path(tmp, "out"))
igv_failed <- run_one("17-150-A-G", "het", file.path(tmp, "igv_failed"), fail_igv = TRUE)
stopifnot(grepl("IGV snapshot: failed simulated IGV failure",
                paste(readLines(file.path(igv_failed, "STATUS.txt")), collapse = ""), fixed = TRUE),
          grepl("simulated IGV failure",
                paste(readLines(file.path(igv_failed, "Report", "index.html"), warn = FALSE), collapse = ""), fixed = TRUE))
vcf_file <- file.path(tmp, "input.vcf")
writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
             "chr17\t150\t.\tA\tG\t.\t.\t."), vcf_file)
vcf_project <- run_one(vcf_file, "het", file.path(tmp, "vcf_input"))
stopifnot(grepl("chr17-150-A-G", paste(readLines(file.path(vcf_project, "variant_input.txt")), collapse = ""), fixed = TRUE),
          file.exists(file.path(vcf_project, "Report", "index.html")))
stopifnot(grepl("counts_available_exploratory", paste(readLines(file.path(project, "ASE", "allele_counts.tsv")), collapse = "")),
          file.exists(file.path(project, "Report", "index.html")),
          file.exists(file.path(project, "Report", "sashimi.pdf")))
junctions <- read.delim(file.path(project, "Junctions", "junction_counts.tsv"), check.names = FALSE)
report <- paste(readLines(file.path(project, "Report", "index.html"), warn = FALSE), collapse = "")
stopifnot(junctions$left_boundary_exon == "exon 1 (exact)",
          junctions$right_boundary_exon == "exon 2 (exact)",
          junctions$junction_id == "J1", junctions$variant_distance_bp == 1,
          junctions$near_variant,
          junctions[["1_nt_per_million_mapped"]] == 10,
          junctions[["2_puro_per_million_mapped"]] == 5,
          file.exists(file.path(project, "Report", "junction_normalized.svg")),
          grepl("exon 1 (exact)", report, fixed = TRUE),
          grepl("exon 2 (exact)", report, fixed = TRUE),
          grepl('<tr class="near-variant">', report, fixed = TRUE))
xml <- xml2::read_xml(file.path(project, "IGV", "paired_IGV.xml"))
tracks <- xml2::xml_find_all(xml, ".//Track")
names <- xml2::xml_attr(tracks, "name")
alignments <- tracks[xml2::xml_attr(tracks, "clazz") == "org.broad.igv.sam.AlignmentTrack"]
stopifnot(length(alignments) == 6L, all(xml2::xml_attr(alignments, "visible") == "false"))
other_tracks <- tracks[xml2::xml_attr(tracks, "clazz") %in%
    c("org.broad.igv.sam.CoverageTrack", "org.broad.igv.sam.SpliceJunctionTrack")]
stopifnot(length(other_tracks) == 12L, all(is.na(xml2::xml_attr(other_tracks, "visible"))))
stopifnot("1 -puro BAM" %in% names, "2 +puro BAM" %in% names,
          "1 -puro + strand" %in% names, "1 -puro - strand" %in% names,
          all(c("1 -puro raw BAM", "2 +puro raw BAM", "Control pool (treatment unknown) raw BAM") %in% names),
          xml2::xml_attr(tracks[names == "Variant"], "color") == "255,0,0",
          xml2::xml_attr(tracks[names == "1 -puro BAM"], "color") != xml2::xml_attr(tracks[names == "2 +puro BAM"], "color"),
          xml2::xml_attr(tracks[names == "1 -puro BAM"], "color") != xml2::xml_attr(tracks[names == "Control pool (treatment unknown) BAM"], "color"))
resources <- normalizePath(xml2::xml_attr(xml2::xml_find_all(xml, ".//Resource"), "path"), winslash = "/", mustWork = TRUE)
expected_bams <- normalizePath(c(file.path(tmp, "bam", c("1.bam", "2.bam")), file.path(tmp, "refs", "pool.bam")),
                               winslash = "/", mustWork = TRUE)
stopifnot(all(expected_bams %in% resources))
indel <- run_one("17-150-AC-A", "het", file.path(tmp, "indel"))
stopifnot(grepl("indel_not_supported", paste(readLines(file.path(indel, "ASE", "allele_counts.tsv")), collapse = "")))
sai <- run_one("17-150-AC-A", "het", file.path(tmp, "sai"), enable_sai = TRUE)
stopifnot(grepl("RefSeqTx-NM_1.1", paste(readLines(file.path(sai, "Report", "index.html")), collapse = ""), fixed = TRUE),
          grepl("score-low", paste(readLines(file.path(sai, "Report", "index.html")), collapse = ""), fixed = TRUE))
stopifnot(grepl("prediction_available_exploratory", paste(readLines(file.path(sai, "STATUS.txt")), collapse = ""), fixed = TRUE),
          grepl("partial_exon_deletion", paste(readLines(file.path(sai, "Report", "index.html")), collapse = ""), fixed = TRUE),
          file.exists(file.path(sai, "SAI10k", "predictions.json")))
sai_missing <- run_one("17-150-A-G", "het", file.path(tmp, "sai_missing"), enable_sai = TRUE, missing_sai = TRUE)
stopifnot(grepl("spliceai_annotation_unavailable", paste(readLines(file.path(sai_missing, "STATUS.txt")), collapse = ""), fixed = TRUE))
err_sai <- tryCatch(run_one("17-150-AC-A", "het", file.path(tmp, "sai_failed"), enable_sai = TRUE, fail_sai = TRUE), error = identity)
stopifnot(inherits(err_sai, "error"), grepl("simulated SAI predictor failure", conditionMessage(err_sai), fixed = TRUE),
          startsWith(readLines(file.path(tmp, "sai_failed", "GENE_1_nt_2_puro", "STATUS.txt"))[1L], "INCOMPLETE"))
error <- tryCatch(run_one("17-150-A-G", "het", file.path(tmp, "failed"), TRUE), error = identity)
failed_status <- file.path(tmp, "failed", "GENE_1_nt_2_puro", "STATUS.txt")
stopifnot(inherits(error, "error"), grepl("Simulated coverage failure", conditionMessage(error), fixed = TRUE),
          grepl("INCOMPLETE", readLines(failed_status)[1L], fixed = TRUE))
error <- tryCatch(run_one("17-150-A-G", "het", file.path(tmp, "failed_plot"), fail_sashimi = TRUE), error = identity)
failed_status <- file.path(tmp, "failed_plot", "GENE_1_nt_2_puro", "STATUS.txt")
stopifnot(inherits(error, "error"), grepl("Simulated sashimipy failure", conditionMessage(error), fixed = TRUE),
          grepl("INCOMPLETE", readLines(failed_status)[1L], fixed = TRUE))
cat("Mocked report pipeline tests passed\n")
