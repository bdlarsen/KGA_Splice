#!/usr/bin/env Rscript
args <- grep("^--file=", commandArgs(), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)))
source(file.path(base, "scripts", "clinical_core.R"))
source(file.path(base, "scripts", "report_core.R"))
source(file.path(base, "scripts", "sai_core.R"))
expect_error <- function(expr, expected) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    if (is.null(error) || !grepl(expected, conditionMessage(error), fixed = TRUE)) stop("Expected error: ", expected)
}
tmp <- tempfile("sai_test_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
for (d in c("snp", "indel", "missing", "bad", "outside", "logs", "check")) dir.create(file.path(tmp, d))
reference <- list(chrom = "chr17", strand = "+", lines = c(
    'chr17\ttest\texon\t101\t150\t.\t+\t.\tgene_id "GENE"; transcript_id "NM_1.1";',
    'chr17\ttest\texon\t201\t250\t.\t+\t.\tgene_id "GENE"; transcript_id "NM_1.1";',
    'chr17\ttest\tCDS\t111\t150\t.\t+\t0\tgene_id "GENE"; transcript_id "NM_1.1";',
    'chr17\ttest\tCDS\t201\t240\t.\t+\t0\tgene_id "GENE"; transcript_id "NM_1.1";'))
module <- file.path(tmp, "module.py")
writeLines("# mock module", module)
inputs <- prepare_sai_transcript(reference, "NM_1.1", 150L, file.path(tmp, "check"))
annotation <- read.delim(inputs$annotation, check.names = FALSE)
transcript <- read.delim(inputs$transcript)
stopifnot(annotation[["#NAME"]] == "RefSeqTx-NM_1.1", annotation$TX_START == 100,
          transcript$exon_starts == "101,201", transcript$cds_start == 111)
reverse_reference <- reference
reverse_reference$strand <- "-"
reverse_reference$lines <- sub("\t\\+\t", "\t-\t", reference$lines)
reverse_inputs <- prepare_sai_transcript(reverse_reference, "NM_1.1", 150L, file.path(tmp, "check"))
reverse_transcript <- read.delim(reverse_inputs$transcript)
stopifnot(reverse_transcript$strand == "-", reverse_transcript$exon_starts == "101,201",
          reverse_transcript$cds_start == 111, reverse_transcript$cds_end == 240)
bad_reference <- reference
bad_reference$lines <- c(reference$lines, 'chr17\ttest\texon\t140\t160\t.\t+\t.\ttranscript_id "NM_1.1";')
expect_error(prepare_sai_transcript(bad_reference, "NM_1.1", 150L, file.path(tmp, "check")), "overlapping")
bad_reference$lines <- reference$lines[1:2]
expect_error(prepare_sai_transcript(bad_reference, "NM_1.1", 150L, file.path(tmp, "check")), "no CDS")
stopifnot(is.null(prepare_sai_transcript(reference, "NM_1.1", 250L, file.path(tmp, "check"), 251L)))
calls <- character()
missing <- FALSE
bad <- FALSE
run_checked <- function(command, args, stdout, stderr) {
    calls <<- c(calls, args[5L])
    writeLines(character(), stdout)
    writeLines(character(), stderr)
    if (args[5L] == "spliceai") {
        stopifnot(args[which(args == "-D") + 1L] == "4999", args[which(args == "-M") + 1L] == "0")
        input <- readLines(args[which(args == "-I") + 1L])
        row <- strsplit(tail(input, 1L), "\t", fixed = TRUE)[[1L]]
        row[8L] <- if (missing) "." else paste0("SpliceAI=", row[5L], "|RefSeqTx-NM_1.1|0.3|0.0|0.0|0.0|1|0|0|0")
        writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", paste(row, collapse = "\t")),
                   args[which(args == "-O") + 1L])
    } else {
        stopifnot(args[5L] == "python")
        data <- data.frame(transcript = if (bad) "NM_wrong" else "NM_1.1", status = "predicted",
                           effect = "partial_exon_deletion", description = "donor shift", consequence = "exon deletion", size_bp = "1",
                           frameshift = "True", affects_coding = "True")
        write.table(data, args[which(args == "--output") + 1L], sep = "\t", quote = FALSE, row.names = FALSE)
        writeLines("{}", args[which(args == "--json") + 1L])
    }
}
invoke <- function(project, var) {
    vcf <- file.path(tmp, paste0(project, ".vcf"))
    writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                 paste("chr17", var$pos, ".", var$ref, var$alt, ".", ".", ".", sep = "\t")), vcf)
    run_sai(module, "KGA_Splice", vcf, "genome.fa", "NM_1.1", reference, var,
            file.path(tmp, project), file.path(tmp, "logs"), "adapter.py")
}
snp <- list(pos = 150L, ref = "A", alt = "G")
result <- invoke("snp", snp)
stopifnot(result$status == "prediction_available_exploratory", result$data$effect == "partial_exon_deletion",
          result$scores$DS_AG == "0.3",
          identical(calls, c("spliceai", "python")))
calls <- character()
result <- invoke("indel", list(pos = 150L, ref = "AC", alt = "A"))
stopifnot(result$status == "prediction_available_exploratory", identical(calls, c("spliceai", "python")))
calls <- character()
missing <- TRUE
result <- invoke("missing", snp)
stopifnot(result$status == "spliceai_annotation_unavailable", identical(calls, "spliceai"))
missing <- FALSE
bad <- TRUE
expect_error(invoke("bad", snp), "Invalid SAI predictor output")
calls <- character()
result <- invoke("outside", list(pos = 251L, ref = "A", alt = "G"))
stopifnot(result$status == "variant_outside_selected_transcript", !length(calls))
cat("SAI integration tests passed\n")
