# Beregn transkriptspecifikke scores lokalt og send dem til Broad-adapteren.
prepare_sai_transcript <- function(reference, tx, variant_pos, outdir, variant_end = variant_pos) {
    records <- strsplit(reference$lines, "\t", fixed = TRUE)
    selected <- Filter(function(x) {
        m <- regmatches(x[9L], regexec('(?:^|;[[:space:]]*)transcript_id "([^"]+)"', x[9L], perl = TRUE))[[1L]]
        length(m) > 1L && identical(m[2L], tx)
    }, records)
    exons <- Filter(function(x) x[3L] == "exon", selected)
    if (!length(exons)) stop_input("SAI: selected transcript has no GTF exons.")
    ex <- unique(data.frame(start = as.integer(vapply(exons, `[`, "", 4L)),
                            end = as.integer(vapply(exons, `[`, "", 5L))))
    ex <- ex[order(ex$start, ex$end), , drop = FALSE]
    if (anyNA(ex) || any(ex$start < 1L | ex$start > ex$end) ||
        (nrow(ex) > 1L && any(ex$start[-1L] <= ex$end[-nrow(ex)]))) stop_input("SAI: invalid or overlapping transcript exons.")
    if (variant_pos < ex$start[1L] || variant_end > ex$end[nrow(ex)]) return(NULL)
    cds <- Filter(function(x) x[3L] %in% c("CDS", "stop_codon"), selected)
    if (!any(vapply(cds, function(x) x[3L] == "CDS", logical(1)))) {
        stop_input("SAI: selected transcript has no CDS annotation; coding status cannot be established.")
    }
    cds_start <- cds_end <- NA_integer_
    if (length(cds)) {
        cds_starts <- as.integer(vapply(cds, `[`, "", 4L))
        cds_ends <- as.integer(vapply(cds, `[`, "", 5L))
        if (anyNA(cds_starts) || anyNA(cds_ends) || any(cds_starts > cds_ends) ||
            any(!vapply(seq_along(cds_starts), function(i) any(ex$start <= cds_starts[i] & ex$end >= cds_ends[i]), logical(1)))) {
            stop_input("SAI: invalid CDS coordinates in GTF.")
        }
        cds_start <- min(cds_starts)
        cds_end <- max(cds_ends)
    }
    annot <- file.path(outdir, "spliceai_transcript.tsv")
    tab <- data.frame("#NAME" = paste0("RefSeqTx-", tx), CHROM = reference$chrom,
        STRAND = reference$strand, TX_START = ex$start[1L] - 1L, TX_END = ex$end[nrow(ex)],
        EXON_START = paste0(paste(ex$start - 1L, collapse = ","), ","),
        EXON_END = paste0(paste(ex$end, collapse = ","), ","), check.names = FALSE)
    write.table(tab, annot, sep = "\t", quote = FALSE, row.names = FALSE)
    transcript <- file.path(outdir, "transcript.tsv")
    write.table(data.frame(transcript = tx, strand = reference$strand,
                           exon_starts = paste(ex$start, collapse = ","),
                           exon_ends = paste(ex$end, collapse = ","),
                           cds_start = cds_start, cds_end = cds_end),
                transcript, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    list(annotation = annot, transcript = transcript)
}

run_sai <- function(module, env, vcf, fasta, tx, reference, variant, project, logs, adapter) {
    outdir <- file.path(project, "SAI10k")
    if (!dir.create(outdir)) stop_input("Cannot create SAI10k output directory.")
    writeLines(c(paste("Broad module:", normalizePath(module, mustWork = TRUE)),
                 paste("Broad module MD5:", unname(tools::md5sum(module))),
                 paste("GTF transcript:", tx), paste("Conda environment:", env),
                 "Genome: GRCh38; SpliceAI distance: 4999; mask: 0",
                 "Standard SpliceAI delta scores only; no raw REF/ALT or per-position scores.",
                 "No FASTA passed to Broad predictor: premature-stop and amino acid sequence not assessed."),
               file.path(outdir, "provenance.txt"))
    if (nchar(variant$ref) > 1L && nchar(variant$alt) > 1L)
        return(list(status = "complex_variant_not_supported_by_spliceai", data = data.frame()))
    inputs <- prepare_sai_transcript(reference, tx, variant$pos, outdir,
                                     variant$pos + nchar(variant$ref) - 1L)
    if (is.null(inputs)) return(list(status = "variant_outside_selected_transcript", data = data.frame()))
    scored <- file.path(outdir, "spliceai_transcript.vcf")
    run_checked("conda", c("run", "--no-capture-output", "-n", "spliceai_env", "spliceai",
                           "-I", vcf, "-O", scored, "-R", fasta, "-A", inputs$annotation, "-D", "4999", "-M", "0"),
                file.path(logs, "sai_spliceai.out"), file.path(logs, "sai_spliceai.err"))
    fields <- validate_annotated_vcf(scored, variant, reference$chrom)
    if (!has_info_annotation(fields[8L], "SpliceAI")) return(list(status = "spliceai_annotation_unavailable", data = data.frame()))
    scores <- parse_spliceai(fields[8L], paste0("RefSeqTx-", tx))
    if ("status" %in% names(scores) || nrow(scores) != 1L || scores$ALLELE != variant$alt) {
        return(list(status = "transcript_scores_unavailable", data = data.frame()))
    }
    cols <- c("DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")
    values <- suppressWarnings(as.numeric(unlist(scores[1L, cols], use.names = FALSE)))
    if (any(!is.finite(values)) || any(values[1:4] < 0 | values[1:4] > 1) ||
        any(values[5:8] != floor(values[5:8])) || any(abs(values[5:8]) > 4999)) {
        stop_input("SAI: invalid transcript SpliceAI scores or positions.")
    }
    score_file <- file.path(outdir, "scores.tsv")
    write.table(scores, score_file, sep = "\t", quote = FALSE, row.names = FALSE)
    output <- file.path(outdir, "predictions.tsv")
    json <- file.path(outdir, "predictions.json")
    run_checked("conda", c("run", "--no-capture-output", "-n", env, "python", adapter,
                           "--module", module, "--transcript", inputs$transcript, "--scores", score_file,
                           "--chrom", reference$chrom, "--pos", as.character(variant$pos),
                           "--ref", variant$ref, "--alt", variant$alt, "--output", output, "--json", json),
                file.path(logs, "sai_predictor.out"), file.path(logs, "sai_predictor.err"))
    if (!file.exists(output) || !file.exists(json)) stop_input("SAI output missing; inspect predictor logs.")
    result <- read.delim(output, check.names = FALSE, colClasses = "character", quote = "", comment.char = "")
    required <- c("transcript", "status", "effect", "description", "consequence", "size_bp", "frameshift", "affects_coding")
    if (!identical(names(result), required) || !nrow(result) || anyNA(result$transcript) ||
        any(result$transcript != tx) || anyNA(result$status) || anyNA(result$effect) ||
        any(!nzchar(result$effect)) || any(!result$status %in% c("predicted", "no_effect_predicted")) ||
        (any(result$status == "predicted") && any(result$effect == "none")) ||
        (any(result$status == "no_effect_predicted") && (nrow(result) != 1L || result$effect[1L] != "none"))) {
        stop_input("Invalid SAI predictor output: ", output)
    }
    list(status = if (result$status[1L] == "predicted") "prediction_available_exploratory" else "no_effect_predicted_exploratory",
         data = result, scores = scores)
}
