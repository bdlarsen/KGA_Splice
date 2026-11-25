#!/usr/bin/env Rscript
# Individuel junctionbank: init, add og compare.
if (sys.nframe() == 0L) {
    script <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(script) != 1L) stop("Run with Rscript.")
    script_dir <- dirname(normalizePath(sub("^--file=", "", script), mustWork = TRUE))
} else {
    script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = TRUE))
}
source(file.path(script_dir, "clinical_core.R"))
source(file.path(script_dir, "reference_core.R"))

bank_regtools <- function() Sys.which("regtools")

bank_samples <- function(bank) {
    dirs <- list.dirs(file.path(bank, "samples"), full.names = TRUE, recursive = FALSE)
    if (!length(dirs)) return(data.frame(sample_id = character(), person_id = character(), condition = character(),
                                         bam = character(), size = numeric(), modified = numeric()))
    rows <- lapply(dirs, function(d) {
        path <- file.path(d, "metadata.rds")
        if (!file.exists(path)) stop_reference("Incomplete published sample: ", d)
        readRDS(path)
    })
    samples <- do.call(rbind, rows)
    samples$condition <- reference_condition(samples$condition)
    samples
}

# Opret banken og fastlås regtools-indstillinger.
init_bank <- function(opt) {
    out <- opt[["--bank"]]
    if (file.exists(out)) stop_reference("Bank already exists: ", out)
    ref <- reference_id(opt[["--reference-id"]], "reference-id")
    tool <- bank_regtools()
    if (!nzchar(tool)) stop_reference("regtools not found on PATH.")
    if (!dir.create(out)) stop_reference("Cannot create bank: ", out)
    if (!dir.create(file.path(out, "samples"))) stop_reference("Cannot create samples directory.")
    if (!dir.create(file.path(out, "staging"))) stop_reference("Cannot create staging directory.")
    settings <- list(schema = 1L, reference_id = ref, regtools_md5 = unname(tools::md5sum(tool)),
                     regtools_path = normalizePath(tool), regtools_args = c("-s", "XS", "-a", "8", "-m", "50", "-M", "500000"))
    saveRDS(settings, file.path(out, "settings.rds"))
    writeLines(c("ClinV1 individual junction bank; research/validation only.",
                 paste("Reference identifier:", ref),
                 paste("Regtools binary MD5:", settings$regtools_md5),
                 "Samples are append-only. A missing junction does not prove biological absence."),
               file.path(out, "README.txt"))
    cat("Created bank: ", out, "\n", sep = "")
}

# Tilføj godkendte prøver enkeltvis; eksisterende data overskrives ikke.
add_samples <- function(opt) {
    bank <- normalizePath(opt[["--bank"]], mustWork = TRUE)
    settings <- readRDS(file.path(bank, "settings.rds"))
    tool <- bank_regtools()
    if (!nzchar(tool) || unname(tools::md5sum(tool)) != settings$regtools_md5) stop_reference("regtools differs from the bank's frozen binary.")
    manifest <- read_bank_manifest(opt[["--manifest"]])
    manifest$bam <- vapply(manifest$bam, reference_bam, character(1))
    if (anyDuplicated(manifest$bam)) stop_reference("Duplicate BAM in manifest.")
    lock <- file.path(bank, ".add.lock")
    if (!dir.create(lock, showWarnings = FALSE)) stop_reference("Bank is locked: ", lock, ". Verify no writer is active before removing a stale lock.")
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    existing <- bank_samples(bank)
    for (i in seq_len(nrow(manifest))) {
        x <- manifest[i, ]
        info <- file.info(x$bam)
        found <- existing[existing$sample_id == x$sample_id, , drop = FALSE]
        if (nrow(found)) {
            if (found$person_id != x$person_id || found$condition != x$condition || found$bam != x$bam ||
                found$size != info$size || found$modified != as.numeric(info$mtime)) {
                stop_reference("Existing sample has changed: ", x$sample_id, ". Create a new bank after review.")
            }
            cat("Already present: ", x$sample_id, "\n", sep = "")
            next
        }
        if (any(existing$person_id == x$person_id & existing$condition == x$condition)) stop_reference("Person/condition already exists: ", x$person_id)
        if (x$bam %in% existing$bam) stop_reference("BAM already belongs to another sample: ", x$bam)
        stage <- tempfile(pattern = ".stage_", tmpdir = file.path(bank, "staging"))
        if (!dir.create(stage)) stop_reference("Cannot create staging directory.")
        stage_done <- FALSE
        tryCatch({
            prefix <- file.path(stage, "input")
            sq <- reference_header(x$bam, prefix)
            dict_path <- file.path(bank, "sequence_dictionary.rds")
            if (file.exists(dict_path) && !identical(readRDS(dict_path), sq)) stop_reference("BAM sequence dictionary differs from bank: ", x$bam)
            bed <- file.path(stage, "junctions.bed")
            reference_run("regtools", c("junctions", "extract", settings$regtools_args, x$bam),
                          bed, file.path(stage, "regtools.err"))
            after <- file.info(x$bam)
            if (is.na(after$size) || after$size != info$size || as.numeric(after$mtime) != as.numeric(info$mtime)) {
                stop_reference("BAM changed while processing: ", x$bam)
            }
            junctions <- read_reference_bed(bed)
            contigs <- sub("^SN:([^\t]+)\t.*$", "\\1", sq)
            if (any(!junctions$chrom %in% contigs)) stop_reference("Junction contig missing from BAM header.")
            inventory <- data.frame(chrom = character(), filename = character(), md5 = character())
            for (chrom in unique(junctions$chrom)) {
                index <- match(chrom, contigs)
                filename <- sprintf("chrom_%05d.rds", index)
                saveRDS(junctions[junctions$chrom == chrom, , drop = FALSE],
                        file.path(stage, filename), compress = TRUE)
                inventory <- rbind(inventory, data.frame(chrom = chrom, filename = filename,
                    md5 = unname(tools::md5sum(file.path(stage, filename)))))
            }
            saveRDS(inventory, file.path(stage, "inventory.rds"))
            metadata <- data.frame(sample_id = x$sample_id, person_id = x$person_id, condition = x$condition,
                                   bam = x$bam, size = info$size, modified = as.numeric(info$mtime), stringsAsFactors = FALSE)
            saveRDS(metadata, file.path(stage, "metadata.rds"))
            write.table(metadata, file.path(stage, "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
            if (!file.exists(dict_path)) saveRDS(sq, dict_path)
            target <- file.path(bank, "samples", x$sample_id)
            if (file.exists(target) || !file.rename(stage, target)) stop_reference("Cannot publish sample: ", target)
            stage_done <- TRUE
            existing <- rbind(existing, metadata)
            cat("Added ", x$sample_id, " (", nrow(junctions), " junctions)\n", sep = "")
        }, finally = {
            if (!stage_done && dir.exists(stage)) {
                cat("Failed sample staging retained for inspection: ", stage, "\n", sep = "")
            }
        })
    }
}

# Sammenlign én patientregion med eksplicit valgte, behandlingsmatchede prøver.
compare_samples <- function(opt) {
    bank <- normalizePath(opt[["--bank"]], mustWork = TRUE)
    condition <- reference_condition(opt[["--condition"]])
    person <- reference_id(opt[["--patient"]], "patient")
    samples <- bank_samples(bank)
    selected <- trimws(readLines(opt[["--controls"]], warn = FALSE))
    selected <- selected[nzchar(selected) & !startsWith(selected, "#")]
    if (!length(selected) || anyDuplicated(selected)) stop_reference("Control list must contain distinct sample IDs.")
    if (any(!selected %in% samples$sample_id)) stop_reference("Control list contains unknown sample ID.")
    samples <- samples[match(selected, samples$sample_id), , drop = FALSE]
    if (any(samples$condition != condition)) stop_reference("Control condition does not match patient condition.")
    if (any(samples$person_id == person)) stop_reference("Patient cannot be their own control.")
    if (anyDuplicated(samples$person_id)) stop_reference("Controls must be independent persons.")
    region <- regexec("^([^:]+):([0-9]+)-([0-9]+)$", opt[["--region"]])
    fields <- regmatches(opt[["--region"]], region)[[1L]]
    if (length(fields) != 4L) stop_reference("Region must be chromosome:start-end, 1-based inclusive.")
    chrom <- fields[2L]
    bounds <- suppressWarnings(as.integer(fields[3:4]))
    if (anyNA(bounds) || bounds[1L] < 1L || bounds[2L] < bounds[1L]) stop_reference("Invalid region bounds.")
    target <- read_reference_bed(opt[["--junctions"]])
    dict <- readRDS(file.path(bank, "sequence_dictionary.rds"))
    contigs <- sub("^SN:([^\t]+)\t.*$", "\\1", dict)
    if (!chrom %in% contigs || any(!target$chrom %in% contigs)) stop_reference("Contig is absent from bank reference.")
    within <- function(x) x[x$chrom == chrom & x$start < bounds[2L] & x$end > bounds[1L] - 1L, , drop = FALSE]
    target <- within(target)
    control_data <- lapply(samples$sample_id, function(id) {
        sample_dir <- file.path(bank, "samples", id)
        inventory_file <- file.path(sample_dir, "inventory.rds")
        if (!file.exists(inventory_file)) stop_reference("Missing sample inventory: ", id)
        inventory <- readRDS(inventory_file)
        entry <- inventory[inventory$chrom == chrom, , drop = FALSE]
        if (nrow(entry) > 1L) stop_reference("Duplicate chromosome inventory: ", id)
        if (nrow(entry) == 0L) return(target[FALSE, , drop = FALSE])
        rds <- file.path(sample_dir, entry$filename)
        if (!file.exists(rds) || is.na(tools::md5sum(rds)) || unname(tools::md5sum(rds)) != entry$md5) {
            stop_reference("Missing or damaged chromosome file for sample: ", id)
        }
        within(readRDS(rds))
    })
    keys <- c("chrom", "start", "end", "strand")
    all_keys <- unique(do.call(rbind, c(list(target[keys]), lapply(control_data, function(x) x[keys]))))
    patient <- target
    names(patient)[names(patient) == "reads"] <- "patient_reads"
    patient <- merge(all_keys, patient, by = keys, all.x = TRUE, sort = FALSE)
    out <- opt[["--output-dir"]]
    if (file.exists(out)) stop_reference("Output directory exists: ", out)
    if (!dir.create(out)) stop_reference("Cannot create output directory: ", out)
    writeLines("INCOMPLETE", file.path(out, "STATUS.txt"))
    write.table(samples, file.path(out, "controls_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(target, file.path(out, "patient_junctions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    audit <- lapply(samples$sample_id, function(id) {
        inv <- readRDS(file.path(bank, "samples", id, "inventory.rds"))
        row <- inv[inv$chrom == chrom, , drop = FALSE]
        if (!nrow(row)) row <- data.frame(chrom = chrom, filename = "", md5 = "")
        data.frame(sample_id = id, row, stringsAsFactors = FALSE)
    })
    write.table(do.call(rbind, audit), file.path(out, "control_inventory.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    output <- file(file.path(out, "comparison.tsv"), open = "wt")
    on.exit(close(output), add = TRUE)
    writeLines("chrom\tintron_start_0\tintron_end_0\tstrand\tpatient_reads\tpatient_status\tcontrol_sample\tcontrol_person\tcontrol_reads\tcontrol_status", output)
    for (i in seq_len(nrow(samples))) {
        control <- samples[i, ]
        ref <- control_data[[i]]
        names(ref)[names(ref) == "reads"] <- "control_reads"
        matched <- merge(patient, ref, by = keys, all.x = TRUE, sort = FALSE)
        matched$control_sample <- rep(control$sample_id, nrow(matched))
        matched$control_person <- rep(control$person_id, nrow(matched))
        matched$patient_status <- ifelse(is.na(matched$patient_reads), "not_observed_coverage_unassessed", "observed")
        matched$control_status <- ifelse(is.na(matched$control_reads), "not_observed_coverage_unassessed", "observed")
        result <- matched[c(keys, "patient_reads", "patient_status", "control_sample", "control_person", "control_reads", "control_status")]
        write.table(result, output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "NA")
    }
    close(output)
    on.exit(NULL, add = FALSE)
    writeLines(c("COMPLETE - descriptive research/validation output, not a clinical classification.",
                 paste("Bank:", bank), paste("Reference identifier:", readRDS(file.path(bank, "settings.rds"))$reference_id),
                 paste("regtools executable MD5:", readRDS(file.path(bank, "settings.rds"))$regtools_md5),
                 paste("Patient BED MD5:", unname(tools::md5sum(opt[["--junctions"]]))),
                 paste("Patient:", person), paste("Condition:", condition),
                 paste("Region:", opt[["--region"]]),
                 "Junctions overlapping the region from patient or selected controls are included.",
                 "Not observed means unknown coverage, not an assessed biological zero.",
                 "No normalized splice usage, retention, statistical test, or RNA effect classification is calculated.",
                 "Controls were explicitly selected by the operator; gene relevance and QC are not verified by software."),
               file.path(out, "STATUS.txt"))
    cat("Comparison ready: ", out, "\n", sep = "")
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if (!length(args)) stop_reference("Use init, add, or compare; see README.md.")
    action <- args[1L]
    if (action == "init") {
        init_bank(reference_options(args[-1L], c("--bank", "--reference-id"), c("--bank", "--reference-id")))
    } else if (action == "add") {
        add_samples(reference_options(args[-1L], c("--bank", "--manifest"), c("--bank", "--manifest")))
    } else if (action == "compare") {
        compare_samples(reference_options(args[-1L],
            c("--bank", "--condition", "--patient", "--junctions", "--controls", "--region", "--output-dir"),
            c("--bank", "--condition", "--patient", "--junctions", "--controls", "--region", "--output-dir")))
    } else stop_reference("Unknown action: ", action)
}
