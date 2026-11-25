# Fælles validering for pool og individuel referencebank.
stop_reference <- function(...) stop(..., call. = FALSE)

# Gamle betingelsesnavne normaliseres uden at ændre eksisterende bankfiler.
reference_condition <- function(value) {
    if (anyNA(value) || any(!value %in% c("nt", "puro", "minus", "plus"))) {
        stop_reference("Condition must be nt or puro (legacy aliases: minus or plus).")
    }
    value[value == "minus"] <- "nt"
    value[value == "plus"] <- "puro"
    value
}

reference_options <- function(args, allowed, required) {
    if (length(args) %% 2L) stop_reference("Every option requires a value.")
    result <- list()
    if (length(args)) for (i in seq.int(1L, length(args), by = 2L)) {
        key <- args[i]
        if (!(key %in% allowed)) stop_reference("Unknown option: ", key)
        if (!is.null(result[[key]])) stop_reference("Duplicate option: ", key)
        if (!nzchar(args[i + 1L])) stop_reference("Empty value: ", key)
        result[[key]] <- args[i + 1L]
    }
    missing <- required[!vapply(required, function(x) !is.null(result[[x]]), logical(1))]
    if (length(missing)) stop_reference("Missing options: ", paste(missing, collapse = ", "))
    result
}

reference_run <- function(command, args, stdout, stderr) {
    cat("Running:", command, paste(shQuote(args), collapse = " "), "\n")
    status <- system2(command, shQuote(args), stdout = stdout, stderr = stderr)
    if (!identical(status, 0L)) stop_reference(command, " failed (", status, "); see ", stderr)
}

reference_id <- function(value, label) {
    if (length(value) != 1L || is.na(value) || !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", value)) {
        stop_reference("Invalid ", label, ": ", value)
    }
    value
}

reference_bam <- function(path) {
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) stop_reference("BAM paths must be absolute: ", path)
    if (grepl("[\t\r\n]", path)) stop_reference("BAM paths cannot contain tabs or newlines.")
    if (!file.exists(path) || dir.exists(path)) stop_reference("BAM not found: ", path)
    if (!grepl("\\.bam$", path, ignore.case = TRUE)) stop_reference("Expected a BAM file: ", path)
    bam <- normalizePath(path, mustWork = TRUE)
    indexes <- c(paste0(bam, ".bai"), sub("\\.bam$", ".bai", bam, ignore.case = TRUE),
                 paste0(bam, ".csi"), sub("\\.bam$", ".csi", bam, ignore.case = TRUE))
    if (!any(file.exists(indexes))) stop_reference("Missing BAM index: ", bam)
    bam
}

# Kontroller BAM-sortering og sekvensdictionary før pooling eller banktilføjelse.
reference_header <- function(bam, prefix) {
    output <- paste0(prefix, ".header.sam")
    reference_run("samtools", c("quickcheck", bam), paste0(prefix, ".quickcheck.out"), paste0(prefix, ".quickcheck.err"))
    reference_run("samtools", c("view", "-H", bam), output, paste0(prefix, ".header.err"))
    lines <- readLines(output, warn = FALSE)
    hd <- lines[startsWith(lines, "@HD\t")]
    if (length(hd) != 1L || !("SO:coordinate" %in% strsplit(hd, "\t", fixed = TRUE)[[1L]])) {
        stop_reference("BAM header must declare coordinate sorting: ", bam)
    }
    sq <- lines[startsWith(lines, "@SQ\t")]
    if (!length(sq)) stop_reference("BAM has no sequence dictionary: ", bam)
    canonical <- vapply(strsplit(sq, "\t", fixed = TRUE), function(parts) {
        sn <- parts[startsWith(parts, "SN:")]
        ln <- parts[startsWith(parts, "LN:")]
        if (length(sn) != 1L || length(ln) != 1L || !grepl("^LN:[1-9][0-9]*$", ln)) stop_reference("Invalid @SQ in ", bam)
        paste(sn, ln, sep = "\t")
    }, character(1))
    if (anyDuplicated(sub("\t.*", "", canonical))) stop_reference("Duplicate contigs in BAM: ", bam)
    canonical
}

read_reference_bed <- function(path) {
    data <- read_regtools_junctions(path, "sample")
    data$sample <- NULL
    names(data)[names(data) == "intron_start_0"] <- "start"
    names(data)[names(data) == "intron_end_0"] <- "end"
    if (anyDuplicated(data[c("chrom", "start", "end", "strand")])) {
        data <- aggregate(reads ~ chrom + start + end + strand, data, sum)
    }
    data
}

# Manifestet kræver eksplicit godkendt QC og én prøve pr. person/betingelse.
read_bank_manifest <- function(path) {
    dat <- read.delim(path, colClasses = "character", check.names = FALSE, quote = "", comment.char = "", na.strings = character())
    required <- c("sample_id", "person_id", "condition", "bam", "qc")
    if (!identical(names(dat), required) || !nrow(dat)) stop_reference("Manifest must have columns sample_id, person_id, condition, bam, qc, in that order.")
    if (anyNA(dat) || any(!nzchar(as.matrix(dat)))) stop_reference("Manifest has empty values.")
    invisible(lapply(dat$sample_id, reference_id, "sample_id"))
    invisible(lapply(dat$person_id, reference_id, "person_id"))
    dat$condition <- reference_condition(dat$condition)
    if (any(dat$qc != "pass")) stop_reference("Manifest may contain only explicitly QC-approved samples: qc=pass.")
    if (anyDuplicated(dat$sample_id) || anyDuplicated(dat[c("person_id", "condition")])) {
        stop_reference("Duplicate sample or person/condition in manifest.")
    }
    dat
}
