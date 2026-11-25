# Fælles inputkontrol og hjælpefunktioner til ClinV1.
stop_input <- function(...) stop(..., call. = FALSE)

write_vswarehouse_igv_session <- function(input, output) {
    session <- xml2::read_xml(input)
    if (xml2::xml_name(xml2::xml_root(session)) != "Session")
        stop_input("IGV session must have a Session root: ", input)
    resources <- xml2::xml_find_all(session, ".//Resource[@path]")
    if (!length(resources)) stop_input("IGV session has no resource paths: ", input)
    paths <- xml2::xml_attr(resources, "path")
    mapped <- paths == "/drift" | startsWith(paths, "/drift/")
    if (any(!mapped)) warning(
        "VSwarehouse IGV session: resource paths outside /drift were not translated; ",
        "verify access in VSwarehouse: ", paste(unique(paths[!mapped]), collapse = ", "),
        call. = FALSE)
    translate <- function(value) {
        mapped <- value == "/drift" | startsWith(value, "/drift/")
        value[mapped] <- paste0("/home/ghuser/Workspace", value[mapped])
        value
    }
    xml2::xml_set_attr(resources, "path", translate(paths))
    tracks <- xml2::xml_find_all(session, ".//Track[@id]")
    xml2::xml_set_attr(tracks, "id", translate(xml2::xml_attr(tracks, "id")))
    xml2::write_xml(session, output)
    invisible(output)
}

clinv1_help <- function() {
    paste(c(
        "ClinV1: one variant and one sample or nt/puro sample pair",
        "Usage: Rscript ClinV1.R -variant VARIANT -gene GENE -tx TRANSCRIPT",
        "       (-p OUTPUT_PARENT | -project-dir RESERVED_DIR) -bam BAM_DIR -s ID or ID_NT,ID_PURO [options]",
        "",
        "Variant: chromosome-position-REF-ALT (e.g. 22-28725241-AC-A)",
        "         or a VCF/VCF.GZ file containing exactly one biallelic variant.",
        "Samples: one ID (treatment unknown), or nt,puro IDs (without/with puromycin).",
        "Optional: -gt het|hmz|hemi, -refs CONTROL_BAM_DIR,",
        "          -sai-module PATH, -sai-env ENV, -g FASTA, -gtf GTF, -vep CACHE.",
        "          -igv-snapshot PATH_TO_IGV_SH (requires xvfb-run and network isolation).",
        "VCF genotype is not used; provide independently verified -gt when appropriate.",
        "A reserved project directory must exist and contain only an optional tmp directory.",
        "For full documentation see Pipeline/ClinV1/README.md."
    ), collapse = "\n")
}

check_required_commands <- function(commands, find_command = Sys.which) {
    found <- find_command(commands)
    missing <- commands[is.na(found) | !nzchar(found)]
    if (length(missing)) {
        stop_input("Required commands not found on PATH: ", paste(missing, collapse = ", "),
                   ". Install or expose them in the current environment before running ClinV1.")
    }
    invisible(found)
}

parse_options <- function(args) {
    defaults <- list(
        "-variant" = NULL, "-gene" = NULL, "-tx" = NULL, "-p" = NULL, "-project-dir" = NULL,
        "-bam" = NULL, "-s" = NULL, "-g" = "/drift/RNA/TotalRNA_Splice/Refs/references/genome.fa",
        "-gtf" = "/drift/RNA/TotalRNA_Splice/GRCh38_Reference.gtf",
        "-vep" = "/drift/RNA/TotalRNA_Splice/VEP_Refs", "-refs" = NULL, "-gt" = NULL,
        "-sai-module" = NULL, "-sai-env" = NULL, "-igv-snapshot" = NULL
    )
    if (length(args) %% 2L) stop_input("Each option must have a value.")
    if (!length(args)) stop_input("Missing arguments. See README.md.")
    for (i in seq.int(1L, length(args), by = 2L)) {
        key <- args[i]
        if (!(key %in% names(defaults))) stop_input("Unknown option: ", key)
        if (i > 1L && key %in% args[seq.int(1L, i - 2L, by = 2L)]) stop_input("Duplicate option: ", key)
        defaults[key] <- list(args[i + 1L])
    }
    for (name in c("-variant", "-gene", "-tx", "-bam", "-s")) {
        if (is.null(defaults[[name]]) || !nzchar(defaults[[name]])) stop_input("Missing option: ", name)
    }
    if (sum(vapply(defaults[c("-p", "-project-dir")], function(x) !is.null(x) && nzchar(x), logical(1))) != 1L)
        stop_input("Specify exactly one of -p or -project-dir.")
    defaults
}

prepare_project <- function(parent, reserved, project_name) {
    if (is.null(reserved)) {
        if (!dir.exists(parent)) stop_input("Output parent directory not found.")
        project <- file.path(normalizePath(parent, mustWork = TRUE), project_name)
        if (file.exists(project)) stop_input("Project already exists; choose a fresh output directory: ", project)
        if (!dir.create(project)) stop_input("Could not create project: ", project)
        return(project)
    }
    if (!dir.exists(reserved) || nzchar(Sys.readlink(reserved)))
        stop_input("Reserved project directory must exist and must not be a symlink: ", reserved)
    project <- normalizePath(reserved, mustWork = TRUE)
    entries <- list.files(project, all.files = TRUE, no.. = TRUE)
    if (any(entries != "tmp")) stop_input("Reserved project directory must be empty apart from tmp: ", project)
    if ("tmp" %in% entries && (!dir.exists(file.path(project, "tmp")) || nzchar(Sys.readlink(file.path(project, "tmp")))))
        stop_input("Project tmp must be a directory, not a symlink: ", project)
    project
}

parse_variant <- function(value) {
    fields <- strsplit(value, "-", fixed = TRUE)[[1]]
    if (length(fields) != 4L || !grepl("^(chr)?([0-9]+|X|Y|M|MT)$", fields[1], ignore.case = TRUE) ||
        !grepl("^[0-9]+$", fields[2]) || !grepl("^[ACGT]+$", fields[3], ignore.case = TRUE) ||
        !grepl("^[ACGT]+$", fields[4], ignore.case = TRUE)) {
        stop_input("Variant must be chromosome-position-REF-ALT with DNA bases, e.g. 22-123-AC-A.")
    }
    position <- suppressWarnings(as.integer(fields[2]))
    if (is.na(position) || position < 1L) stop_input("Invalid variant position.")
    if (toupper(fields[3]) == toupper(fields[4])) stop_input("REF and ALT must differ.")
    list(chrom = fields[1], pos = position, ref = toupper(fields[3]), alt = toupper(fields[4]))
}

read_single_vcf_variant <- function(path) {
    if (!file.exists(path) || dir.exists(path)) stop_input("VCF file not found: ", path)
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con))
    lines <- readLines(con, warn = FALSE)
    headers <- lines[startsWith(lines, "#CHROM\t")]
    rows <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (length(headers) != 1L || length(rows) != 1L)
        stop_input("VCF must contain a #CHROM header and exactly one variant record: ", path)
    header <- strsplit(headers, "\t", fixed = TRUE)[[1L]]
    row <- strsplit(rows, "\t", fixed = TRUE)[[1L]]
    if (length(header) < 8L || !identical(header[1:8], c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")) ||
        length(row) != length(header)) stop_input("Invalid VCF header or record: ", path)
    if (grepl(",", row[5L], fixed = TRUE)) stop_input("VCF must contain one biallelic variant; multiallelic ALT is not supported.")
    parse_variant(paste(row[c(1L, 2L, 4L, 5L)], collapse = "-"))
}

parse_variant_input <- function(value) {
    if (grepl("\\.vcf(\\.gz)?$", value, ignore.case = TRUE)) return(read_single_vcf_variant(value))
    parse_variant(value)
}

parse_sample_ids <- function(value) {
    if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*(,[A-Za-z0-9][A-Za-z0-9._-]*)?$", value))
        stop_input("-s requires one ID or two IDs using letters, digits, dot, underscore or hyphen: nt,puro.")
    ids <- strsplit(value, ",", fixed = TRUE)[[1L]]
    if (length(ids) == 2L && ids[1L] == ids[2L]) stop_input("The two sample IDs must be different.")
    ids
}

find_sample_bam <- function(directory, id) {
    files <- list.files(directory, pattern = "\\.bam$", full.names = TRUE)
    names <- basename(files)
    ids <- unique(c(id, paste0("G150-", id)))
    matched <- vapply(names, function(name) any(vapply(ids, function(candidate)
        name == paste0(candidate, ".bam") || startsWith(name, paste0(candidate, "_")) ||
            startsWith(name, paste0(candidate, "-")), logical(1))), logical(1))
    files <- files[matched]
    if (length(files) != 1L) stop_input("Expected one BAM for ", id, "; found ", length(files), ".")
    normalizePath(files, mustWork = TRUE)
}

check_bam_index <- function(bam) {
    candidates <- c(paste0(bam, ".bai"), sub("\\.bam$", ".bai", bam),
                    paste0(bam, ".csi"), sub("\\.bam$", ".csi", bam))
    if (!any(file.exists(candidates))) stop_input("Missing BAM index: ", bam)
}

choose_contig <- function(requested, contigs) {
    candidates <- unique(c(requested, sub("^chr", "", requested, ignore.case = TRUE),
                           paste0("chr", sub("^chr", "", requested, ignore.case = TRUE))))
    matches <- candidates[candidates %in% contigs]
    if (length(matches) != 1L) stop_input("Chromosome ", requested, " does not map uniquely to reference contigs.")
    matches
}

# Udvælg gen og transkript fra reference-GTF.
read_gene_gtf <- function(path, gene, transcript) {
    con <- file(path, open = "r")
    on.exit(close(con))
    lines <- character()
    repeat {
        chunk <- readLines(con, n = 100000L, warn = FALSE)
        if (!length(chunk)) break
        # Match complete attribute values rather than partial gene names.
        chunk <- chunk[grepl(gene, chunk, fixed = TRUE)]
        keep <- vapply(strsplit(chunk, "\t", fixed = TRUE), function(parts) {
            if (length(parts) != 9L) return(FALSE)
            attrs <- parts[9L]
            matches <- regmatches(attrs, gregexpr('(?:^|;[[:space:]]*)gene_(?:name|id) "[^"]+"', attrs, perl = TRUE))[[1L]]
            values <- sub('.*"([^"]+)"$', "\\1", matches)
            gene %in% values
        }, logical(1))
        lines <- c(lines, chunk[keep])
    }
    if (!length(lines)) stop_input("Gene not found in GTF: ", gene)
    fields <- strsplit(lines, "\t", fixed = TRUE)
    transcript_match <- vapply(fields, function(x) {
        value <- regmatches(x[9L], regexec('(?:^|;[[:space:]]*)transcript_id "([^"]+)"', x[9L], perl = TRUE))[[1L]]
        length(value) > 1L && identical(value[2L], transcript)
    }, logical(1))
    if (!any(transcript_match)) stop_input("Transcript ", transcript, " not found for gene ", gene, ".")
    chromosomes <- unique(vapply(fields, `[`, "", 1L))
    strands <- unique(vapply(fields, `[`, "", 7L))
    if (length(chromosomes) != 1L || length(strands) != 1L || !strands %in% c("+", "-")) {
        stop_input("Gene has ambiguous chromosome or strand in GTF.")
    }
    starts <- as.integer(vapply(fields, `[`, "", 4L))
    ends <- as.integer(vapply(fields, `[`, "", 5L))
    if (anyNA(starts) || anyNA(ends)) stop_input("Invalid gene coordinates in GTF.")
    list(lines = lines, chrom = chromosomes, strand = strands, start = min(starts), end = max(ends))
}

# Læs BED12-junctions og behold usikre strandværdier.
read_regtools_junctions <- function(path, sample) {
    lines <- readLines(path, warn = FALSE)
    result <- data.frame(chrom = character(), intron_start_0 = integer(), intron_end_0 = integer(),
                         strand = character(), sample = character(), reads = integer())
    if (!length(lines)) return(result)
    records <- strsplit(lines, "\t", fixed = TRUE)
    rows <- lapply(seq_along(records), function(i) {
        x <- records[[i]]
        if (length(x) < 12L || x[10L] != "2") stop_input("Invalid regtools BED12 output: ", path, " line ", i)
        blocks <- strsplit(sub(",$", "", x[11L]), ",", fixed = TRUE)[[1L]]
        values <- suppressWarnings(as.integer(c(x[2L], x[3L], x[5L], blocks)))
        if (length(values) != 5L || anyNA(values) || !x[6L] %in% c("+", "-", ".", "*", "?")) {
            stop_input("Invalid junction record: ", path, " line ", i, "; strand=", x[6L])
        }
        start <- values[1L] + values[4L]
        end <- values[2L] - values[5L]
        if (start >= end || values[3L] < 0L) stop_input("Invalid junction coordinates or read count: ", path, " line ", i)
        data.frame(chrom = x[1L], intron_start_0 = start, intron_end_0 = end,
                   strand = x[6L], sample = sample, reads = values[3L])
    })
    do.call(rbind, rows)
}

validate_annotated_vcf <- function(path, variant, chrom) {
    if (!file.exists(path) || file.info(path)$size == 0L) stop_input("Annotation output missing or empty: ", path)
    lines <- readLines(path, warn = FALSE)
    rows <- lines[!startsWith(lines, "#") & nzchar(lines)]
    if (length(rows) != 1L) stop_input("Annotation output must contain exactly one variant: ", path)
    fields <- strsplit(rows, "\t", fixed = TRUE)[[1L]]
    if (length(fields) < 8L || !identical(fields[1L], chrom) ||
        !identical(fields[2L], as.character(variant$pos)) ||
        !identical(toupper(fields[4L]), variant$ref) || !identical(toupper(fields[5L]), variant$alt)) {
        stop_input("Annotation output variant differs from input: ", path)
    }
    fields
}

has_info_annotation <- function(info, key) {
    parts <- strsplit(info, ";", fixed = TRUE)[[1L]]
    matches <- parts[startsWith(parts, paste0(key, "="))]
    if (length(matches) != 1L) return(FALSE)
    value <- substring(matches, nchar(key) + 2L)
    nzchar(value) && value != "."
}

# Sammenlign rå junctiontællinger og annotér eksakte introngrænser.
compare_junctions <- function(tables, samples) {
    all <- do.call(rbind, tables)
    if (is.null(all) || !nrow(all)) {
        result <- data.frame(chrom = character(), intron_start_0 = integer(), intron_end_0 = integer(), strand = character())
        for (sample in samples) result[[sample]] <- integer()
        return(result)
    }
    keys <- c("chrom", "intron_start_0", "intron_end_0", "strand")
    grouped <- aggregate(reads ~ chrom + intron_start_0 + intron_end_0 + strand + sample, all, sum)
    result <- unique(grouped[keys])
    for (sample in samples) {
        subset <- grouped[grouped$sample == sample, c(keys, "reads"), drop = FALSE]
        names(subset)[5L] <- sample
        result <- merge(result, subset, by = keys, all.x = TRUE, sort = FALSE)
        result[[sample]][is.na(result[[sample]])] <- 0L
    }
    result[order(result$chrom, result$intron_start_0, result$intron_end_0, result$strand), , drop = FALSE]
}

annotate_junctions <- function(junctions, gene_gtf, selected_tx) {
    records <- strsplit(gene_gtf$lines, "\t", fixed = TRUE)
    exon_records <- records[vapply(records, function(x) x[3L] == "exon", logical(1))]
    exons <- lapply(exon_records, function(x) {
        m <- regmatches(x[9L], regexec('(?:^|;[[:space:]]*)transcript_id "([^"]+)"', x[9L], perl = TRUE))[[1L]]
        if (length(m) < 2L) stop_input("GTF exon is missing transcript_id.")
        data.frame(tx = m[2L], start = as.integer(x[4L]), end = as.integer(x[5L]))
    })
    if (!length(exons)) stop_input("No exons in GTF for selected gene.")
    exons <- do.call(rbind, exons)
    if (anyNA(exons$start) || anyNA(exons$end) || any(exons$end < exons$start)) stop_input("Invalid exon coordinates.")
    make_introns <- function(e) {
        transcript_name <- e$tx[1L]
        e <- unique(e[order(e$start, e$end), c("start", "end"), drop = FALSE])
        if (nrow(e) < 2L) return(character())
        left <- e$end[-nrow(e)]
        right <- e$start[-1L] - 1L
        if (any(left >= right)) stop_input("Overlapping or invalid exons in transcript ", transcript_name, ".")
        paste(left, right, sep = ":")
    }
    transcript_introns <- lapply(split(exons, exons$tx), make_introns)
    if (!(selected_tx %in% names(transcript_introns))) stop_input("Selected transcript has no GTF exons.")
    selected <- unique(exons[exons$tx == selected_tx, c("start", "end"), drop = FALSE])
    selected <- selected[order(selected$start, selected$end), , drop = FALSE]
    selected$number <- if (gene_gtf$strand == "+") seq_len(nrow(selected)) else rev(seq_len(nrow(selected)))
    nearest_exon <- function(boundary, positions) {
        distances <- abs(positions - boundary)
        closest <- which(distances == min(distances))
        labels <- paste0("exon ", selected$number[closest])
        paste0(paste(labels, collapse = " / "), if (distances[closest[1L]] == 0) " (exact)" else
                   paste0(" (nearest; ", distances[closest[1L]], " bp)"))
    }
    junctions$left_boundary_exon <- vapply(seq_len(nrow(junctions)), function(i) {
        if (junctions$chrom[i] != gene_gtf$chrom) return("not available: other contig")
        nearest_exon(junctions$intron_start_0[i], selected$end)
    }, character(1))
    junctions$right_boundary_exon <- vapply(seq_len(nrow(junctions)), function(i) {
        if (junctions$chrom[i] != gene_gtf$chrom) return("not available: other contig")
        nearest_exon(junctions$intron_end_0[i] + 1L, selected$start)
    }, character(1))
    keys <- paste(junctions$intron_start_0, junctions$intron_end_0, sep = ":")
    all_introns <- unique(unlist(transcript_introns, use.names = FALSE))
    junctions$reference_status <- ifelse(!junctions$strand %in% c("+", "-"), "strand_unknown",
        ifelse(junctions$chrom != gene_gtf$chrom | junctions$strand != gene_gtf$strand, "other_strand_or_contig",
            ifelse(keys %in% transcript_introns[[selected_tx]], "selected_transcript",
                ifelse(keys %in% all_introns, "other_gene_transcript", "not_in_gene_annotation"))))
    junctions
}

# IDs are local to the result set; proximity uses 1-based intronic boundary bases.
annotate_junction_proximity <- function(junctions, variant, chrom, window = 100L) {
    keys <- c("chrom", "intron_start_0", "intron_end_0", "strand")
    if (anyDuplicated(junctions[keys])) stop_input("Duplicate junction coordinates and strand.")
    if (length(window) != 1L || is.na(window) || window < 0 || window != floor(window))
        stop_input("Invalid junction proximity window.")
    ord <- order(junctions$chrom, junctions$intron_start_0, junctions$intron_end_0, junctions$strand)
    ids <- character(nrow(junctions))
    ids[ord] <- paste0("J", seq_len(nrow(junctions)))
    junctions$junction_id <- ids
    variant_end <- variant$pos + nchar(variant$ref) - 1
    distance <- function(base) pmax(variant$pos - base, base - variant_end, 0)
    junctions$variant_distance_bp <- ifelse(junctions$chrom == chrom,
        pmin(distance(junctions$intron_start_0 + 1), distance(junctions$intron_end_0)), NA_real_)
    junctions$near_variant <- !is.na(junctions$variant_distance_bp) & junctions$variant_distance_bp <= window
    junctions[, c("junction_id", setdiff(names(junctions), "junction_id")), drop = FALSE]
}

read_mapped_alignments <- function(path) {
    lines <- readLines(path, warn = FALSE)
    fields <- strsplit(lines, "\t", fixed = TRUE)
    if (!length(fields) || any(lengths(fields) != 4L)) stop_input("Invalid samtools idxstats output: ", path)
    values <- do.call(rbind, fields)
    if (any(!grepl("^[0-9]+$", values[, 2L])) || any(!grepl("^[0-9]+$", values[, 3L])) ||
        any(!grepl("^[0-9]+$", values[, 4L])) || anyDuplicated(values[, 1L]) ||
        sum(values[, 1L] == "*") != 1L) stop_input("Invalid samtools idxstats output: ", path)
    mapped <- as.numeric(values[, 3L])
    if (any(!is.finite(mapped))) stop_input("Invalid mapped alignment count: ", path)
    sum(mapped)
}

normalize_junctions <- function(junctions, labels, mapped) {
    if (length(labels) != length(mapped) || anyDuplicated(labels) ||
        any(!labels %in% names(junctions)) || any(!is.finite(mapped) | mapped < 0))
        stop_input("Invalid junction normalization inputs.")
    for (i in seq_along(labels)) {
        column <- paste0(labels[i], "_per_million_mapped")
        if (column %in% names(junctions)) stop_input("Duplicate normalized junction column: ", column)
        junctions[[column]] <- if (mapped[i] == 0) rep(NA_real_, nrow(junctions)) else
            junctions[[labels[i]]] * 1e6 / mapped[i]
    }
    junctions
}

run_checked <- function(command, args, stdout, stderr) {
    cat("Running:", command, paste(shQuote(args), collapse = " "), "\n")
    status <- system2(command, args = shQuote(args), stdout = stdout, stderr = stderr)
    if (!identical(status, 0L)) stop_input("Command failed (", status, "): ", command, ". See ", stderr)
}
