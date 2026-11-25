#!/usr/bin/env Rscript
# Hovedforløb: målrettet analyse af ét patientpar og eventuelle referencer.

# Indlæs fælles funktioner og kontrollér argumenter og afhængigheder.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript scripts/ClinV1.R.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
source(file.path(script_dir, "clinical_core.R"))
source(file.path(script_dir, "report_core.R"))
source(file.path(script_dir, "sai_core.R"))
source(file.path(script_dir, "igv_snapshot_core.R"))
cli_args <- commandArgs(trailingOnly = TRUE)
if (!length(cli_args) || identical(cli_args, "--help") || identical(cli_args, "-h")) {
    cat(clinv1_help(), "\n")
    quit(save = "no", status = 0L)
}
opt <- parse_options(cli_args)
variant <- parse_variant_input(opt[["-variant"]])
gene <- opt[["-gene"]]
tx <- opt[["-tx"]]
genotype <- if (is.null(opt[["-gt"]])) "unknown" else tolower(opt[["-gt"]])
if (!genotype %in% c("het", "hmz", "hemi", "unknown")) stop_input("-gt must be het, hmz or hemi.")
ase_eligible <- genotype == "het" && nchar(variant$ref) == 1L && nchar(variant$alt) == 1L
sai_enabled <- !is.null(opt[["-sai-module"]])
if (!sai_enabled && !is.null(opt[["-sai-env"]])) stop_input("-sai-env requires -sai-module.")
sai_env <- if (is.null(opt[["-sai-env"]])) "KGA_Splice" else opt[["-sai-env"]]
if (!grepl("^[A-Za-z0-9_.-]+$", sai_env)) stop_input("Invalid -sai-env name.")
if (sai_enabled && (!file.exists(opt[["-sai-module"]]) || dir.exists(opt[["-sai-module"]]) ||
                    file.info(opt[["-sai-module"]])$size == 0)) stop_input("SAI module missing or empty: ", opt[["-sai-module"]])
if (!grepl("^[A-Za-z0-9_.-]+$", gene) || !grepl("^[A-Za-z0-9_.-]+$", tx)) stop_input("Invalid gene or transcript.")
ids <- parse_sample_ids(opt[["-s"]])
if (!dir.exists(opt[["-bam"]])) stop_input("BAM directory not found.")
if (!is.null(opt[["-p"]]) && !dir.exists(opt[["-p"]])) stop_input("Output parent directory not found.")
if (!dir.exists(opt[["-vep"]])) stop_input("VEP cache not found.")
if (!requireNamespace("xml2", quietly = TRUE)) stop_input("R package xml2 is required.")
check_required_commands(c("samtools", "sambamba", "regtools", "conda", "bedtools", "bash", if (ase_eligible) "gatk"))
for (path in c(opt[["-g"]], paste0(opt[["-g"]], ".fai"), opt[["-gtf"]])) {
    if (!file.exists(path)) stop_input("Required reference file not found: ", path)
}

# Find BAM'er og kontrollér, at gen, transkript og reference passer sammen.
patient <- vapply(ids, function(id) find_sample_bam(opt[["-bam"]], id), character(1))
controls <- character()
if (!is.null(opt[["-refs"]])) {
    if (!dir.exists(opt[["-refs"]])) stop_input("Reference BAM directory not found.")
    controls <- normalizePath(list.files(opt[["-refs"]], pattern = "\\.bam$", full.names = TRUE), mustWork = TRUE)
    if (!length(controls)) stop_input("No BAM files found in reference directory.")
}
bams <- c(unname(patient), controls)
if (anyDuplicated(bams)) stop_input("The same BAM cannot be used twice.")
invisible(lapply(bams, check_bam_index))
control_labels <- if (length(controls)) paste0("control_", seq_along(controls)) else character()
single_sample <- length(ids) == 1L
patient_labels <- if (single_sample) paste0(ids, "_treatment_unknown") else paste0(ids, c("_nt", "_puro"))
labels <- c(patient_labels, control_labels)
reference <- read_gene_gtf(opt[["-gtf"]], gene, tx)
fai <- read.delim(paste0(opt[["-g"]], ".fai"), header = FALSE, stringsAsFactors = FALSE)
if (ncol(fai) < 2L) stop_input("Invalid FASTA index.")
chrom <- choose_contig(variant$chrom, fai[[1L]])
if (reference$chrom != chrom) stop_input("GTF and FASTA chromosome names differ: ", reference$chrom, " vs ", chrom)
if (variant$pos < reference$start || variant$pos > reference$end) stop_input("Variant is outside the selected gene coordinates.")
chrom_length <- fai[[2L]][match(chrom, fai[[1L]])]
if (variant$pos + nchar(variant$ref) - 1L > chrom_length) stop_input("Variant exceeds chromosome length.")

# Opret et nyt projekt; STATUS forbliver INCOMPLETE, hvis et trin fejler.
project_name <- if (single_sample) paste0(gene, "_", ids[1L], "_treatment_unknown") else
    paste0(gene, "_", ids[1L], "_nt_", ids[2L], "_puro")
project <- prepare_project(opt[["-p"]], opt[["-project-dir"]], project_name)
writeLines("INCOMPLETE: run did not finish. Inspect console and Logs.", file.path(project, "STATUS.txt"))
for (dir in c("BAM", "Junctions", "IGV", "Logs", "Coverage", "ASE", "Report")) {
    if (!dir.create(file.path(project, dir))) stop_input("Cannot create output directory: ", dir)
}
logs <- file.path(project, "Logs")
run_checked("conda", c("run", "--no-capture-output", "-n", "KGA_Splice", "sashimipy", "--help"),
            file.path(logs, "sashimi_check.out"), file.path(logs, "sashimi_check.err"))
if (sai_enabled) run_checked("conda", c("run", "--no-capture-output", "-n", sai_env, "python",
    file.path(script_dir, "sai_lookup.py"), "--module", opt[["-sai-module"]], "--check"),
    file.path(logs, "sai_check.out"), file.path(logs, "sai_check.err"))
gene_gtf <- file.path(project, paste0(gene, ".reference.gtf"))
writeLines(reference$lines, gene_gtf)
start <- max(1, reference$start - 10000)
end <- min(chrom_length, reference$end + 10000)
region <- paste0(chrom, ":", start, "-", end)
bed <- file.path(project, paste0(gene, ".regions.bed"))
writeLines(paste(chrom, start - 1, end, sep = "\t"), bed)
variant_bed <- file.path(project, paste0(gene, ".variant.bed"))
writeLines(paste(chrom, variant$pos - 1L, variant$pos - 1L + nchar(variant$ref), "variant", sep = "\t"), variant_bed)
sizes <- file.path(project, "Coverage", "chrom.sizes")
write.table(fai[, 1:2], sizes, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Kontrollér REF og lav fælles variantannotering én gang for patientparret.
ref_fa <- file.path(logs, "reference_allele.fa")
run_checked("samtools", c("faidx", opt[["-g"]], paste0(chrom, ":", variant$pos, "-", variant$pos + nchar(variant$ref) - 1L)),
            ref_fa, file.path(logs, "reference_allele.err"))
ref_lines <- readLines(ref_fa, warn = FALSE)
ref_seq <- toupper(paste(ref_lines[!startsWith(ref_lines, ">")], collapse = ""))
if (!identical(ref_seq, variant$ref)) stop_input("REF allele does not match FASTA at ", chrom, ":", variant$pos, ".")
vcf <- file.path(project, paste0(gene, ".vcf"))
resolved_variant <- paste(chrom, variant$pos, variant$ref, variant$alt, sep = "-")
writeLines(c(paste("Input:", opt[["-variant"]]), paste("Resolved variant:", resolved_variant),
             "Genome assembly: GRCh38; input VCF genotype is not used."),
           file.path(project, "variant_input.txt"))
writeLines(c("##fileformat=VCFv4.2", paste0("##contig=<ID=", chrom, ",length=", chrom_length, ">"),
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
             paste(chrom, variant$pos, ".", variant$ref, variant$alt, ".", ".", ".", sep = "\t")), vcf)

vep <- file.path(project, paste0(gene, ".ann.vep.vcf"))
spliceai <- file.path(project, paste0(gene, ".spliceai.vcf"))
run_checked("conda", c("run", "--no-capture-output", "-n", "vep_env", "vep", "-i", vcf, "-o", vep,
                       "--cache", "--offline", "--dir_cache", opt[["-vep"]], "--refseq", "--vcf",
                       "--fasta", opt[["-g"]], "--assembly", "GRCh38", "--everything", "--no_stats"),
            file.path(logs, "vep.out"), file.path(logs, "vep.err"))
run_checked("conda", c("run", "--no-capture-output", "-n", "spliceai_env", "spliceai", "-I", vcf,
                       "-O", spliceai, "-R", opt[["-g"]], "-A", "grch38", "-D", "150"),
            file.path(logs, "spliceai.out"), file.path(logs, "spliceai.err"))
vep_fields <- validate_annotated_vcf(vep, variant, chrom)
ai_fields <- validate_annotated_vcf(spliceai, variant, chrom)
vep_status <- if (has_info_annotation(vep_fields[8L], "CSQ")) "annotation present" else "not annotated; inspect logs and cache"
ai_status <- if (has_info_annotation(ai_fields[8L], "SpliceAI")) "annotation present" else "not annotated; inspect logs and reference compatibility"
sai <- if (sai_enabled) run_sai(opt[["-sai-module"]], sai_env, vcf, opt[["-g"]],
                                tx, reference, variant, project, logs, file.path(script_dir, "sai_lookup.py")) else
    list(status = "not_enabled", data = data.frame())

# ASE kræver en brugerangivet heterozygot SNP.
ase_vcf <- NULL
if (ase_eligible) {
    fasta <- opt[["-g"]]
    dict <- sub("\\.(fa|fasta|fna)$", ".dict", fasta, ignore.case = TRUE)
    if (identical(dict, fasta) || !file.exists(dict)) stop_input("GATK reference dictionary is required for ASE: ", dict)
    ase_vcf <- file.path(project, "ASE", "patient_het.vcf")
    writeLines(c("##fileformat=VCFv4.2", paste0("##contig=<ID=", chrom, ",length=", chrom_length, ">"),
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPATIENT",
                 paste(chrom, variant$pos, ".", variant$ref, variant$alt, ".", ".", ".", "GT", "0/1", sep = "\t")), ase_vcf)
    run_checked("gatk", c("IndexFeatureFile", "-I", ase_vcf), file.path(logs, "ase_vcf_index.out"),
                file.path(logs, "ase_vcf_index.err"))
}
# Udtræk og analysér kun genregionen for hver patient- og reference-BAM.
tables <- vector("list", length(bams))
clipped <- character(length(bams))
depth_tables <- vector("list", length(bams))
region_counts <- numeric(length(bams))
mapped_counts <- numeric(length(bams))
coverage_plus <- coverage_minus <- character(length(bams))
coverage_plus_reads <- coverage_minus_reads <- numeric(length(bams))
ase_results <- vector("list", length(ids))
for (i in seq_along(bams)) {
    label <- labels[i]
    cat("\nProcessing ", label, ": ", bams[i], "\n", sep = "")
    header <- file.path(logs, paste0(label, "_header.sam"))
    run_checked("samtools", c("view", "-H", bams[i]), header, file.path(logs, paste0(label, "_header.err")))
    header_lines <- readLines(header, warn = FALSE)
    sq <- strsplit(header_lines[startsWith(header_lines, "@SQ\t")], "\t", fixed = TRUE)
    matching_sq <- Filter(function(x) paste0("SN:", chrom) %in% x, sq)
    if (length(matching_sq) != 1L) stop_input("BAM ", label, " does not contain exactly one matching chromosome ", chrom)
    ln <- matching_sq[[1L]][startsWith(matching_sq[[1L]], "LN:")]
    if (length(ln) != 1L || is.na(suppressWarnings(as.numeric(sub("^LN:", "", ln)))) ||
        as.numeric(sub("^LN:", "", ln)) != chrom_length) stop_input("BAM ", label, " chromosome length does not match FASTA.")
    run_checked("samtools", c("quickcheck", bams[i]), file.path(logs, paste0(label, "_input_quickcheck.out")),
                file.path(logs, paste0(label, "_input_quickcheck.err")))
    idxstats <- file.path(logs, paste0(label, "_idxstats.tsv"))
    run_checked("samtools", c("idxstats", bams[i]), idxstats, file.path(logs, paste0(label, "_idxstats.err")))
    mapped_counts[i] <- read_mapped_alignments(idxstats)
    clipped[i] <- file.path(project, "BAM", paste0(label, ".", gene, ".bam"))
    run_checked("samtools", c("view", "-b", bams[i], region), clipped[i],
                file.path(logs, paste0(label, "_extract.err")))
    run_checked("samtools", c("quickcheck", clipped[i]), file.path(logs, paste0(label, "_quickcheck.out")),
                file.path(logs, paste0(label, "_quickcheck.err")))
    run_checked("sambamba", c("index", clipped[i]), file.path(logs, paste0(label, "_index.out")),
                file.path(logs, paste0(label, "_index.err")))
    run_checked("sambamba", c("view", "-c", clipped[i]), file.path(logs, paste0(label, "_region_read_count.txt")),
                file.path(logs, paste0(label, "_region_read_count.err")))
    count_text <- readLines(file.path(logs, paste0(label, "_region_read_count.txt")), warn = FALSE)
    if (length(count_text) != 1L || !grepl("^[0-9]+$", count_text)) stop_input("Invalid region read count for ", label)
    region_counts[i] <- as.numeric(count_text)
    for (strand in c("plus", "minus")) {
        count_file <- file.path(logs, paste0(label, "_coverage_", strand, "_reads.txt"))
        args <- if (strand == "plus") c("view", "-c", "-d", "XS:+", "-F", "16", clipped[i]) else
            c("view", "-c", "-d", "XS:-", "-f", "16", clipped[i])
        run_checked("samtools", args, count_file, file.path(logs, paste0(label, "_coverage_", strand, "_reads.err")))
        count <- readLines(count_file, warn = FALSE)
        if (length(count) != 1L || !grepl("^[0-9]+$", count)) stop_input("Invalid coverage read count for ", label, " ", strand)
        if (strand == "plus") coverage_plus_reads[i] <- as.numeric(count) else coverage_minus_reads[i] <- as.numeric(count)
    }
    if (coverage_plus_reads[i] + coverage_minus_reads[i] > region_counts[i]) stop_input("Coverage counts exceed region reads for ", label)
    coverage_plus[i] <- file.path(project, "Coverage", paste0(label, ".plus.bedgraph"))
    coverage_minus[i] <- file.path(project, "Coverage", paste0(label, ".minus.bedgraph"))
    run_checked("bash", c(file.path(script_dir, "stranded_coverage.sh"), clipped[i], sizes, coverage_plus[i], coverage_minus[i]),
                file.path(logs, paste0(label, "_coverage.out")), file.path(logs, paste0(label, "_coverage.err")))
    depth_file <- file.path(logs, paste0(label, "_variant_depth.tsv"))
    locus <- paste0(chrom, ":", variant$pos, "-", variant$pos + nchar(variant$ref) - 1L)
    run_checked("samtools", c("depth", "-aa", "-r", locus, clipped[i]), depth_file,
                file.path(logs, paste0(label, "_variant_depth.err")))
    depth <- read.delim(depth_file, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(depth) != 3L || nrow(depth) != nchar(variant$ref) ||
        any(depth[[1L]] != chrom) || any(depth[[2L]] != seq.int(variant$pos, length.out = nchar(variant$ref))) ||
        anyNA(depth[[3L]]) || any(depth[[3L]] < 0)) stop_input("Invalid variant depth output for ", label)
    depth_tables[[i]] <- data.frame(sample = label, chrom = depth[[1L]], position_1 = depth[[2L]], depth = depth[[3L]])
    junctions <- file.path(project, "Junctions", paste0(label, ".bed"))
    run_checked("regtools", c("junctions", "extract", "-s", "XS", "-a", "8", "-m", "50", "-M", "500000", clipped[i]),
                junctions, file.path(logs, paste0(label, "_regtools.err")))
    tables[[i]] <- read_regtools_junctions(junctions, label)
    if (i <= length(ids)) {
        if (ase_eligible) {
            ase_file <- file.path(project, "ASE", paste0(label, ".ase_counts.csv"))
            run_checked("gatk", c("ASEReadCounter", "-R", opt[["-g"]], "-V", ase_vcf, "-I", clipped[i],
                                   "-O", ase_file, "--output-format", "CSV"),
                        file.path(logs, paste0(label, "_ase.out")), file.path(logs, paste0(label, "_ase.err")))
            counts <- read_ase_counts(ase_file, chrom, variant$pos, variant$ref, variant$alt)
        } else {
            status <- if (genotype == "unknown") "genotype_not_provided" else if (genotype != "het") "not_heterozygous" else "indel_not_supported"
            counts <- list(status = status, ref = NA_integer_, alt = NA_integer_, total = NA_integer_, fraction = NA_real_)
        }
        ase_results[[i]] <- data.frame(sample = label, status = counts$status, ref_count = counts$ref,
                                       alt_count = counts$alt, informative_count = counts$total, alt_fraction = counts$fraction)
    }
}
# Saml junctiontællinger, dækning og alleltællinger.
comparison <- annotate_junction_proximity(
    annotate_junctions(compare_junctions(tables, labels), reference, tx), variant, chrom)
comparison <- normalize_junctions(comparison, labels, mapped_counts)
write.table(comparison, file.path(project, "Junctions", "junction_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(do.call(rbind, depth_tables), file.path(project, "variant_depth.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
ase <- do.call(rbind, ase_results)
write.table(ase, file.path(project, "ASE", "allele_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
samples <- data.frame(label = labels, condition = c(if (single_sample) "treatment_unknown" else c("nt", "puro"),
                                                   rep("control_unspecified", length(controls))),
                      input_bam = bams, region_bam = clipped, region_reads = region_counts,
                      mapped_alignments_idxstats = mapped_counts,
                      xs_plus_flag_plus = coverage_plus_reads, xs_minus_flag_minus = coverage_minus_reads,
                      excluded_from_stranded = region_counts - coverage_plus_reads - coverage_minus_reads)
write.table(samples, file.path(project, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# IGV: regionsspor først og originale input-BAM'er nederst.
session <- xml2::xml_new_root("Session", genome = "hg38", locus = paste0(chrom, ":", max(1, variant$pos - 1000), "-", variant$pos + 1000), version = "8")
resources <- xml2::xml_add_child(session, "Resources")
coverage_paths <- c(coverage_plus, coverage_minus)
coverage_paths <- coverage_paths[file.info(coverage_paths)$size > 0]
for (path in unique(c(variant_bed, gene_gtf, clipped, coverage_paths, bams))) xml2::xml_add_child(resources, "Resource", path = path)
panel <- xml2::xml_add_child(session, "Panel", name = "DataPanel", height = "1000")
xml2::xml_add_child(panel, "Track", clazz = "org.broad.igv.track.FeatureTrack", id = variant_bed, name = "Variant",
                    color = "255,0,0", height = "50", displayMode = "EXPANDED", featureVisibilityWindow = "-1")
xml2::xml_add_child(panel, "Track", clazz = "org.broad.igv.track.FeatureTrack", id = gene_gtf, name = paste("REFERENCE", gene, tx))
track_names <- c(if (single_sample) paste(ids, "(treatment unknown)") else paste(ids, c("-puro", "+puro")),
                 if (length(controls)) paste("Control", tools::file_path_sans_ext(basename(controls)), "(treatment unknown)") else character())
track_colors <- c(if (single_sample) "38,76,151" else c("38,76,151", "94,149,214"),
                  rep("0,128,95", length(controls)))
for (i in seq_along(clipped)) {
    for (strand in c("plus", "minus")) {
        path <- if (strand == "plus") coverage_plus[i] else coverage_minus[i]
        if (file.info(path)$size > 0) xml2::xml_add_child(panel, "Track", clazz = "org.broad.igv.track.DataTrack",
            id = path, name = paste(track_names[i], if (strand == "plus") "+ strand" else "- strand"),
            color = track_colors[i], height = "80", autoScale = "true")
    }
    xml2::xml_add_child(panel, "Track", clazz = "org.broad.igv.sam.CoverageTrack", id = paste0(clipped[i], "_coverage"), name = paste(track_names[i], "coverage"), color = track_colors[i])
    xml2::xml_add_child(panel, "Track", clazz = "org.broad.igv.sam.SpliceJunctionTrack", id = paste0(clipped[i], "_junctions"), name = paste(track_names[i], "junctions"), color = track_colors[i], altColor = track_colors[i])
    read_track <- xml2::xml_add_child(panel, "Track", clazz = "org.broad.igv.sam.AlignmentTrack", id = clipped[i],
                                     name = paste(track_names[i], "BAM"), color = track_colors[i],
                                     altColor = track_colors[i], experimentType = "RNA", visible = "false")
    xml2::xml_add_child(read_track, "RenderOptions")
}
for (i in seq_along(bams)) {
    xml2::xml_add_child(panel, "Track", clazz = "org.broad.igv.sam.CoverageTrack", id = paste0(bams[i], "_coverage"), name = paste(track_names[i], "raw coverage"), color = track_colors[i])
    xml2::xml_add_child(panel, "Track", clazz = "org.broad.igv.sam.SpliceJunctionTrack", id = paste0(bams[i], "_junctions"), name = paste(track_names[i], "raw junctions"), color = track_colors[i], altColor = track_colors[i])
    read_track <- xml2::xml_add_child(panel, "Track", clazz = "org.broad.igv.sam.AlignmentTrack", id = bams[i],
                                     name = paste(track_names[i], "raw BAM"), color = track_colors[i],
                                     altColor = track_colors[i], experimentType = "RNA", visible = "false")
    xml2::xml_add_child(read_track, "RenderOptions")
}
xml2::write_xml(session, file.path(project, "IGV", "paired_IGV.xml"))
write_vswarehouse_igv_session(file.path(project, "IGV", "paired_IGV.xml"),
                              file.path(project, "IGV", "paired_IGV_vswarehouse.xml"))
snapshot_locus <- paste0(chrom, ":", max(1, variant$pos - 1000), "-", min(chrom_length, variant$pos + 1000))
igv_snapshot <- render_igv_snapshot(opt[["-igv-snapshot"]], file.path(project, "IGV", "paired_IGV.xml"),
                                    opt[["-g"]], project, snapshot_locus)
# Generér rapport og afslut først status efter vellykket visualisering.
plot_junctions(file.path(project, "Report", "junction_counts.svg"), comparison, labels)
plot_junctions(file.path(project, "Report", "junction_normalized.svg"), comparison, labels, normalized = TRUE)
run_sashimi(file.path(project, "Report", "sashimi.pdf"), clipped, labels, gene_gtf, reference, tx,
            variant$pos, start, end, logs)
write_report(file.path(project, "Report", "index.html"), gene, opt[["-variant"]], tx, genotype,
             parse_vep(vep, gene, tx), parse_spliceai(ai_fields[8L], gene), samples,
             do.call(rbind, depth_tables), ase, comparison, sai, igv_snapshot)
unknown_strand <- sum(!comparison$strand %in% c("+", "-"))
writeLines(c("ClinV1 research/validation candidate - not approved for clinical use",
             paste("Variant:", opt[["-variant"]]), paste("Gene:", gene), paste("Transcript requested:", tx),
             paste("Resolved genomic variant:", resolved_variant),
             if (single_sample) "Single patient sample: treatment unknown; no nt/puro comparison." else
                 "Patient samples: first nt, second puro; verify treatment metadata.",
             paste("Region:", region), paste("VEP:", vep_status), paste("SpliceAI:", ai_status),
             paste("SAI-10k-calc:", sai$status),
             paste("IGV snapshot:", igv_snapshot$status, igv_snapshot$detail),
             paste("Samples with zero region reads:", if (any(region_counts == 0)) paste(labels[region_counts == 0], collapse = ", ") else "none"),
             paste("Samples with zero depth at any variant REF base:",
                   if (any(vapply(depth_tables, function(x) any(x$depth == 0), logical(1))))
                       paste(labels[vapply(depth_tables, function(x) any(x$depth == 0), logical(1))], collapse = ", ") else "none"),
             paste("Distinct junctions:", nrow(comparison)), paste("Junctions with unknown strand:", unknown_strand),
             paste("Samples with zero mapped alignments (normalization unavailable):",
                   if (any(mapped_counts == 0)) paste(labels[mapped_counts == 0], collapse = ", ") else "none"),
             "Junctions per million mapped alignments are descriptive, not PSI or a statistical comparison.",
             "Junction counts are descriptive, not PSI or a statistical comparison.",
             "Junction reference status uses exact boundaries and strand in the supplied GTF only.",
             "Zero junction reads do not exclude an event; inspect region counts, local coverage and assay sensitivity.",
             "Input quickcheck checks header and EOF only; it does not scan for internal BAM corruption.",
             "Variant depth uses samtools depth defaults; it does not establish allele-specific coverage or assay adequacy.",
             "Controls are not assigned a treatment or assumed independent; verify metadata.",
             "Stranded bedgraphs use XS plus matching alignment strand as in v6.4; they exclude other reads and are not validated transcript strand.",
             "ASE is exploratory, target heterozygous SNP only; indels and other loci are not assessed.",
             "No automatic exon skipping, retention, pseudoexon, frame, NMD or clinical classification is performed.",
             "RNA sample identity, transcript relevance and genotype must be independently verified.",
             "Regtools uses XS strand tags; missing/incorrect tags can invalidate strand-dependent interpretation.",
             "Sashimipy assumes fr-firststrand and XS tags, as in v6.4; junction plot threshold is 20 reads.",
             "IGV paths are container paths and may require translation on the IGV host."),
           file.path(project, "STATUS.txt"))
cat("\nClinV1 candidate results: ", project, "\nReview STATUS.txt before interpretation.\n", sep = "")
