#!/usr/bin/env Rscript
args <- grep("^--file=", commandArgs(), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)))
source(file.path(base, "scripts", "clinical_core.R"))
source(file.path(base, "scripts", "report_core.R"))
expect_error <- function(expr, expected) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    if (is.null(error) || !grepl(expected, conditionMessage(error), fixed = TRUE)) stop("Expected error: ", expected)
}
tmp <- tempfile("report_test_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
vep <- file.path(tmp, "vep.vcf")
writeLines(c('##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: SYMBOL|Consequence|Feature|HGVSc|IMPACT">',
             '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
             'chr17\t150\t.\tA\tG\t.\t.\tCSQ=GENE|other|NM_2.1|x|LOW,GENE|splice_donor_variant|NM_1.1|c.4+1A>G|HIGH'), vep)
v <- parse_vep(vep, "GENE", "NM_1.1")
stopifnot(nrow(v) == 1L, v$Consequence == "splice_donor_variant", v$HGVSc == "c.4+1A>G",
          parse_vep(vep, "GENE", "NM_missing")$status == "Selected transcript not annotated by VEP")
ai <- parse_spliceai("SpliceAI=G|GENE|0.02|0|1|0.5|-12|0|3|2", "GENE")
stopifnot(ai$DS_DG == "1", ai$DP_DG == "3",
          parse_spliceai(".", "GENE")$status == "SpliceAI annotation unavailable")
expect_error(parse_spliceai("SpliceAI=G|GENE|0.1", "GENE"), "Invalid SpliceAI")
asefile <- file.path(tmp, "counts.csv")
writeLines(c("contig,position,refAllele,altAllele,refCount,altCount",
             "chr17,150,A,G,3,2"), asefile)
ase <- read_ase_counts(asefile, "chr17", 150L, "A", "G")
stopifnot(ase$status == "low_depth_exploratory", ase$total == 5, ase$fraction == 0.4)
writeLines(c("contig,position,refAllele,altAllele,refCount,altCount",
             "chr19,51024634,C,T,12,8"), asefile)
ase_t <- read_ase_counts(asefile, "chr19", 51024634L, "C", "T")
stopifnot(ase_t$status == "counts_available_exploratory", ase_t$ref == 12,
          ase_t$alt == 8, ase_t$fraction == 0.4)
writeLines(c("contig,position,refAllele,altAllele,refCount,altCount",
             "chr17,150,A,G,3,2"), asefile)
expect_error(read_ase_counts(asefile, "chr17", 151L, "A", "G"), "variant mismatch")
writeLines(c("contig,position,refAllele,altAllele,refCount,altCount", "chr17,150,A,G,0,0"), asefile)
stopifnot(read_ase_counts(asefile, "chr17", 150L, "A", "G")$status == "zero_informative_reads")
for (value in c("Inf", "NA", "-1", "1.5")) {
    writeLines(c("contig,position,refAllele,altAllele,refCount,altCount",
                 paste0("chr17,150,A,G,", value, ",1")), asefile)
    expect_error(read_ase_counts(asefile, "chr17", 150L, "A", "G"), "Invalid ASE counts")
}
writeLines(c("contig,position,refAllele,altAllele,refCount,altCount", "chr17,NA,A,G,1,1"), asefile)
expect_error(read_ase_counts(asefile, "chr17", 150L, "A", "G"), "variant mismatch")
writeLines(c("contig,position,refAllele,altAllele,refCount,altCount",
             "chr17,150,A,G,1,1", "chr17,150,A,G,1,1"), asefile)
expect_error(read_ase_counts(asefile, "chr17", 150L, "A", "G"), "multiple records")
writeLines("contig,position,refAllele,altAllele,refCount,altCount", asefile)
stopifnot(read_ase_counts(asefile, "chr17", 150L, "A", "G")$status == "not_reported_by_ASEReadCounter")
stopifnot(html_escape("<&\"") == "&lt;&amp;&quot;")
stopifnot(identical(unname(vapply(c("0", "0.199", "0.2", "0.499", "0.5", "0.799", "0.8", "1", ".", "1.1"),
                           score_band, "")), c("low", "low", "recall", "recall", "recommended", "recommended",
                                             "precision", "precision", "unavailable", "unavailable")))
samples <- data.frame(label = c("minus", "plus"), region_reads = c(10, 20))
junctions <- data.frame(chrom = "chr17", intron_start_0 = 150L, intron_end_0 = 200L, strand = "+", minus = 2L, plus = 3L,
                        left_boundary_exon = "exon 1 (exact)", right_boundary_exon = "exon 2 (nearest; 4 bp)")
junctions <- annotate_junction_proximity(junctions, parse_variant("17-150-A-G"), "chr17")
stopifnot(junctions$junction_id == "J1", junctions$variant_distance_bp == 1, junctions$near_variant,
          plot_junction_labels(junctions) == "J1  E1 - E2 (~4bp) [+] *")
depth <- data.frame(sample = "minus", position_1 = 150L, depth = 1L)
ase_table <- data.frame(sample = "minus", status = "not_reported")
plot_junctions(file.path(tmp, "junction_counts.svg"), junctions, samples$label)
write_report(file.path(tmp, "index.html"), "GENE", "17-150-A-G", "NM_1.1", "het",
             v, ai, samples, depth, ase_table, junctions)
html <- paste(readLines(file.path(tmp, "index.html"), warn = FALSE), collapse = "")
stopifnot(grepl("splice_donor_variant", html, fixed = TRUE), grepl("c.4+1A&gt;G", html, fixed = TRUE),
          grepl("junction_counts.svg", html, fixed = TRUE), grepl("sashimi.pdf", html, fixed = TRUE),
          grepl("junction_normalized.svg", html, fixed = TRUE),
          grepl('class="score-low"', html, fixed = TRUE), grepl('class="score-precision"', html, fixed = TRUE),
          grepl("left_boundary_exon", html, fixed = TRUE), grepl("exon 1 (exact)", html, fixed = TRUE),
          grepl("exon 2 (nearest; 4 bp)", html, fixed = TRUE),
          grepl('<tr class="near-variant">', html, fixed = TRUE),
          grepl("J1", html, fixed = TRUE), grepl("variant_distance_bp", html, fixed = TRUE),
          file.info(file.path(tmp, "junction_counts.svg"))$size > 0)
write_report(file.path(tmp, "igv_report.html"), "GENE", "17-150-A-G", "NM_1.1", "het",
             v, ai, samples, depth, ase_table, junctions,
             igv_snapshot = list(status = "available", detail = "chr17:1-100"))
igv_html <- paste(readLines(file.path(tmp, "igv_report.html"), warn = FALSE), collapse = "")
stopifnot(grepl('src="igv_snapshot.svg"', igv_html, fixed = TRUE),
          grepl("chr17:1-100", igv_html, fixed = TRUE))
large <- junctions[rep(1L, 45L), ]
large$junction_id <- paste0("J", seq_len(nrow(large)))
large$minus <- 45:1
large$plus <- 0L
large$near_variant <- seq_len(nrow(large)) > 20L
large$variant_distance_bp <- seq_len(nrow(large))
selected <- select_junctions(large, samples$label, 20L, 20L)
stopifnot(nrow(selected) == 40L, "J40" %in% selected$junction_id,
          !("J41" %in% selected$junction_id), nrow(select_junctions(large, samples$label, 30L)) == 45L)
plot_junctions(file.path(tmp, "many.svg"), large, samples$label)
plot_junctions(file.path(tmp, "empty.svg"), junctions[FALSE, ], samples$label)
large <- normalize_junctions(large, samples$label, c(0, 2000000))
plot_junctions(file.path(tmp, "normalized.svg"), large, samples$label, normalized = TRUE)
large <- normalize_junctions(large[, setdiff(names(large), paste0(samples$label, "_per_million_mapped"))],
                             samples$label, c(0, 0))
plot_junctions(file.path(tmp, "unavailable.svg"), large, samples$label, normalized = TRUE)
stopifnot(file.info(file.path(tmp, "many.svg"))$size > 0, file.info(file.path(tmp, "empty.svg"))$size > 0,
          file.info(file.path(tmp, "normalized.svg"))$size > 0, file.info(file.path(tmp, "unavailable.svg"))$size > 0)
cat("ClinV1 report tests passed\n")
