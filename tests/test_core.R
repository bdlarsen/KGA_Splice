#!/usr/bin/env Rscript
args <- grep("^--file=", commandArgs(), value = TRUE)
base <- dirname(dirname(normalizePath(sub("^--file=", "", args[1L]), mustWork = TRUE)))
source(file.path(base, "scripts", "clinical_core.R"))
expect_error <- function(expr, expected) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    if (is.null(error) || !grepl(expected, conditionMessage(error), fixed = TRUE)) {
        stop("Expected error containing: ", expected)
    }
}
check_required_commands(c("samtools", "regtools"), function(x) setNames(c("/usr/bin/samtools", "/usr/bin/regtools"), x))
expect_error(check_required_commands(c("samtools", "regtools"), function(x) setNames(c("/usr/bin/samtools", ""), x)),
             "Required commands not found on PATH: regtools")
valid_args <- c("-variant", "17-100-A-G", "-gene", "GENE", "-tx", "NM_1.1",
                "-p", "out", "-bam", "bam", "-s", "1,2")
stopifnot(parse_options(valid_args)[["-s"]] == "1,2")
project_args <- valid_args
project_args[project_args == "-p"] <- "-project-dir"
stopifnot(parse_options(project_args)[["-project-dir"]] == "out")
expect_error(parse_options(c(valid_args, "-project-dir", "out")), "exactly one")
expect_error(parse_options(valid_args[-c(7L, 8L)]), "exactly one")
expect_error(parse_options(c(valid_args, "-gene", "OTHER")), "Duplicate option")
expect_error(parse_options(c(valid_args, "-bad", "x")), "Unknown option")
expect_error(parse_options(c(valid_args, "-gene")), "Each option")
stopifnot(identical(parse_sample_ids("040966-2498-B,040966-2498-A"), c("040966-2498-B", "040966-2498-A")),
          identical(parse_sample_ids("1,2"), c("1", "2")),
          identical(parse_sample_ids("040966-2498"), "040966-2498"))
expect_error(parse_sample_ids("a,a"), "must be different")
expect_error(parse_sample_ids("../a,b"), "-s requires one ID")
expect_error(parse_sample_ids("a,b,c"), "-s requires one ID")
expect_error(parse_sample_ids("a, b"), "-s requires one ID")
stopifnot(parse_variant("chr17-100-ac-a")$ref == "AC")
expect_error(parse_variant("17-0-A-G"), "Invalid variant position")
expect_error(parse_variant("17-100-A-A"), "must differ")
expect_error(parse_variant("17-100-A-<DEL>"), "Variant must be")
stopifnot(choose_contig("17", c("chr17", "chr1")) == "chr17")
expect_error(choose_contig("17", c("chr17", "17")), "does not map uniquely")

tmp <- tempfile("clin_test_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE))
reserved <- file.path(tmp, "reserved")
dir.create(reserved)
dir.create(file.path(reserved, "tmp"))
writeLines("temporary", file.path(reserved, "tmp", "work"))
stopifnot(prepare_project(NULL, reserved, "unused") == normalizePath(reserved))
writeLines("existing", file.path(reserved, ".existing"))
expect_error(prepare_project(NULL, reserved, "unused"), "empty apart from tmp")
expect_error(prepare_project(NULL, file.path(tmp, "missing"), "unused"), "must exist")
standard <- prepare_project(tmp, NULL, "standard")
stopifnot(dir.exists(standard))
expect_error(prepare_project(tmp, NULL, "standard"), "Project already exists")
vcf_input <- file.path(tmp, "input.vcf")
header <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
record <- "chr17\t100\t.\tAC\tA\t.\t.\t."
writeLines(c("##fileformat=VCFv4.2", header, record), vcf_input)
stopifnot(identical(parse_variant_input(vcf_input), parse_variant("chr17-100-AC-A")))
gzip_input <- file.path(tmp, "input.vcf.gz")
con <- gzfile(gzip_input, "wt")
writeLines(c("##fileformat=VCFv4.2", header, record), con)
close(con)
stopifnot(identical(parse_variant_input(gzip_input), parse_variant("chr17-100-AC-A")))
writeLines(c(header, record, record), vcf_input)
expect_error(parse_variant_input(vcf_input), "exactly one variant")
writeLines(c(header, "chr17\t100\t.\tA\tG,T\t.\t.\t."), vcf_input)
expect_error(parse_variant_input(vcf_input), "multiallelic")
writeLines(c(header, "chr17\t100\t.\tA\t<DEL>\t.\t.\t."), vcf_input)
expect_error(parse_variant_input(vcf_input), "Variant must be")
writeLines(c(record), vcf_input)
expect_error(parse_variant_input(vcf_input), "#CHROM header")
expect_error(parse_variant_input(file.path(tmp, "missing.vcf")), "VCF file not found")
stopifnot(grepl("VCF/VCF.GZ", clinv1_help(), fixed = TRUE))
file.create(file.path(tmp, "1_x.bam"), file.path(tmp, "12_x.bam"))
stopifnot(basename(find_sample_bam(tmp, "1")) == "1_x.bam")
expect_error(find_sample_bam(tmp, "2"), "found 0")
file.create(file.path(tmp, "1.bam"))
expect_error(find_sample_bam(tmp, "1"), "found 2")
expect_error(check_bam_index(file.path(tmp, "1.bam")), "Missing BAM index")
file.create(file.path(tmp, "1.bam.bai"))
check_bam_index(file.path(tmp, "1.bam"))
file.create(file.path(tmp, "040966-2498-B_illumina.bam"), file.path(tmp, "040966-2498-A_illumina.bam"),
            file.path(tmp, "040966-2498-B2_illumina.bam"))
stopifnot(basename(find_sample_bam(tmp, "040966-2498-B")) == "040966-2498-B_illumina.bam",
          basename(find_sample_bam(tmp, "040966-2498-A")) == "040966-2498-A_illumina.bam")
file.create(file.path(tmp, "040966-2498-B_other.bam"))
expect_error(find_sample_bam(tmp, "040966-2498-B"), "found 2")
file.create(file.path(tmp, "G150-113248405352-B-PTEN_target.bam"),
            file.path(tmp, "G150-113248405352-A-PTEN_target.bam"),
            file.path(tmp, "G150-113248405352-B2_target.bam"))
stopifnot(basename(find_sample_bam(tmp, "113248405352-B")) == "G150-113248405352-B-PTEN_target.bam",
          basename(find_sample_bam(tmp, "113248405352-A")) == "G150-113248405352-A-PTEN_target.bam")
file.create(file.path(tmp, "113248405352-B.bam"))
expect_error(find_sample_bam(tmp, "113248405352-B"), "found 2")

gtf <- file.path(tmp, "test.gtf")
writeLines(c(
    'chr17\ttest\texon\t101\t150\t.\t+\t.\tgene_id "G001"; gene_name "GENE"; transcript_id "NM_1.1";',
    'chr17\ttest\texon\t201\t250\t.\t+\t.\tgene_id "G001"; gene_name "GENE"; transcript_id "NM_1.1";',
    'chr17\ttest\texon\t101\t150\t.\t+\t.\tgene_id "G001"; gene_name "GENE"; transcript_id "NM_2.1";',
    'chr17\ttest\texon\t301\t350\t.\t+\t.\tgene_id "G001"; gene_name "GENE"; transcript_id "NM_2.1";',
    'chr17\ttest\texon\t101\t150\t.\t+\t.\tgene_id "GENE2"; transcript_id "NO";'), gtf)
gene <- read_gene_gtf(gtf, "GENE", "NM_1.1")
stopifnot(length(gene$lines) == 4L, gene$start == 101L, gene$end == 350L)
expect_error(read_gene_gtf(gtf, "GENE", "NM_missing"), "not found")
bed1 <- file.path(tmp, "one.bed")
bed2 <- file.path(tmp, "two.bed")
writeLines(c("chr17\t140\t210\tJ1\t10\t+\t140\t210\t0\t2\t10,10\t0,60",
             "chr17\t140\t310\tJ2\t4\t+\t140\t310\t0\t2\t10,10\t0,160",
             "chr17\t140\t230\tJ5\t2\t+\t140\t230\t0\t2\t10,10\t0,80",
             "chr17\t140\t220\tJ3\t5\t.\t140\t220\t0\t2\t10,10\t0,70",
             "chr17\t140\t240\tJ6\t6\t*\t140\t240\t0\t2\t10,10\t0,90",
             "chr17\t140\t250\tJ7\t2\t?\t140\t250\t0\t2\t10,10\t0,100",
             "chr17\t140\t210\tJ4\t3\t-\t140\t210\t0\t2\t10,10\t0,60"), bed1)
writeLines("chr17\t140\t210\tJ1\t7\t+\t140\t210\t0\t2\t10,10\t0,60", bed2)
one <- read_regtools_junctions(bed1, "minus")
two <- read_regtools_junctions(bed2, "plus")
stopifnot(one$intron_start_0[1] == 150L, one$intron_end_0[1] == 200L)
combined <- annotate_junctions(compare_junctions(list(one, two), c("minus", "plus")), gene, "NM_1.1")
known <- combined[combined$intron_end_0 == 200 & combined$strand == "+", ]
stopifnot(known$minus == 10L, known$plus == 7L, known$reference_status == "selected_transcript")
stopifnot(known$left_boundary_exon == "exon 1 (exact)",
          known$right_boundary_exon == "exon 2 (exact)")
stopifnot(combined$reference_status[combined$intron_end_0 == 300] == "other_gene_transcript")
stopifnot(combined$reference_status[combined$intron_end_0 == 220] == "not_in_gene_annotation")
stopifnot(combined$right_boundary_exon[combined$intron_end_0 == 220] == "exon 2 (nearest; 20 bp)",
          combined$left_boundary_exon[combined$strand == "."] == "exon 1 (exact)")
stopifnot(combined$reference_status[combined$strand == "."] == "strand_unknown")
stopifnot(all(combined$reference_status[combined$strand %in% c("*", "?")] == "strand_unknown"),
          combined$minus[combined$strand == "*"] == 6L)
stopifnot(combined$reference_status[combined$strand == "-"] == "other_strand_or_contig")
minus_gene <- gene
minus_gene$strand <- "-"
minus <- annotate_junctions(combined[, setdiff(names(combined), "reference_status")], minus_gene, "NM_1.1")
stopifnot(minus$reference_status[minus$strand == "-" & minus$intron_end_0 == 200] == "selected_transcript")
reverse_known <- minus[minus$strand == "-" & minus$intron_end_0 == 200, ]
stopifnot(reverse_known$left_boundary_exon == "exon 2 (exact)",
          reverse_known$right_boundary_exon == "exon 1 (exact)")
extra <- data.frame(chrom = c("chr17", "chr17", "chr18"), intron_start_0 = c(200L, 150L, 150L),
                    intron_end_0 = c(250L, 300L, 200L), strand = c("+", "+", "+"))
three_exon_gene <- gene
three_exon_gene$lines <- c(gene$lines,
    'chr17\ttest\texon\t301\t350\t.\t+\t.\tgene_id "G001"; gene_name "GENE"; transcript_id "NM_1.1";')
extra <- annotate_junctions(extra, three_exon_gene, "NM_1.1")
stopifnot(extra$left_boundary_exon[1] == "exon 1 / exon 2 (nearest; 50 bp)",
          extra$right_boundary_exon[1] == "exon 2 / exon 3 (nearest; 50 bp)",
          extra$left_boundary_exon[2] == "exon 1 (exact)",
          extra$right_boundary_exon[2] == "exon 3 (exact)",
          extra$reference_status[2] == "other_gene_transcript",
          extra$left_boundary_exon[3] == "not available: other contig",
          extra$right_boundary_exon[3] == "not available: other contig")
empty <- file.path(tmp, "empty.bed")
file.create(empty)
zero <- read_regtools_junctions(empty, "zero")
empty_combined <- annotate_junctions(compare_junctions(list(zero, zero), c("minus", "plus")), gene, "NM_1.1")
stopifnot(nrow(empty_combined) == 0L,
          all(c("left_boundary_exon", "right_boundary_exon") %in% names(empty_combined)))
proximity <- data.frame(chrom = c("chr17", "chr17", "chr17", "chr18"),
                        intron_start_0 = c(250L, 251L, 149L, 150L),
                        intron_end_0 = c(400L, 401L, 200L, 200L), strand = "+")
proximity <- annotate_junction_proximity(proximity, parse_variant("17-150-AC-A"), "chr17")
stopifnot(identical(proximity$variant_distance_bp[1:3], c(100, 101, 0)),
          identical(proximity$near_variant, c(TRUE, FALSE, TRUE, FALSE)),
          is.na(proximity$variant_distance_bp[4]),
          identical(proximity$junction_id, c("J2", "J3", "J1", "J4")))
stopifnot(nrow(annotate_junction_proximity(empty_combined, parse_variant("17-150-A-G"), "chr17")) == 0L)
expect_error(annotate_junction_proximity(rbind(proximity, proximity[1, ]), parse_variant("17-150-A-G"), "chr17"),
             "Duplicate junction")
idx <- file.path(tmp, "idxstats")
writeLines(c("chr17\t1000\t800000\t2", "chr18\t500\t200000\t0", "*\t0\t0\t10"), idx)
stopifnot(read_mapped_alignments(idx) == 1000000)
normalized <- normalize_junctions(data.frame(a = c(10L, 0L), b = c(20L, 2L), c = c(1L, 0L)),
                                  c("a", "b", "c"), c(1000000, 2000000, 0))
stopifnot(identical(normalized$a_per_million_mapped, c(10, 0)),
          identical(normalized$b_per_million_mapped, c(10, 1)),
          all(is.na(normalized$c_per_million_mapped)))
writeLines("chr17\t1000\tbroken\t0", idx)
expect_error(read_mapped_alignments(idx), "Invalid samtools idxstats")
writeLines("chr17\t1\t2\tbad", empty)
expect_error(read_regtools_junctions(empty, "bad"), "Invalid regtools")
writeLines("chr17\t140\t210\tJ1\t10\t!\t140\t210\t0\t2\t10,10\t0,60", empty)
expect_error(read_regtools_junctions(empty, "bad"), "line 1; strand=!")
vcf <- file.path(tmp, "test.vcf")
writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
             "chr17\t100\t.\tA\tG\t.\t.\t."), vcf)
stopifnot(validate_annotated_vcf(vcf, parse_variant("17-100-A-G"), "chr17")[1] == "chr17")
expect_error(validate_annotated_vcf(vcf, parse_variant("17-101-A-G"), "chr17"), "differs from input")
stopifnot(has_info_annotation("X=1;SpliceAI=G|GENE|0.1;Y=2", "SpliceAI"),
          !has_info_annotation("SpliceAI=.", "SpliceAI"),
          !has_info_annotation("SpliceAI=", "SpliceAI"),
          !has_info_annotation("OtherSpliceAI=1", "SpliceAI"))
cat("ClinV1 core tests passed\n")
