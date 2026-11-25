# Læs annotationer og alleltællinger til rapporten.
html_escape <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
}

parse_vep <- function(path, gene, tx) {
    lines <- readLines(path, warn = FALSE)
    header <- lines[startsWith(lines, "##INFO=<ID=CSQ,")]
    if (length(header) == 0L) return(data.frame(status = "VEP CSQ header missing"))
    if (length(header) != 1L || !grepl("Format: ", header, fixed = TRUE)) stop_input("Invalid VEP CSQ header: ", path)
    fieldtext <- sub('".*$', "", sub("^.*Format: ", "", header))
    fields <- strsplit(fieldtext, "|", fixed = TRUE)[[1L]]
    if (!all(c("SYMBOL", "Feature", "Consequence") %in% fields) || anyDuplicated(fields)) stop_input("VEP CSQ field schema incomplete.")
    row <- lines[!startsWith(lines, "#") & nzchar(lines)]
    if (length(row) != 1L) stop_input("Invalid VEP data rows.")
    info <- strsplit(strsplit(row, "\t", fixed = TRUE)[[1L]][8L], ";", fixed = TRUE)[[1L]]
    csq <- info[startsWith(info, "CSQ=")]
    if (length(csq) == 0L) return(data.frame(status = "VEP CSQ not present"))
    if (length(csq) != 1L) stop_input("Multiple CSQ entries in VEP INFO.")
    annotations <- strsplit(sub("^CSQ=", "", csq), ",", fixed = TRUE)[[1L]]
    records <- lapply(annotations, function(a) {
        v <- strsplit(paste0(a, "|END"), "|", fixed = TRUE)[[1L]]
        v <- v[-length(v)]
        if (length(v) > length(fields)) stop_input("VEP CSQ has more fields than header.")
        length(v) <- length(fields)
        setNames(as.list(v), fields)
    })
    records <- Filter(function(x) identical(x$SYMBOL, gene) && identical(x$Feature, tx), records)
    if (!length(records)) return(data.frame(status = "Selected transcript not annotated by VEP"))
    show <- intersect(c("Feature", "Consequence", "IMPACT", "HGVSc", "HGVSp", "EXON", "INTRON", "CLIN_SIG"), fields)
    result <- do.call(rbind, lapply(records, function(x) as.data.frame(x[show], stringsAsFactors = FALSE)))
    result[is.na(result) | result == ""] <- "not provided"
    result
}

parse_spliceai <- function(info, gene) {
    tag <- strsplit(info, ";", fixed = TRUE)[[1L]]
    tag <- tag[startsWith(tag, "SpliceAI=")]
    if (!length(tag) || tag[1L] %in% c("SpliceAI=", "SpliceAI=.")) return(data.frame(status = "SpliceAI annotation unavailable"))
    if (length(tag) != 1L) stop_input("Multiple SpliceAI INFO entries.")
    entries <- strsplit(sub("^SpliceAI=", "", tag), ",", fixed = TRUE)[[1L]]
    vals <- lapply(entries, function(x) strsplit(x, "|", fixed = TRUE)[[1L]])
    if (any(lengths(vals) != 10L)) stop_input("Invalid SpliceAI annotation format.")
    vals <- Filter(function(x) x[2L] == gene, vals)
    if (!length(vals)) return(data.frame(status = "SpliceAI annotation unavailable for selected gene"))
    names <- c("ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")
    result <- as.data.frame(do.call(rbind, vals), stringsAsFactors = FALSE)
    colnames(result) <- names
    result
}

read_ase_counts <- function(path, chrom, pos, ref, alt) {
    data <- read.csv(path, colClasses = "character", comment.char = "#", check.names = FALSE)
    needed <- c("contig", "position", "refAllele", "altAllele", "refCount", "altCount")
    if (!all(needed %in% names(data))) stop_input("ASEReadCounter output missing required columns: ", path)
    if (nrow(data) > 1L) stop_input("ASEReadCounter returned multiple records: ", path)
    if (!nrow(data)) return(list(status = "not_reported_by_ASEReadCounter", ref = NA_integer_, alt = NA_integer_, total = NA_integer_, fraction = NA_real_))
    if (anyNA(data[1L, c("contig", "position", "refAllele", "altAllele")]) ||
        data$contig != chrom || data$position != pos || data$refAllele != ref || data$altAllele != alt) stop_input("ASEReadCounter variant mismatch: ", path)
    counts <- suppressWarnings(as.numeric(c(data$refCount, data$altCount)))
    if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) stop_input("Invalid ASE counts: ", path)
    total <- sum(counts)
    list(status = if (total == 0) "zero_informative_reads" else if (total < 10) "low_depth_exploratory" else "counts_available_exploratory",
         ref = counts[1L], alt = counts[2L], total = total, fraction = if (total) counts[2L] / total else NA_real_)
}

score_band <- function(value) {
    x <- suppressWarnings(as.numeric(as.character(value)))
    if (length(x) != 1L || !is.finite(x) || x < 0 || x > 1) return("unavailable")
    if (x < 0.2) return("low")
    if (x < 0.5) return("recall")
    if (x < 0.8) return("recommended")
    "precision"
}

html_table <- function(data, highlight = NULL, score_columns = character()) {
    if (!nrow(data)) return("<p>Ingen data.</p>")
    cells <- function(row, tag) paste0("<", tag, ">", html_escape(ifelse(is.na(row), "NA", row)), "</", tag, ">", collapse = "")
    row_cells <- function(i) {
        paste(vapply(seq_along(data), function(j) {
            value <- data[[j]][i]
            text <- if (is.na(value)) "NA" else as.character(value)
            if (names(data)[j] %in% score_columns) {
                band <- score_band(value)
                return(paste0('<td class="score-', band, '" title="SpliceAI: ', band, '">',
                              html_escape(text), ' <small>(', band, ')</small></td>'))
            }
            paste0("<td>", html_escape(text), "</td>")
        }, ""), collapse = "")
    }
    paste0("<table><thead><tr>", cells(names(data), "th"), "</tr></thead><tbody>",
           paste(vapply(seq_len(nrow(data)), function(i) paste0(if (!is.null(highlight) && isTRUE(highlight[i])) '<tr class="near-variant">' else "<tr>", row_cells(i), "</tr>"), ""),
                 collapse = ""), "</tbody></table>")
}

select_junctions <- function(junctions, sample_cols, top, near_limit = Inf) {
    if (!nrow(junctions)) return(junctions)
    counts <- rowSums(junctions[sample_cols])
    ranked <- order(-counts, junctions$junction_id)
    near <- which(junctions$near_variant)
    near <- near[order(junctions$variant_distance_bp[near], -counts[near], junctions$junction_id[near])]
    keep <- unique(c(head(ranked, top), head(near, near_limit)))
    junctions[keep, , drop = FALSE]
}

plot_junction_labels <- function(data) {
    compact <- function(x) {
        x <- gsub("exon ", "E", x, fixed = TRUE)
        x <- gsub(" (exact)", "", x, fixed = TRUE)
        x <- sub(" \\(nearest; ([0-9]+) bp\\)", " (~\\1bp)", x)
        x
    }
    paste0(data$junction_id, "  ", compact(data$left_boundary_exon), " - ",
           compact(data$right_boundary_exon), " [", data$strand, "]",
           ifelse(data$near_variant, " *", ""))
}

# Skriv HTML med beskrivende resultater og tydelige begrænsninger.
write_report <- function(path, gene, variant, tx, genotype, vep, spliceai, samples, depth, ase, junctions,
                         sai = list(status = "not_enabled", data = data.frame()),
                         igv_snapshot = list(status = "not_enabled", detail = "")) {
    sample_cols <- intersect(samples$label, names(junctions))
    sorted <- select_junctions(junctions, sample_cols, 30L)
    doc <- c('<!doctype html><html lang="da"><head><meta charset="utf-8"><title>ClinV1 report</title>',
        '<style>body{font:15px sans-serif;max-width:1250px;margin:2em auto;padding:0 1em;color:#20232a}table{border-collapse:collapse;display:block;overflow-x:auto;margin:1em 0}td,th{border:1px solid #bbb;padding:.4em;text-align:left}th{background:#edf1f5}tr.near-variant{background:#fff1d2}img{max-width:100%}.warning{background:#fff1d2;padding:1em}.score-low{background:#dcecff}.score-recall{background:#b9e2e8}.score-recommended{background:#ffdf9b}.score-precision{background:#f7b2b2}.score-unavailable{background:#e5e5e5}</style></head><body>',
        "<h1>ClinV1: observationsrapport</h1>",
        '<p class="warning">Forsknings-/valideringskandidat. Ingen automatisk klinisk klassifikation. nt og puro er ikke biologiske replikater; referencepool er ikke en uafhængig kontrolkohorte.</p>',
        if ("treatment_unknown" %in% samples$condition)
            '<p class="warning">Enkelt patientprøve med ukendt behandling; ingen nt/puro-sammenligning. Kontrollér prøvens metadata.</p>',
        paste0("<p>Gen: ", html_escape(gene), "; variant: ", html_escape(variant), "; transkript: ", html_escape(tx),
               "; angivet zygositet: ", html_escape(genotype), " (ikke bekræftet af RNA).</p>"),
        "<h2>VEP: valgt transkript</h2><p>VEP annoterer varianten i forhold til transkriptet og tilgængelige databaser; konsekvensen er ikke en observeret RNA-effekt.</p>", html_table(vep),
        '<h2>SpliceAI</h2><p>DS_AG/AL og DS_DG/DL beskriver prædikteret acceptor- eller donor-gain/loss (0-1); DP angiver position relativt til varianten. Standardkørslen bruger indbygget genannotation og D=150, ikke nødvendigvis det valgte transkript. Farverne gælder kun DS-værdier, ikke DP. For variantens samlede score bruges normalt maksimum af de fire DS.</p>',
        '<p>Vejledende SpliceAI-grænser: <span class="score-low">&lt;0,2 low</span>; <span class="score-recall">≥0,2 recall</span> (høj følsomhed); <span class="score-recommended">≥0,5 recommended</span> (forfatterens anbefalede balance); <span class="score-precision">≥0,8 precision</span> (høj præcision). Grå unavailable betyder manglende/ugyldig score. Farve er ikke en klinisk klassifikation; lav score udelukker ikke splice-effekt. Tærskler afhænger af kontekst og validering. Kilde: <a href="https://github.com/Illumina/SpliceAI">SpliceAI dokumentation</a>.</p>',
        html_table(spliceai, score_columns = c("DS_AG", "DS_AL", "DS_DG", "DS_DL")),
        "<h2>SAI-10k: prædikteret effekt</h2><p>Broad Lookup-algoritmen bruger lokale, transkriptspecifikke SpliceAI-deltascores. Separate REF/ALT- og alle positionsscores anvendes ikke, og stopkodon/proteinsekvens vurderes ikke. Resultatet kan afvige fra websiden. Hypoteserne er ikke observeret RNA eller en klinisk klassifikation; fravær af prædiktion udelukker ikke effekt.</p>",
        paste0("<p>Status: ", html_escape(sai$status), "</p>"), html_table(sai$data),
        "<p>SAI-effekten er kategorisk og har ikke en selvstændig 0-1-score. Tabellen nedenfor viser de transkriptspecifikke SpliceAI-input (D=4999, raw) på samme farveskala som ovenfor; grænserne er ikke validerede SAI-10k-effekttærskler. De kan afvige fra standardkørslens scores.</p>",
        if (is.null(sai$scores)) '<p class="score-unavailable">Transkriptspecifikke scorer er ikke tilgængelige.</p>' else
            html_table(sai$scores, score_columns = c("DS_AG", "DS_AL", "DS_DG", "DS_DL")),
        "<h2>Prøver og regionsstatistik</h2><p>Antal reads gælder BAM-udsnittet, ikke normaliseret ekspression. XS/strandfiltrerede coverage-tracks udelader andre reads.</p>", html_table(samples),
        "<h2>Dækning ved variant</h2><p>Samtools depth-standardfiltre; nul udelukker ikke splejsningseffekt, især ikke ved introniske positioner.</p>", html_table(depth),
        "<h2>Alleludtryk (ASE/MAE)</h2><p>Deskriptive alleltællinger for angivet heterozygot SNP. Lav dækning, mappingsbias, PCR-duplikater, splice-mønster og manglende fase kan påvirke fortolkningen. Ingen automatisk MAE- eller klinisk konklusion.</p>", html_table(ase),
        "<h2>Observerede junctions</h2><p>Coverage og junctions vises med sashimipy som i v6.4 (fr-firststrand og XS-tag antages; kontrollér biblioteket). Plotgrænsen er 20 junction-reads, så fravær i figuren er ikke bevis for fravær. Figuren bruger BAM-reads og kan adskille sig fra regtools-tabellen. Søjlefiguren viser top 20 efter samlede reads samt op til 20 ekstra nærliggende junctions, udvalgt efter afstand og derefter antal reads. Tabellen viser top 30 og alle nærliggende junctions; alle junctions findes i junction_counts.tsv. J-ID er fælles for denne kørsels figur, tabel og TSV, men kan ændres mellem kørsler. Det er ikke statistiske test; kontrollér dækning og reads i IGV.</p>",
        "<p>Normaliseret værdi = junction-reads × 1.000.000 / mappede alignments i hele input-BAM ifølge samtools idxstats. Det er ikke PSI eller en test for differential splejsning og korrigerer ikke genekspression, library-præparation, mapping eller kvalitet. idxstats kan tælle sekundære/supplementære alignments; regtools og nævneren kan anvende forskellige filtre. En pooled reference er ikke uafhængige prøver. Nul mappede alignments giver NA, ikke nul; sådanne prøver udelades fra den normaliserede figur. Rå tællinger bevares, og begge figurer viser samme junctionudvalg.</p>",
        "<p>Stjerne og rød figurtekst samt farvet tabelrække markerer højst 100 bp mellem en af junctionens første/sidste intronbaser og nærmeste base i variantens REF-interval (1-baseret). Det er en praktisk visuel grænse, ikke en biologisk eller klinisk tærskel. Ved insertion anvendes VCF-ankerbasen. Nærhed beviser ikke sammenhæng, og fjernere junctions kan være relevante, eksempelvis ved exon skipping eller dybt introniske effekter. Ukendt eller modsat strand markeres også positionsmæssigt. E1/E2 er exoner i valgt transkript; ~ angiver nærmeste grænse med afstand, og rækkefølgen i figuren er genomisk venstre til højre.</p>",
        "<p>Exonkolonnerne sammenligner hver genomisk venstre/højre introngrænse med henholdsvis exon-slut/start i det valgte transkript. Exoner nummereres i transkriptets 5'-til-3'-retning ud fra GTF-koordinater; på minus-strand går numrene derfor fra højre mod venstre. Exact betyder eksakt grænse; nearest viser den nærmeste tilsvarende grænse og absolut afstand i bp, ikke et bekræftet splice-site. Lige nære exoner vises alle. Disse positionsbaserede sammenligninger vises også ved ukendt eller modsat strand og bekræfter ikke junctionens strand eller transkript. Reference-status kræver derimod eksakte introngrænser og korrekt strand; andre transkripter kan have andre exoner.</p>",
        '<p><a href="sashimi.pdf">Åbn sashimi-plot (PDF)</a></p><object data="sashimi.pdf" type="application/pdf" width="100%" height="700"><p>Åbn PDF-linket ovenfor.</p></object>',
        "<h2>IGV-snapshot</h2>",
        if (identical(igv_snapshot$status, "available"))
            paste0('<p>Statisk visning af ', html_escape(igv_snapshot$detail),
                   '. Alignment-tracks er skjult. Kontrollér den interaktive IGV-session; billedet kan være ufuldstændigt og er ikke en automatisk fortolkning.</p>',
                   '<p><a href="igv_snapshot.svg">Åbn SVG</a></p><img src="igv_snapshot.svg" alt="IGV-snapshot">')
        else paste0('<p class="warning">IGV-snapshot: ', html_escape(igv_snapshot$status), ". ",
                    html_escape(igv_snapshot$detail), "</p>"),
        '<h3>Rå junctiontællinger</h3><img src="junction_counts.svg" alt="Rå junctiontællinger">',
        '<h3>Junctions per million mappede alignments</h3><img src="junction_normalized.svg" alt="Dybdenormaliserede junctiontællinger">',
        html_table(sorted, sorted$near_variant),
        '<p>Se også IGV-session, rå VCF, BED-filer og logs. Relative figurfiler skal forblive sammen med rapporten.</p></body></html>')
    writeLines(doc, path, useBytes = TRUE)
}

# Afgræns omkring transkriptets exons og kald v6.4's sashimipy.
run_sashimi <- function(path, bams, labels, gtf, reference, transcript, variant_pos, region_start, region_end, logs) {
    records <- strsplit(reference$lines, "\t", fixed = TRUE)
    exon_rows <- Filter(function(x) {
        x[3L] == "exon" && {
            match <- regmatches(x[9L], regexec('(?:^|;[[:space:]]*)transcript_id "([^"]+)"', x[9L], perl = TRUE))[[1L]]
            length(match) > 1L && identical(match[2L], transcript)
        }
    }, records)
    if (!length(exon_rows) || !reference$strand %in% c("+", "-")) stop_input("Sashimi requires selected-transcript exons and a known strand.")
    exons <- unique(data.frame(start = as.integer(vapply(exon_rows, `[`, "", 4L)),
                               end = as.integer(vapply(exon_rows, `[`, "", 5L))))
    exons <- exons[order(exons$start, exons$end), , drop = FALSE]
    distance <- pmax(exons$start - variant_pos, variant_pos - exons$end, 0)
    nearest <- which.min(distance)
    left <- if (nearest > 2L) exons$start[nearest - 2L] else exons$start[nearest] - 2500L
    right <- if (nearest + 2L <= nrow(exons)) exons$end[nearest + 2L] else exons$end[nearest] + 2500L
    left <- max(region_start, min(left, variant_pos), 1L)
    right <- min(region_end, max(right, variant_pos))
    if (left > variant_pos || right < variant_pos || left >= right) stop_input("Invalid sashimi plot window.")
    density <- file.path(logs, "sashimi_density.tsv")
    colors <- rep(c("#264C97", "#5E95D6", "#00805F"), length.out = length(bams))
    if (length(bams) > 2L) colors[3:length(bams)] <- "#00805F"
    writeLines(paste(bams, "bam", labels, colors, "frf", sep = "\t"), density)
    sites <- sort(unique(c(exons$start, exons$end)))
    sites <- sites[sites >= left & sites <= right]
    args <- c("run", "--no-capture-output", "-n", "KGA_Splice", "sashimipy",
              "-e", paste0(reference$chrom, ":", left, "-", right, ":", reference$strand),
              "--density", density, "-r", gtf, "-o", path, "-d", "300", "-t", "20",
              "--show-junction-num", "--density-by-strand", "--same-y", "--intron-scale", "0.3",
              "--reverse-minus", "--rs", paste0("XS:A:", reference$strand),
              if (length(sites)) c("--sites", paste(sites, collapse = ",")),
              "--stroke", paste0(variant_pos, "-", variant_pos, ":", variant_pos, "-", variant_pos, "@red"), "--show-id")
    run_checked("conda", args, file.path(logs, "sashimi.out"), file.path(logs, "sashimi.err"))
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < 5L ||
        !identical(rawToChar(readBin(path, "raw", n = 5L)), "%PDF-")) stop_input("sashimipy did not create a valid PDF header: ", path)
}

plot_junctions <- function(path, junctions, sample_cols, normalized = FALSE) {
    data <- select_junctions(junctions, sample_cols, 20L, 20L)
    grDevices::svg(path, width = 14, height = max(7, 3 + nrow(data) * 0.42))
    on.exit(grDevices::dev.off())
    if (!nrow(junctions)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No junctions observed; coverage not assessed")
        return(invisible(NULL))
    }
    data <- data[nrow(data):1L, , drop = FALSE]
    columns <- if (normalized) paste0(sample_cols, "_per_million_mapped") else sample_cols
    if (any(!columns %in% names(data))) stop_input("Missing junction plot columns.")
    available <- vapply(data[columns], function(x) all(is.finite(x)), logical(1))
    if (!any(available)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Normalization unavailable: no valid mapped alignment counts")
        return(invisible(NULL))
    }
    counts <- t(as.matrix(data[columns[available]]))
    labels <- plot_junction_labels(data)
    colors <- grDevices::hcl.colors(length(sample_cols), "Dark 3")[available]
    graphics::par(mar = c(5, 1, 4, 1), xpd = NA)
    positions <- graphics::barplot(counts, beside = TRUE, horiz = TRUE, axes = FALSE,
        col = colors, xlab = if (normalized) "Junction reads per million mapped alignments" else "Junction reads (raw)",
        main = if (normalized) "Depth-normalized junctions" else "Top and variant-near junctions",
        xlim = c(-max(counts, 1) * 1.2, max(counts, 1) * 1.05), names.arg = rep("", nrow(data)))
    graphics::axis(1, at = pretty(c(0, max(counts, 1)))[pretty(c(0, max(counts, 1))) >= 0])
    centers <- if (is.matrix(positions)) colMeans(positions) else positions
    graphics::text(0, centers, labels = labels, pos = 2, cex = 0.72,
                   col = ifelse(data$near_variant, "#B00020", "#20232a"))
    graphics::legend("topright", legend = sample_cols[available], fill = colors, cex = 0.7)
}
