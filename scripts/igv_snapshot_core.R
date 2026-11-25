render_igv_snapshot <- function(igv, session_path, fasta, project, locus,
                                find_command = Sys.which, runner = system2) {
    if (is.null(igv)) return(list(status = "not_enabled", detail = "Aktivér med -igv-snapshot."))
    logs <- file.path(project, "Logs")
    fail <- function(message) {
        message("IGV snapshot failed: ", message)
        list(status = "failed", detail = message)
    }
    if (!file.exists(igv) || dir.exists(igv) || file.access(igv, 1L) != 0L)
        return(fail(paste("IGV launcher missing or not executable:", igv)))
    xvfb <- unname(find_command("xvfb-run"))
    if (length(xvfb) != 1L || is.na(xvfb) || !nzchar(xvfb))
        return(fail("xvfb-run not found on PATH."))
    paths <- c(session_path, fasta, project, igv)
    if (any(grepl("[[:space:]]", paths)))
        return(fail("IGV batch paths must not contain whitespace."))
    tryCatch({
        session <- xml2::read_xml(session_path)
        xml2::xml_set_attr(xml2::xml_root(session), "genome", normalizePath(fasta, mustWork = TRUE))
        local_session <- file.path(logs, "igv_snapshot_session.xml")
        xml2::write_xml(session, local_session)
        work <- file.path(logs, "igv_runtime")
        if (!dir.exists(work) && !dir.create(work)) stop_input("Cannot create IGV runtime directory.")
        batch <- file.path(logs, "igv_snapshot.batch")
        output <- file.path(project, "Report", "igv_snapshot.svg")
        writeLines(c(paste("genome", fasta), paste("load", local_session), paste("goto", locus),
                     "maxPanelHeight 3000", "scrollToTop",
                     paste("snapshotDirectory", file.path(project, "Report")),
                     "snapshot igv_snapshot.svg", "exit"), batch)
        args <- c("-a", igv, "-g", fasta, "--igvDirectory", work, "-b", batch)
        stdout <- file.path(logs, "igv_snapshot.out")
        stderr <- file.path(logs, "igv_snapshot.err")
        status <- runner(xvfb, args = shQuote(args), stdout = stdout, stderr = stderr, timeout = 180L)
        if (!identical(status, 0L)) stop_input("IGV command failed or timed out (", status, "); see Logs/igv_snapshot.err and .out.")
        if (!file.exists(stdout) || !file.exists(stderr)) stop_input("IGV logs missing; snapshot cannot be verified.")
        log_lines <- c(readLines(stdout, warn = FALSE), readLines(stderr, warn = FALSE))
        if (any(grepl("(^|[[:space:]])(SEVERE|ERROR)([[:space:]]|:|\\[|$)", log_lines, ignore.case = TRUE, perl = TRUE)))
            stop_input("IGV reported an error; see Logs/igv_snapshot.err and .out. SVG may be incomplete.")
        if (!file.exists(output) || file.info(output)$size == 0L) stop_input("IGV did not create an SVG; inspect IGV logs.")
        svg <- xml2::read_xml(output)
        if (xml2::xml_name(xml2::xml_root(svg)) != "svg") stop_input("IGV output is not an SVG document.")
        list(status = "available", detail = locus)
    }, error = function(e) fail(conditionMessage(e)))
}
