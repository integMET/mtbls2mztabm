## Retrieval of ISA-TAB files. By default files are read from the study's
## MetaboLights FTP folder through MsBackendMetaboLights; with `path` set to
## a local directory holding the ISA-TAB files the same files are read from
## disk (used for the bundled MTBLS103 and MTBLS40 example data and for
## offline work).

.check_path <- function(path) {
    if (!is.null(path)) {
        if (!is.character(path) || length(path) != 1L || !dir.exists(path))
            stop("'path' must be an existing directory with ISA-TAB files")
    }
    invisible(path)
}

.mtbls_files <- function(mtblsId, path = NULL, pattern = NULL) {
    if (is.null(path)) {
        retry(mtbls_list_files(mtblsId, pattern = pattern),
              ntimes = 5L, sleep_mult = 4)
    } else {
        list.files(path, pattern = pattern)
    }
}

.mtbls_file_path <- function(mtblsId, file, path = NULL) {
    if (is.null(path))
        paste0(mtbls_ftp_path(mtblsId, mustWork = FALSE), file)
    else file.path(path, file)
}

.mtbls_read_table <- function(mtblsId, file, path = NULL) {
    f <- .mtbls_file_path(mtblsId, file, path)
    x <- retry(read.table(f, header = TRUE, sep = "\t", quote = "",
                          comment.char = "", check.names = FALSE,
                          stringsAsFactors = FALSE, encoding = "UTF-8"),
               ntimes = 5L, sleep_mult = 4)
    .clean_table(x)
}

.mtbls_read_lines <- function(mtblsId, file, path = NULL) {
    f <- .mtbls_file_path(mtblsId, file, path)
    retry(readLines(f, warn = FALSE, encoding = "UTF-8"),
          ntimes = 5L, sleep_mult = 4)
}

## The sample table (s_*.txt).
.mtbls_sample_table <- function(mtblsId, path = NULL) {
    if (is.null(path)) {
        x <- retry(mtbls_sample_data(mtblsId), ntimes = 5L, sleep_mult = 4)
        return(.clean_table(x))
    }
    sfile <- .mtbls_files(mtblsId, path, pattern = "^s_.*\\.txt$")
    if (length(sfile) != 1L)
        stop("Expected exactly one sample file (s_*.txt) in '", path, "'")
    .mtbls_read_table(mtblsId, sfile, path)
}

## Assay table joined to the sample table (all columns, duplicated names
## such as "Term Source REF" kept in place). This mirrors
## `MsBackendMetaboLights::mtbls_metadata()`, which is used for remote data.
.mtbls_isa_table <- function(mtblsId, assayName, path = NULL) {
    if (is.null(path)) {
        x <- retry(mtbls_metadata(mtblsId, assayName = assayName,
                                  keepOntology = TRUE, keepProtocol = TRUE,
                                  simplify = FALSE),
                   ntimes = 5L, sleep_mult = 4)
        return(.clean_table(x))
    }
    a <- .mtbls_read_table(mtblsId, assayName, path)
    s <- .mtbls_sample_table(mtblsId, path)
    o <- match(a[["Sample Name"]], s[["Sample Name"]])
    cbind(a, s[o, , drop = FALSE])
}

## Pick the assay file to work with. An explicit `assayName` is validated;
## otherwise a single assay is used as is and, among several, the first mass
## spectrometry assay (by the investigation's technology type) is chosen.
.mtbls_select_assay <- function(mtblsId, assayName = character(),
                                path = NULL) {
    afiles <- .mtbls_files(mtblsId, path, pattern = "^a_.*\\.txt$")
    if (!length(afiles))
        stop("No assay file (a_*.txt) found for ", mtblsId)
    if (length(assayName)) {
        assayName <- as.character(assayName)[1L]
        if (!assayName %in% afiles)
            stop("Assay \"", assayName, "\" not found. Available assay ",
                 "files:\n  ", paste(afiles, collapse = "\n  "))
        return(assayName)
    }
    if (length(afiles) == 1L) return(afiles)
    ms <- character()
    inv <- tryCatch(mtbls_investigation(mtblsId, path),
                    error = function(e) NULL)
    if (!is.null(inv)) {
        fn <- .inv_field(inv, "Study Assay File Name")
        tt <- .inv_field(inv, "Study Assay Technology Type")
        if (length(fn) && length(fn) == length(tt))
            ms <- fn[grepl("mass spectrometry", tt, ignore.case = TRUE)]
    }
    ms <- intersect(ms, afiles)
    if (!length(ms))
        ms <- grep("mass.?spec|LC-?MS|GC-?MS|MS", afiles, value = TRUE)
    sel <- if (length(ms)) ms[1L] else afiles[1L]
    message("Study ", mtblsId, " has ", length(afiles), " assays; using \"",
            sel, "\".\nUse 'assayName' to select another one:\n  ",
            paste(afiles, collapse = "\n  "))
    sel
}

## The assay files in the order to export them: those of the investigation
## first, without the assays of another technology such as NMR (all assays if
## none is left, so that the export reports what is wrong with them).
.mtbls_export_assays <- function(mtblsId, afiles, path = NULL) {
    inv <- tryCatch(mtbls_investigation(mtblsId, path),
                    error = function(e) NULL)
    fn <- .inv_field(inv, "Study Assay File Name")
    tt <- .inv_field(inv, "Study Assay Technology Type")
    if (!length(fn) || length(fn) != length(tt))
        return(afiles)
    other <- fn[!is.na(tt) & nzchar(tt) &
                !grepl("mass spectrometry", tt, ignore.case = TRUE)]
    ms <- setdiff(afiles, other)
    if (!length(ms))
        return(afiles)
    c(intersect(fn, ms), setdiff(ms, fn))
}

## The metabolite assignment files (MAFs) of a study to export: the assay of
## each MAF, named by the MAF. A MAF shared by several assays gets the first of
## them. A MAF without assay can not be exported, as there are no runs to
## describe in the metadata section.
.mtbls_maf_assays <- function(mtblsId, path = NULL) {
    .check_path(path)
    afiles <- .mtbls_files(mtblsId, path, pattern = "^a_.*\\.txt$")
    mafs <- lapply(setNames(afiles, afiles), function(a) {
        m <- .col_or_na(.mtbls_read_table(mtblsId, a, path),
                        "Metabolite Assignment File")
        unique(m[!is.na(m) & nzchar(m)])
    })
    noassay <- setdiff(.mtbls_files(mtblsId, path, pattern = "^m_.*\\.tsv$"),
                       unlist(mafs))
    if (length(noassay))
        warning("No mzTab-M file for the metabolite assignment file(s) ",
                paste0("\"", noassay, "\"", collapse = ", "), " of ", mtblsId,
                ": no assay refers to them.", call. = FALSE)
    use <- .mtbls_export_assays(mtblsId, afiles, path)
    skip <- setdiff(names(mafs)[lengths(mafs) > 0L], use)
    if (length(skip))
        message("Not exporting the assay(s) of another technology: ",
                paste(skip, collapse = ", "))
    nomaf <- use[!lengths(mafs[use])]
    if (length(nomaf))
        message("Not exporting the assay(s) without metabolite assignment ",
                "file: ", paste(nomaf, collapse = ", "))
    mafs <- mafs[use]
    assay <- rep(names(mafs), lengths(mafs))
    maf <- unlist(mafs, use.names = FALSE)
    first <- !duplicated(maf)
    for (m in unique(maf[!first]))
        message("Metabolite assignment file \"", m, "\" is shared by the ",
                "assays ", paste(assay[maf == m], collapse = ", "),
                ": its mzTab-M file describes the runs of the first. Use ",
                "'assayName' to select another one.")
    setNames(assay[first], maf[first])
}

## The part of the name of a MAF that tells it from the others of its study.
.maf_stems <- function(maf, mtblsId) {
    s <- sub("(_v[0-9]+)?_maf\\.tsv$", "", sub("^m_", "", maf))
    s <- sub(paste0("^", mtblsId, "_+"), "", s)
    make.unique(gsub("[^[:alnum:]._-]+", "_", s), sep = "_")
}

## Everything needed to describe one assay: its name, the joined assay and
## sample table, and the metabolite assignment file the assay refers to
## (name and content).
.mtbls_fetch <- function(mtblsId, assayName = character(), path = NULL) {
    .check_path(path)
    assayName <- .mtbls_select_assay(mtblsId, assayName, path)
    isa <- .mtbls_isa_table(mtblsId, assayName, path)
    if (!nrow(isa))
        stop("Assay \"", assayName, "\" has no samples")
    mafName <- unique(isa[["Metabolite Assignment File"]])
    mafName <- mafName[!is.na(mafName)]
    if (length(mafName) != 1L)
        stop("Assay \"", assayName, "\" must reference exactly one ",
             "metabolite assignment file (MAF), but references ",
             length(mafName), ". The MAF provides the small molecules of ",
             "the mzTab-M file.")
    maf <- mtbls_maf(mtblsId, mafName, path)
    list(assayName = assayName, isa = isa, mafName = mafName, maf = maf)
}

## The spectral data file of each row of an assay table: the derived file, or
## the raw file where no derived file is deposited (GC-MS studies often list
## only their raw netCDF files).
.data_file <- function(isa) {
    derived <- .col_or_na(isa, "Derived Spectral Data File")
    raw <- .col_or_na(isa, "Raw Spectral Data File")
    use_raw <- is.na(derived) & !is.na(raw)
    if (any(use_raw))
        message("No derived spectral data file for ", sum(use_raw),
                " sample(s); using the raw spectral data file instead.")
    ifelse(use_raw, raw, derived)
}
