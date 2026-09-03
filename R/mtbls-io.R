## Retrieval of ISA-TAB files. By default files are read from the study's
## MetaboLights FTP folder through MsBackendMetaboLights; with `path` set to
## a local directory holding the ISA-TAB files the same files are read from
## disk (used for the bundled MTBLS103 example data and for offline work).

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

## Everything needed to build the SummarizedExperiment of one assay.
.mtbls_fetch <- function(mtblsId, assayName = character(), path = NULL) {
    .check_path(path)
    assayName <- .mtbls_select_assay(mtblsId, assayName, path)
    isa <- .mtbls_isa_table(mtblsId, assayName, path)
    if (!nrow(isa))
        stop("Assay \"", assayName, "\" has no samples")
    mafName <- unique(isa[["Metabolite Assignment File"]])
    mafName <- mafName[!is.na(mafName)]
    if (length(mafName) != 1L)
        stop("The assay must reference exactly one metabolite assignment ",
             "file; found ", length(mafName))
    maf <- mtbls_maf(mtblsId, mafName, path)
    list(assayName = assayName, isa = isa, maf = maf)
}
