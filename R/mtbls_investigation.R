#' Read the ISA-TAB investigation file of a MetaboLights study
#'
#' The investigation file (`i_Investigation.txt`) holds the study-level
#' metadata of a MetaboLights data set: title, description, publications,
#' contacts, protocols, factors and the list of assay files. It is a
#' key-value file with one tab-separated record per line. This function
#' returns it as a named `list` in which every entry is the `character`
#' vector of values of one record (multi-valued records such as
#' `"Study Assay File Name"` have one element per assay).
#'
#' @param mtblsId `character(1)` with the MetaboLights study identifier,
#'     e.g. `"MTBLS103"`.
#'
#' @param path optional `character(1)` with a local directory holding the
#'     study's ISA-TAB files. If `NULL` (the default) the files are read from
#'     the MetaboLights FTP server through the *MsBackendMetaboLights*
#'     package.
#'
#' @return Named `list` of `character` vectors, one per record of the
#'     investigation file. Section headers (such as `"STUDY FACTORS"`) are
#'     dropped. Records with the same key (which can happen across sections)
#'     are all kept; `x[["key"]]` returns the first.
#'
#' @author Kozo Nishida
#'
#' @export
#'
#' @examples
#' ## The ISA-TAB files of MTBLS103 shipped with the package
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' inv <- mtbls_investigation("MTBLS103", path = ext)
#' inv[["Study Title"]]
#' inv[["Study Assay File Name"]]
#' inv[["Study Publication DOI"]]
#'
#' ## The same from the MetaboLights FTP server (needs internet access):
#' ## inv <- mtbls_investigation("MTBLS103")
mtbls_investigation <- function(mtblsId, path = NULL) {
    .check_path(path)
    ifile <- .mtbls_files(mtblsId, path, pattern = "^i_.*\\.txt$")
    if (length(ifile) != 1L)
        stop("Expected exactly one investigation file (i_*.txt) for ",
             mtblsId, ", found ", length(ifile))
    lns <- .mtbls_read_lines(mtblsId, ifile, path)
    lns <- lns[nzchar(trimws(lns))]
    parts <- strsplit(lns, "\t", fixed = TRUE)
    keys <- vapply(parts, `[`, character(1), 1L)
    vals <- lapply(parts, function(p) .unquote(p[-1L]))
    ## section headers are all upper case and carry no value
    is_section <- !grepl("[a-z]", keys) & lengths(vals) == 0L
    res <- vals[!is_section]
    names(res) <- keys[!is_section]
    res
}

## Values of one investigation record, `character()` if absent.
.inv_field <- function(inv, key) {
    v <- inv[[key]]
    if (is.null(v)) character() else v
}
