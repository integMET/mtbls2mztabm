#' Study factors of a MetaboLights study
#'
#' ISA-TAB keeps the experimental factors of a study in two places: their
#' *definitions* in the investigation file (`STUDY FACTORS` block) and their
#' per-sample *values* in the sample file, as one `"Factor Value[<name>]"`
#' column per factor, each optionally followed by `"Unit"`,
#' `"Term Source REF"` and `"Term Accession Number"` qualifier columns.
#'
#' - `mtbls_factors()` returns the factor definitions.
#' - `mtbls_factor_values()` returns the factor values per sample or, when
#'   `assayName` is given, per data file of that assay (using the assay
#'   table to link samples to data files).
#'
#' @inheritParams mtbls_investigation
#'
#' @param assayName optional `character(1)` with the name of an assay file
#'     (`"a_*.txt"`) of the study. If given, factor values are reported per
#'     row of that assay table (i.e. per MS data file) and the result gains a
#'     `data_file` column.
#'
#' @return
#' For `mtbls_factors()`: a `data.frame` with one row per factor and the
#' columns `name`, `type`, `type_term_accession` and `type_term_source`.
#'
#' For `mtbls_factor_values()`: a `data.frame` with `sample_name`
#' (and `data_file` for an assay) followed by one column per factor, named
#' after the factor in lower case with non-alphanumeric characters replaced
#' by `_`. Where the ISA-TAB file provides them, `<factor>_unit`,
#' `<factor>_unit_term_source`, `<factor>_unit_term_accession`,
#' `<factor>_term_source` and `<factor>_term_accession` columns follow.
#' Empty cells are reported as `NA`.
#'
#' @author Kozo Nishida
#'
#' @export
#'
#' @examples
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#'
#' ## Definitions: MTBLS103 has a single factor, "HIAE"
#' mtbls_factors("MTBLS103", path = ext)
#'
#' ## Values per sample
#' fv <- mtbls_factor_values("MTBLS103", path = ext)
#' head(fv)
#' table(fv$hiae)
#'
#' ## Values per data file of the HILIC assay
#' fa <- mtbls_factor_values(
#'     "MTBLS103", path = ext,
#'     assayName = "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt")
#' head(fa)
#'
#' ## From the MetaboLights FTP server (needs internet access):
#' ## mtbls_factors("MTBLS103")
mtbls_factors <- function(mtblsId, path = NULL) {
    inv <- mtbls_investigation(mtblsId, path)
    nm <- .inv_field(inv, "Study Factor Name")
    nm <- nm[nzchar(nm)]
    pad <- function(x) {
        x <- as.character(x)
        length(x) <- length(nm)
        x[is.na(x) | !nzchar(x)] <- NA_character_
        x
    }
    data.frame(
        name = nm,
        type = pad(.inv_field(inv, "Study Factor Type")),
        type_term_accession = pad(
            .inv_field(inv, "Study Factor Type Term Accession Number")),
        type_term_source = pad(
            .inv_field(inv, "Study Factor Type Term Source REF")),
        stringsAsFactors = FALSE)
}

#' @rdname mtbls_factors
#'
#' @export
mtbls_factor_values <- function(mtblsId, assayName = character(),
                                path = NULL) {
    .check_path(path)
    if (length(assayName)) {
        assayName <- .mtbls_select_assay(mtblsId, assayName, path)
        x <- .mtbls_isa_table(mtblsId, assayName, path)
        res <- data.frame(sample_name = x[["Sample Name"]],
                          data_file = .col_or_na(x,
                                                 "Derived Spectral Data File"),
                          stringsAsFactors = FALSE)
    } else {
        x <- .mtbls_sample_table(mtblsId, path)
        res <- data.frame(sample_name = x[["Sample Name"]],
                          stringsAsFactors = FALSE)
    }
    .factor_values(x, res)
}

## Extract "Factor Value[...]" columns (and their qualifiers) from an ISA-TAB
## sample or assay+sample table `x`, appending them to `res`.
.factor_values <- function(x, res = data.frame(
                               sample_name = x[["Sample Name"]],
                               stringsAsFactors = FALSE)) {
    fcols <- grep("^Factor Value\\[", names(x))
    if (!length(fcols))
        warning("No 'Factor Value[...]' column found in the sample data")
    quals <- c("Unit", "Term Source REF", "Term Accession Number")
    for (i in fcols) {
        fname <- sub("^Factor Value\\[(.*)\\]$", "\\1", names(x)[i])
        key <- .key(fname)
        res[[key]] <- x[[i]]
        j <- i + 1L
        unit_seen <- FALSE
        while (j <= ncol(x) && names(x)[j] %in% quals) {
            suffix <- switch(names(x)[j],
                             Unit = "_unit",
                             `Term Source REF` = if (unit_seen)
                                 "_unit_term_source" else "_term_source",
                             `Term Accession Number` = if (unit_seen)
                                 "_unit_term_accession" else "_term_accession")
            if (names(x)[j] == "Unit") unit_seen <- TRUE
            if (!all(is.na(x[[j]])))
                res[[paste0(key, suffix)]] <- x[[j]]
            j <- j + 1L
        }
    }
    res
}
