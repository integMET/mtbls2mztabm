#' Metabolite assignment file of a MetaboLights study
#'
#' The metabolite assignment file (MAF, `m_*_maf.tsv`) of a MetaboLights
#' assay lists the identified/annotated metabolites or features (one per row)
#' with their annotation columns (`database_identifier`, `chemical_formula`,
#' `smiles`, `inchi`, `metabolite_identification`, `mass_to_charge`,
#' `retention_time`, ...) followed by one abundance column per sample.
#'
#' @inheritParams mtbls_investigation
#'
#' @param mafName optional `character(1)` with the name of the MAF to read.
#'     Required if the study has more than one MAF.
#'
#' @return `data.frame` with the content of the MAF; empty cells are `NA`.
#'
#' @author Kozo Nishida
#'
#' @export
#'
#' @examples
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' maf <- mtbls_maf("MTBLS103", path = ext)
#' dim(maf)
#' maf[, c("database_identifier", "metabolite_identification")]
#' ## abundance columns are named by sample
#' head(colnames(maf), 25)
mtbls_maf <- function(mtblsId, mafName = character(), path = NULL) {
    .check_path(path)
    if (!length(mafName)) {
        mafs <- .mtbls_files(mtblsId, path, pattern = "^m_.*\\.tsv$")
        if (length(mafs) != 1L)
            stop("Study ", mtblsId, " has ", length(mafs), " metabolite ",
                 "assignment files; specify 'mafName':\n  ",
                 paste(mafs, collapse = "\n  "))
        mafName <- mafs
    }
    .mtbls_read_table(mtblsId, mafName, path)
}

## Abundance matrix (features x assay rows) from a MAF. Abundance columns are
## matched to the assay rows by sample name, MS assay name or data file name,
## whichever matches most rows; unmatched rows get NA columns.
.maf_abundance <- function(maf, isa) {
    files <- .col_or_na(isa, "Derived Spectral Data File")
    cands <- list(
        `Sample Name` = .col_or_na(isa, "Sample Name"),
        `MS Assay Name` = .col_or_na(isa, "MS Assay Name"),
        `Derived Spectral Data File` = basename(files),
        `Derived Spectral Data File (without extension)` =
            file_path_sans_ext(basename(files)))
    hits <- vapply(cands, function(v)
        sum(!is.na(v) & v %in% colnames(maf)), integer(1))
    if (!any(hits > 0L))
        stop("None of the sample names, MS assay names or data file names ",
             "of the assay matches an abundance column of the MAF")
    best <- names(cands)[which.max(hits)]
    key <- cands[[best]]
    ok <- !is.na(key) & key %in% colnames(maf)
    m <- matrix(NA_real_, nrow = nrow(maf), ncol = length(key))
    vals <- lapply(maf[key[ok]], .as_numeric)
    m[, ok] <- do.call(cbind, vals)
    if (any(!ok))
        message(sum(!ok), " of ", length(key), " assay rows have no ",
                "abundance column in the MAF (matched on '", best,
                "'); their abundances are NA")
    m
}
