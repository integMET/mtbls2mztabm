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
#' @seealso [mtbls_maf_to_parquet()] to save the MAF as a Parquet file.
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

#' Save the metabolite assignment file of a MetaboLights assay as Parquet
#'
#' Reads the metabolite assignment file (MAF, `m_*_maf.tsv`) of an assay of a
#' MetaboLights study and saves it as an
#' [Apache Parquet](https://parquet.apache.org) file with the
#' *nanoparquet* package. The Parquet file is the second of the two
#' intermediate files of the export workflow (the first one, written by
#' [mtbls_factors_to_toml()], describes the study); [mtbls_parquet_to_sml()]
#' turns it into the mzTab-M `SML` section.
#'
#' The content of the MAF is kept as is, with the column names unchanged, but
#' the column types are made stable: the annotation columns of the MAF
#' (`database_identifier`, `smiles`, `reliability`, ...) are text and
#' `mass_to_charge`, `retention_time` and the other numeric standard columns
#' are numbers also when a study leaves them empty, so that the file has the
#' same schema for every study. Empty cells are `NA` (Parquet null). The
#' study identifier and the name of the MAF are stored as key-value metadata
#' of the file (`mtbls_id` and `maf_file`).
#'
#' @inheritParams mtbls_factors_to_toml
#'
#' @param file `character(1)` with the path of the Parquet file to write.
#'     Default: `"<mtblsId>.parquet"` in the working directory.
#'
#' @param assayName optional `character(1)` with the assay file (`"a_*.txt"`)
#'     whose MAF to save. With a single assay it is used as is; with several
#'     and no `assayName` the first mass spectrometry assay is chosen and a
#'     message names it. Ignored if `mafName` is given.
#'
#' @param mafName optional `character(1)` with the name of the MAF to save,
#'     to select it directly instead of through its assay.
#'
#' @return Invisibly, the path of the Parquet file.
#'
#' @author Kozo Nishida
#'
#' @seealso [mtbls_maf()] to read the MAF into R,
#'     [mtbls_parquet_to_sml()] to read the file and [mtbls_to_mztabm()] for
#'     the complete workflow.
#'
#' @export
#'
#' @examples
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' f <- mtbls_maf_to_parquet(
#'     "MTBLS103", file = tempfile(fileext = ".parquet"), path = ext,
#'     assayName = "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt")
#' maf <- nanoparquet::read_parquet(f)
#' dim(maf)
#' maf[, c("database_identifier", "metabolite_identification")]
#'
#' ## From the MetaboLights FTP server (needs internet access):
#' ## mtbls_maf_to_parquet("MTBLS103")
mtbls_maf_to_parquet <- function(mtblsId, file = paste0(mtblsId, ".parquet"),
                                 assayName = character(),
                                 mafName = character(), path = NULL) {
    stopifnot(is.character(mtblsId), length(mtblsId) == 1L, nzchar(mtblsId),
              is.character(file), length(file) == 1L)
    if (length(mafName)) {
        mafName <- as.character(mafName)[1L]
        maf <- mtbls_maf(mtblsId, mafName, path)
    } else {
        dat <- .mtbls_fetch(mtblsId, assayName, path)
        mafName <- dat$mafName
        maf <- dat$maf
    }
    .write_maf_parquet(maf, file, mtblsId, mafName)
}

## Standard annotation columns of a MAF that hold text and numbers.
.maf_text_columns <- c("database_identifier", "chemical_formula", "smiles",
                       "inchi", "metabolite_identification", "fragmentation",
                       "modifications", "taxid", "species", "database",
                       "database_version", "reliability", "uri",
                       "search_engine")
.maf_numeric_columns <- c("mass_to_charge", "charge", "retention_time",
                          "search_engine_score",
                          "smallmolecule_abundance_sub",
                          "smallmolecule_abundance_stdev_sub",
                          "smallmolecule_abundance_std_error_sub")

## Give the columns of a MAF the types they have in every study: standard
## text columns are `character` and numbers are `double` (also if a column
## holds integers only, or if all values are missing, which `read.table()`
## reads as `logical`). Standard numeric columns with non-numeric entries and
## all other columns keep their values.
.maf_for_parquet <- function(maf) {
    nm <- names(maf)
    for (i in seq_along(maf)) {
        v <- maf[[i]]
        if (nm[i] %in% .maf_text_columns) {
            v <- as.character(v)
        } else if (is.logical(v)) {
            v <- if (all(is.na(v))) as.numeric(v) else as.character(v)
        } else if (is.integer(v)) {
            v <- as.numeric(v)
        } else if (nm[i] %in% .maf_numeric_columns && !is.numeric(v)) {
            n <- .as_numeric(v)
            v <- if (identical(is.na(n), is.na(v))) n else as.character(v)
        }
        maf[[i]] <- v
    }
    maf
}

## Write the MAF `maf` of study `mtblsId` to the Parquet file `file`.
.write_maf_parquet <- function(maf, file, mtblsId, mafName) {
    .ensure_dir(file)
    ver <- tryCatch(as.character(packageVersion("mtbls2mztabm")),
                    error = function(e) "unknown")
    write_parquet(.maf_for_parquet(maf), file,
                  metadata = c(mtbls_id = mtblsId, maf_file = mafName,
                               generator = paste("mtbls2mztabm", ver)))
    invisible(file)
}

## Read a Parquet file written by `mtbls_maf_to_parquet()`.
.read_maf_parquet <- function(file) {
    if (!is.character(file) || length(file) != 1L || !file.exists(file))
        stop("Parquet file '", paste(file, collapse = ", "),
             "' does not exist")
    as.data.frame(read_parquet(file), stringsAsFactors = FALSE)
}

## Columns of a MAF that are annotations rather than sample abundances.
.maf_annotation_columns <- c(.maf_text_columns, .maf_numeric_columns)

## Name of the MAF column holding the abundances of each row of the assay
## table `isa` (NA if the MAF has none). The columns are matched on the
## sample name, the MS assay name or the data file name (with or without
## extension), whichever matches most rows. If none of the columns found holds
## a number (`filled`, the columns with numbers, see `.maf_filled()`), the
## columns of `.maf_fallback()` are used if there is a one-to-one alternative.
.maf_columns <- function(isa, maf_names, filled = maf_names) {
    maf_names <- maf_names[!maf_names %in% .maf_annotation_columns]
    repeated <- unique(maf_names[duplicated(maf_names)])
    if (length(repeated))
        message("The MAF repeats ", length(repeated), " abundance column ",
                "names (", .some(repeated), "); only the first column of ",
                "each name can be used")
    derived <- .col_or_na(isa, "Derived Spectral Data File")
    raw <- .col_or_na(isa, "Raw Spectral Data File")
    cands <- list(
        `Sample Name` = .col_or_na(isa, "Sample Name"),
        `MS Assay Name` = .col_or_na(isa, "MS Assay Name"),
        `Derived Spectral Data File` = basename(derived),
        `Derived Spectral Data File (without extension)` =
            file_path_sans_ext(basename(derived)),
        `Raw Spectral Data File` = basename(raw),
        `Raw Spectral Data File (without extension)` =
            file_path_sans_ext(basename(raw)))
    hits <- vapply(cands, function(v)
        sum(!is.na(v) & v %in% maf_names), integer(1))
    best <- names(cands)[which.max(hits)]
    key <- cands[[best]]
    res <- ifelse(!is.na(key) & key %in% maf_names, key, NA_character_)
    if (!any(res %in% filled)) {
        alt <- .maf_fallback(cands, setdiff(intersect(filled, maf_names),
                                            repeated))
        if (!is.null(alt)) {
            message("No abundance column matched on the names of the assay ",
                    "holds a number; using the ", length(alt), " columns ",
                    "with numbers that correspond one-to-one to the runs ",
                    "(e.g. '", alt[1L], "' for the run '", key[1L], "')")
            res <- alt
        } else if (any(hits > 0L)) {
            warning("None of the ", sum(!is.na(res)), " abundance columns ",
                    "of the MAF matched on '", best, "' holds a number and ",
                    "no other columns correspond one-to-one to the runs",
                    call. = FALSE)
        }
    }
    if (all(is.na(res))) {
        warning("None of the sample names, MS assay names or data file ",
                "names of the assay matches an abundance column of the MAF",
                call. = FALSE)
        return(res)
    }
    if (anyNA(res))
        message(sum(is.na(res)), " of ", length(res), " assay rows have no ",
                "abundance column in the MAF (matched on '", best,
                "'); their abundances are NA")
    used <- res[!is.na(res)]
    shared <- unique(used[duplicated(used)])
    if (length(shared))
        warning(sum(used %in% shared), " assay rows share ", length(shared),
                " abundance columns of the MAF (", .some(shared), "): the ",
                "values are repeated for each run, which gives these columns ",
                "more weight in the abundance of their study variables and ",
                "understates the variation", call. = FALSE)
    res
}

## The first values of `x` for a message.
.some <- function(x, n = 3L) {
    paste0(paste(x[seq_len(min(n, length(x)))], collapse = ", "),
           if (length(x) > n) ", ...")
}

## Columns of the MAF that correspond one-to-one to the runs, for an assay
## whose matched columns hold no number. `cands` are the names of the runs and
## `pool` the columns with numbers whose name is not repeated in the MAF. A
## column corresponds to a run if the names are equal or else if one name is
## the other one plus a suffix after a separator ("DFCI-01583-GC-01" and
## "DFCI-01583"). All runs need a column of their own and all names of the runs
## must lead to the same columns; otherwise `NULL`, nothing is guessed.
.maf_fallback <- function(cands, pool) {
    if (!length(pool))
        return(NULL)
    unique_map <- function(maps) {
        maps <- unique(Filter(Negate(is.null), maps))
        if (length(maps) == 1L) maps[[1L]]
    }
    one_to_one <- function(key, related) {
        cols <- lapply(key, function(k)
            if (is.na(k)) character() else pool[related(k, pool)])
        if (all(lengths(cols) == 1L) && !anyDuplicated(unlist(cols)))
            unlist(cols)
    }
    equal <- unique_map(lapply(cands, one_to_one, `==`))
    if (!is.null(equal))
        return(equal)
    unique_map(lapply(cands, one_to_one, function(k, pool)
        .has_suffix(k, pool) | .has_suffix(pool, k)))
}

## Is `x` the name `prefix` followed by a separator and a suffix?
.has_suffix <- function(x, prefix) {
    startsWith(x, prefix) & nchar(x) > nchar(prefix) &
        substring(x, nchar(prefix) + 1L, nchar(prefix) + 1L) %in%
        c("-", "_", ".", " ")
}

## A number with thousands separators.
.maf_grouped_number <- "^[-+]?[1-9][0-9]{0,2}(,[0-9]{3})+(\\.[0-9]*)?$"

## The abundance columns of the MAF `maf` as a list of numeric vectors named
## after the columns. Thousands separators ("1,157,497") are removed if all
## cells with a comma of the MAF are written that way and at least one of them
## cannot be a number with a decimal comma (two commas, or a comma and a point).
## Otherwise the cells with a comma are NA and their number is the attribute
## "ambiguous" of the result.
.maf_numbers <- function(maf) {
    cols <- as.list(maf)[!names(maf) %in% .maf_annotation_columns]
    cells <- trimws(unlist(lapply(cols, as.character), use.names = FALSE))
    commas <- cells[!is.na(cells) & grepl(",", cells, fixed = TRUE)]
    grouped <- length(commas) && all(grepl(.maf_grouped_number, commas)) &&
        any(grepl(",.*,|\\.", commas))
    if (grouped)
        cols <- lapply(cols, function(v) gsub(",", "", v, fixed = TRUE))
    res <- lapply(cols, .as_numeric)
    attr(res, "ambiguous") <- if (grouped) 0L else length(commas)
    res
}

## Names of the abundance columns of the MAF `maf` holding at least one number.
## With a repeated name, this is the first column, which is the one used.
.maf_filled <- function(maf) {
    nums <- .maf_numbers(maf)
    nums <- nums[!duplicated(names(nums))]
    names(nums)[vapply(nums, function(v) any(!is.na(v)), logical(1))]
}

## Abundance matrix (MAF rows x runs) from the MAF columns of the runs
## (`columns`, NA for runs without one, whose abundances are NA).
.maf_abundance <- function(maf, columns) {
    nums <- .maf_numbers(maf)
    ok <- !is.na(columns) & columns %in% names(nums)
    if (!any(ok))
        stop("None of the runs of the study description has an abundance ",
             "column in the metabolite assignment file")
    if (attr(nums, "ambiguous") > 0L)
        warning(attr(nums, "ambiguous"), " abundances of the metabolite ",
                "assignment file have a comma that is not clearly a ",
                "thousands separator; they are reported as missing",
                call. = FALSE)
    m <- matrix(NA_real_, nrow = nrow(maf), ncol = length(columns))
    m[, ok] <- do.call(cbind, nums[columns[ok]])
    m
}
