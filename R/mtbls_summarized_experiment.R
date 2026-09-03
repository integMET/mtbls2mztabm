#' SummarizedExperiment of a MetaboLights assay from its ISA-TAB files
#'
#' Builds a [SummarizedExperiment::SummarizedExperiment()] for one assay of a
#' MetaboLights study from its ISA-TAB files: the abundance matrix of the
#' metabolite assignment file (MAF) becomes the assay data, the sample and
#' assay tables become the `colData()` and the MAF annotation columns the
#' `rowData()`.
#'
#' The layout follows the `SummarizedExperiment` produced by the *xcms*
#' preprocessing workflow for MetaboLights data (as used for MTBLS8735) so
#' that both can be handled by the same downstream code, in particular
#' [se_to_mztabm()]:
#'
#' - `assays`: `"raw"` and `"raw_filled"`, both holding the MAF abundances
#'   (a MAF reports a single, final abundance table, so there is nothing to
#'   gap-fill and the two assays are identical).
#' - `colData()` columns, in this order: `sample_name`,
#'   `derived_spectra_data_file` (the *raw* spectral data file where the
#'   assay table lists no derived one, as in GC-MS studies deposited as
#'   netCDF files), `metabolite_asssignment_file` (spelling kept from the
#'   reference workflow), `source_name`, `organism`,
#'   `blood_sample_type` (`"Characteristics[Organism part]"`),
#'   `sample_type`, `age` (`integer`, from an *Age* factor or
#'   characteristic if present), `unit` (of age), `phenotype` (the study
#'   factor, see `factor`), `injection_index` (the row order of the assay
#'   table), `species` and `tissue` (mzTab-M CV parameters built from the
#'   ontology terms of organism and organism part), `polarity`
#'   (`"Parameter Value[Scan polarity]"`) and `instrument` (the mzTab-M
#'   instrument index, `"1"`). Fields the study does not deposit are `NA`.
#' - `rowData()` columns: the *xcms* feature columns `mzmed` (MAF
#'   `mass_to_charge`), `mzmin`, `mzmax`, `rtmed` (MAF `retention_time`),
#'   `rtmin`, `rtmax`, `npeaks` (number of samples with a reported
#'   abundance), one count column per phenotype level, `ms_level` and
#'   `feature_id`, followed by the MAF annotations `database_identifier`,
#'   `chemical_formula`, `smiles`, `inchi` and `metabolite_identification`
#'   where present.
#' - Row names are *xcms*-style feature IDs (`"FT01"`, ...), column names
#'   the base names of the spectral data files.
#'
#' @inheritParams mtbls_factor_values
#'
#' @param assayName optional `character(1)` with the assay file to use. With
#'     a single assay it is used as is; with several and no `assayName` the
#'     first mass spectrometry assay is chosen and a message names it.
#'
#' @param factor optional `character(1)` with the name of the study factor
#'     that becomes the `phenotype` column (the ISA-TAB column
#'     `"Factor Value[<factor>]"`). By default the first factor that is not
#'     *Age* is used.
#'
#' @return A `SummarizedExperiment` with features in rows and the data files
#'     of the assay in columns.
#'
#' @author Kozo Nishida
#'
#' @seealso [se_to_mztabm()] to export the object to mzTab-M and
#'     [mtbls_to_mztabm()] to do both steps at once.
#'
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' se <- mtbls_summarized_experiment("MTBLS103", path = ext)
#' se
#' colData(se)[1:3, c("sample_name", "phenotype", "species", "tissue")]
#' rowData(se)[, c("feature_id", "npeaks", "CTR", "HIAE", "TREAT",
#'                 "metabolite_identification")]
#' assay(se, "raw")[1:3, 1:3]
#'
#' ## MTBLS40 (GC-TOF-MS, Arabidopsis thaliana) deposits raw netCDF files
#' ## only, so these are used as the spectral data files
#' ext40 <- system.file("extdata", "MTBLS40", package = "mtbls2mztabm")
#' se40 <- mtbls_summarized_experiment("MTBLS40", path = ext40)
#' se40
#' head(se40$derived_spectra_data_file)
#'
#' ## From the MetaboLights FTP server, choosing the assay explicitly
#' ## (needs internet access):
#' ## se <- mtbls_summarized_experiment(
#' ##     "MTBLS103",
#' ##     assayName = "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt")
mtbls_summarized_experiment <- function(mtblsId, assayName = character(),
                                        factor = character(), path = NULL) {
    dat <- .mtbls_fetch(mtblsId, assayName, path)
    .build_se(dat$isa, dat$maf, factor = factor)
}

## Assemble the SummarizedExperiment from the joined assay+sample table and
## the MAF.
.build_se <- function(isa, maf, factor = character()) {
    stopifnot(is.data.frame(isa), is.data.frame(maf))
    raw <- .maf_abundance(maf, isa)
    fid <- .feature_ids(nrow(raw))
    files <- .data_file(isa)
    dimnames(raw) <- list(fid, .col_names(isa, files))
    cd <- .col_data(isa, factor = factor, files = files)
    rownames(cd) <- colnames(raw)
    rd <- .row_data(maf, raw, cd$phenotype, fid)
    SummarizedExperiment(assays = list(raw = raw, raw_filled = raw),
                         rowData = rd, colData = cd)
}

## The spectral data file of each sample: the derived file, or the raw file
## where no derived file is deposited (GC-MS studies often list only their
## raw netCDF files).
.data_file <- function(isa) {
    derived <- .col_or_na(isa, "Derived Spectral Data File")
    raw <- .col_or_na(isa, "Raw Spectral Data File")
    use_raw <- is.na(derived) & !is.na(raw)
    if (any(use_raw))
        message("No derived spectral data file for ", sum(use_raw),
                " sample(s); using the raw spectral data file instead.")
    ifelse(use_raw, raw, derived)
}

.col_names <- function(isa, files = .data_file(isa)) {
    cn <- ifelse(is.na(files), .col_or_na(isa, "Sample Name"),
                 basename(files))
    make.unique(cn)
}

## Age (integer) and its unit from an Age factor or characteristic.
.age_unit <- function(isa) {
    n <- nrow(isa)
    fcols <- grep("^Factor Value\\[", names(isa), value = TRUE)
    fnames <- sub("^Factor Value\\[(.*)\\]$", "\\1", fcols)
    age_col <- c(fcols[tolower(fnames) == "age"],
                 grep("^Characteristics\\[Age\\]$", names(isa),
                      ignore.case = TRUE, value = TRUE))[1L]
    if (is.na(age_col))
        return(list(age = rep(NA_integer_, n), unit = rep(NA_character_, n)))
    list(age = as.integer(round(.as_numeric(isa[[age_col]]))),
         unit = .qualifier(isa, age_col, "Unit"))
}

## The "Factor Value[...]" column that becomes `phenotype`: the requested
## factor, or the first factor that is not Age; NA if there is none.
.phenotype_column <- function(isa, factor = character()) {
    fcols <- grep("^Factor Value\\[", names(isa), value = TRUE)
    fnames <- sub("^Factor Value\\[(.*)\\]$", "\\1", fcols)
    if (length(factor)) {
        pcol <- paste0("Factor Value[", factor[1L], "]")
        if (!pcol %in% names(isa))
            stop("Factor \"", factor[1L], "\" not found. Available ",
                 "factors: ", paste(fnames, collapse = ", "))
        return(pcol)
    }
    cand <- fcols[tolower(fnames) != "age"]
    if (length(cand)) cand[1L] else NA_character_
}

## The 15 colData columns of the reference layout.
.col_data <- function(isa, factor = character(), files = .data_file(isa)) {
    n <- nrow(isa)
    au <- .age_unit(isa)
    pcol <- .phenotype_column(isa, factor)
    phenotype <- if (is.na(pcol)) rep(NA_character_, n)
                 else as.character(isa[[pcol]])
    organism <- .col_or_na(isa, "Characteristics[Organism]")
    org_part <- .col_or_na(isa, "Characteristics[Organism part]")
    term <- function(nm, what) .qualifier(isa, nm, what)
    DataFrame(
        sample_name = .col_or_na(isa, "Sample Name"),
        derived_spectra_data_file = files,
        metabolite_asssignment_file = .col_or_na(isa,
                                                 "Metabolite Assignment File"),
        source_name = .col_or_na(isa, "Source Name"),
        organism = organism,
        blood_sample_type = org_part,
        sample_type = .col_or_na(isa, "Characteristics[Sample type]"),
        age = au$age,
        unit = au$unit,
        phenotype = phenotype,
        injection_index = seq_len(n),
        species = .cv_param(
            term("Characteristics[Organism]", "Term Accession Number"),
            organism,
            term("Characteristics[Organism]", "Term Source REF")),
        tissue = .cv_param(
            term("Characteristics[Organism part]", "Term Accession Number"),
            org_part,
            term("Characteristics[Organism part]", "Term Source REF")),
        polarity = .col_or_na(isa, "Parameter Value[Scan polarity]"),
        instrument = rep("1", n))
}

## The xcms-style rowData columns followed by the MAF annotations.
.row_data <- function(maf, raw, phenotype, fid) {
    n <- nrow(raw)
    named_na <- setNames(rep(NA_real_, n), fid)
    num <- function(col) {
        if (col %in% names(maf)) .as_numeric(maf[[col]])
        else rep(NA_real_, n)
    }
    rd <- DataFrame(mzmed = num("mass_to_charge"),
                    mzmin = named_na, mzmax = named_na,
                    rtmed = num("retention_time"),
                    rtmin = named_na, rtmax = named_na,
                    npeaks = as.numeric(rowSums(!is.na(raw))),
                    row.names = fid)
    groups <- sort(unique(phenotype[!is.na(phenotype)]))
    for (g in groups)
        rd[[g]] <- as.numeric(rowSums(
            !is.na(raw[, which(phenotype == g), drop = FALSE])))
    rd$ms_level <- rep(1L, n)
    rd$feature_id <- fid
    ann <- intersect(c("database_identifier", "chemical_formula", "smiles",
                       "inchi", "metabolite_identification"), names(maf))
    for (cn in ann)
        rd[[cn]] <- as.character(maf[[cn]])
    rd
}
