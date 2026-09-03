#' Export a SummarizedExperiment of MetaboLights data to mzTab-M
#'
#' `se_to_mztabm()` maps a `SummarizedExperiment` in the layout produced by
#' [mtbls_summarized_experiment()] onto an mzTab-M 2.1 object of the
#' *RmzTabM* package. The `colData()` columns `sample_name`, `species` and
#' `tissue` fill the mzTab-M *sample* section, `derived_spectra_data_file`
#' and `polarity` the *ms_run* and *assay* sections, `phenotype` becomes the
#' study variable, and the `rowData()` columns `mzmed`, `rtmed`,
#' `feature_id` and the MAF annotations fill the small molecule feature
#' (SMF) table together with the `"raw_filled"` (or first) assay.
#'
#' Study-level metadata is filled in from the ISA-TAB files when they are
#' provided: the instrument, ion source and mass analyzer from the assay
#' table (`isa`), and the title, description, publication and MetaboLights
#' URI from the investigation (`investigation`). The compound database is
#' derived from the identifiers in `rowData()$database_identifier` (ChEBI
#' identifiers are recognised). Contacts are only written when the
#' investigation carries ORCIDs, which mzTab-M requires for a contact.
#'
#' @param se `SummarizedExperiment` as returned by
#'     [mtbls_summarized_experiment()]. Required `colData()` columns:
#'     `sample_name`, `derived_spectra_data_file`, `species`, `tissue`,
#'     `polarity`, `phenotype`; required `rowData()` column: `mzmed`.
#'
#' @param id `character(1)` with the mzTab-M identifier of the file,
#'     usually the MetaboLights study ID.
#'
#' @param isa optional `data.frame` with the joined assay and sample table
#'     of the study (as used internally by [mtbls_summarized_experiment()])
#'     from which instrument metadata is taken.
#'
#' @param investigation optional named `list` as returned by
#'     [mtbls_investigation()] from which title, description, publication
#'     and contacts are taken.
#'
#' @param ftp_path optional `character(1)` with the study's FTP folder; if
#'     given, the ms_run locations are written as full URIs.
#'
#' @return An `RmzTabM::MzTabM` object. Write it with
#'     [RmzTabM::writeMzTabM()].
#'
#' @author Kozo Nishida
#'
#' @seealso [mtbls_to_mztabm()] for the one-call export of a study.
#'
#' @export
#'
#' @examples
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' se <- mtbls_summarized_experiment("MTBLS103", path = ext)
#' inv <- mtbls_investigation("MTBLS103", path = ext)
#' m <- se_to_mztabm(se, id = "MTBLS103", investigation = inv)
#' m
#' head(RmzTabM::mtd(m), 12)
#' ## the sample section built from colData(se)
#' RmzTabM::mtd(m)[grepl("^sample\\[1\\]", RmzTabM::mtd(m)[, 1]), ]
se_to_mztabm <- function(se, id, isa = NULL, investigation = NULL,
                         ftp_path = NULL) {
    stopifnot(is(se, "SummarizedExperiment"),
              is.character(id), length(id) == 1L, nzchar(id))
    req <- c("sample_name", "derived_spectra_data_file", "species",
             "tissue", "polarity", "phenotype")
    miss <- setdiff(req, colnames(colData(se)))
    if (length(miss))
        stop("colData(se) lacks the column(s): ", paste(miss, collapse = ", "))
    if (!"mzmed" %in% colnames(rowData(se)))
        stop("rowData(se) lacks the mandatory column 'mzmed'")
    if (!length(assayNames(se)))
        stop("'se' has no assay")
    rowData(se)$feature_id <- rownames(se)
    prep <- .export_columns(se, ftp_path)
    se <- prep$se
    assay_name <- if ("raw_filled" %in% assayNames(se)) "raw_filled"
                  else assayNames(se)[1L]
    m <- MzTabM(se, id = id,
                sampleCols. = sampleCols(sample = "sample_name",
                                         species = "species",
                                         tissue = "tissue"),
                msRunCols. = do.call(msRunCols, prep$ms_args),
                assayCols. = assayCols(assay = "derived_spectra_data_file"),
                groups = "phenotype", assayName = assay_name,
                smfCols. = do.call(smfCols, .smf_args(se)))
    ver <- tryCatch(as.character(packageVersion("mtbls2mztabm")),
                    error = function(e) "unknown")
    m <- setMtdField(m, "software",
                     sprintf("[,, mtbls2mztabm, %s]", ver))
    if (!is.null(isa))
        m <- .set_instrument(m, isa)
    if ("database_identifier" %in% colnames(rowData(se)))
        m <- .set_database(m, rowData(se)$database_identifier)
    if (!is.null(investigation))
        m <- .set_study_info(m, investigation, id)
    m
}

## Export-only colData columns (full file URIs, CV parameters accepted by
## RmzTabM) and the matching msRunCols() arguments.
.export_columns <- function(se, ftp_path = NULL) {
    cd <- colData(se)
    loc <- cd$derived_spectra_data_file
    if (length(ftp_path) && !is.na(ftp_path[1L]))
        loc <- ifelse(is.na(loc), NA_character_, paste0(ftp_path[1L], loc))
    colData(se)$ms_run_location <- loc
    colData(se)$ms_run_format <- .file_format_cv(cd$derived_spectra_data_file)
    colData(se)$ms_run_id_format <-
        .id_format_cv(cd$derived_spectra_data_file)
    colData(se)$scan_polarity_cv <- .polarity_plain(cd$polarity)
    ms_args <- list(location = "ms_run_location",
                    scan_polarity = "scan_polarity_cv")
    ## mzTab-M wants format and id_format together, or neither
    if (!anyNA(colData(se)$ms_run_format) &&
        !anyNA(colData(se)$ms_run_id_format)) {
        ms_args$format <- "ms_run_format"
        ms_args$id_format <- "ms_run_id_format"
    }
    list(se = se, ms_args = ms_args)
}

## smfCols() arguments: the mandatory m/z, then the annotation present.
.smf_args <- function(se) {
    args <- list(exp_mass_to_charge = "mzmed", feature_id = "feature_id")
    if ("rtmed" %in% colnames(rowData(se)))
        args$retention_time_in_seconds <- "rtmed"
    extra <- intersect(c("metabolite_identification", "database_identifier",
                         "chemical_formula"), colnames(rowData(se)))
    c(args, setNames(as.list(extra), extra))
}

## Instrument section from the ISA-TAB assay table.
.set_instrument <- function(m, isa) {
    val <- function(nm) {
        v <- .col_or_na(isa, nm)
        v[!is.na(v)][1L]
    }
    acc <- function(nm) {
        v <- .qualifier(isa, nm, "Term Accession Number")
        v[!is.na(v)][1L]
    }
    inst <- val("Parameter Value[Instrument]")
    src <- val("Parameter Value[Ion source]")
    ana <- val("Parameter Value[Mass analyzer]")
    det <- val("Parameter Value[Detector]")
    if (all(is.na(c(inst, src, ana, det))))
        return(m)
    setMtdInstrument(
        m,
        name = .ms_cv(acc("Parameter Value[Instrument]"), inst,
                      "MS:1000031", "instrument model"),
        source = .ms_cv(acc("Parameter Value[Ion source]"), src,
                        "MS:1000008", "ionization type"),
        analyzer = c(`analyzer[1]` = .ms_cv(
            acc("Parameter Value[Mass analyzer]"), ana,
            "MS:1000443", "mass analyzer type")),
        detector = .ms_cv(acc("Parameter Value[Detector]"), det,
                          "MS:1000026", "detector type"))
}

## Compound database from the identifier prefixes of the MAF.
.set_database <- function(m, ids) {
    ids <- ids[!is.na(ids)]
    prefix <- unique(toupper(sub(":.*$", "", ids[grepl(":", ids)])))
    if (length(prefix) != 1L)
        return(m)
    if (prefix == "CHEBI") {
        setMtdDatabase(m, name = "[MIRIAM, MIR:00000002, ChEBI, ]",
                       prefix = "CHEBI", version = "unknown",
                       uri = "https://www.ebi.ac.uk/chebi/", replace = TRUE)
    } else {
        setMtdDatabase(m, name = sprintf("[,, %s, ]", prefix),
                       prefix = prefix, version = "unknown", uri = "null",
                       replace = TRUE)
    }
}

## Title, description, publication, URI and (with ORCIDs) contacts.
.set_study_info <- function(m, inv, id) {
    first <- function(key) {
        v <- .inv_field(inv, key)
        v <- v[!is.na(v) & nzchar(v)]
        if (length(v)) v[1L] else NA_character_
    }
    title <- first("Study Title")
    if (!is.na(title))
        m <- setMtdField(m, "title", title)
    descr <- first("Study Description")
    if (!is.na(descr))
        m <- setMtdField(m, "description", descr)
    doi <- .inv_field(inv, "Study Publication DOI")
    pmid <- .inv_field(inv, "Study PubMed ID")
    pub <- character()
    for (i in seq_len(max(length(doi), length(pmid)))) {
        p <- c(if (i <= length(doi) && nzchar(doi[i]))
                   paste0("doi:", sub("^doi:", "", doi[i])),
               if (i <= length(pmid) && nzchar(pmid[i]))
                   paste0("pubmed:", pmid[i]))
        if (length(p)) pub <- c(pub, paste(p, collapse = " | "))
    }
    for (p in pub)
        m <- setMtdField(m, "publication", p)
    m <- setMtdField(m, "uri",
                     paste0("https://www.ebi.ac.uk/metabolights/", id))
    .set_contacts(m, inv)
}

## Contacts need an ORCID in mzTab-M; ISA-TAB investigation files rarely
## provide one, so contacts are only written when they do.
.set_contacts <- function(m, inv) {
    last <- .inv_field(inv, "Study Person Last Name")
    firstn <- .inv_field(inv, "Study Person First Name")
    aff <- .inv_field(inv, "Study Person Affiliation")
    email <- .inv_field(inv, "Study Person Email")
    okey <- grep("ORCID", names(inv), ignore.case = TRUE, value = TRUE)
    orcid <- if (length(okey)) .inv_field(inv, okey[1L]) else character()
    for (i in seq_along(last)) {
        oi <- if (i <= length(orcid)) orcid[i] else ""
        oi <- sub("^https?://orcid.org/", "", oi)
        if (!grepl("^\\d{4}-\\d{4}-\\d{4}-\\d{3}[0-9X]$", oi))
            next
        m <- setMtdContact(
            m, name = trimws(paste(firstn[i], last[i])),
            affiliation = if (i <= length(aff) && nzchar(aff[i])) aff[i]
                          else "null",
            email = if (i <= length(email) && nzchar(email[i])) email[i]
                    else "null",
            orcid = oi)
    }
    m
}

#' Export a MetaboLights study to an mzTab-M file
#'
#' The main function of the package: given a MetaboLights study identifier
#' it retrieves the study's ISA-TAB files (through the
#' *MsBackendMetaboLights* package), builds the `SummarizedExperiment` of
#' one of its mass spectrometry assays with
#' [mtbls_summarized_experiment()], maps it to mzTab-M with
#' [se_to_mztabm()] and writes the file with [RmzTabM::writeMzTabM()].
#'
#' After writing, the file is read back to count its `MTD` and `SMF` rows,
#' and a message reports these together with any metadata still at a
#' placeholder value (no instrument, no compound database, no contact).
#'
#' @inheritParams mtbls_summarized_experiment
#'
#' @param file `character(1)` with the path of the mzTab-M file to write.
#'     Default: `"<mtblsId>.mztab"` in the working directory. An extension
#'     other than `.mztab`/`.mzTab` gets `.mztab` appended by *RmzTabM*.
#'
#' @param verbose `logical(1)` whether to report the written sections.
#'
#' @return Invisibly, the path of the written file. The `SummarizedExperiment`
#'     and the `MzTabM` object are attached as attributes `"se"` and
#'     `"mztabm"`.
#'
#' @author Kozo Nishida
#'
#' @export
#'
#' @examples
#' ## From the ISA-TAB files of MTBLS103 shipped with the package
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' out <- mtbls_to_mztabm("MTBLS103", file = tempfile(fileext = ".mztab"),
#'                        path = ext)
#' readLines(out, n = 8)
#' RmzTabM::readMzTabM(out)
#'
#' ## MTBLS40, a GC-TOF-MS study of Arabidopsis thaliana genotypes whose
#' ## assay table lists raw netCDF files only: these become the ms_runs, the
#' ## factor Genotype the study variables
#' ext40 <- system.file("extdata", "MTBLS40", package = "mtbls2mztabm")
#' out40 <- mtbls_to_mztabm("MTBLS40", file = tempfile(fileext = ".mztab"),
#'                          path = ext40)
#' lns <- readLines(out40)
#' lns[grepl("^MTD\\s(ms_run\\[1\\]|instrument|study_variable\\[[1-4]\\]\\s)",
#'           lns)]
#' se40 <- attr(out40, "se")
#' se40
#' head(SummarizedExperiment::rowData(se40)[, c("feature_id", "rtmed",
#'     "database_identifier", "metabolite_identification")])
#'
#' ## The complete study from the MetaboLights FTP server (needs internet
#' ## access); with several assays the first MS assay is used unless
#' ## 'assayName' is given:
#' ## mtbls_to_mztabm("MTBLS103")
#' ## mtbls_to_mztabm("MTBLS103", file = "MTBLS103_hilic.mztab",
#' ##     assayName = "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt")
#' ## mtbls_to_mztabm("MTBLS40")
mtbls_to_mztabm <- function(mtblsId, file = paste0(mtblsId, ".mztab"),
                            assayName = character(), factor = character(),
                            path = NULL, verbose = TRUE) {
    stopifnot(is.character(mtblsId), length(mtblsId) == 1L, nzchar(mtblsId),
              is.character(file), length(file) == 1L)
    dat <- .mtbls_fetch(mtblsId, assayName, path)
    se <- .build_se(dat$isa, dat$maf, factor = factor)
    inv <- tryCatch(mtbls_investigation(mtblsId, path), error = function(e) {
        warning("Could not read the investigation file: ",
                conditionMessage(e), call. = FALSE)
        NULL
    })
    ftp <- if (is.null(path)) mtbls_ftp_path(mtblsId, mustWork = FALSE)
           else NULL
    m <- se_to_mztabm(se, id = mtblsId, isa = dat$isa, investigation = inv,
                      ftp_path = ftp)
    writeMzTabM(m, file)
    out <- if (grepl("\\.mz[tT]ab$", file)) file else paste0(file, ".mztab")
    if (!file.exists(out))
        stop("writing ", out, " failed")
    if (verbose) {
        lns <- readLines(out, warn = FALSE)
        mtd <- RmzTabM::mtd(m)
        placeholders <- c(
            if (!any(grepl("^instrument\\[", mtd[, 1L]))) "instrument",
            if (any(grepl("no database", mtd[mtd[, 1L] == "database[1]", 2L])))
                "compound database",
            if (!any(grepl("^contact\\[", mtd[, 1L]))) "contact (no ORCID)")
        message("Wrote ", out, "\n  assay: ", dat$assayName,
                "\n  MTD ", sum(startsWith(lns, "MTD\t")), " rows, SMF ",
                sum(startsWith(lns, "SMF\t")), " rows (", nrow(se),
                " features x ", ncol(se), " assays)",
                if (length(placeholders))
                    paste0("\n  not available from the ISA-TAB files: ",
                           paste(placeholders, collapse = ", ")))
    }
    attr(out, "se") <- se
    attr(out, "mztabm") <- m
    invisible(out)
}
