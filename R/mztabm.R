#' Build the mzTab-M metadata (MTD) section from a TOML study description
#'
#' `mtbls_toml_to_mtd()` reads a study description written by
#' [mtbls_factors_to_toml()] and builds the *metadata* (`MTD`) section of an
#' mzTab-M 2.1 file with the *RmzTabM* package:
#'
#' - the `[[sample]]` and `[[run]]` entries become the mzTab-M *samples*,
#'   *ms_runs* and *assays* (through [RmzTabM::mtdFromSampleData()]); the
#'   organism and organism part of a sample are written as the sample's
#'   species and tissue. The format and native identifier format of an MS
#'   run are derived from the extension of its data file (mzML, mzXML,
#'   mzData and netCDF files);
#' - the `[[factor]]` entries become the *study variable groups*, the groups
#'   (levels) of a factor the *study variables*, each referring to the assays
#'   of the samples in the group;
#' - title, description, publications and the MetaboLights URI of the study
#'   as well as the instrument, ion source, mass analyzer and detector and
#'   contacts (with ORCID) are added where the description provides them. All
#'   MS runs refer to the instrument.
#'
#' @param toml `character(1)` with the path of a TOML file written by
#'     [mtbls_factors_to_toml()].
#'
#' @param factor optional `character()` with the names of the study factors
#'     that become study variable groups. By default the first factor that is
#'     not *Age* is used (an experiment with a continuous factor such as the
#'     age would get one study variable per value); without any factor a
#'     single study variable `"undefined"` is defined. Use
#'     `factor = c("Genotype", "Treatment")` to use several factors.
#'
#' @param id optional `character(1)` with the mzTab-M identifier of the file.
#'     By default the study identifier of the TOML file.
#'
#' @return A two-column `character` matrix with the MTD content, as accepted
#'     by [RmzTabM::MzTabM()].
#'
#' @author Kozo Nishida
#'
#' @seealso [mtbls_factors_to_toml()] for the TOML file,
#'     [mtbls_parquet_to_sml()] for the second section and
#'     [mtd_sml_to_mztabm()] to combine them.
#'
#' @export
#'
#' @examples
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' f <- mtbls_factors_to_toml(
#'     "MTBLS103", file = tempfile(fileext = ".toml"), path = ext,
#'     assayName = "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt")
#' mtd <- mtbls_toml_to_mtd(f)
#' head(mtd, 12)
#' ## the study variables are the groups of the study factor
#' mtd[grepl("^study_variable(_group)?\\[[0-9]\\]$", mtd[, 1L]), ]
mtbls_toml_to_mtd <- function(toml, factor = character(), id = NULL) {
    x <- .read_study_toml(toml)
    if (is.null(id))
        id <- x$study$identifier
    stopifnot(is.character(id), length(id) == 1L, nzchar(id))
    groups <- .study_variable_groups(x, factor)
    sd <- .sample_data(x, groups)
    ver <- function(pkg) tryCatch(as.character(packageVersion(pkg)),
                                  error = function(e) "unknown")
    mtd <- mtdSkeleton(id = id, software = sprintf("[,, mtbls2mztabm, %s]",
                                                   ver("mtbls2mztabm")))
    mtd <- setMtdField(mtd, "software",
                       sprintf("[,, RmzTabM, %s]", ver("RmzTabM")))
    msrun <- c(location = ".mtbls_location", scan_polarity = ".mtbls_polarity")
    ## the validator rejects an instrument no MS run refers to
    if (!is.null(x$instrument))
        msrun <- c(msrun, instrument_ref = ".mtbls_instrument_ref")
    ## mzTab-M wants format and id_format together, or neither
    if (!anyNA(sd$.mtbls_format) && !anyNA(sd$.mtbls_id_format))
        msrun <- c(msrun, format = ".mtbls_format",
                   id_format = ".mtbls_id_format")
    runs <- mtdFromSampleData(
        sd, sampleCols. = c(sample = ".mtbls_sample"), msRunCols. = msrun,
        assayCols. = c(assay = ".mtbls_assay"),
        groups = .cv_name(groups),
        group_description = .group_descriptions(x, groups))
    ## The samples are written from the first run of each sample, in the
    ## order `mtdFromSampleData()` refers to them from the assays. Species
    ## and tissue are left out for samples without ontology term.
    runs <- runs[!grepl("^sample\\[[0-9]+\\]", runs[, 1L]), , drop = FALSE]
    first <- match(unique(sd$.mtbls_sample), sd$.mtbls_sample)
    terms <- function(v) lapply(v[first], function(z)
        if (is.na(z)) NULL else z)
    samples <- mtdSample(sample = sd$.mtbls_sample[first],
                         species = terms(sd$.mtbls_species),
                         tissue = terms(sd$.mtbls_tissue))
    mtd <- .mtd_study_info(rbind(mtd, samples, runs), x)
    mtdSort(mtd)
}

## Read a study description written by `mtbls_factors_to_toml()` into
## data.frames: `study`, `instrument`, `publications` and `contacts` (lists
## and data.frames), `factors`, `samples` and `runs`.
.read_study_toml <- function(file) {
    x <- .toml_read(file)
    if (is.null(x$schema_version) || is.null(x$study$identifier) ||
        is.null(x$run))
        stop("'", file, "' is not a study description written by ",
             "mtbls_factors_to_toml(): the 'schema_version', [study] or ",
             "[[run]] entries are missing")
    if (x$schema_version > .toml_schema_version)
        stop("'", file, "' has the layout version ", x$schema_version,
             ", but this version of mtbls2mztabm reads version ",
             .toml_schema_version, " at most")
    recs <- function(key, cols) {
        r <- x[[key]]
        res <- lapply(setNames(cols, cols), function(cl)
            vapply(r, function(e) {
                v <- e[[cl]]
                if (is.null(v) || !length(v)) NA_character_
                else as.character(v)[1L]
            }, character(1)))
        as.data.frame(res, stringsAsFactors = FALSE)
    }
    samples <- recs("sample", c("name", "organism", "organism_accession",
                                "organism_source", "organism_part",
                                "organism_part_accession",
                                "organism_part_source"))
    samples$factors <- lapply(x$sample, function(e) e$factors)
    factors <- recs("factor", c("name", "type", "type_term_accession",
                                "type_term_source"))
    factors$levels <- lapply(x$factor, function(e)
        as.character(unlist(e$levels)))
    runs <- recs("run", c("assay", "sample", "ms_assay_name", "data_file",
                          "polarity", "maf_column"))
    if (!nrow(runs))
        stop("The study description '", file, "' lists no runs")
    if (anyNA(runs$assay) || anyDuplicated(runs$assay))
        stop("The [[run]] entries of '", file, "' need unique 'assay' names")
    list(study = x$study, source = x$source, instrument = x$instrument,
         publications = recs("publication", c("doi", "pubmed_id")),
         contacts = recs("contact", c("name", "affiliation", "email",
                                      "orcid")),
         factors = factors, samples = samples, runs = runs)
}

## Names of the factors that become study variable groups (see the argument
## `factor` of `mtbls_toml_to_mtd()`).
.study_variable_groups <- function(x, factor = character()) {
    f <- x$factors
    has_values <- vapply(f$levels, length, integer(1)) > 0L
    if (length(factor)) {
        factor <- as.character(factor)
        miss <- setdiff(factor, f$name)
        if (length(miss))
            stop("Factor \"", miss[1L], "\" not found. Available factors: ",
                 paste(f$name, collapse = ", "))
        empty <- setdiff(factor, f$name[has_values])
        if (length(empty))
            stop("Factor \"", empty[1L], "\" has no values in the assay")
        return(unique(factor))
    }
    cand <- f$name[has_values & tolower(f$name) != "age"]
    if (length(cand)) cand[1L] else character()
}

## Description of each study variable group from the factor definition.
.group_descriptions <- function(x, groups) {
    f <- x$factors
    if (!length(groups))
        return(character())
    vapply(groups, function(g) {
        type <- f$type[match(g, f$name)]
        paste0("Study factor ", g, if (!is.na(type)) paste0(" (", type, ")"))
    }, character(1), USE.NAMES = FALSE)
}

## Sample data of the runs (one row per run) in the layout expected by
## `RmzTabM::mtdFromSampleData()`: the technical columns carry a ".mtbls_"
## prefix, one column per study variable group is named after its factor.
.sample_data <- function(x, groups = character()) {
    runs <- x$runs
    smp <- x$samples
    idx <- match(runs$sample, smp$name)
    if (anyNA(idx))
        warning("No [[sample]] entry for the sample(s) ",
                paste(unique(runs$sample[is.na(idx)]), collapse = ", "),
                call. = FALSE)
    pol <- .polarity_plain(runs$polarity)
    if (anyNA(pol))
        stop("No scan polarity for ", sum(is.na(pol)), " of ", nrow(runs),
             " runs. mzTab-M describes mass spectrometry runs, each with a ",
             "scan polarity: add 'polarity = \"positive\"' (or ",
             "\"negative\") to the [[run]] entries of the study ",
             "description, or write it again with the argument 'polarity'. ",
             "(An assay that is not a mass spectrometry assay, such as NMR, ",
             "cannot be exported.)")
    uri <- x$study$files_uri
    ## the study description keeps the deposited file names, a location has to
    ## be a valid URI (a name that is already encoded is not encoded again)
    loc <- paste0(if (length(uri)) uri, URLencode(runs$data_file))
    ## mzTab-M needs a location for every run and the validator does not
    ## accept "null": use the location of the study for runs without file
    nofile <- is.na(runs$data_file)
    if (any(nofile)) {
        study_uri <- c(uri, x$study$uri, "null")[1L]
        message("No data file for ", sum(nofile), " run(s): the study ",
                "location (", study_uri, ") is used as their ms_run location")
        loc[nofile] <- study_uri
    }
    sd <- data.frame(
        .mtbls_assay = runs$assay,
        .mtbls_sample = runs$sample,
        .mtbls_location = loc,
        .mtbls_format = .file_format_cv(runs$data_file),
        .mtbls_id_format = .id_format_cv(runs$data_file),
        .mtbls_polarity = pol,
        .mtbls_instrument_ref = 1L,
        .mtbls_species = .cv_param(smp$organism_accession[idx],
                                   smp$organism[idx],
                                   smp$organism_source[idx]),
        .mtbls_tissue = .cv_param(smp$organism_part_accession[idx],
                                  smp$organism_part[idx],
                                  smp$organism_part_source[idx]),
        stringsAsFactors = FALSE)
    for (g in groups) {
        v <- vapply(idx, function(i) {
            x <- if (is.na(i)) NULL else smp$factors[[i]][[g]]
            if (is.null(x)) NA_character_ else as.character(x)[1L]
        }, character(1))
        ## RmzTabM makes the runs without a value a study variable of their own
        if (anyNA(v))
            message("Factor \"", g, "\" has no value for ", sum(is.na(v)),
                    " of ", length(v), " runs (for example quality control ",
                    "samples): they form the study variable \"NA\". Add ",
                    "the value to the [sample.factors] of the study ",
                    "description to change this.")
        sd[[.cv_name(g)]] <- v
    }
    sd
}

## A name as it is written into a CV parameter "[,, name, ]": a name with a
## comma has to be quoted.
.cv_name <- function(x) {
    x <- as.character(x)
    quote <- grepl(",", x, fixed = TRUE)
    x[quote] <- paste0("\"", x[quote], "\"")
    x
}

## Title, description, publications, URI, instrument and contacts.
.mtd_study_info <- function(mtd, x) {
    st <- x$study
    for (key in c("title", "description"))
        if (length(st[[key]]) && !is.na(st[[key]]))
            mtd <- setMtdField(mtd, key, st[[key]])
    pub <- x$publications
    for (i in seq_len(nrow(pub))) {
        p <- c(if (!is.na(pub$doi[i])) paste0("doi:", pub$doi[i]),
               if (!is.na(pub$pubmed_id[i]))
                   paste0("pubmed:", pub$pubmed_id[i]))
        if (length(p))
            mtd <- setMtdField(mtd, "publication", paste(p, collapse = " | "))
    }
    if (length(st$uri))
        mtd <- setMtdField(mtd, "uri", st$uri)
    mtd <- .mtd_instrument(mtd, x$instrument)
    ct <- x$contacts
    for (i in seq_len(nrow(ct)))
        mtd <- setMtdContact(
            mtd, name = ct$name[i],
            affiliation = if (is.na(ct$affiliation[i])) "null"
                          else ct$affiliation[i],
            email = if (is.na(ct$email[i])) "null" else ct$email[i],
            orcid = ct$orcid[i])
    mtd
}

## Instrument section from the `[instrument]` table of the description.
.mtd_instrument <- function(mtd, inst) {
    if (is.null(inst))
        return(mtd)
    val <- function(k) {
        v <- inst[[k]]
        if (is.null(v)) NA_character_ else as.character(v)[1L]
    }
    cv <- function(k, parent_acc, parent_name)
        .ms_cv(val(paste0(k, "_accession")), val(k), parent_acc, parent_name)
    setMtdInstrument(
        mtd,
        name = cv("name", "MS:1000031", "instrument model"),
        source = cv("source", "MS:1000008", "ionization type"),
        analyzer = c(`analyzer[1]` = cv("analyzer", "MS:1000443",
                                        "mass analyzer type")),
        detector = cv("detector", "MS:1000026", "detector type"))
}

#' Build the mzTab-M small molecule (SML) section from a Parquet file
#'
#' `mtbls_parquet_to_sml()` reads the metabolite assignment file saved by
#' [mtbls_maf_to_parquet()] and builds the *small molecule summary* (`SML`)
#' section of an mzTab-M 2.1 file with the *RmzTabM* package
#' ([RmzTabM::smlCreate()]). Each row of the file becomes one small
#' molecule:
#'
#' - the abundances are those of the MAF column of each run, in the order of
#'   the assays of the metadata section (the columns are looked up through
#'   the `maf_column` of the `[[run]]` entries of the `toml` file); the
#'   average and variation of every study variable are added with
#'   [RmzTabM::smlAddStudyVariableColumns()];
#' - `database_identifier`, `chemical_formula`, `smiles`, `inchi` and `uri`
#'   are taken over from the MAF and `metabolite_identification` becomes the
#'   `chemical_name`. Placeholders such as `"unknown"` are reported as
#'   `null`. Where a molecule has several identifiers (separated by `|`) a
#'   single formula, structure or name is repeated for each of them. The
#'   annotations must be as many as the identifiers (a molecule without
#'   identifier has none): those that are not, are reported as `null` and
#'   kept in an optional column named after the MAF column, e.g.
#'   `opt_global_metabolite_identification` for the name of a peak that
#'   was not identified;
#' - `reliability` is mandatory in mzTab-M. The value of the MAF is used if
#'   it is an identification level of the Metabolomics Standards Initiative
#'   (`0` to `4`), which is the scheme declared in the metadata section;
#'   otherwise the level given with the argument `reliability`;
#' - the retention time and the m/z of the MAF, which the `SML` section has
#'   no column for, are kept as optional columns `opt_global_retention_time`
#'   and `opt_global_mass_to_charge` if the MAF reports any.
#'
#' Abundances are numbers; non-numeric entries and missing values are
#' reported as `null`. The abundance of a study variable is the mean over its
#' assays without the missing values and is `null` only if all of them are
#' missing. Its variation is calculated by `RmzTabM`, also without the missing
#' values, and is `null` if fewer than two abundances are available.
#'
#' @param parquet `character(1)` with the path of a Parquet file written by
#'     [mtbls_maf_to_parquet()].
#'
#' @param toml `character(1)` with the path of the TOML study description of
#'     the same assay, written by [mtbls_factors_to_toml()].
#'
#' @param mtd optional two-column matrix with the metadata section built from
#'     `toml` by [mtbls_toml_to_mtd()], from which the study variables are
#'     taken. By default it is built with the study factors selected by
#'     `factor`.
#'
#' @param factor optional `character()` with the study factors that become
#'     study variables, see [mtbls_toml_to_mtd()]. Ignored if `mtd` is given.
#'
#' @param reliability `character(1)`, a level from `"0"` to `"4"`, reported
#'     for the molecules for which the MAF gives no identification level. The
#'     default, `"4"`, is the lowest level (*unknown compound*) and thus
#'     makes no claim about the confidence of the identification.
#'
#' @return A `data.frame` with the SML section, as accepted by
#'     [RmzTabM::MzTabM()].
#'
#' @author Kozo Nishida
#'
#' @seealso [mtbls_maf_to_parquet()] for the Parquet file,
#'     [mtbls_toml_to_mtd()] for the metadata section and
#'     [mtd_sml_to_mztabm()] to combine them.
#'
#' @export
#'
#' @examples
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' assay <- "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt"
#' toml <- mtbls_factors_to_toml("MTBLS103", tempfile(fileext = ".toml"),
#'                               assayName = assay, path = ext)
#' parquet <- mtbls_maf_to_parquet("MTBLS103", tempfile(fileext = ".parquet"),
#'                                 assayName = assay, path = ext)
#' sml <- mtbls_parquet_to_sml(parquet, toml)
#' dim(sml)
#' sml[, c("SML_ID", "database_identifier", "chemical_name", "reliability",
#'         "abundance_assay[1]", "abundance_study_variable[1]")]
mtbls_parquet_to_sml <- function(parquet, toml, mtd = NULL,
                                 factor = character(), reliability = "4") {
    reliability <- as.character(reliability)
    if (length(reliability) != 1L || !grepl("^[0-4]$", reliability))
        stop("'reliability' must be one of \"0\", \"1\", \"2\", \"3\", \"4\"")
    maf <- .read_maf_parquet(parquet)
    if (!nrow(maf))
        stop("The metabolite assignment file '", parquet, "' has no rows")
    x <- .read_study_toml(toml)
    if (is.null(mtd))
        mtd <- mtbls_toml_to_mtd(toml, factor = factor)
    n_assay <- sum(grepl("^assay\\[[0-9]+\\]$", mtd[, 1L]))
    if (n_assay != nrow(x$runs))
        stop("The metadata section has ", n_assay, " assays, but '", toml,
             "' lists ", nrow(x$runs), " runs. Build both from the same ",
             "study description.")
    amat <- .maf_abundance(maf, x$runs$maf_column)
    ann <- .sml_annotations(maf, reliability)
    sml <- do.call(smlCreate, c(list(x = amat), ann))
    sml <- .sml_study_variable_means(smlAddStudyVariableColumns(sml, mtd), mtd)
    ## a list of numbers cannot be "null|null" (RmzTabM repeats a missing
    ## value for each identifier of a molecule)
    sml$theoretical_neutral_mass[
        grepl("^null(\\|null)+$", sml$theoretical_neutral_mass)] <- "null"
    sml
}

## `RmzTabM::smlAddStudyVariableColumns()` reports the mean of a study variable
## as missing if the abundance of any of its assays is: calculate the study
## variables with the average function "mean" from the available abundances.
.sml_study_variable_means <- function(sml, mtd) {
    key <- mtd[, 1L]
    avg <- grep("^study_variable\\[[0-9]+\\]-average_function$", key)
    for (i in avg[grepl("MS:1002962", mtd[avg, 2L], fixed = TRUE)]) {
        sv <- sub("-average_function$", "", key[i])
        refs <- strsplit(mtd[key == paste0(sv, "-assay_refs"), 2L], "|",
                         fixed = TRUE)[[1L]]
        ab <- as.matrix(sml[, paste0("abundance_", refs), drop = FALSE])
        m <- rowMeans(ab, na.rm = TRUE)
        ## no abundance at all: the mean is NaN, reported as null
        m[is.nan(m)] <- NA_real_
        sml[[paste0("abundance_", sv)]] <- unname(m)
    }
    sml
}

## Identifiers that mean "no identification".
.id_placeholders <- c("unknown", "na", "n/a", "nan", "null", "none", "-", "?")

## Annotation arguments of `smlCreate()` from the MAF `maf`.
.sml_annotations <- function(maf, reliability = "4") {
    n <- nrow(maf)
    col <- function(nm) {
        v <- if (nm %in% names(maf)) trimws(as.character(maf[[nm]]))
             else rep(NA_character_, n)
        v[!is.na(v) & !nzchar(v)] <- NA_character_
        v
    }
    ids <- col("database_identifier")
    ids[tolower(ids) %in% .id_placeholders] <- NA_character_
    n_id <- ifelse(is.na(ids), 1L, lengths(strsplit(ids, "|", fixed = TRUE)))
    ## The annotations are lists with one entry per identifier: a single value
    ## is repeated for each identifier of a molecule, a value with another
    ## number of entries cannot be reported.
    align <- function(v, name = FALSE) {
        k <- ifelse(is.na(v), NA_integer_,
                    lengths(strsplit(v, "|", fixed = TRUE)))
        if (name) {
            ## a name is free text: a "|" in it is not a separator
            fix <- !is.na(k) & k != 1L & k != n_id
            v[fix] <- gsub("|", "/", v[fix], fixed = TRUE)
            k[fix] <- 1L
        }
        one <- !is.na(k) & k == 1L & n_id > 1L
        v[one] <- vapply(which(one), function(i)
            paste(rep(v[i], n_id[i]), collapse = "|"), character(1))
        v[!is.na(k) & k != 1L & k != n_id] <- NA_character_
        v
    }
    ## MAF column of each annotation and the optional column keeping the
    ## values that cannot be reported in the standard one
    from <- c(chemical_formula = "chemical_formula", smiles = "smiles",
              inchi = "inchi", chemical_name = "metabolite_identification",
              uri = "uri")
    raw <- lapply(from, col)
    std <- Map(function(v, nm) align(v, name = nm == "chemical_name"), raw,
               names(raw))
    ## without an identifier there is nothing to annotate (the validator
    ## requires the lists to be as long as the identifier list)
    std <- lapply(std, function(v) {
        v[is.na(ids)] <- NA_character_
        v
    })
    res <- c(list(database_identifier = ids), std)
    lost <- Map(function(r, s) !is.na(r) & is.na(s), raw, std)
    for (k in names(from)[vapply(lost, any, logical(1))])
        res[[paste0("opt_global_", from[[k]])]] <-
            ifelse(lost[[k]], raw[[k]], NA_character_)
    if (any(unlist(lost)))
        message(sum(unlist(lost)), " annotation value(s) of molecules ",
                "without an identifier or with an unmatched number of ",
                "\"|\"-separated entries are kept in optional columns ",
                "(opt_global_*)")
    rel <- col("reliability")
    ok <- !is.na(rel) & grepl("^[0-4]$", rel)
    if (any(!ok))
        message("No identification level (0 to 4) in the MAF for ",
                sum(!ok), " of ", n, " molecules: reported as level ",
                reliability, " (argument 'reliability')")
    rel[!ok] <- reliability
    res$reliability <- rel
    ## retention time and m/z have no SML column: optional columns
    for (nm in c("retention_time", "mass_to_charge"))
        if (nm %in% names(maf) && any(!is.na(.as_numeric(maf[[nm]]))))
            res[[paste0("opt_global_", nm)]] <- .as_numeric(maf[[nm]])
    res
}

#' Combine metadata and small molecule sections into an mzTab-M object
#'
#' `mtd_sml_to_mztabm()` merges the `MTD` section built by
#' [mtbls_toml_to_mtd()] and the `SML` section built by
#' [mtbls_parquet_to_sml()] into an `RmzTabM::MzTabM` object, which
#' [RmzTabM::writeMzTabM()] writes as an mzTab-M 2.1 file.
#'
#' The sections are checked against each other (the abundance columns of the
#' `SML` must match the assays and the study variables of the `MTD`) and the
#' metadata is completed with what depends on both: the compound databases
#' are derived from the prefixes of the `database_identifier` column (ChEBI
#' and HMDB are recognised, other prefixes are reported by name) unless the
#' `MTD` defines databases already, and the unit of the small molecule
#' *feature* quantification is left out because the file has no feature
#' section. The `mzTab-profile` is set by `RmzTabM` according to the sections
#' present.
#'
#' @param mtd two-column `matrix` with the `MTD` section, see
#'     [mtbls_toml_to_mtd()].
#'
#' @param sml `data.frame` with the `SML` section, see
#'     [mtbls_parquet_to_sml()].
#'
#' @param file optional `character(1)` with the path of an mzTab-M file to
#'     write. An extension other than `.mztab`/`.mzTab` gets `.mztab`
#'     appended. By default no file is written.
#'
#' @param validate `logical(1)` whether to validate the written file with the
#'     official mzTab-M validator, see [RmzTabM::writeMzTabM()]. Only used if
#'     `file` is given.
#'
#' @return An `RmzTabM::MzTabM` object, invisibly if it was written to `file`.
#'
#' @author Kozo Nishida
#'
#' @seealso [mtbls_to_mztabm()] for the complete workflow.
#'
#' @export
#'
#' @examples
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' assay <- "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt"
#' toml <- mtbls_factors_to_toml("MTBLS103", tempfile(fileext = ".toml"),
#'                               assayName = assay, path = ext)
#' parquet <- mtbls_maf_to_parquet("MTBLS103", tempfile(fileext = ".parquet"),
#'                                 assayName = assay, path = ext)
#' mtd <- mtbls_toml_to_mtd(toml)
#' sml <- mtbls_parquet_to_sml(parquet, toml, mtd)
#' m <- mtd_sml_to_mztabm(mtd, sml)
#' m
#'
#' ## ... and write the mzTab-M file (validate = TRUE checks it with the
#' ## official validator, which RmzTabM downloads on first use)
#' f <- tempfile(fileext = ".mztab")
#' mtd_sml_to_mztabm(mtd, sml, file = f, validate = FALSE)
#' readLines(f, n = 6)
mtd_sml_to_mztabm <- function(mtd, sml, file = NULL, validate = TRUE) {
    stopifnot(is.matrix(mtd) || is.data.frame(mtd), ncol(mtd) == 2L,
              is.matrix(sml) || is.data.frame(sml),
              is.null(file) || (is.character(file) && length(file) == 1L))
    count <- function(pattern, x) sum(grepl(pattern, x))
    n_assay <- count("^assay\\[[0-9]+\\]$", mtd[, 1L])
    n_abund <- count("^abundance_assay\\[[0-9]+\\]$", colnames(sml))
    if (n_assay != n_abund)
        stop("The MTD section has ", n_assay, " assays, but the SML section ",
             "has ", n_abund, " abundance columns")
    n_var <- count("^study_variable\\[[0-9]+\\]$", mtd[, 1L])
    n_avg <- count("^abundance_study_variable\\[[0-9]+\\]$", colnames(sml))
    if (n_var != n_avg)
        stop("The MTD section has ", n_var, " study variables, but the SML ",
             "section has ", n_avg, " abundance_study_variable columns; ",
             "see RmzTabM::smlAddStudyVariableColumns()")
    ## no feature section: no feature quantification unit
    mtd <- mtd[mtd[, 1L] != "small_molecule_feature-quantification_unit", ,
               drop = FALSE]
    ## compound databases, unless the metadata defines them already
    no_db <- all(mtd[mtd[, 1L] == "database[1]-prefix", 2L] %in% "null")
    if (no_db && "database_identifier" %in% colnames(sml))
        mtd <- .mtd_databases(mtd, sml[["database_identifier"]])
    m <- MzTabM(mtd = mtdSort(mtd), sml = sml)
    if (is.null(file))
        return(m)
    if (!grepl("\\.mz[tT]ab$", file))
        file <- paste0(file, ".mztab")
    writeMzTabM(m, file, validate = validate)
    invisible(m)
}

## Compound databases from the prefixes of the identifiers in the SML
## ("PREFIX:id", several identifiers separated by "|").
.mtd_databases <- function(mtd, ids) {
    ids <- unlist(strsplit(as.character(ids), "|", fixed = TRUE))
    ids <- ids[grepl(":", ids, fixed = TRUE)]
    prefix <- unique(sub(":.*$", "", trimws(ids)))
    prefix <- prefix[!duplicated(toupper(prefix))]
    if (!length(prefix))
        return(mtd)
    known <- list(
        CHEBI = c("[MIRIAM, MIR:00000002, ChEBI, ]",
                  "https://www.ebi.ac.uk/chebi/"),
        HMDB = c("[MIRIAM, MIR:00100079, HMDB, ]", "https://hmdb.ca/"))
    ## the validator does not accept "null" as the URI of a database: point
    ## to the registry entry of an unknown prefix
    info <- lapply(prefix, function(p) {
        k <- known[[toupper(p)]]
        if (is.null(k))
            c(sprintf("[,, %s, ]", .cv_name(p)),
              paste0("https://registry.identifiers.org/registry?query=",
                     URLencode(p, reserved = TRUE)))
        else k
    })
    setMtdDatabase(mtd, name = vapply(info, `[`, character(1), 1L),
                   prefix = prefix, version = rep("unknown", length(prefix)),
                   uri = vapply(info, `[`, character(1), 2L), replace = TRUE)
}

#' Export a MetaboLights study to an mzTab-M file
#'
#' @description
#' The main function of the package: given a MetaboLights study identifier
#' it creates one mzTab-M 2.1 file for each metabolite assignment file (MAF)
#' of the study's mass spectrometry assays in the following steps, each of
#' which is also available as a function of its own:
#'
#' 1. the study factors (sample groups), the samples, runs and study metadata
#'    are saved to a TOML file with [mtbls_factors_to_toml()];
#' 2. the metabolite assignment file (MAF) of the assay is saved to a Parquet
#'    file with [mtbls_maf_to_parquet()];
#' 3. the metadata (`MTD`) section is built from the TOML file with
#'    [mtbls_toml_to_mtd()];
#' 4. the small molecule (`SML`) section is built from the Parquet file with
#'    [mtbls_parquet_to_sml()];
#' 5. both are combined with [mtd_sml_to_mztabm()] and written with
#'    [RmzTabM::writeMzTabM()].
#'
#' @details
#' The ISA-TAB files are retrieved through the *MsBackendMetaboLights*
#' package, or read from `path`. The two intermediate files are kept, next to
#' the mzTab-M file by default, so that they can be inspected or edited and
#' the file regenerated from them.
#'
#' After writing, the file is read back to count its `MTD` and `SML` rows,
#' and a message reports these together with any metadata still at a
#' placeholder value (no instrument, no compound database, no contact).
#'
#' A study with several MAFs gets one set of files per MAF, with the part of
#' the MAF name that tells them apart added to the paths: the MAF
#' `m_MTBLS1066_NEG_RP_LCMS_v2_maf.tsv` is written to
#' `MTBLS1066_NEG_RP_LCMS.mztab`. A MAF that fails gives a warning and the
#' others are still written; it is an error only if all of them fail. Some
#' MAFs need attention:
#'
#' - a MAF shared by several assays gets a single file, which describes the
#'   runs of the first of these assays (a message names them);
#' - a MAF that no assay refers to has no runs to describe and gets no file (a
#'   warning names it);
#' - the MAFs of assays of another technology, such as NMR, get no file unless
#'   the study has no mass spectrometry assay, which is then an error.
#'
#' @inheritParams mtbls_factors_to_toml
#' @inheritParams mtbls_toml_to_mtd
#' @inheritParams mtbls_parquet_to_sml
#'
#' @param assayName optional `character(1)` with an assay file (`"a_*.txt"`) to
#'     write a single mzTab-M file for its MAF, at the path `file`. By default
#'     a file is written for each MAF of the study.
#'
#' @param file `character(1)` with the path of the mzTab-M file to write.
#'     Default: `"<mtblsId>.mztab"` in the working directory. An extension
#'     other than `.mztab`/`.mzTab` gets `.mztab` appended. See the details for
#'     a study with several MAFs.
#'
#' @param toml,parquet `character(1)` with the paths of the TOML and Parquet
#'     files to write. By default they are written next to `file` with the
#'     extensions `.toml` and `.parquet`.
#'
#' @param validate `logical(1)` whether to validate the written file with the
#'     official mzTab-M validator, see [RmzTabM::writeMzTabM()]. The
#'     validator is downloaded by *RmzTabM* on first use, which needs
#'     internet access; the file is removed again if it is not valid.
#'
#' @param verbose `logical(1)` whether to report what is done: the written
#'     sections and the substitutions made for information the study does not
#'     provide (messages). With `FALSE` no message is shown.
#'
#' @return Invisibly, the path of the mzTab-M file, with the attributes
#'     `"toml"` and `"parquet"` (the paths of the intermediate files) and
#'     `"mztabm"` (the `MzTabM` object). With several MAFs, the paths of the
#'     files written, named by their MAF, with the same attributes in the same
#'     order: `"mztabm"` is then a `list` of `MzTabM` objects.
#'
#' @author Kozo Nishida
#'
#' @export
#'
#' @examples
#' ## From the ISA-TAB files of MTBLS103 shipped with the package
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' out <- mtbls_to_mztabm("MTBLS103", file = tempfile(fileext = ".mztab"),
#'                        path = ext, validate = FALSE)
#' readLines(out, n = 8)
#' RmzTabM::readMzTabM(out, validate = FALSE)
#'
#' ## the intermediate files the file was built from
#' readLines(attr(out, "toml"), n = 12)
#' nanoparquet::read_parquet(attr(out, "parquet"))[1:3, 1:6]
#'
#' ## MTBLS40, a GC-TOF-MS study of Arabidopsis thaliana genotypes whose
#' ## assay table lists raw netCDF files only: these become the ms_runs, the
#' ## factor Genotype the study variables
#' ext40 <- system.file("extdata", "MTBLS40", package = "mtbls2mztabm")
#' out40 <- mtbls_to_mztabm("MTBLS40", file = tempfile(fileext = ".mztab"),
#'                          path = ext40, validate = FALSE)
#' lns <- readLines(out40)
#' lns[grepl("^MTD\\s(ms_run\\[1\\]|instrument|study_variable\\[[1-4]\\]\\s)",
#'           lns)]
#'
#' ## The complete study from the MetaboLights FTP server (needs internet
#' ## access), with one file per MAF of its MS assays unless 'assayName' is
#' ## given. The files are checked with the mzTab-M validator:
#' ## mtbls_to_mztabm("MTBLS103")
#' ## mtbls_to_mztabm("MTBLS103", file = "MTBLS103_hilic.mztab",
#' ##     assayName = "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt")
#' ## mtbls_to_mztabm("MTBLS40")
mtbls_to_mztabm <- function(mtblsId, file = paste0(mtblsId, ".mztab"),
                            assayName = character(), factor = character(),
                            polarity = character(), path = NULL,
                            toml = NULL, parquet = NULL, reliability = "4",
                            validate = TRUE, verbose = TRUE) {
    stopifnot(is.character(mtblsId), length(mtblsId) == 1L, nzchar(mtblsId),
              is.character(file), length(file) == 1L)
    if (!grepl("\\.mz[tT]ab$", file))
        file <- paste0(file, ".mztab")
    base <- sub("\\.mz[tT]ab$", "", file)
    if (is.null(toml))
        toml <- paste0(base, ".toml")
    if (is.null(parquet))
        parquet <- paste0(base, ".parquet")
    one <- function(assay, file, toml, parquet)
        .export_mztabm(mtblsId, file, toml, parquet, assay, factor, polarity,
                       path, reliability, validate)
    export <- function() {
        assays <- if (length(assayName)) assayName
                  else .mtbls_maf_assays(mtblsId, path)
        if (length(assays) < 2L)
            return(one(assays, file, toml, parquet))
        .export_mafs(mtblsId, assays, one, file, toml, parquet)
    }
    invisible(if (verbose) export() else suppressMessages(export()))
}

## One mzTab-M file for each metabolite assignment file (the names of
## `assays`), with the part of the MAF name that tells them apart added to the
## file names. A MAF that fails does not stop the others.
.export_mafs <- function(mtblsId, assays, one, file, toml, parquet) {
    message("Study ", mtblsId, " has ", length(assays), " metabolite ",
            "assignment files to export: writing an mzTab-M file for each.")
    named <- function(p, stem) {
        ext <- file_ext(p)
        paste0(file_path_sans_ext(p), "_", stem,
               if (nzchar(ext)) paste0(".", ext))
    }
    res <- Map(function(assay, stem)
        tryCatch(one(assay, named(file, stem), named(toml, stem),
                     named(parquet, stem)),
                 error = function(e) e),
        assays, .maf_stems(names(assays), mtblsId))
    failed <- vapply(res, inherits, logical(1), "error")
    why <- vapply(names(res)[failed], function(m)
        paste0("\"", m, "\": ", conditionMessage(res[[m]])), character(1))
    if (all(failed))
        stop("No mzTab-M file could be written for ", mtblsId, ":\n",
             paste(why, collapse = "\n"), call. = FALSE)
    for (w in why)
        warning("No mzTab-M file for the metabolite assignment file ", w,
                call. = FALSE)
    res <- res[!failed]
    out <- vapply(res, as.character, character(1))
    attr(out, "toml") <- unname(vapply(res, attr, character(1), "toml"))
    attr(out, "parquet") <- unname(vapply(res, attr, character(1), "parquet"))
    attr(out, "mztabm") <- lapply(res, attr, "mztabm")
    out
}

## The steps of `mtbls_to_mztabm()`; returns the path of the mzTab-M file
## with the intermediate files and the MzTabM object as attributes.
.export_mztabm <- function(mtblsId, file, toml, parquet, assayName, factor,
                           polarity, path, reliability, validate) {
    dat <- .mtbls_fetch(mtblsId, assayName, path)
    .write_study_toml(mtblsId, dat, toml, polarity = polarity, path = path)
    .write_maf_parquet(dat$maf, parquet, mtblsId, dat$mafName)
    mtd <- mtbls_toml_to_mtd(toml, factor = factor)
    sml <- mtbls_parquet_to_sml(parquet, toml, mtd = mtd,
                                reliability = reliability)
    m <- mtd_sml_to_mztabm(mtd, sml, file = file, validate = validate)
    if (!file.exists(file))
        stop("writing ", file, " failed")
    lns <- readLines(file, warn = FALSE)
    mt <- RmzTabM::mtd(m)
    placeholders <- c(
        if (!any(grepl("^instrument\\[", mt[, 1L]))) "instrument",
        if (any(grepl("no database", mt[mt[, 1L] == "database[1]", 2L])))
            "compound database",
        if (!any(grepl("^contact\\[", mt[, 1L]))) "contact (no ORCID)")
    message("Wrote ", file, "\n  assay: ", dat$assayName,
            "\n  MAF: ", dat$mafName, "\n  intermediate files: ", toml,
            ", ", parquet,
            "\n  MTD ", sum(startsWith(lns, "MTD\t")), " rows, SML ",
            sum(startsWith(lns, "SML\t")), " rows (", nrow(sml),
            " small molecules x ", nrow(.mtd_assays(m)), " assays)",
            if (length(placeholders))
                paste0("\n  not available from the ISA-TAB files: ",
                       paste(placeholders, collapse = ", ")))
    attr(file, "toml") <- toml
    attr(file, "parquet") <- parquet
    attr(file, "mztabm") <- m
    file
}

## The assays defined in the MTD section of an MzTabM object.
.mtd_assays <- function(m) {
    mt <- RmzTabM::mtd(m)
    mt[grepl("^assay\\[[0-9]+\\]$", mt[, 1L]), , drop = FALSE]
}
