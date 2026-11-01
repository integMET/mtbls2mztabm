#' Save the factors and study metadata of a MetaboLights assay as TOML
#'
#' Reads the ISA-TAB files of a MetaboLights study and writes everything
#' needed to describe one of its mass spectrometry assays in the *metadata*
#' section of an mzTab-M file to a [TOML](https://toml.io) file. The file is
#' the first of the two intermediate files of the export workflow (the
#' second one, written by [mtbls_maf_to_parquet()], holds the metabolite
#' assignment file); [mtbls_toml_to_mtd()] turns it into the mzTab-M `MTD`
#' section.
#'
#' The file is human readable and can be edited, for example to add a scan
#' polarity or to correct an ontology term before the metadata section is
#' generated. It has the following layout (see `vignette("mtbls2mztabm")` for
#' an example):
#'
#' - top level: `schema_version`;
#' - `[study]`: `identifier`, `title`, `description`, `uri` and, for studies
#'   read from the MetaboLights FTP server, `files_uri` (the folder holding
#'   the data files);
#' - `[source]`: the `assay_file` and the metabolite assignment file
#'   (`maf_file`) the description belongs to;
#' - `[instrument]` (if the assay table lists one): `name`, `source`,
#'   `analyzer` and `detector`, each optionally with its ontology term as
#'   `name_accession`, `source_accession`, `analyzer_accession` and
#'   `detector_accession`;
#' - `[[publication]]` (`doi`, `pubmed_id`) and `[[contact]]` (`name`,
#'   `affiliation`, `email`, `orcid`; only persons with an ORCID, which
#'   mzTab-M requires for a contact);
#' - `[[factor]]`: the *study factors*, i.e. the definitions of the sample
#'   groups: `name`, `type`, `type_term_accession`, `type_term_source` and
#'   the `levels` (groups) found in the assay;
#' - `[[sample]]`: one entry per sample of the assay with `name`, the
#'   organism (`organism`, `organism_accession`, `organism_source`) and
#'   organism part (`organism_part`, `organism_part_accession`,
#'   `organism_part_source`) and, in the sub-table `[sample.factors]`, the
#'   value of each study factor, that is the sample group the sample belongs
#'   to;
#' - `[[run]]`: one entry per row of the assay table, i.e. per MS run and
#'   mzTab-M assay: `assay` (its unique name), `sample`, `ms_assay_name`,
#'   `data_file`, `polarity` and `maf_column`, the name of the column of the
#'   metabolite assignment file holding the abundances of the run.
#'
#' Fields the study does not provide are left out.
#'
#' The `maf_column` is the column named after the sample, the MS assay or the
#' data file of the run, whichever matches most runs. If none of these columns
#' holds a number, other columns are used only if they correspond one-to-one to
#' the runs, by equal names or by a name that is the other one plus a suffix
#' after a separator (`DFCI-01583` for the run `DFCI-01583-GC-01`); a message
#' reports this. Otherwise a warning says that the run has no abundances. Runs
#' sharing a column also give a warning: such a column has more weight in the
#' abundance of its study variable and the variation is understated.
#'
#' @inheritParams mtbls_investigation
#'
#' @param file `character(1)` with the path of the TOML file to write.
#'     Default: `"<mtblsId>.toml"` in the working directory.
#'
#' @param assayName optional `character(1)` with the assay file (`"a_*.txt"`)
#'     to describe. With a single assay it is used as is; with several and no
#'     `assayName` the first mass spectrometry assay is chosen and a message
#'     names it.
#'
#' @param polarity optional `character(1)`, `"positive"` or `"negative"`, used
#'     for the runs whose assay table lists no scan polarity. Without it, such
#'     runs are reported as positive, with a message, if their ion source is
#'     electron ionization.
#'
#' @return Invisibly, the path of the TOML file.
#'
#' @author Kozo Nishida
#'
#' @seealso [mtbls_maf_to_parquet()] for the second intermediate file,
#'     [mtbls_toml_to_mtd()] to read the file and [mtbls_to_mztabm()] for the
#'     complete workflow.
#'
#' @export
#'
#' @examples
#' ext <- system.file("extdata", "MTBLS103", package = "mtbls2mztabm")
#' f <- mtbls_factors_to_toml(
#'     "MTBLS103", file = tempfile(fileext = ".toml"), path = ext,
#'     assayName = "a_MTBLS103_ibanez_02_HILIC_mass_spectrometry.txt")
#' lns <- readLines(f)
#' ## study metadata and the study factor with its sample groups
#' head(lns, 22)
#' lns[grep("^\\[\\[factor\\]\\]", lns) + 0:4]
#'
#' ## From the MetaboLights FTP server (needs internet access):
#' ## mtbls_factors_to_toml("MTBLS103")
mtbls_factors_to_toml <- function(mtblsId, file = paste0(mtblsId, ".toml"),
                                  assayName = character(),
                                  polarity = character(), path = NULL) {
    stopifnot(is.character(mtblsId), length(mtblsId) == 1L, nzchar(mtblsId),
              is.character(file), length(file) == 1L)
    dat <- .mtbls_fetch(mtblsId, assayName, path)
    .write_study_toml(mtblsId, dat, file, polarity = polarity, path = path)
}

## Version of the TOML layout read by mtbls_toml_to_mtd().
.toml_schema_version <- 1L

## Write the study description of the assay `dat` (as returned by
## `.mtbls_fetch()`) to `file`.
.write_study_toml <- function(mtblsId, dat, file, polarity = character(),
                              path = NULL) {
    inv <- tryCatch(mtbls_investigation(mtblsId, path), error = function(e) {
        warning("Could not read the investigation file: ",
                conditionMessage(e), call. = FALSE)
        NULL
    })
    files_uri <- if (is.null(path)) mtbls_ftp_path(mtblsId, mustWork = FALSE)
                 else NULL
    x <- .study_description(mtblsId, dat, inv, files_uri, polarity)
    ver <- tryCatch(as.character(packageVersion("mtbls2mztabm")),
                    error = function(e) "unknown")
    .toml_write(x, file, comment = c(
        paste("Description of MetaboLights study", mtblsId, "for mzTab-M"),
        paste0("Written by mtbls2mztabm ", ver, "; read with ",
               "mtbls_toml_to_mtd()")))
    invisible(file)
}

## The study description as a nested list ready for `.toml_write()`.
.study_description <- function(mtblsId, dat, inv, files_uri = NULL,
                               polarity = character()) {
    isa <- dat$isa
    runs <- .toml_runs(isa, colnames(dat$maf), polarity,
                       filled = .maf_filled(dat$maf))
    list(schema_version = .toml_schema_version,
         study = list(identifier = mtblsId,
                      title = .inv_first(inv, "Study Title"),
                      description = .inv_first(inv, "Study Description"),
                      uri = paste0("https://www.ebi.ac.uk/metabolights/",
                                   mtblsId),
                      files_uri = files_uri),
         source = list(assay_file = dat$assayName, maf_file = dat$mafName),
         instrument = .toml_instrument(isa),
         publication = .toml_publications(inv),
         contact = .toml_contacts(inv),
         factor = .toml_factors(isa, inv),
         sample = .toml_samples(isa),
         run = runs)
}

## First non-missing, non-empty value of an investigation record.
.inv_first <- function(inv, key) {
    v <- .inv_field(inv, key)
    v <- v[!is.na(v) & nzchar(v)]
    if (length(v)) v[1L] else NA_character_
}

## `v[i]`, NA if there is no (non-empty) i-th value.
.nth_or_na <- function(v, i) {
    if (is.na(i) || i > length(v) || is.na(v[i]) || !nzchar(v[i]))
        NA_character_
    else v[i]
}

## Publications of the investigation (DOI and PubMed ID, matched by position).
.toml_publications <- function(inv) {
    doi <- .inv_field(inv, "Study Publication DOI")
    pmid <- .inv_field(inv, "Study PubMed ID")
    res <- lapply(seq_len(max(length(doi), length(pmid))), function(i) {
        d <- sub("^(doi:|https?://(dx\\.)?doi\\.org/)", "",
                 .nth_or_na(doi, i), ignore.case = TRUE)
        list(doi = d, pubmed_id = .nth_or_na(pmid, i))
    })
    res[!vapply(res, function(p) all(is.na(unlist(p))), logical(1))]
}

## Contacts with an ORCID (mzTab-M needs one; ISA-TAB investigation files
## rarely provide one).
.toml_contacts <- function(inv) {
    last <- .inv_field(inv, "Study Person Last Name")
    firstn <- .inv_field(inv, "Study Person First Name")
    aff <- .inv_field(inv, "Study Person Affiliation")
    email <- .inv_field(inv, "Study Person Email")
    okey <- grep("ORCID", names(inv), ignore.case = TRUE, value = TRUE)
    orcid <- if (length(okey)) .inv_field(inv, okey[1L]) else character()
    res <- lapply(seq_along(last), function(i) {
        oi <- sub("^https?://orcid.org/", "", .nth_or_na(orcid, i))
        if (is.na(oi) || !grepl("^\\d{4}-\\d{4}-\\d{4}-\\d{3}[0-9X]$", oi))
            return(NULL)
        nm <- c(.nth_or_na(firstn, i), .nth_or_na(last, i))
        nm <- trimws(paste(nm[!is.na(nm)], collapse = " "))
        if (!nzchar(nm))
            return(NULL)
        list(name = nm, affiliation = .nth_or_na(aff, i),
             email = .nth_or_na(email, i), orcid = oi)
    })
    res[!vapply(res, is.null, logical(1))]
}

## Instrument, ion source, mass analyzer and detector of the assay table
## (first value listed) with their ontology accessions.
.toml_instrument <- function(isa) {
    parts <- c(name = "Instrument", source = "Ion source",
               analyzer = "Mass analyzer", detector = "Detector")
    res <- list()
    for (k in names(parts)) {
        col <- paste0("Parameter Value[", parts[[k]], "]")
        v <- .col_or_na(isa, col)
        a <- .qualifier(isa, col, "Term Accession Number")
        res[[k]] <- v[!is.na(v)][1L]
        res[[paste0(k, "_accession")]] <- a[!is.na(a)][1L]
    }
    if (all(is.na(unlist(res)))) NULL else res
}

## The "Factor Value[...]" columns of the assay+sample table.
.factor_columns <- function(isa) {
    fcols <- grep("^Factor Value\\[", names(isa), value = TRUE)
    setNames(fcols, sub("^Factor Value\\[(.*)\\]$", "\\1", fcols))
}

## Study factors: the definition of the investigation file (matched by
## name) and the groups (levels) present in the assay.
.toml_factors <- function(isa, inv) {
    fcols <- .factor_columns(isa)
    def <- .inv_field(inv, "Study Factor Name")
    pick <- function(key, nm) {
        v <- .inv_field(inv, key)
        i <- match(tolower(nm), tolower(def))
        .nth_or_na(v, i)
    }
    lapply(seq_along(fcols), function(i) {
        nm <- names(fcols)[i]
        lev <- unique(as.character(isa[[fcols[i]]]))
        list(name = nm,
             type = pick("Study Factor Type", nm),
             type_term_accession = pick(
                 "Study Factor Type Term Accession Number", nm),
             type_term_source = pick("Study Factor Type Term Source REF", nm),
             levels = I(lev[!is.na(lev)]))
    })
}

## One entry per sample of the assay: organism, organism part and the value
## of every study factor.
.toml_samples <- function(isa) {
    fcols <- .factor_columns(isa)
    nm <- .col_or_na(isa, "Sample Name")
    first <- !duplicated(nm)
    term <- function(col, what) .qualifier(isa, col, what)
    org <- "Characteristics[Organism]"
    part <- "Characteristics[Organism part]"
    d <- data.frame(
        name = nm,
        organism = .col_or_na(isa, org),
        organism_accession = term(org, "Term Accession Number"),
        organism_source = term(org, "Term Source REF"),
        organism_part = .col_or_na(isa, part),
        organism_part_accession = term(part, "Term Accession Number"),
        organism_part_source = term(part, "Term Source REF"),
        stringsAsFactors = FALSE)[first, , drop = FALSE]
    fv <- lapply(fcols, function(col) as.character(isa[[col]])[first])
    lapply(seq_len(nrow(d)), function(i) {
        s <- as.list(d[i, , drop = FALSE])
        f <- lapply(fv, `[`, i)
        f <- f[!vapply(f, is.na, logical(1))]
        if (length(f)) s$factors <- f
        s
    })
}

## One entry per row of the assay table (MS run / mzTab-M assay).
.toml_runs <- function(isa, maf_names, polarity = character(),
                       filled = maf_names) {
    n <- nrow(isa)
    sample <- .col_or_na(isa, "Sample Name")
    if (anyNA(sample))
        stop("The assay table has rows without a 'Sample Name'")
    files <- .data_file(isa)
    ms_name <- .col_or_na(isa, "MS Assay Name")
    ## the assay is named after its data file, else after the MS assay
    assay <- make.unique(ifelse(!is.na(files), files,
                                ifelse(!is.na(ms_name), ms_name, sample)))
    pol <- .polarity_plain(.col_or_na(isa, "Parameter Value[Scan polarity]"))
    if (length(polarity)) {
        pol[is.na(pol)] <- .polarity_plain(polarity[1L])
    }
    ## electron ionization measures positive ions
    ei <- is.na(pol) &
        .is_electron_ionization(.col_or_na(isa, "Parameter Value[Ion source]"))
    if (any(ei)) {
        message("No scan polarity for ", sum(ei), " of ", n, " runs: reported ",
                "as positive because their ion source is electron ionization ",
                "(argument 'polarity')")
        pol[ei] <- "positive"
    }
    cols <- .maf_columns(isa, maf_names, filled)
    lapply(seq_len(n), function(i)
        list(assay = assay[i], sample = sample[i], ms_assay_name = ms_name[i],
             data_file = files[i], polarity = pol[i],
             maf_column = cols[i]))
}
