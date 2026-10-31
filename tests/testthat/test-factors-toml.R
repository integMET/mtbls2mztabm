test_that("mtbls_factors_to_toml saves the study description of MTBLS103", {
    f <- tempfile(fileext = ".toml")
    expect_silent(out <- mtbls_factors_to_toml("MTBLS103", f, path = ext))
    expect_identical(out, f)
    expect_true(file.exists(f))
    expect_match(readLines(f, n = 1L), "MTBLS103")
    x <- .toml_read(f)
    expect_identical(x$schema_version, 1L)
    expect_identical(x$study$identifier, "MTBLS103")
    expect_match(x$study$title, "HDL particles")
    expect_identical(x$study$uri, "https://www.ebi.ac.uk/metabolights/MTBLS103")
    ## local files: no location of the data files
    expect_null(x$study$files_uri)
    expect_identical(x$source$assay_file, hilic)
    expect_identical(x$source$maf_file,
                     "m_MTBLS103_ibanez_02_HILIC_mass_spectrometry_v2_maf.tsv")
    ## instrument, ion source and mass analyzer with their ontology terms
    expect_identical(x$instrument$name, "Agilent 6550 iFunnel Q-TOF")
    expect_identical(x$instrument$source, "electrospray ionization")
    expect_identical(x$instrument$source_accession,
                     "http://purl.obolibrary.org/obo/MS_1000073")
    expect_identical(x$instrument$analyzer, "quadrupole time-of-flight")
    expect_null(x$instrument$detector)
    expect_identical(x$publication,
                     list(list(doi = "10.1038/srep11496",
                               pubmed_id = "26099471")))
    ## no ORCID in the investigation file: no contact
    expect_null(x$contact)
})

test_that("the study factors and their sample groups are saved", {
    x <- .toml_read(study_files("MTBLS103")$toml)
    expect_length(x$factor, 1L)
    fac <- x$factor[[1L]]
    expect_identical(fac$name, "HIAE")
    expect_identical(fac$type, "Group")
    expect_match(fac$type_term_accession, "NCIT_C43359$")
    expect_identical(fac$type_term_source, "NCIT")
    expect_identical(fac$levels, c("CTR", "HIAE", "TREAT"))
    ## every sample carries its group
    expect_length(x$sample, 32L)
    groups <- vapply(x$sample, function(s) s$factors$HIAE, character(1))
    expect_identical(as.vector(table(groups)), c(14L, 12L, 6L))
    s1 <- x$sample[[1L]]
    expect_identical(s1$name, "CON_BASA_567795")
    expect_identical(s1$organism, "Homo sapiens")
    expect_identical(s1$organism_accession,
                     "http://purl.obolibrary.org/obo/NCBITaxon_9606")
    expect_identical(s1$organism_source, "NCBITAXON")
    expect_identical(s1$organism_part, "blood serum")
    expect_identical(s1$organism_part_source, "BTO")
    ## the values agree with the factor values of the sample file
    fv <- mtbls_factor_values("MTBLS103", path = ext)
    expect_identical(groups, fv$hiae[match(vapply(x$sample, `[[`, "", "name"),
                                           fv$sample_name)])
})

test_that("the runs map assays to samples, data files and MAF columns", {
    x <- .toml_read(study_files("MTBLS103")$toml)
    expect_length(x$run, 32L)
    r1 <- x$run[[1L]]
    expect_identical(r1$assay, "FILES/HILIC_CON_BASA_567795.mzXML")
    expect_identical(r1$sample, "CON_BASA_567795")
    expect_identical(r1$ms_assay_name, "CON_BASA_567795")
    expect_identical(r1$data_file, "FILES/HILIC_CON_BASA_567795.mzXML")
    expect_identical(r1$polarity, "positive")
    expect_identical(r1$maf_column, "CON_BASA_567795")
    maf <- mtbls_maf("MTBLS103", path = ext)
    cols <- vapply(x$run, `[[`, "", "maf_column")
    expect_true(all(cols %in% colnames(maf)))
    expect_false(anyDuplicated(cols) > 0L)
    assay <- vapply(x$run, `[[`, "", "assay")
    expect_false(anyDuplicated(assay) > 0L)
})

test_that("MTBLS40 (GC-MS, raw netCDF files) uses the raw files as data files", {
    f <- tempfile(fileext = ".toml")
    expect_message(mtbls_factors_to_toml("MTBLS40", f, path = ext40),
                   "raw spectral data file")
    x <- .toml_read(f)
    expect_identical(x$study$identifier, "MTBLS40")
    expect_length(x$run, 51L)
    expect_identical(x$run[[1L]]$data_file, "FILES/L_01_1.cdf")
    expect_identical(x$run[[1L]]$maf_column, "L_01_1")
    expect_true(all(endsWith(vapply(x$run, `[[`, "", "data_file"), ".cdf")))
    expect_identical(x$factor[[1L]]$name, "Genotype")
    expect_setequal(x$factor[[1L]]$levels, c("alkane", "Col-0", "mto1", "tt4"))
    expect_identical(x$sample[[1L]]$organism, "Arabidopsis thaliana")
    expect_identical(x$sample[[1L]]$organism_part_accession,
                     "http://purl.obolibrary.org/obo/BTO_0001658")
    expect_identical(x$instrument$name, "LECO Pegasus III")
    expect_length(x$publication, 2L)
    expect_identical(x$publication[[2L]]$pubmed_id, "17556050")
})

test_that("mtbls_factors_to_toml validates its input", {
    f <- tempfile(fileext = ".toml")
    expect_error(mtbls_factors_to_toml(c("a", "b"), f, path = ext))
    expect_error(mtbls_factors_to_toml("MTBLS103", f, path = "/no/dir"),
                 "existing directory")
    expect_error(mtbls_factors_to_toml("MTBLS103", f, path = ext,
                                       assayName = "a_none.txt"), "not found")
    expect_false(file.exists(f))
    expect_error(mtbls_factors_to_toml("MTBLS103", f, path = ext,
                                       polarity = "sideways"), "Unsupported")
})

test_that("a missing scan polarity can be filled in", {
    isa <- data.frame(`Sample Name` = c("A", "B"), check.names = FALSE,
                      `Parameter Value[Scan polarity]` = c("negative", NA),
                      stringsAsFactors = FALSE)
    r <- .toml_runs(isa, c("A", "B"))
    expect_identical(vapply(r, function(e) e$polarity, ""), c("negative", NA))
    r <- .toml_runs(isa, c("A", "B"), polarity = "Positive")
    expect_identical(vapply(r, `[[`, "", "polarity"),
                     c("negative", "positive"))
})

test_that("runs are named after the data file, the MS assay or the sample", {
    isa <- data.frame(
        `Sample Name` = c("S1", "S2", "S3", "S3"),
        `MS Assay Name` = c("M1", "M2", NA, NA),
        `Derived Spectral Data File` = c("f1.mzML", NA, NA, NA),
        check.names = FALSE, stringsAsFactors = FALSE)
    expect_warning(r <- .toml_runs(isa, c("S1", "S2", "S3")), "share 1")
    expect_identical(vapply(r, `[[`, "", "assay"),
                     c("f1.mzML", "M2", "S3", "S3.1"))
    expect_identical(vapply(r, `[[`, "", "sample"), c("S1", "S2", "S3", "S3"))
    isa$`Sample Name`[2L] <- NA
    expect_error(.toml_runs(isa, c("S1", "S2", "S3")),
                 "without a 'Sample Name'")
})

test_that("publications and contacts come from the investigation", {
    inv <- list(`Study Publication DOI` = c("doi:10.1/a", "", "https://doi.org/10.2/b"),
                `Study PubMed ID` = c("1", "2", ""))
    p <- .toml_publications(inv)
    expect_identical(p, list(list(doi = "10.1/a", pubmed_id = "1"),
                             list(doi = NA_character_, pubmed_id = "2"),
                             list(doi = "10.2/b", pubmed_id = NA_character_)))
    expect_identical(.toml_publications(list()), list())

    inv <- list(`Study Person Last Name` = c("Doe", "Roe", "Poe"),
                `Study Person First Name` = c("Jane", "", "Edgar"),
                `Study Person Affiliation` = c("RIKEN", "", ""),
                `Study Person Email` = c("j@x.org", "", ""),
                `Study Person ORCID` = c("https://orcid.org/0000-0002-1825-0097",
                                         "0000-0002-1825-0098", "none"))
    ct <- .toml_contacts(inv)
    expect_length(ct, 2L)
    expect_identical(ct[[1L]], list(name = "Jane Doe", affiliation = "RIKEN",
                                    email = "j@x.org",
                                    orcid = "0000-0002-1825-0097"))
    expect_identical(ct[[2L]]$name, "Roe")
    expect_true(is.na(ct[[2L]]$email))
    expect_identical(.toml_contacts(list(`Study Person Last Name` = "Doe")),
                     list())
})

test_that("factor definitions are matched by name and levels are listed", {
    isa <- data.frame(`Sample Name` = c("a", "b", "c"),
                      `Factor Value[Dose]` = c("low", NA, "high"),
                      `Factor Value[Empty]` = NA_character_,
                      `Factor Value[Undefined]` = c("x", "x", "y"),
                      check.names = FALSE, stringsAsFactors = FALSE)
    inv <- list(`Study Factor Name` = c("dose", "Other"),
                `Study Factor Type` = c("compound", "x"))
    f <- .toml_factors(isa, inv)
    expect_identical(vapply(f, `[[`, "", "name"),
                     c("Dose", "Empty", "Undefined"))
    expect_identical(f[[1L]]$type, "compound")
    expect_identical(unclass(f[[1L]]$levels), c("low", "high"))
    expect_length(f[[2L]]$levels, 0L)
    expect_true(is.na(f[[3L]]$type))
    expect_identical(.toml_factors(isa[, 1L, drop = FALSE], inv), list())
    ## samples: only the values a sample has
    s <- .toml_samples(isa)
    expect_identical(s[[1L]]$factors, list(Dose = "low", Undefined = "x"))
    expect_identical(s[[2L]]$factors, list(Undefined = "x"))
})

test_that("an assay without a metabolite assignment file is rejected", {
    d <- tempfile()
    dir.create(d)
    file.copy(list.files(ext, full.names = TRUE), d)
    a <- file.path(d, hilic)
    tab <- read.delim(a, check.names = FALSE, quote = "", comment.char = "")
    tab$`Metabolite Assignment File` <- NA
    write.table(tab, a, sep = "\t", quote = FALSE, row.names = FALSE,
                na = "")
    expect_error(mtbls_factors_to_toml("MTBLS103", tempfile(), path = d),
                 "must reference exactly one metabolite assignment file")
})
