test_that("mtbls_toml_to_mtd builds the metadata section of MTBLS103", {
    mtd <- mtbls_toml_to_mtd(study_files("MTBLS103")$toml)
    expect_true(is.matrix(mtd))
    expect_identical(ncol(mtd), 2L)
    expect_identical(mtd_value(mtd, "mzTab-version"), "2.1.0-M")
    expect_identical(mtd_value(mtd, "mzTab-ID"), "MTBLS103")
    expect_match(mtd_value(mtd, "title"), "HDL")
    expect_identical(mtd_value(mtd, "publication[1]"),
                     "doi:10.1038/srep11496 | pubmed:26099471")
    expect_identical(mtd_value(mtd, "uri[1]"),
                     "https://www.ebi.ac.uk/metabolights/MTBLS103")
    ## software: this package and RmzTabM
    expect_match(mtd_value(mtd, "software[1]"), "^\\[,, mtbls2mztabm, ")
    expect_match(mtd_value(mtd, "software[2]"), "^\\[,, RmzTabM, ")
    ## samples, MS runs and assays
    expect_identical(sum(grepl("^sample\\[[0-9]+\\]$", mtd[, 1L])), 32L)
    expect_identical(sum(grepl("^ms_run\\[[0-9]+\\]-location$", mtd[, 1L])), 32L)
    expect_identical(sum(grepl("^assay\\[[0-9]+\\]$", mtd[, 1L])), 32L)
    expect_identical(mtd_value(mtd, "sample[1]"), "CON_BASA_567795")
    expect_identical(mtd_value(mtd, "sample[1]-species[1]"),
                     "[NCBITaxon, NCBITaxon:9606, Homo sapiens, ]")
    expect_identical(mtd_value(mtd, "sample[1]-tissue[1]"),
                     "[BTO, BTO:0000133, blood serum, ]")
    expect_identical(mtd_value(mtd, "ms_run[1]-location"),
                     "FILES/HILIC_CON_BASA_567795.mzXML")
    expect_identical(mtd_value(mtd, "assay[1]"),
                     "FILES/HILIC_CON_BASA_567795.mzXML")
    expect_identical(mtd_value(mtd, "assay[1]-sample_ref"), "sample[1]")
    expect_identical(mtd_value(mtd, "assay[1]-ms_run_ref"), "ms_run[1]")
    expect_identical(unique(unname(mtd[grepl("-scan_polarity\\[1\\]$", mtd[, 1L]), 2L])),
                     "[MS, MS:1000130, positive scan, ]")
    expect_identical(unique(unname(mtd[grepl("-format$", mtd[, 1L]), 2L])),
                     "[MS, MS:1000566, ISB mzXML format, ]")
    expect_identical(unique(unname(mtd[grepl("-id_format$", mtd[, 1L]), 2L])),
                     "[MS, MS:1000776, scan number only nativeID format, ]")
    ## instrument of the assay table
    expect_identical(
        mtd_value(mtd, "instrument[1]-name"),
        "[MS, MS:1000031, instrument model, Agilent 6550 iFunnel Q-TOF]")
    expect_identical(mtd_value(mtd, "instrument[1]-source"),
                     "[MS, MS:1000073, electrospray ionization, ]")
    expect_identical(
        mtd_value(mtd, "instrument[1]-analyzer[1]"),
        "[MS, MS:1000443, mass analyzer type, quadrupole time-of-flight]")
    expect_identical(mtd_value(mtd, "instrument[1]-detector"),
                     "[MS, MS:1000026, detector type, ]")
    ## every MS run refers to the instrument
    refs <- mtd[grepl("^ms_run\\[[0-9]+\\]-instrument_ref$", mtd[, 1L]), 2L]
    expect_length(refs, 32L)
    expect_identical(unique(unname(refs)), "instrument[1]")
    ## no ORCID, no contact; databases are added when the SML is known
    expect_false(any(grepl("^contact", mtd[, 1L])))
})

test_that("MS runs refer to no instrument if the study has none", {
    f <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        x$instrument <- NULL
        x
    })
    mtd <- mtbls_toml_to_mtd(f)
    expect_false(any(grepl("^instrument\\[", mtd[, 1L])))
    expect_false(any(grepl("-instrument_ref$", mtd[, 1L])))
})

test_that("the factors become study variable groups with their sample groups", {
    mtd <- mtbls_toml_to_mtd(study_files("MTBLS103")$toml)
    expect_identical(mtd_value(mtd, "study_variable_group[1]"), "[,,HIAE,]")
    expect_identical(mtd_value(mtd, "study_variable_group[1]-description"),
                     "Study factor HIAE (Group)")
    expect_identical(
        mtd_value(mtd, "study_variable_group[1]-study_variable_refs"),
        "study_variable[1]|study_variable[2]|study_variable[3]")
    expect_identical(sum(grepl("^study_variable\\[[0-9]+\\]$", mtd[, 1L])), 3L)
    expect_identical(unname(mtd[grepl("^study_variable\\[[0-9]+\\]$", mtd[, 1L]), 2L]),
                     c("CTR", "HIAE", "TREAT"))
    ## the study variable refers to the assays of its sample group
    refs <- strsplit(mtd_value(mtd, "study_variable[3]-assay_refs"), "|",
                     fixed = TRUE)[[1L]]
    expect_identical(refs, paste0("assay[", 27:32, "]"))
    expect_identical(length(strsplit(
        mtd_value(mtd, "study_variable[1]-assay_refs"), "|", fixed = TRUE)[[1L]]),
        14L)
    expect_identical(mtd_value(mtd, "study_variable[1]-average_function"),
                     "[MS, MS:1002962, mean, ]")
})

test_that("mtbls_toml_to_mtd builds the metadata section of MTBLS40", {
    mtd <- mtbls_toml_to_mtd(study_files("MTBLS40")$toml)
    expect_identical(sum(grepl("^ms_run\\[[0-9]+\\]-location$", mtd[, 1L])), 51L)
    expect_identical(mtd_value(mtd, "ms_run[1]-location"), "FILES/L_01_1.cdf")
    expect_identical(mtd_value(mtd, "ms_run[1]-format"),
                     "[MS, MS:1000825, Andi-MS format, ]")
    expect_identical(mtd_value(mtd, "ms_run[1]-id_format"),
                     "[MS, MS:1000776, scan number only nativeID format, ]")
    expect_identical(
        mtd_value(mtd, "instrument[1]-name"),
        "[MS, MS:1000031, instrument model, LECO Pegasus III]")
    expect_identical(mtd_value(mtd, "instrument[1]-analyzer[1]"),
                     "[MS, MS:1000084, time-of-flight, ]")
    expect_identical(sum(grepl("^study_variable\\[[0-9]+\\]$", mtd[, 1L])), 4L)
    expect_setequal(unname(mtd[grepl("^study_variable\\[[0-9]+\\]$", mtd[, 1L]), 2L]),
                    c("alkane", "Col-0", "mto1", "tt4"))
    expect_identical(mtd_value(mtd, "study_variable_group[1]"), "[,,Genotype,]")
    expect_identical(mtd_value(mtd, "publication[2]"),
                     "doi:10.1016/j.jchromb.2007.05.002 | pubmed:17556050")
    expect_identical(mtd_value(mtd, "sample[1]-species[1]"),
                     "[NCBITAXON, NCBITAXON:3702, Arabidopsis thaliana, ]")
    expect_identical(mtd_value(mtd, "sample[1]-tissue[1]"),
                     "[BTO, BTO:0001658, aerial part, ]")
})

test_that("mtbls_toml_to_mtd selects the study factors", {
    toml <- study_files("MTBLS103")$toml
    ## a second factor, one with a comma in its name
    f <- modified_toml(toml, function(x) {
        x$factor <- c(x$factor, list(list(
            name = "Dose, mg", type = "dose", levels = I(c("low", "high")))))
        for (i in seq_along(x$sample))
            x$sample[[i]]$factors[["Dose, mg"]] <- c("low", "high")[1 + i %% 2]
        x
    })
    ## default: the first factor
    mtd <- mtbls_toml_to_mtd(f)
    expect_identical(sum(grepl("^study_variable_group\\[[0-9]+\\]$", mtd[, 1L])),
                     1L)
    ## several factors; a comma in the name is quoted in the CV parameter
    mtd <- mtbls_toml_to_mtd(f, factor = c("Dose, mg", "HIAE"))
    expect_identical(mtd_value(mtd, "study_variable_group[1]"), "[,,\"Dose, mg\",]")
    expect_identical(mtd_value(mtd, "study_variable_group[2]"), "[,,HIAE,]")
    expect_identical(sum(grepl("^study_variable\\[[0-9]+\\]$", mtd[, 1L])), 5L)
    expect_error(mtbls_toml_to_mtd(f, factor = "nope"),
                 "Factor \"nope\" not found. Available factors: HIAE, Dose, mg")
    ## no factor: a single undefined study variable
    f0 <- modified_toml(toml, function(x) {
        x$factor <- NULL
        for (i in seq_along(x$sample)) x$sample[[i]]$factors <- NULL
        x
    })
    mtd <- mtbls_toml_to_mtd(f0)
    expect_identical(mtd_value(mtd, "study_variable[1]"), "undefined")
    expect_error(mtbls_toml_to_mtd(f0, factor = "HIAE"), "not found")
})

test_that("runs without a value of the study factor are reported", {
    f <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        for (i in 1:3) x$sample[[i]]$factors <- NULL
        x
    })
    expect_message(mtd <- mtbls_toml_to_mtd(f),
                   "Factor \"HIAE\" has no value for 3 of 32 runs")
    ## they form a study variable of their own, as RmzTabM makes it
    expect_identical(sum(grepl("^study_variable\\[[0-9]+\\]$", mtd[, 1L])), 4L)
    expect_identical(
        mtd_value(mtd, "study_variable[1]-assay_refs"), "assay[1]|assay[2]|assay[3]")
    ## no message if no factor is used
    f0 <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        x$factor <- NULL
        x
    })
    expect_silent(mtbls_toml_to_mtd(f0))
})

test_that("Age is not a study variable group unless requested", {
    f <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        x$factor <- c(list(list(name = "Age", levels = I(c("12", "13")))),
                      x$factor)
        for (i in seq_along(x$sample))
            x$sample[[i]]$factors$Age <- c("12", "13")[1 + i %% 2]
        x
    })
    mtd <- mtbls_toml_to_mtd(f)
    expect_identical(mtd_value(mtd, "study_variable_group[1]"), "[,,HIAE,]")
    mtd <- mtbls_toml_to_mtd(f, factor = "Age")
    expect_identical(mtd_value(mtd, "study_variable_group[1]"), "[,,Age,]")
    expect_identical(sum(grepl("^study_variable\\[[0-9]+\\]$", mtd[, 1L])), 2L)
    ## a factor without values cannot be used
    f2 <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        x$factor[[1L]]$levels <- NULL
        x
    })
    expect_error(mtbls_toml_to_mtd(f2, factor = "HIAE"), "has no values")
    expect_identical(mtd_value(mtbls_toml_to_mtd(f2), "study_variable[1]"),
                     "undefined")
})

test_that("mtbls_toml_to_mtd uses the location of remote data files", {
    f <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        x$study$files_uri <- "ftp://ftp.example.org/MTBLS103/"
        x
    })
    mtd <- mtbls_toml_to_mtd(f)
    expect_identical(mtd_value(mtd, "ms_run[1]-location"),
                     "ftp://ftp.example.org/MTBLS103/FILES/HILIC_CON_BASA_567795.mzXML")
    ## the assay keeps the name of the data file
    expect_identical(mtd_value(mtd, "assay[1]"),
                     "FILES/HILIC_CON_BASA_567795.mzXML")
})

test_that("runs without data file get the study location", {
    f <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        for (i in 1:2) {
            x$run[[i]]$data_file <- NULL
            x$run[[i]]$assay <- x$run[[i]]$ms_assay_name
        }
        x
    })
    expect_message(mtd <- mtbls_toml_to_mtd(f), "No data file for 2 run")
    expect_identical(mtd_value(mtd, "ms_run[1]-location"),
                     "https://www.ebi.ac.uk/metabolights/MTBLS103")
    expect_identical(mtd_value(mtd, "assay[1]"), "CON_BASA_567795")
    expect_identical(mtd_value(mtd, "ms_run[3]-location"),
                     "FILES/HILIC_CON_BASA_580552.mzXML")
    ## formats are only reported if every run has one
    expect_false(any(grepl("-format$", mtd[, 1L])))
    expect_false(any(grepl("-id_format$", mtd[, 1L])))
})

test_that("samples without ontology terms get no species or tissue", {
    f <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        for (i in 1:3)
            x$sample[[i]][c("organism", "organism_accession", "organism_source",
                            "organism_part", "organism_part_accession",
                            "organism_part_source")] <- NULL
        x
    })
    mtd <- mtbls_toml_to_mtd(f)
    expect_false(any(grepl("^sample\\[[1-3]\\]-(species|tissue)", mtd[, 1L])))
    expect_identical(mtd_value(mtd, "sample[4]-species[1]"),
                     "[NCBITaxon, NCBITaxon:9606, Homo sapiens, ]")
    expect_identical(sum(grepl("^sample\\[[0-9]+\\]-species\\[1\\]$", mtd[, 1L])),
                     29L)
    f2 <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        for (i in seq_along(x$sample))
            x$sample[[i]][grepl("^organism", names(x$sample[[i]]))] <- NULL
        x
    })
    mtd <- mtbls_toml_to_mtd(f2)
    expect_false(any(grepl("-(species|tissue)\\[", mtd[, 1L])))
    expect_identical(sum(grepl("^sample\\[[0-9]+\\]$", mtd[, 1L])), 32L)
})

test_that("a sample measured twice is one sample with two assays", {
    f <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        r <- x$run[[1L]]
        r$assay <- "FILES/replicate.mzXML"
        r$data_file <- "FILES/replicate.mzXML"
        x$run <- c(x$run, list(r))
        x
    })
    mtd <- mtbls_toml_to_mtd(f)
    expect_identical(sum(grepl("^sample\\[[0-9]+\\]$", mtd[, 1L])), 32L)
    expect_identical(sum(grepl("^assay\\[[0-9]+\\]$", mtd[, 1L])), 33L)
    expect_identical(mtd_value(mtd, "assay[33]-sample_ref"), "sample[1]")
    expect_identical(mtd_value(mtd, "assay[33]-ms_run_ref"), "ms_run[33]")
    expect_match(mtd_value(mtd, "study_variable[1]-assay_refs"),
                 "assay[33]", fixed = TRUE)
})

test_that("the scan polarity of every run is needed", {
    f <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        for (i in 1:2) x$run[[i]]$polarity <- NULL
        x
    })
    expect_error(mtbls_toml_to_mtd(f), "No scan polarity for 2 of 32 runs")
    f2 <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        for (i in 1:2) x$run[[i]]$polarity <- "Negative"
        x
    })
    mtd <- mtbls_toml_to_mtd(f2)
    expect_identical(mtd_value(mtd, "ms_run[1]-scan_polarity[1]"),
                     "[MS, MS:1000129, negative scan, ]")
    expect_identical(mtd_value(mtd, "ms_run[3]-scan_polarity[1]"),
                     "[MS, MS:1000130, positive scan, ]")
})

test_that("contacts with an ORCID and publications are written", {
    f <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        x$contact <- list(
            list(name = "Jane Doe", affiliation = "RIKEN", email = "j@x.org",
                 orcid = "0000-0002-1825-0097"),
            list(name = "Roe", orcid = "0000-0002-1825-0098"))
        x$publication <- list(list(doi = "10.1/a"), list(pubmed_id = "2"))
        x
    })
    mtd <- mtbls_toml_to_mtd(f)
    expect_identical(mtd_value(mtd, "contact[1]-name"), "Jane Doe")
    expect_identical(mtd_value(mtd, "contact[1]-affiliation"), "RIKEN")
    expect_identical(mtd_value(mtd, "contact[1]-orcid"), "0000-0002-1825-0097")
    expect_identical(mtd_value(mtd, "contact[2]-email"), "null")
    expect_identical(mtd_value(mtd, "publication[1]"), "doi:10.1/a")
    expect_identical(mtd_value(mtd, "publication[2]"), "pubmed:2")
})

test_that("mtbls_toml_to_mtd checks the study description", {
    expect_error(mtbls_toml_to_mtd(tempfile()), "does not exist")
    f <- tempfile(fileext = ".toml")
    writeLines("a = 1", f)
    expect_error(mtbls_toml_to_mtd(f), "not a study description")
    newer <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        x$schema_version <- 99L
        x
    })
    expect_error(mtbls_toml_to_mtd(newer), "layout version 99")
    dup <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        x$run[[2L]]$assay <- x$run[[1L]]$assay
        x
    })
    expect_error(mtbls_toml_to_mtd(dup), "unique 'assay' names")
    noruns <- modified_toml(study_files("MTBLS103")$toml, function(x) {
        x$run <- NULL
        x
    })
    expect_error(mtbls_toml_to_mtd(noruns), "not a study description")
    expect_error(mtbls_toml_to_mtd(study_files("MTBLS103")$toml, id = c("a", "b")))
    ## the identifier can be changed
    mtd <- mtbls_toml_to_mtd(study_files("MTBLS103")$toml, id = "MY-ID")
    expect_identical(mtd_value(mtd, "mzTab-ID"), "MY-ID")
})
