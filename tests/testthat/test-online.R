test_that("the MetaboLights FTP server is used by default", {
    skip_if_offline("ftp.ebi.ac.uk")
    skip_on_cran()
    f <- tryCatch(mtbls_factors("MTBLS103"), error = function(e) e)
    skip_if(inherits(f, "error"), "MetaboLights FTP not reachable")
    expect_identical(f$name, "HIAE")
    fv <- mtbls_factor_values("MTBLS103")
    expect_identical(nrow(fv), 32L)
})

test_that("the study description is saved from the FTP server", {
    skip_if_offline("ftp.ebi.ac.uk")
    skip_on_cran()
    toml <- tempfile(fileext = ".toml")
    res <- tryCatch(mtbls_factors_to_toml("MTBLS103", toml, assayName = hilic),
                    error = function(e) e)
    skip_if(inherits(res, "error"), "MetaboLights FTP not reachable")
    x <- .toml_read(toml)
    ## the data files are located in the study folder of the server
    expect_match(x$study$files_uri, "ftp.*MTBLS103/$")
    expect_length(x$run, 32L)
    mtd <- mtbls_toml_to_mtd(toml)
    expect_match(mtd_value(mtd, "ms_run[1]-location"),
                 "^ftp.*MTBLS103/FILES/HILIC_CON_BASA_567795.mzXML$")
    ## the MAF of the assay was fetched as well
    expect_identical(x$source$maf_file,
                     "m_MTBLS103_ibanez_02_HILIC_mass_spectrometry_v2_maf.tsv")
    expect_true(all(vapply(x$run, function(r) !is.null(r$maf_column), NA)))
})
