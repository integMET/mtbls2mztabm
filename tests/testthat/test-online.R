test_that("the MetaboLights FTP server is used by default", {
    skip_if_offline("ftp.ebi.ac.uk")
    skip_on_cran()
    f <- tryCatch(mtbls_factors("MTBLS103"), error = function(e) e)
    skip_if(inherits(f, "error"), "MetaboLights FTP not reachable")
    expect_identical(f$name, "HIAE")
    fv <- mtbls_factor_values("MTBLS103")
    expect_identical(nrow(fv), 32L)
})
