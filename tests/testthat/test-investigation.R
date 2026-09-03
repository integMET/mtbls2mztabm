test_that("mtbls_investigation parses the investigation file", {
    inv <- mtbls_investigation("MTBLS103", path = ext)
    expect_type(inv, "list")
    expect_identical(inv[["Study Identifier"]], "MTBLS103")
    expect_match(inv[["Study Title"]], "HDL particles")
    expect_length(inv[["Study Assay File Name"]], 4L)
    expect_identical(inv[["Study Publication DOI"]], "10.1038/srep11496")
    expect_false(any(names(inv) == "STUDY FACTORS"))
})

test_that("mtbls_factors returns the factor definitions", {
    f <- mtbls_factors("MTBLS103", path = ext)
    expect_s3_class(f, "data.frame")
    expect_identical(nrow(f), 1L)
    expect_identical(f$name, "HIAE")
    expect_identical(f$type, "Group")
    expect_match(f$type_term_accession, "NCIT_C43359$")
    expect_identical(f$type_term_source, "NCIT")
})
