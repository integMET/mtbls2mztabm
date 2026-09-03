test_that("mtbls_factor_values per sample", {
    fv <- mtbls_factor_values("MTBLS103", path = ext)
    expect_identical(colnames(fv), c("sample_name", "hiae"))
    expect_identical(nrow(fv), 32L)
    expect_identical(as.vector(table(fv$hiae)), c(14L, 12L, 6L))
    expect_identical(sort(unique(fv$hiae)), c("CTR", "HIAE", "TREAT"))
})

test_that("mtbls_factor_values per data file of an assay", {
    fa <- mtbls_factor_values("MTBLS103", assayName = hilic, path = ext)
    expect_identical(colnames(fa), c("sample_name", "data_file", "hiae"))
    expect_identical(nrow(fa), 32L)
    expect_true(all(startsWith(fa$data_file, "FILES/HILIC_")))
    expect_error(mtbls_factor_values("MTBLS103", assayName = "a_none.txt",
                                     path = ext), "not found")
})

test_that(".factor_values handles units and ontology qualifiers", {
    x <- data.frame(
        `Sample Name` = c("A", "B"),
        `Factor Value[Age]` = c(53, NA),
        `Unit` = c("year", "year"),
        `Term Source REF` = c("UO", "UO"),
        `Term Accession Number` = c("UO_0000036", "UO_0000036"),
        `Factor Value[Phenotype]` = c("CVD", NA),
        `Term Source REF` = c(NA, NA),
        `Term Accession Number` = c(NA, NA),
        check.names = FALSE, stringsAsFactors = FALSE)
    r <- .factor_values(x)
    expect_identical(colnames(r), c("sample_name", "age", "age_unit",
                                    "age_unit_term_source",
                                    "age_unit_term_accession", "phenotype"))
    expect_identical(r$age, c(53, NA))
    expect_identical(r$phenotype, c("CVD", NA))
    expect_warning(.factor_values(data.frame(`Sample Name` = "A",
                                             check.names = FALSE)),
                   "No 'Factor Value")
})
