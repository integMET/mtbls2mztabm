test_that("mtbls_maf reads the MAF", {
    maf <- mtbls_maf("MTBLS103", path = ext)
    expect_identical(dim(maf), c(13L, 53L))
    expect_identical(maf$metabolite_identification[1L], "L-glutamine")
    expect_true(all(is.na(maf$mass_to_charge)))
    expect_type(maf$CON_BASA_567795, "double")
})

test_that("mtbls_summarized_experiment matches the reference layout", {
    se <- mtbls_summarized_experiment("MTBLS103", path = ext)
    expect_s4_class(se, "SummarizedExperiment")
    expect_identical(dim(se), c(13L, 32L))
    expect_identical(assayNames(se), c("raw", "raw_filled"))
    expect_identical(assay(se, "raw"), assay(se, "raw_filled"))
    expect_false(anyNA(assay(se, "raw")))
    expect_identical(rownames(se)[c(1L, 13L)], c("FT01", "FT13"))
    expect_identical(colnames(se)[1L], "HILIC_CON_BASA_567795.mzXML")
    expect_equal(assay(se, "raw")[1L, 1L], 981143.812)

    cd <- colData(se)
    expect_identical(colnames(cd), ref_coldata)
    expect_type(cd$age, "integer")
    expect_type(cd$injection_index, "integer")
    expect_true(all(is.na(cd$age)))
    expect_identical(unique(cd$species),
                     "[NCBITaxon, NCBITaxon:9606, Homo sapiens, ]")
    expect_identical(unique(cd$tissue), "[BTO, BTO:0000133, blood serum, ]")
    expect_identical(unique(cd$polarity), "positive")
    expect_identical(as.vector(table(cd$phenotype)), c(14L, 12L, 6L))
    expect_identical(unique(cd$blood_sample_type), "blood serum")

    rd <- rowData(se)
    expect_identical(colnames(rd)[seq_len(12)],
                     c(ref_rowdata, "CTR", "HIAE", "TREAT", "ms_level",
                       "feature_id"))
    expect_identical(names(rd$mzmin), rownames(se))
    expect_true(all(rd$npeaks == 32))
    expect_true(all(rd$CTR == 14 & rd$HIAE == 12 & rd$TREAT == 6))
    expect_type(rd$ms_level, "integer")
    expect_identical(rd$metabolite_identification[8L], "choline")
    expect_identical(rd$database_identifier[1L], "CHEBI:18050")
    expect_identical(S4Vectors::metadata(se), list())
})

test_that("factor selection and errors", {
    se <- mtbls_summarized_experiment("MTBLS103", factor = "HIAE", path = ext)
    expect_identical(colData(se)$phenotype[1L], "CTR")
    expect_error(mtbls_summarized_experiment("MTBLS103", factor = "nope",
                                             path = ext), "not found")
    expect_error(mtbls_summarized_experiment("MTBLS103", path = "/no/dir"),
                 "existing directory")
})

test_that(".maf_abundance matches by sample name and reports missing ones", {
    maf <- data.frame(metabolite_identification = c("a", "b"),
                      S1 = c(1, 2), S2 = c(3, 4), stringsAsFactors = FALSE)
    isa <- data.frame(`Sample Name` = c("S1", "S2", "POOL"),
                      `Derived Spectral Data File` = c("f1", "f2", "f3"),
                      check.names = FALSE, stringsAsFactors = FALSE)
    expect_message(m <- .maf_abundance(maf, isa), "1 of 3 assay rows")
    expect_identical(dim(m), c(2L, 3L))
    expect_identical(m[, 1L], c(1, 2))
    expect_true(all(is.na(m[, 3L])))
    isa$`Sample Name` <- c("x", "y", "z")
    expect_error(.maf_abundance(maf, isa), "matches an abundance column")
})

test_that("a single-row MAF gives a one-row matrix", {
    maf <- data.frame(metabolite_identification = "a", S1 = 1, S2 = 3,
                      stringsAsFactors = FALSE)
    isa <- data.frame(`Sample Name` = c("S1", "S2"), check.names = FALSE,
                      stringsAsFactors = FALSE)
    expect_identical(dim(.maf_abundance(maf, isa)), c(1L, 2L))
})
