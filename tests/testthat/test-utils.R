test_that(".cv_param builds mzTab-M CV parameters", {
    expect_identical(
        .cv_param("http://purl.obolibrary.org/obo/BTO_0000133", "blood serum"),
        "[BTO, BTO:0000133, blood serum, ]")
    expect_identical(
        .cv_param("http://purl.obolibrary.org/obo/NCBITaxon_9606",
                  "Homo sapiens"),
        "[NCBITaxon, NCBITaxon:9606, Homo sapiens, ]")
    ## non-OBO style accession: label from the URL, or the Term Source REF
    expect_identical(
        .cv_param("http://purl.bioontology.org/ontology/NCBITAXON/9606",
                  "Homo sapiens"),
        "[NCBITAXON, NCBITAXON:9606, Homo sapiens, ]")
    expect_identical(.cv_param("CHEBI:1234", "x"), "[CHEBI, CHEBI:1234, x, ]")
    ## no accession: name only; no name: NA
    expect_identical(.cv_param(NA, "blood serum"), "[,, blood serum, ]")
    expect_identical(.cv_param(NA, NA), NA_character_)
    expect_length(.cv_param(c(NA, NA), c("a", "b")), 2L)
})

test_that("polarity and file format CV parameters", {
    expect_identical(.polarity_plain(c("positive", " Negative", "POS")),
                     c("positive", "negative", "positive"))
    expect_identical(.polarity_plain(c(NA, NA)), c(NA_character_, NA))
    ## missing values next to valid ones stay missing
    expect_identical(.polarity_plain(c("positive", NA)),
                     c("positive", NA_character_))
    expect_error(.polarity_plain(c("positive", "other")), "Unsupported")
    expect_error(.polarity_plain(c(NA, "other", "x")), "other, x")
    expect_identical(.file_format_cv(c("FILES/a.mzXML", "b.mzML", "c.txt")),
                     c("[MS, MS:1000566, ISB mzXML format, ]",
                       "[MS, MS:1000584, mzML format, ]", NA))
    expect_identical(.id_format_cv(c("FILES/a.mzXML", "b.mzML")),
                     c("[MS, MS:1000776, scan number only nativeID format, ]",
                       NA))
})

test_that(".ms_cv uses the ISA-TAB accession when present", {
    expect_identical(
        .ms_cv("http://purl.obolibrary.org/obo/MS_1000073",
               "electrospray ionization", "MS:1000008", "ionization type"),
        "[MS, MS:1000073, electrospray ionization, ]")
    expect_identical(
        .ms_cv(NA, "Agilent 6550", "MS:1000031", "instrument model"),
        "[MS, MS:1000031, instrument model, Agilent 6550]")
    expect_identical(.ms_cv(NA, NA, "MS:1000026", "detector type"),
                     "[MS, MS:1000026, detector type, ]")
})

test_that(".qualifier finds only directly following qualifier columns", {
    x <- data.frame(a = 1, `Unit` = "year", `Term Source REF` = "UO",
                    b = 2, `Term Source REF` = "X", check.names = FALSE)
    expect_identical(.qualifier(x, "a", "Unit"), "year")
    expect_identical(.qualifier(x, "a", "Term Source REF"), "UO")
    expect_identical(.qualifier(x, "b", "Term Source REF"), "X")
    expect_identical(.qualifier(x, "a", "Term Accession Number"),
                     NA_character_)
    expect_identical(.qualifier(x, "zz", "Unit"), NA_character_)
})

test_that(".check_path rejects non-directories", {
    expect_error(.check_path("/no/such/dir"), "existing directory")
    expect_null(.check_path(NULL))
})

test_that(".as_numeric coerces without warnings", {
    expect_no_warning(r <- .as_numeric(c("1.5", "abc", NA, "2e3", " 7 ", "")))
    expect_identical(r, c(1.5, NA, NA, 2000, 7, NA))
    expect_identical(.as_numeric(c(1, 2.5)), c(1, 2.5))
    expect_identical(.as_numeric(character()), numeric())
})
