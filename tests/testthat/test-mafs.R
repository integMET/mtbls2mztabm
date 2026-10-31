## A copy of the ISA-TAB files of MTBLS103.
study_copy <- function() {
    d <- tempfile()
    dir.create(d)
    file.copy(list.files(ext, full.names = TRUE), d)
    d
}

hilic_maf <- "m_MTBLS103_ibanez_02_HILIC_mass_spectrometry_v2_maf.tsv"
rp <- "a_MTBLS103_ibanez_02_RP_LC_mass_spectrometry.txt"
rp_maf <- "m_MTBLS103_ibanez_02_RP_LC_mass_spectrometry_v2_maf.tsv"

## Add an assay referring to the MAF `maf`, which is `nmaf` rows of the MAF of
## the bundled assay unless it exists already.
add_assay <- function(d, assay, maf, nmaf = 5L, polarity = TRUE) {
    lns <- readLines(file.path(d, hilic))
    lns <- gsub(hilic_maf, maf, lns, fixed = TRUE)
    writeLines(lns, file.path(d, assay))
    if (nzchar(maf) && !file.exists(file.path(d, maf)))
        writeLines(readLines(file.path(d, hilic_maf))[seq_len(nmaf + 1L)],
                   file.path(d, maf))
    if (!polarity)
        drop_polarity(d, assay)
    invisible(d)
}

drop_polarity <- function(d, assay) {
    lns <- readLines(file.path(d, assay))
    writeLines(gsub("\tpositive\t", "\t\t", lns, fixed = TRUE),
               file.path(d, assay))
}

test_that("a study with several MAFs gets an mzTab-M file for each", {
    d <- add_assay(study_copy(), rp, rp_maf)
    out_dir <- tempfile()
    dir.create(out_dir)
    f <- file.path(out_dir, "S.mztab")
    msgs <- capture_messages(
        out <- mtbls_to_mztabm("MTBLS103", file = f, path = d,
                               validate = FALSE))
    expect_true(any(grepl("has 2 metabolite assignment files", msgs)))
    ## the order of the investigation file: RP before HILIC
    expect_identical(names(out), c(rp_maf, hilic_maf))
    stems <- c("ibanez_02_RP_LC_mass_spectrometry",
               "ibanez_02_HILIC_mass_spectrometry")
    expect_identical(unname(as.character(out)),
                     file.path(out_dir, paste0("S_", stems, ".mztab")))
    expect_true(all(file.exists(out)))
    expect_identical(attr(out, "toml"),
                     file.path(out_dir, paste0("S_", stems, ".toml")))
    expect_identical(attr(out, "parquet"),
                     file.path(out_dir, paste0("S_", stems, ".parquet")))
    ## every file has the content of its own MAF and assay
    n <- nrow(mtbls_maf("MTBLS103", hilic_maf, path = d))
    expect_identical(nrow(nanoparquet::read_parquet(attr(out, "parquet")[1L])),
                     5L)
    expect_identical(nrow(nanoparquet::read_parquet(attr(out, "parquet")[2L])),
                     n)
    expect_identical(.toml_read(attr(out, "toml")[1L])$source,
                     list(assay_file = rp, maf_file = rp_maf))
    expect_identical(.toml_read(attr(out, "toml")[2L])$source,
                     list(assay_file = hilic, maf_file = hilic_maf))
    expect_identical(names(attr(out, "mztabm")), names(out))
    expect_s4_class(attr(out, "mztabm")[[1L]], "MzTabM")
    expect_identical(sum(startsWith(readLines(out[1L]), "SML\t")), 5L)
    expect_identical(sum(startsWith(readLines(out[2L]), "SML\t")), n)
})

test_that("verbose = FALSE silences the export of several MAFs", {
    d <- add_assay(study_copy(), rp, rp_maf)
    expect_silent(mtbls_to_mztabm("MTBLS103", path = d, validate = FALSE,
                                  file = tempfile(fileext = ".mztab"),
                                  verbose = FALSE))
})

test_that("the files of several MAFs pass the official validator", {
    skip_if_no_validator()
    d <- add_assay(study_copy(), rp, rp_maf)
    out <- suppressMessages(
        mtbls_to_mztabm("MTBLS103", file = tempfile(fileext = ".mztab"),
                        path = d, validate = TRUE))
    expect_length(out, 2L)
    expect_true(all(file.exists(out)))
})

test_that("explicit paths of the intermediate files get the MAF name as well", {
    d <- add_assay(study_copy(), rp, rp_maf)
    o <- tempfile()
    dir.create(o)
    out <- mtbls_to_mztabm("MTBLS103", file = file.path(o, "f.mzTab"), path = d,
                           toml = file.path(o, "t.toml"),
                           parquet = file.path(o, "p"), validate = FALSE,
                           verbose = FALSE)
    expect_identical(attr(out, "toml")[1L],
                     file.path(o, "t_ibanez_02_RP_LC_mass_spectrometry.toml"))
    expect_identical(attr(out, "parquet")[1L],
                     file.path(o, "p_ibanez_02_RP_LC_mass_spectrometry"))
    expect_identical(unname(as.character(out))[1L],
                     file.path(o, "f_ibanez_02_RP_LC_mass_spectrometry.mzTab"))
})

test_that("assayName writes a single file in a study with several MAFs", {
    d <- add_assay(study_copy(), rp, rp_maf)
    f <- tempfile(fileext = ".mztab")
    out <- mtbls_to_mztabm("MTBLS103", file = f, path = d, assayName = hilic,
                           validate = FALSE, verbose = FALSE)
    expect_identical(as.character(out), f)
    expect_null(names(out))
    expect_s4_class(attr(out, "mztabm"), "MzTabM")
    expect_identical(.toml_read(attr(out, "toml"))$source$maf_file, hilic_maf)
})

test_that("a MAF shared by several assays gets one file from the first assay", {
    d <- add_assay(study_copy(), rp, hilic_maf)
    f <- tempfile(fileext = ".mztab")
    msgs <- capture_messages(
        out <- mtbls_to_mztabm("MTBLS103", file = f, path = d,
                               validate = FALSE))
    expect_true(any(grepl("is shared by the assays", msgs)))
    expect_identical(as.character(out), f)
    expect_identical(.toml_read(attr(out, "toml"))$source$assay_file, rp)
    ## the other assay is available through assayName
    out <- mtbls_to_mztabm("MTBLS103", file = f, path = d, assayName = hilic,
                           validate = FALSE, verbose = FALSE)
    expect_identical(.toml_read(attr(out, "toml"))$source$assay_file, hilic)
})

test_that("a MAF no assay refers to gives a warning and no file", {
    d <- study_copy()
    file.copy(file.path(d, hilic_maf), file.path(d, "m_other_v2_maf.tsv"))
    f <- tempfile(fileext = ".mztab")
    expect_warning(
        out <- mtbls_to_mztabm("MTBLS103", file = f, path = d,
                               validate = FALSE, verbose = FALSE),
        "m_other_v2_maf.tsv.*no assay refers")
    expect_identical(as.character(out), f)
    expect_true(file.exists(f))
})

test_that("assays of another technology and without MAF are not exported", {
    d <- study_copy()
    nmr <- "a_MTBLS103_ibanez_02_NMR_spectroscopy.txt"
    add_assay(d, nmr, "m_nmr_v2_maf.tsv")
    f <- tempfile(fileext = ".mztab")
    msgs <- capture_messages(
        out <- mtbls_to_mztabm("MTBLS103", file = f, path = d,
                               validate = FALSE))
    expect_true(any(grepl(paste0("another technology: ", nmr), msgs)))
    expect_identical(as.character(out), f)
    expect_identical(.toml_read(attr(out, "toml"))$source$assay_file, hilic)
    ## an assay without MAF
    d <- add_assay(study_copy(), rp, "")
    msgs <- capture_messages(
        out <- mtbls_to_mztabm("MTBLS103", file = f, path = d,
                               validate = FALSE))
    expect_true(any(grepl(paste0("without metabolite assignment file: ", rp),
                          msgs)))
    expect_identical(.toml_read(attr(out, "toml"))$source$assay_file, hilic)
})

test_that("a MAF that fails does not stop the others", {
    d <- add_assay(study_copy(), rp, rp_maf, polarity = FALSE)
    expect_warning(
        out <- mtbls_to_mztabm("MTBLS103", file = tempfile(fileext = ".mztab"),
                               path = d, validate = FALSE, verbose = FALSE),
        paste0(rp_maf, ".*No scan polarity"))
    expect_identical(names(out), hilic_maf)
    expect_true(file.exists(out))
    expect_length(attr(out, "toml"), 1L)
    expect_length(attr(out, "mztabm"), 1L)
})

test_that("it is an error if all MAFs fail", {
    d <- add_assay(study_copy(), rp, rp_maf, polarity = FALSE)
    drop_polarity(d, hilic)
    expect_error(
        mtbls_to_mztabm("MTBLS103", file = tempfile(fileext = ".mztab"),
                        path = d, validate = FALSE, verbose = FALSE),
        "No mzTab-M file could be written for MTBLS103:.*No scan polarity")
})

test_that("a study of only another technology still gives the export's error", {
    d <- study_copy()
    nmr <- "a_MTBLS103_ibanez_02_NMR_spectroscopy.txt"
    file.rename(file.path(d, hilic), file.path(d, nmr))
    drop_polarity(d, nmr)
    expect_error(
        mtbls_to_mztabm("MTBLS103", file = tempfile(fileext = ".mztab"),
                        path = d, validate = FALSE, verbose = FALSE),
        "No scan polarity")
})

test_that(".maf_stems gives file name parts that tell the MAFs apart", {
    expect_identical(
        .maf_stems(c("m_MTBLS1_a b_v2_maf.tsv", "m_MTBLS1__c_maf.tsv",
                     "m_other-1_v2_maf.tsv", "m_MTBLS1_a_b_maf.tsv"),
                   "MTBLS1"),
        c("a_b", "c", "other-1", "a_b_1"))
})
