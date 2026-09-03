test_that("se_to_mztabm builds a complete MzTabM object", {
    se <- mtbls_summarized_experiment("MTBLS103", path = ext)
    inv <- mtbls_investigation("MTBLS103", path = ext)
    m <- se_to_mztabm(se, id = "MTBLS103", investigation = inv)
    expect_s4_class(m, "MzTabM")
    mtd <- RmzTabM::mtd(m)
    val <- function(key) unname(mtd[mtd[, 1L] == key, 2L])
    expect_identical(val("mzTab-ID"), "MTBLS103")
    expect_true(any(grepl("mtbls2mztabm", mtd[, 2L])))
    expect_match(val("title"), "HDL")
    expect_identical(val("publication[1]"),
                     "doi:10.1038/srep11496 | pubmed:26099471")
    expect_identical(val("database[1]-prefix"), "CHEBI")
    expect_identical(sum(grepl("^sample\\[[0-9]+\\]$", mtd[, 1L])), 32L)
    expect_identical(sum(grepl("^ms_run\\[[0-9]+\\]-location$", mtd[, 1L])),
                     32L)
    expect_identical(unique(unname(mtd[grepl("scan_polarity", mtd[, 1L]), 2L])),
                     "[MS, MS:1000130, positive scan, ]")
    expect_identical(unique(unname(mtd[grepl("-format$", mtd[, 1L]), 2L])),
                     "[MS, MS:1000566, ISB mzXML format, ]")
    expect_identical(unique(unname(mtd[grepl("-id_format$", mtd[, 1L]), 2L])),
                     "[MS, MS:1000776, scan number only nativeID format, ]")
    expect_identical(sum(grepl("^study_variable\\[[0-9]+\\]$", mtd[, 1L])),
                     3L)
    expect_identical(nrow(RmzTabM::smf(m)), 13L)
    ## no instrument without the assay table, no contact without ORCID
    expect_false(any(grepl("^instrument", mtd[, 1L])))
    expect_false(any(grepl("^contact", mtd[, 1L])))
})

test_that("se_to_mztabm validates its input", {
    se <- mtbls_summarized_experiment("MTBLS103", path = ext)
    expect_error(se_to_mztabm(se[, 1:3][, 0], "x"))
    bad <- se
    colData(bad)$species <- NULL
    expect_error(se_to_mztabm(bad, "MTBLS103"), "species")
    bad <- se
    rowData(bad)$mzmed <- NULL
    expect_error(se_to_mztabm(bad, "MTBLS103"), "mzmed")
    expect_error(se_to_mztabm(se, character()), "id")
})

test_that("mtbls_to_mztabm writes a valid file", {
    f <- tempfile(fileext = ".mztab")
    expect_message(out <- mtbls_to_mztabm("MTBLS103", file = f, path = ext),
                   "SMF 13 rows")
    expect_identical(as.character(out), f)
    expect_true(file.exists(f))
    lns <- readLines(f)
    expect_identical(lns[1L], "MTD\tmzTab-version\t2.1.0-M")
    expect_gt(sum(startsWith(lns, "MTD\t")), 100L)
    expect_identical(sum(startsWith(lns, "SMF\t")), 13L)
    expect_identical(sum(startsWith(lns, "SFH\t")), 1L)
    ## instrument comes from the assay table
    inst_name <- paste0("instrument[1]-name\t[MS, MS:1000031, ",
                        "instrument model, Agilent 6550 iFunnel Q-TOF]")
    expect_true(any(grepl(inst_name, lns, fixed = TRUE)))
    inst_src <- paste0("instrument[1]-source\t[MS, MS:1000073, ",
                       "electrospray ionization, ]")
    expect_true(any(grepl(inst_src, lns, fixed = TRUE)))
    ## round trip
    back <- RmzTabM::readMzTabM(f)
    expect_s4_class(back, "MzTabM")
    expect_s4_class(attr(out, "se"), "SummarizedExperiment")
    expect_s4_class(attr(out, "mztabm"), "MzTabM")
    ## extension is appended when missing
    f2 <- tempfile()
    out2 <- mtbls_to_mztabm("MTBLS103", file = f2, path = ext,
                            verbose = FALSE)
    expect_identical(as.character(out2), paste0(f2, ".mztab"))
    expect_true(file.exists(out2))
})
