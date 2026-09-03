ext40 <- system.file("extdata", "MTBLS40", package = "mtbls2mztabm")

test_that("MTBLS40 (GC-MS, raw netCDF files) becomes a SummarizedExperiment", {
    expect_message(
        se <- mtbls_summarized_experiment("MTBLS40", path = ext40),
        "raw spectral data file")
    expect_identical(dim(se), c(171L, 51L))
    expect_identical(colnames(colData(se)), ref_coldata)
    expect_true(all(endsWith(se$derived_spectra_data_file, ".cdf")))
    expect_identical(colnames(se)[1L], "L_01_1.cdf")
    expect_setequal(unique(se$phenotype), c("alkane", "Col-0", "mto1", "tt4"))
    expect_identical(
        se$species[1L],
        "[NCBITAXON, NCBITAXON:3702, Arabidopsis thaliana, ]")
    expect_identical(se$tissue[1L], "[BTO, BTO:0001658, aerial part, ]")
    expect_true(all(is.na(rowData(se)$mzmed)))
    expect_false(anyNA(rowData(se)$rtmed))
    expect_identical(rownames(se)[c(1L, 171L)], c("FT001", "FT171"))
})

test_that("mtbls_to_mztabm exports MTBLS40 with the raw files as ms_runs", {
    f <- tempfile(fileext = ".mztab")
    out <- suppressMessages(mtbls_to_mztabm("MTBLS40", file = f, path = ext40))
    lns <- readLines(f)
    has <- function(x) expect_true(any(grepl(x, lns, fixed = TRUE)))
    expect_identical(sum(startsWith(lns, "SMF\t")), 171L)
    expect_identical(sum(grepl("^MTD\tms_run\\[[0-9]+\\]-location\t", lns)),
                     51L)
    has("ms_run[1]-location\tFILES/L_01_1.cdf")
    has("ms_run[1]-format\t[MS, MS:1000825, Andi-MS format, ]")
    has(paste0("ms_run[1]-id_format\t[MS, MS:1000776, ",
               "scan number only nativeID format, ]"))
    has(paste0("instrument[1]-name\t[MS, MS:1000031, instrument model, ",
               "LECO Pegasus III]"))
    has("instrument[1]-analyzer[1]\t[MS, MS:1000084, time-of-flight, ]")
    expect_identical(sum(grepl("^MTD\tstudy_variable\\[[0-9]+\\]\t", lns)), 4L)
    has("publication[2]\tdoi:10.1016/j.jchromb.2007.05.002 | pubmed:17556050")
    has("database[1]-prefix\tCHEBI")
    expect_s4_class(RmzTabM::readMzTabM(f), "MzTabM")
})
