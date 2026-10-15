test_that("MTBLS40 (GC-MS, raw netCDF files) is exported with the raw files as ms_runs", {
    f <- tempfile(fileext = ".mztab")
    msgs <- capture_messages(
        out <- mtbls_to_mztabm("MTBLS40", file = f, path = ext40,
                               validate = FALSE))
    expect_match(msgs, "raw spectral data file", all = FALSE)
    lns <- readLines(f)
    has <- function(x) expect_true(any(grepl(x, lns, fixed = TRUE)))
    expect_identical(sum(startsWith(lns, "SML\t")), 171L)
    expect_identical(sum(startsWith(lns, "SMF\t")), 0L)
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
    has("study_variable_group[1]\t[,,Genotype,]")
    has("publication[2]\tdoi:10.1016/j.jchromb.2007.05.002 | pubmed:17556050")
    has("database[1]-prefix\tCHEBI")
    expect_s4_class(RmzTabM::readMzTabM(f, validate = FALSE), "MzTabM")
})

test_that("the MAF of MTBLS40 keeps its retention times and unidentified peaks", {
    f <- tempfile(fileext = ".mztab")
    suppressMessages(mtbls_to_mztabm("MTBLS40", file = f, path = ext40,
                                     validate = FALSE))
    lns <- readLines(f)
    smh <- strsplit(lns[startsWith(lns, "SMH\t")], "\t")[[1L]]
    sml <- do.call(rbind, strsplit(lns[startsWith(lns, "SML\t")], "\t"))
    expect_identical(nrow(sml), 171L)
    expect_true(all(c("opt_global_retention_time",
                      "opt_global_metabolite_identification") %in% smh))
    rt <- sml[, match("opt_global_retention_time", smh)]
    expect_identical(rt[1:3], c("1046", "1047", "1066"))
    name <- sml[, match("opt_global_metabolite_identification", smh)]
    id <- sml[, match("database_identifier", smh)]
    expect_true(all(name[id == "null"] != "null"))
    expect_true(all(name[id != "null"] == "null"))
    ## missing abundances are null
    ab <- sml[, grepl("^abundance_assay", smh)]
    expect_true(any(ab == "null"))
    expect_false(any(ab == "NA"))
})
