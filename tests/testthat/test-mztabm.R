test_that("mtd_sml_to_mztabm combines the sections into an MzTabM object", {
    files <- study_files("MTBLS103")
    mtd <- mtbls_toml_to_mtd(files$toml)
    sml <- suppressMessages(mtbls_parquet_to_sml(files$parquet, files$toml, mtd))
    m <- mtd_sml_to_mztabm(mtd, sml)
    expect_s4_class(m, "MzTabM")
    mt <- RmzTabM::mtd(m)
    ## the profile follows the sections present
    expect_identical(mtd_value(mt, "mzTab-profile"), "M+S")
    expect_identical(nrow(RmzTabM::sml(m)), 13L)
    expect_identical(nrow(RmzTabM::smf(m)), 0L)
    expect_identical(nrow(RmzTabM::sme(m)), 0L)
    expect_identical(RmzTabM::sml(m), sml)
    ## no feature section, no unit for its quantification
    expect_false("small_molecule_feature-quantification_unit" %in% mt[, 1L])
    expect_true("small_molecule-quantification_unit" %in% mt[, 1L])
    ## the compound database follows from the identifiers of the SML
    expect_identical(mtd_value(mt, "database[1]"),
                     "[MIRIAM, MIR:00000002, ChEBI, ]")
    expect_identical(mtd_value(mt, "database[1]-prefix"), "CHEBI")
    expect_identical(mtd_value(mt, "database[1]-uri"),
                     "https://www.ebi.ac.uk/chebi/")
    expect_false(any(grepl("no database", mt[, 2L])))
    ## the rest of the metadata is untouched
    expect_identical(mtd_value(mt, "mzTab-ID"), "MTBLS103")
    expect_identical(mt[!grepl("^database|feature-quant|mzTab-profile", mt[, 1L]), ],
                     mtd[!grepl("^database|feature-quant|mzTab-profile", mtd[, 1L]), ])
})

test_that(".mtd_databases derives the databases from the identifiers", {
    mtd <- RmzTabM::mtdSkeleton(id = "x", software = "[,, a, 1]")
    expect_identical(mtd_value(mtd, "database[1]-prefix"), "null")
    ## no identifiers: the placeholder stays
    expect_identical(.mtd_databases(mtd, c("null", NA, "unknown")), mtd)
    ## several prefixes, several identifiers per molecule, different case
    ids <- c("CHEBI:1", "HMDB:HMDB0000122|KEGG:C00031", "ChEBI:2", "null",
             "HMDB:HMDB1")
    d <- .mtd_databases(mtd, ids)
    expect_identical(mtd_value(d, "database[1]-prefix"), "CHEBI")
    expect_identical(mtd_value(d, "database[2]-prefix"), "HMDB")
    expect_identical(mtd_value(d, "database[3]-prefix"), "KEGG")
    expect_false("database[4]-prefix" %in% d[, 1L])
    expect_identical(mtd_value(d, "database[2]"),
                     "[MIRIAM, MIR:00100079, HMDB, ]")
    ## unknown databases point to the registry entry of their prefix
    expect_identical(mtd_value(d, "database[3]"), "[,, KEGG, ]")
    expect_identical(mtd_value(d, "database[3]-uri"),
                     "https://registry.identifiers.org/registry?query=KEGG")
    expect_identical(mtd_value(d, "database[3]-version"), "unknown")
    ## names with a comma are quoted
    d2 <- .mtd_databases(mtd, "A,B:1")
    expect_identical(mtd_value(d2, "database[1]"), "[,, \"A,B\", ]")
})

test_that("mtd_sml_to_mztabm writes the file if asked to", {
    files <- study_files("MTBLS103")
    mtd <- mtbls_toml_to_mtd(files$toml)
    sml <- suppressMessages(mtbls_parquet_to_sml(files$parquet, files$toml, mtd))
    f <- tempfile()
    expect_invisible(m <- mtd_sml_to_mztabm(mtd, sml, file = f,
                                            validate = FALSE))
    expect_s4_class(m, "MzTabM")
    expect_true(file.exists(paste0(f, ".mztab")))
    expect_false(file.exists(f))
    lns <- readLines(paste0(f, ".mztab"))
    expect_identical(sum(startsWith(lns, "SML\t")), 13L)
    expect_identical(lns[1L], "MTD\tmzTab-version\t2.1.0-M")
    ## the same as writing the object
    f2 <- tempfile(fileext = ".mzTab")
    RmzTabM::writeMzTabM(m, f2, validate = FALSE)
    expect_identical(readLines(f2), lns)
    ## without file nothing is written and the object is returned visibly
    expect_visible(mtd_sml_to_mztabm(mtd, sml))
    expect_error(mtd_sml_to_mztabm(mtd, sml, file = c("a", "b")))
})

test_that("mtd_sml_to_mztabm keeps databases defined in the metadata", {
    files <- study_files("MTBLS103")
    mtd <- mtbls_toml_to_mtd(files$toml)
    sml <- suppressMessages(mtbls_parquet_to_sml(files$parquet, files$toml, mtd))
    mtd <- RmzTabM::setMtdDatabase(mtd, name = "[,, my database, ]",
                                   prefix = "CHEBI", version = "2024",
                                   uri = "https://example.org/", replace = TRUE)
    mt <- RmzTabM::mtd(mtd_sml_to_mztabm(mtd, sml))
    expect_identical(mtd_value(mt, "database[1]"), "[,, my database, ]")
    expect_identical(mtd_value(mt, "database[1]-version"), "2024")
})

test_that("mtd_sml_to_mztabm checks that the sections belong together", {
    files <- study_files("MTBLS103")
    mtd <- mtbls_toml_to_mtd(files$toml)
    sml <- suppressMessages(mtbls_parquet_to_sml(files$parquet, files$toml, mtd))
    expect_error(mtd_sml_to_mztabm(mtd, sml[, colnames(sml) != "abundance_assay[32]"]),
                 "32 assays, but the SML section has 31 abundance columns")
    ## the study variables must match as well
    no_sv <- sml[, !grepl("study_variable", colnames(sml))]
    expect_error(mtd_sml_to_mztabm(mtd, no_sv),
                 "3 study variables, but the SML section has 0")
    expect_error(mtd_sml_to_mztabm("not a matrix", sml))
    expect_error(mtd_sml_to_mztabm(mtd, "not a data frame"))
    expect_error(mtd_sml_to_mztabm(mtd[, 1L, drop = FALSE], sml))
})

test_that("mtbls_to_mztabm writes the mzTab-M file from the intermediate files", {
    f <- tempfile(fileext = ".mztab")
    msgs <- capture_messages(
        out <- mtbls_to_mztabm("MTBLS103", file = f, path = ext,
                               validate = FALSE))
    expect_match(msgs, "SML 13 rows", all = FALSE)
    expect_identical(as.character(out), f)
    expect_true(file.exists(f))
    lns <- readLines(f)
    expect_identical(lns[1L], "MTD\tmzTab-version\t2.1.0-M")
    expect_gt(sum(startsWith(lns, "MTD\t")), 100L)
    expect_identical(sum(startsWith(lns, "SMH\t")), 1L)
    expect_identical(sum(startsWith(lns, "SML\t")), 13L)
    expect_identical(sum(startsWith(lns, "SMF\t")), 0L)
    expect_true("MTD\tmzTab-profile\tM+S" %in% lns)
    expect_true(any(grepl("^MTD\tinstrument\\[1\\]-name\t\\[MS, MS:1000031, instrument model, Agilent 6550 iFunnel Q-TOF\\]$", lns)))
    ## the intermediate files are kept next to the mzTab-M file
    expect_identical(attr(out, "toml"), sub("mztab$", "toml", f))
    expect_identical(attr(out, "parquet"), sub("mztab$", "parquet", f))
    expect_true(file.exists(attr(out, "toml")))
    expect_true(file.exists(attr(out, "parquet")))
    expect_identical(.toml_read(attr(out, "toml"))$study$identifier, "MTBLS103")
    expect_identical(dim(nanoparquet::read_parquet(attr(out, "parquet"))),
                     c(13L, 53L))
    ## the file reads back
    back <- RmzTabM::readMzTabM(f, validate = FALSE)
    expect_s4_class(back, "MzTabM")
    expect_identical(nrow(RmzTabM::sml(back)), 13L)
    expect_s4_class(attr(out, "mztabm"), "MzTabM")
    expect_setequal(RmzTabM::mtd(attr(out, "mztabm"))[, 1L],
                    RmzTabM::mtd(back)[, 1L])
})

test_that("mtbls_to_mztabm builds the file from the saved intermediate files", {
    f <- tempfile(fileext = ".mztab")
    out <- mtbls_to_mztabm("MTBLS103", file = f, path = ext, validate = FALSE,
                           verbose = FALSE)
    ## the same file is obtained from the intermediate files alone ...
    mtd <- mtbls_toml_to_mtd(attr(out, "toml"))
    sml <- suppressMessages(mtbls_parquet_to_sml(attr(out, "parquet"),
                                                 attr(out, "toml"), mtd))
    f2 <- tempfile(fileext = ".mztab")
    RmzTabM::writeMzTabM(mtd_sml_to_mztabm(mtd, sml), f2, validate = FALSE)
    expect_identical(readLines(f2), readLines(f))
    ## ... and a change of the study description changes the metadata
    toml <- .toml_read(attr(out, "toml"))
    toml$study$title <- "A new title"
    .toml_write(toml, attr(out, "toml"))
    mtd2 <- mtbls_toml_to_mtd(attr(out, "toml"))
    expect_identical(mtd_value(mtd2, "title"), "A new title")
    expect_false(identical(mtd_value(mtd, "title"), "A new title"))
})

test_that("mtbls_to_mztabm names the files and passes its arguments on", {
    ## an extension is appended when missing
    f <- tempfile()
    out <- mtbls_to_mztabm("MTBLS103", file = f, path = ext, validate = FALSE,
                           verbose = FALSE)
    expect_identical(as.character(out), paste0(f, ".mztab"))
    expect_true(file.exists(out))
    expect_identical(attr(out, "toml"), paste0(f, ".toml"))
    ## intermediate files can be put elsewhere
    d <- file.path(tempfile(), "steps")
    out <- mtbls_to_mztabm("MTBLS103", file = tempfile(fileext = ".mzTab"),
                           path = ext, toml = file.path(d, "s.toml"),
                           parquet = file.path(d, "m.parquet"),
                           validate = FALSE, verbose = FALSE)
    expect_identical(sort(list.files(d)), c("m.parquet", "s.toml"))
    expect_identical(attr(out, "toml"), file.path(d, "s.toml"))
    ## assay, factor and reliability reach the steps
    f2 <- tempfile(fileext = ".mztab")
    suppressMessages(mtbls_to_mztabm("MTBLS103", file = f2, path = ext,
                                     assayName = hilic, factor = "HIAE",
                                     reliability = "2", validate = FALSE))
    lns <- readLines(f2)
    smh <- strsplit(lns[startsWith(lns, "SMH\t")], "\t")[[1L]]
    sml <- do.call(rbind, strsplit(lns[startsWith(lns, "SML\t")], "\t"))
    expect_true(all(sml[, match("reliability", smh)] == "2"))
    expect_error(mtbls_to_mztabm("MTBLS103", file = f2, path = ext,
                                 factor = "nope", validate = FALSE),
                 "not found")
    expect_error(mtbls_to_mztabm("MTBLS103", file = f2, path = ext,
                                 assayName = "a_none.txt", validate = FALSE),
                 "not found")
    expect_error(mtbls_to_mztabm(c("MTBLS103", "MTBLS40"), path = ext))
    expect_error(mtbls_to_mztabm("MTBLS103", file = c("a", "b"), path = ext))
    expect_error(mtbls_to_mztabm("MTBLS103", path = "/no/dir"),
                 "existing directory")
})

test_that("mtbls_to_mztabm reports what is written", {
    msgs <- capture_messages(
        mtbls_to_mztabm("MTBLS103", file = tempfile(fileext = ".mztab"),
                        path = ext, validate = FALSE))
    expect_match(msgs, "13 small molecules x 32 assays", all = FALSE)
    expect_match(msgs, "not available from the ISA-TAB files: contact \\(no ORCID\\)",
                 all = FALSE)
    ## the substitutions made are reported too
    expect_match(msgs, "No identification level", all = FALSE)
    ## verbose = FALSE silences everything
    expect_silent(
        out <- mtbls_to_mztabm("MTBLS103", file = tempfile(fileext = ".mztab"),
                               path = ext, validate = FALSE, verbose = FALSE))
    expect_true(file.exists(out))
    expect_silent(mtbls_to_mztabm("MTBLS40", file = tempfile(fileext = ".mztab"),
                                  path = ext40, validate = FALSE,
                                  verbose = FALSE))
})

test_that("the mzTab-M files of the bundled studies pass the official validator", {
    skip_if_no_validator()
    for (study in c("MTBLS103", "MTBLS40")) {
        path <- system.file("extdata", study, package = "mtbls2mztabm")
        ## validate = TRUE removes the file and stops if it is not valid
        out <- suppressMessages(
            mtbls_to_mztabm(study, file = tempfile(fileext = ".mztab"),
                            path = path, validate = TRUE))
        expect_true(file.exists(out))
        res <- suppressMessages(RmzTabM::mzTabMValidator(out))
        expect_equal(res$n_errors, 0, info = study)
        expect_equal(res$n_warnings, 0, info = study)
        expect_equal(res$total_messages, 0, info = study)
    }
})

test_that("files with unusual studies and MAFs pass the official validator", {
    skip_if_no_validator()
    files <- study_files("MTBLS103")
    ## a study description with replicates, partial organism information, a
    ## second factor named with a comma and runs without data file ...
    toml <- modified_toml(files$toml, function(x) {
        for (i in 1:3)
            x$sample[[i]][grepl("^organism", names(x$sample[[i]]))] <- NULL
        r <- x$run[[5L]]
        r$assay <- "replicate"
        r$data_file <- NULL
        x$run <- c(x$run, list(r))
        x$factor <- c(x$factor, list(list(name = "Dose, mg", levels = I(c("a", "b")))))
        for (i in seq_along(x$sample))
            x$sample[[i]]$factors[["Dose, mg"]] <- c("a", "b")[1 + i %% 2]
        x
    })
    ## ... and a MAF with several identifiers, placeholders, missing values
    ## and a retention time
    maf <- mtbls_maf("MTBLS103", path = ext)
    maf$database_identifier[1:4] <- c("CHEBI:18050|HMDB:HMDB0000641", "unknown",
                                      "KEGG:C00064", NA)
    maf$chemical_formula[1L] <- "C5H10N2O3|C5H10N2O3"
    maf$metabolite_identification[2L] <- "a|b"
    maf$retention_time <- seq(10, by = 1.5, length.out = nrow(maf))
    maf$reliability[5:6] <- c("2", "1")
    maf[3L, 22:24] <- NA
    pq <- tempfile(fileext = ".parquet")
    .write_maf_parquet(maf, pq, "MTBLS103", "x")
    mtd <- suppressMessages(mtbls_toml_to_mtd(toml, factor = c("HIAE", "Dose, mg")))
    sml <- suppressMessages(mtbls_parquet_to_sml(pq, toml, mtd))
    f <- tempfile(fileext = ".mztab")
    RmzTabM::writeMzTabM(mtd_sml_to_mztabm(mtd, sml), f, validate = FALSE)
    res <- suppressMessages(RmzTabM::mzTabMValidator(f))
    expect_equal(res$n_errors, 0)
    expect_equal(res$total_messages, 0)
})
