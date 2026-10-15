test_that("mtbls_maf reads the MAF", {
    maf <- mtbls_maf("MTBLS103", path = ext)
    expect_identical(dim(maf), c(13L, 53L))
    expect_identical(maf$metabolite_identification[1L], "L-glutamine")
    expect_true(all(is.na(maf$mass_to_charge)))
    expect_type(maf$CON_BASA_567795, "double")
    expect_error(mtbls_maf("MTBLS103", path = "/no/dir"), "existing directory")
})

test_that("mtbls_maf_to_parquet saves the MAF", {
    f <- file.path(tempfile(), "maf.parquet")
    expect_silent(out <- mtbls_maf_to_parquet("MTBLS103", f, path = ext))
    expect_identical(out, f)
    expect_true(file.exists(f))
    back <- nanoparquet::read_parquet(f)
    maf <- mtbls_maf("MTBLS103", path = ext)
    expect_identical(dim(back), dim(maf))
    expect_identical(colnames(back), colnames(maf))
    expect_identical(back$metabolite_identification,
                     maf$metabolite_identification)
    expect_identical(back$database_identifier[c(1L, 8L)],
                     c("CHEBI:18050", "CHEBI:15354"))
    ## abundances are kept exactly
    expect_identical(back$CON_BASA_567795, maf$CON_BASA_567795)
    expect_equal(back$CON_BASA_567795[1L], 981143.812)
    ## provenance
    kv <- nanoparquet::read_parquet_metadata(f)$file_meta_data$key_value_metadata[[1L]]
    expect_identical(kv$value[kv$key == "mtbls_id"], "MTBLS103")
    expect_identical(
        kv$value[kv$key == "maf_file"],
        "m_MTBLS103_ibanez_02_HILIC_mass_spectrometry_v2_maf.tsv")
})

test_that("the Parquet file has the same column types for every study", {
    f <- study_files("MTBLS103")$parquet
    sch <- nanoparquet::read_parquet_schema(f)
    sch <- sch[sch$name != "schema", ]
    types <- setNames(sch$r_type, sch$name)
    ## columns that are completely empty in MTBLS103 keep their type
    for (nm in c("database_identifier", "smiles", "reliability", "uri",
                 "fragmentation", "taxid"))
        expect_identical(unname(types[nm]), "character", info = nm)
    for (nm in c("mass_to_charge", "retention_time", "charge",
                 "search_engine_score", "smallmolecule_abundance_sub"))
        expect_identical(unname(types[nm]), "double", info = nm)
    expect_true(all(types[-seq_len(21L)] == "double"))
    ## MTBLS40 reports retention times, but no m/z
    sch40 <- nanoparquet::read_parquet_schema(study_files("MTBLS40")$parquet)
    t40 <- setNames(sch40$r_type, sch40$name)
    expect_identical(unname(t40["mass_to_charge"]), "double")
    expect_identical(unname(t40["retention_time"]), "double")
    expect_identical(unname(t40["reliability"]), "character")
})

test_that(".maf_for_parquet fixes the types of standard columns only", {
    maf <- data.frame(
        database_identifier = c(NA, NA), reliability = c(1, 2),
        mass_to_charge = c("100.5", NA), retention_time = c("1", "n/a"),
        charge = c(NA, NA), S1 = c(NA, NA), S2 = c("1", "2"),
        S3 = c(TRUE, NA), S4 = 1:2, stringsAsFactors = FALSE)
    r <- .maf_for_parquet(maf)
    expect_identical(r$database_identifier, c(NA_character_, NA_character_))
    expect_identical(r$reliability, c("1", "2"))
    expect_identical(r$mass_to_charge, c(100.5, NA))
    ## non-numeric entries are not lost
    expect_identical(r$retention_time, c("1", "n/a"))
    expect_identical(r$charge, c(NA_real_, NA_real_))
    ## all-missing sample columns are numbers, other columns are kept
    expect_identical(r$S1, c(NA_real_, NA_real_))
    expect_identical(r$S2, c("1", "2"))
    expect_identical(r$S3, c("TRUE", NA))
    ## integers are numbers like all other abundances
    expect_identical(r$S4, c(1, 2))
    expect_identical(.maf_for_parquet(data.frame(retention_time = 1:2))$retention_time,
                     c(1, 2))
})

test_that("mtbls_maf_to_parquet selects the MAF through the assay or by name", {
    d <- tempfile()
    dir.create(d)
    file.copy(list.files(ext, full.names = TRUE), d)
    maf2 <- "m_other_maf.tsv"
    file.copy(file.path(ext, "m_MTBLS103_ibanez_02_HILIC_mass_spectrometry_v2_maf.tsv"),
              file.path(d, maf2))
    ## two MAFs: the one of the assay is saved ...
    f <- tempfile(fileext = ".parquet")
    expect_silent(mtbls_maf_to_parquet("MTBLS103", f, assayName = hilic,
                                       path = d))
    kv <- nanoparquet::read_parquet_metadata(f)$file_meta_data$key_value_metadata[[1L]]
    expect_match(kv$value[kv$key == "maf_file"], "HILIC_mass_spectrometry_v2")
    ## ... or the one named explicitly
    f2 <- tempfile(fileext = ".parquet")
    mtbls_maf_to_parquet("MTBLS103", f2, mafName = maf2, path = d)
    kv2 <- nanoparquet::read_parquet_metadata(f2)$file_meta_data$key_value_metadata[[1L]]
    expect_identical(kv2$value[kv2$key == "maf_file"], maf2)
    expect_error(mtbls_maf_to_parquet(c("a", "b"), tempfile(), path = ext))
})

test_that(".read_maf_parquet reports a missing file", {
    expect_error(.read_maf_parquet(tempfile()), "does not exist")
    x <- .read_maf_parquet(study_files("MTBLS103")$parquet)
    expect_identical(class(x), "data.frame")
    expect_identical(dim(x), c(13L, 53L))
})

test_that(".maf_columns finds the abundance column of each assay row", {
    maf_names <- c("metabolite_identification", "S1", "S2", "f2")
    isa <- data.frame(`Sample Name` = c("S1", "S2", "POOL"),
                      `Derived Spectral Data File` = c("x/f1.mzML", "x/f2.mzML",
                                                       "x/f3.mzML"),
                      check.names = FALSE, stringsAsFactors = FALSE)
    expect_message(cols <- .maf_columns(isa, maf_names), "1 of 3 assay rows")
    expect_identical(cols, c("S1", "S2", NA))
    ## the key that matches most rows wins: here the data file names
    isa$`Sample Name` <- c("x", "y", "z")
    maf_names <- c("metabolite_identification", "f1", "f2", "f3")
    expect_silent(cols <- .maf_columns(isa, maf_names))
    expect_identical(cols, c("f1", "f2", "f3"))
    ## raw files are used where the assay lists no derived files
    isa2 <- data.frame(`Sample Name` = c("x", "y"),
                       `Raw Spectral Data File` = c("FILES/r1.cdf", "FILES/r2.cdf"),
                       check.names = FALSE, stringsAsFactors = FALSE)
    expect_identical(.maf_columns(isa2, c("r1", "r2")), c("r1", "r2"))
    expect_identical(.maf_columns(isa2, c("r1.cdf", "r2.cdf")),
                     c("r1.cdf", "r2.cdf"))
    ## no match at all
    expect_warning(cols <- .maf_columns(isa2, c("a", "b")),
                   "matches an abundance column")
    expect_identical(cols, c(NA_character_, NA_character_))
    ## annotation columns are never abundance columns
    isa3 <- data.frame(`Sample Name` = "species", check.names = FALSE,
                       stringsAsFactors = FALSE)
    expect_warning(.maf_columns(isa3, c("species", "smiles")), "matches")
})

test_that(".maf_abundance builds the numeric abundance matrix", {
    maf <- data.frame(metabolite_identification = c("a", "b"),
                      S1 = c(1, 2), S2 = c("3", "n/d"), stringsAsFactors = FALSE)
    m <- .maf_abundance(maf, c("S1", NA, "S2", "S1"))
    expect_identical(dim(m), c(2L, 4L))
    expect_identical(m[, 1L], c(1, 2))
    expect_true(all(is.na(m[, 2L])))
    expect_identical(m[, 3L], c(3, NA))
    expect_identical(m[, 4L], c(1, 2))
    expect_error(.maf_abundance(maf, c(NA, "other")), "abundance column")
})

test_that("a single-row MAF gives a one-row matrix", {
    maf <- data.frame(metabolite_identification = "a", S1 = 1, S2 = 3,
                      stringsAsFactors = FALSE)
    expect_identical(dim(.maf_abundance(maf, c("S1", "S2"))), c(1L, 2L))
})
