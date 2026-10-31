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

test_that(".maf_columns falls back to corresponding columns with numbers", {
    runs <- paste0("D-0", 1:3, "-GC-01")
    isa <- data.frame(`Sample Name` = runs, `MS Assay Name` = runs,
                      check.names = FALSE, stringsAsFactors = FALSE)
    short <- paste0("D-0", 1:3)
    ## the matched columns are empty, a suffix leads to the columns with values
    expect_message(cols <- .maf_columns(isa, c("species", runs, short), short),
                   "'D-01' for the run 'D-01-GC-01'")
    expect_identical(cols, short)
    ## also if no column is matched at all
    expect_message(cols <- .maf_columns(isa, c("species", short)), "one-to-one")
    expect_identical(cols, short)
    ## the suffix can be on the name of the column
    isa2 <- data.frame(`Sample Name` = short, check.names = FALSE,
                       stringsAsFactors = FALSE)
    expect_message(cols <- .maf_columns(isa2, c(short, runs), runs),
                   "one-to-one")
    expect_identical(cols, runs)
    ## the other names of the runs are tried first
    isa3 <- data.frame(`Sample Name` = short, `MS Assay Name` = c("a", "b", "c"),
                       check.names = FALSE, stringsAsFactors = FALSE)
    expect_message(cols <- .maf_columns(isa3, c(short, letters[1:3]),
                                        letters[1:3]), "one-to-one")
    expect_identical(cols, letters[1:3])
    ## nothing changes if a matched column holds numbers
    expect_silent(cols <- .maf_columns(isa, c(runs, short), c(runs[1L], short)))
    expect_identical(cols, runs)
    ## a separator is needed: "D-0" is not the sample of "D-01-GC-01"
    isa4 <- data.frame(`Sample Name` = "D-01", check.names = FALSE,
                       stringsAsFactors = FALSE)
    expect_warning(cols <- .maf_columns(isa4, "D-0"), "matches")
    expect_identical(cols, NA_character_)
})

test_that(".maf_columns does not guess if the mapping is ambiguous", {
    runs <- paste0("D-0", 1:3, "-GC-01")
    isa <- data.frame(`Sample Name` = runs, check.names = FALSE,
                      stringsAsFactors = FALSE)
    empty <- "holds a number"
    ## two columns for a run
    filled <- c("D-01", "D-01-GC", "D-02", "D-03")
    expect_warning(cols <- .maf_columns(isa, c(runs, filled), filled), empty)
    expect_identical(cols, runs)
    ## one column for several runs
    isa2 <- isa
    isa2$`Sample Name` <- c("D-01-a", "D-01-b", "D-02-a")
    filled <- c("D-01", "D-02")
    expect_warning(cols <- .maf_columns(isa2, c(isa2$`Sample Name`, filled),
                                        filled), empty)
    expect_identical(cols, isa2$`Sample Name`)
    ## not all runs have a column
    filled <- c("D-01", "D-02")
    expect_warning(cols <- .maf_columns(isa, c(runs, filled), filled), empty)
    expect_identical(cols, runs)
    ## a repeated name: which of the columns belongs to the run is not known
    filled <- c("D-01", "D-02", "D-03")
    expect_warning(
        expect_message(cols <- .maf_columns(isa, c(runs, filled, "D-03"),
                                            filled), "repeats 1"), empty)
    expect_identical(cols, runs)
    ## the names of the runs lead to different columns
    isa3 <- data.frame(`Sample Name` = paste0("S", 1:3),
                       `MS Assay Name` = paste0("M", 1:3),
                       check.names = FALSE, stringsAsFactors = FALSE)
    filled <- c("S1-x", "S2-x", "S3-x", "M1-y", "M2-y", "M3-y")
    expect_warning(cols <- .maf_columns(isa3, c(paste0("S", 1:3), filled),
                                        filled), empty)
    expect_identical(cols, paste0("S", 1:3))
    ## no columns with numbers at all, as for an empty deposit
    expect_warning(cols <- .maf_columns(isa, runs, character()), empty)
    expect_identical(cols, runs)
})

test_that(".maf_columns warns about columns used by several runs", {
    isa <- data.frame(`Sample Name` = c("B1", "B2", "S1", "S2"),
                      `MS Assay Name` = c("B", "B", "S1", "S2"),
                      check.names = FALSE, stringsAsFactors = FALSE)
    expect_warning(cols <- .maf_columns(isa, c("B", "S1", "S2")),
                   "2 assay rows share 1 abundance columns of the MAF \\(B\\)")
    expect_identical(cols, c("B", "B", "S1", "S2"))
})

test_that(".maf_filled finds the columns with numbers", {
    maf <- data.frame(species = c("a", "b"), S1 = c(NA, NA), S2 = c("1,234,567", NA),
                      S3 = c("n/d", NA), S4 = c(0, NA), stringsAsFactors = FALSE)
    expect_identical(.maf_filled(maf), c("S2", "S4"))
    ## the first column of a repeated name is the one that is used
    names(maf)[2:3] <- "S2"
    expect_identical(.maf_filled(maf), "S4")
})

test_that("thousands separators in abundances are recognized", {
    num <- function(...) {
        .maf_numbers(data.frame(mass_to_charge = "1,5", S1 = c(...),
                                stringsAsFactors = FALSE))
    }
    r <- num("600,461", "1,157,497", "512", "2,000.5", NA, "-1,000,000", "n/d")
    expect_identical(r$S1, c(600461, 1157497, 512, 2000.5, NA, -1e6, NA))
    expect_identical(attr(r, "ambiguous"), 0L)
    ## each of these may be a number with a decimal comma
    for (x in list(c("1,234", "12,345"), c("1,157,497", "1,5"),
                   c("1,157,497", "0,123"), c("1,157,497", "1234,567"))) {
        r <- num(x)
        expect_identical(r$S1[grepl(",", x)], rep(NA_real_, sum(grepl(",", x))),
                         info = x)
        expect_identical(attr(r, "ambiguous"), length(x), info = x)
    }
    maf <- data.frame(S1 = c("1,234", "12,345"), S2 = c(1, 2),
                      stringsAsFactors = FALSE)
    expect_warning(m <- .maf_abundance(maf, c("S1", "S2")), "2 abundances")
    expect_identical(m, cbind(c(NA_real_, NA_real_), c(1, 2)))
    ## other columns than abundances stay as they are
    expect_identical(.as_numeric("1,234"), NA_real_)
})

test_that("abundances of empty matched columns are taken from their samples", {
    ## the layout of MTBLS135: the assay names the runs "<sample>-GC-01", the
    ## columns of this name are empty and the values, with thousands
    ## separators, are in the columns named after the sample
    assay <- hilic
    abundant <- function(maf) {
        s <- setdiff(names(maf), .maf_annotation_columns)
        maf[s] <- lapply(maf[s], function(v)
            ifelse(is.na(v), NA_character_, formatC(v, format = "f", digits = 4L,
                                                    big.mark = ",")))
        maf[paste0(s, "-GC-01")] <- NA_real_
        maf
    }
    d <- modified_study("MTBLS103", function(x) paste0(x, "-GC-01"), abundant)
    expect_message(
        toml <- mtbls_factors_to_toml("MTBLS103", tempfile(fileext = ".toml"),
                                      assayName = assay, path = d),
        "one-to-one")
    x <- .toml_read(toml)
    expect_identical(x$run[[1L]]$sample, "CON_BASA_567795-GC-01")
    expect_identical(x$run[[1L]]$maf_column, "CON_BASA_567795")
    ## the abundances are the ones of the MAF, in the order of the runs
    maf <- mtbls_maf("MTBLS103", path = ext)
    parquet <- mtbls_maf_to_parquet("MTBLS103", tempfile(fileext = ".parquet"),
                                    assayName = assay, path = d)
    sml <- suppressMessages(mtbls_parquet_to_sml(parquet, toml))
    runs <- vapply(x$run, `[[`, "", "maf_column")
    expect_length(runs, 32L)
    expect_false(anyDuplicated(runs) > 0L)
    for (i in seq_along(runs))
        expect_equal(sml[[paste0("abundance_assay[", i, "]")]],
                     maf[[runs[i]]], tolerance = 1e-6, info = runs[i])
})

test_that("runs that the MAF does not tell apart stay without abundances", {
    ## the layout of MTBLS737: the MAF repeats a column name, which is no reason
    ## to give a run without a column of its own one of these columns
    samples <- function(x) sub("^CON_BASA_567795$", "CON_BASA_567795.1", x)
    twice <- function(maf) {
        names(maf)[names(maf) == "CON_BASA_574488"] <- "CON_BASA_567795"
        maf
    }
    d <- modified_study("MTBLS103", samples, twice)
    expect_message(
        toml <- mtbls_factors_to_toml("MTBLS103", tempfile(fileext = ".toml"),
                                      assayName = hilic, path = d),
        "repeats 1 abundance column names \\(CON_BASA_567795\\)")
    cols <- vapply(.toml_read(toml)$run, function(r)
        if (is.null(r$maf_column)) NA_character_ else r$maf_column, "")
    expect_length(cols, 32L)
    expect_identical(sum(is.na(cols)), 2L)
    expect_true(is.na(cols[1L]))
    expect_false(anyDuplicated(cols[!is.na(cols)]) > 0L)
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
