test_that("mtbls_parquet_to_sml builds the small molecule section of MTBLS103", {
    files <- study_files("MTBLS103")
    expect_message(sml <- mtbls_parquet_to_sml(files$parquet, files$toml),
                   "reported as level 4")
    expect_s3_class(sml, "data.frame")
    expect_identical(dim(sml), c(13L, 52L))
    expect_identical(
        colnames(sml)[1:14],
        c("SMH", "SML_ID", "SMF_ID_REFS", "database_identifier",
          "chemical_formula", "smiles", "inchi", "chemical_name", "uri",
          "theoretical_neutral_mass", "adduct_ions", "reliability",
          "best_id_confidence_measure", "best_id_confidence_value"))
    expect_identical(
        colnames(sml)[15:ncol(sml)],
        c(paste0("abundance_assay[", 1:32, "]"),
          paste0("abundance_study_variable[", 1:3, "]"),
          paste0("abundance_variation_study_variable[", 1:3, "]")))
    expect_true(all(sml$SMH == "SML"))
    expect_identical(sml$SML_ID, 1:13)
    ## annotations of the MAF
    expect_identical(sml$database_identifier[c(1L, 8L)],
                     c("CHEBI:18050", "CHEBI:15354"))
    expect_identical(sml$chemical_name[c(1L, 8L)], c("L-glutamine", "choline"))
    expect_identical(sml$chemical_formula[1L], "C5H10N2O3")
    expect_match(sml$inchi[1L], "^InChI=1S/C5H10N2O3")
    expect_identical(sml$smiles[1L], "C(CC(=O)N)C(C(=O)O)N")
    ## what the MAF does not provide is null
    expect_true(all(sml$uri == "null"))
    expect_true(all(sml$theoretical_neutral_mass == "null"))
    expect_true(all(sml$adduct_ions == "null"))
    expect_true(all(sml$SMF_ID_REFS == "null"))
    expect_false(any(grepl("^opt_", colnames(sml))))
    ## reliability is mandatory: the lowest level where the MAF has none
    expect_true(all(sml$reliability == "4"))
})

test_that("the abundances of the SML are those of the MAF", {
    files <- study_files("MTBLS103")
    sml <- suppressMessages(mtbls_parquet_to_sml(files$parquet, files$toml))
    maf <- mtbls_maf("MTBLS103", path = ext)
    expect_equal(sml[["abundance_assay[1]"]], maf$CON_BASA_567795)
    col32 <- .toml_read(files$toml)$run[[32L]]$maf_column
    expect_equal(sml[["abundance_assay[32]"]], maf[[col32]])
    expect_equal(sml[1L, "abundance_assay[1]"], 981143.812)
    ## study variables: mean and coefficient of variation per sample group
    grp <- .toml_read(files$toml)$sample
    grp <- vapply(grp, function(s) s$factors$HIAE, "")
    ab <- as.matrix(sml[, paste0("abundance_assay[", 1:32, "]")])
    expect_equal(sml[["abundance_study_variable[1]"]],
                 unname(rowMeans(ab[, grp == "CTR"])))
    expect_equal(sml[["abundance_study_variable[3]"]],
                 unname(rowMeans(ab[, grp == "TREAT"])))
    cv <- apply(ab[, grp == "HIAE"], 1L, function(v) sd(v) / mean(v))
    expect_equal(sml[["abundance_variation_study_variable[2]"]], unname(cv))
})

test_that("MTBLS40 gives optional columns for what the SML cannot hold", {
    files <- study_files("MTBLS40")
    expect_message(
        expect_message(sml <- mtbls_parquet_to_sml(files$parquet, files$toml),
                       "kept in optional columns"),
        "reported as level 4")
    expect_identical(dim(sml), c(171L, 14L + 51L + 8L + 2L))
    ## the retention time of the MAF
    expect_identical(sml$opt_global_retention_time[1:3], c("1046", "1047", "1066"))
    expect_false("opt_global_mass_to_charge" %in% colnames(sml))
    ## the label of molecules without identifier
    expect_identical(sml$database_identifier[1:2], c("null", "CHEBI:24996"))
    expect_identical(sml$chemical_name[1:2], c("null", "Lactate"))
    expect_identical(sml$opt_global_metabolite_identification[1:2],
                     c("PR-MST", "null"))
    expect_identical(sum(sml$database_identifier == "null"),
                     sum(sml$opt_global_metabolite_identification != "null"))
    ## missing abundances are NA (written as null)
    expect_true(anyNA(sml[, paste0("abundance_assay[", 1:51, "]")]))
    expect_false(anyNA(sml$SML_ID))
})

test_that(".sml_annotations repeats single values and rejects mismatches", {
    maf <- data.frame(
        database_identifier = c("CHEBI:1", "CHEBI:1|CHEBI:2", "CHEBI:1|CHEBI:2",
                                "CHEBI:1|CHEBI:2|CHEBI:3", "unknown", "N/A",
                                NA, "HMDB0000123"),
        chemical_formula = c("C1", "C1", "C1|C2", "C1|C2", "C9", NA, NA, ""),
        metabolite_identification = c("a", "b", "c|d", "e", "f", "g|h", "i",
                                      "j|k"),  # "|" is no separator with one identifier
        smiles = NA_character_, inchi = NA_character_, reliability = "1",
        stringsAsFactors = FALSE)
    expect_message(a <- .sml_annotations(maf), "kept in optional columns")
    expect_identical(a$database_identifier,
                     c("CHEBI:1", "CHEBI:1|CHEBI:2", "CHEBI:1|CHEBI:2",
                       "CHEBI:1|CHEBI:2|CHEBI:3", NA, NA, NA, "HMDB0000123"))
    ## a single formula is repeated for each identifier, as is a single name
    expect_identical(a$chemical_formula,
                     c("C1", "C1|C1", "C1|C2", NA, NA, NA, NA, NA))
    expect_identical(a$chemical_name,
                     c("a", "b|b", "c|d", "e|e|e", NA, NA, NA, "j/k"))
    ## values that cannot be reported are kept in an optional column
    expect_identical(a$opt_global_chemical_formula,
                     c(NA, NA, NA, "C1|C2", "C9", NA, NA, NA))
    expect_identical(a$opt_global_metabolite_identification,
                     c(NA, NA, NA, NA, "f", "g|h", "i", NA))
    expect_false("opt_global_smiles" %in% names(a))
    expect_false("opt_global_inchi" %in% names(a))
})

test_that(".sml_annotations handles reliability and optional columns", {
    maf <- data.frame(
        database_identifier = "CHEBI:1",
        reliability = c("0", "1", "4", "5", "2a", " 3 ", NA, "x"),
        retention_time = c(1.5, NA, 3, 4, 5, 6, 7, 8),
        mass_to_charge = NA_real_, stringsAsFactors = FALSE)
    ## levels 0 to 4 are kept, anything else gets the given level
    expect_message(a <- .sml_annotations(maf, "2"), "4 of 8 molecules")
    expect_identical(a$reliability, c("0", "1", "4", "2", "2", "3", "2", "2"))
    expect_identical(a$opt_global_retention_time, c(1.5, NA, 3, 4, 5, 6, 7, 8))
    expect_false("opt_global_mass_to_charge" %in% names(a))
    maf$mass_to_charge <- c(100, rep(NA, 7))
    a <- suppressMessages(.sml_annotations(maf))
    expect_identical(a$opt_global_mass_to_charge, c(100, rep(NA, 7)))
    ## no message if the MAF gives all levels
    maf$reliability <- "1"
    expect_silent(.sml_annotations(maf))
    ## a MAF without annotation columns
    expect_message(a <- .sml_annotations(data.frame(S1 = 1:2)), "level 4")
    expect_identical(a$database_identifier, c(NA_character_, NA_character_))
    expect_identical(a$reliability, c("4", "4"))
})

test_that("a molecule with several identifiers gives a valid null mass", {
    files <- study_files("MTBLS103")
    maf <- mtbls_maf("MTBLS103", path = ext)
    maf$database_identifier[1L] <- "CHEBI:18050|CHEBI:1"
    maf$chemical_formula[1L] <- "C5H10N2O3|C1"
    pq <- tempfile(fileext = ".parquet")
    .write_maf_parquet(maf, pq, "MTBLS103", "x")
    sml <- suppressMessages(mtbls_parquet_to_sml(pq, files$toml))
    expect_identical(sml$database_identifier[1L], "CHEBI:18050|CHEBI:1")
    expect_identical(sml$chemical_formula[1L], "C5H10N2O3|C1")
    expect_identical(sml$chemical_name[1L], "L-glutamine|L-glutamine")
    ## "null|null" is not a list of numbers
    expect_identical(sml$theoretical_neutral_mass[1L], "null")
})

test_that("mtbls_parquet_to_sml takes the study variables from the metadata", {
    files <- study_files("MTBLS103")
    toml <- modified_toml(files$toml, function(x) {
        x$factor <- c(x$factor, list(list(name = "Sex", levels = I(c("f", "m")))))
        for (i in seq_along(x$sample))
            x$sample[[i]]$factors$Sex <- c("f", "m")[1 + i %% 2]
        x
    })
    sml1 <- suppressMessages(mtbls_parquet_to_sml(files$parquet, toml))
    sml2 <- suppressMessages(mtbls_parquet_to_sml(files$parquet, toml,
                                                  factor = c("HIAE", "Sex")))
    expect_identical(sum(grepl("^abundance_study_variable", colnames(sml1))), 3L)
    expect_identical(sum(grepl("^abundance_study_variable", colnames(sml2))), 5L)
    mtd <- mtbls_toml_to_mtd(toml, factor = "Sex")
    sml3 <- suppressMessages(mtbls_parquet_to_sml(files$parquet, toml, mtd))
    expect_identical(sum(grepl("^abundance_study_variable", colnames(sml3))), 2L)
    ## an MTD of another description does not fit
    other <- mtbls_toml_to_mtd(study_files("MTBLS40")$toml)
    expect_error(mtbls_parquet_to_sml(files$parquet, files$toml, other),
                 "metadata section has 51 assays, but .* lists 32 runs")
})

test_that("mtbls_parquet_to_sml validates its input", {
    files <- study_files("MTBLS103")
    expect_error(mtbls_parquet_to_sml(files$parquet, files$toml,
                                      reliability = "5"), "'reliability'")
    expect_error(mtbls_parquet_to_sml(files$parquet, files$toml,
                                      reliability = c("1", "2")),
                 "'reliability'")
    expect_error(mtbls_parquet_to_sml(tempfile(), files$toml), "does not exist")
    expect_error(mtbls_parquet_to_sml(files$parquet, tempfile()),
                 "does not exist")
    ## no rows
    empty <- tempfile(fileext = ".parquet")
    nanoparquet::write_parquet(mtbls_maf("MTBLS103", path = ext)[0L, ], empty)
    expect_error(mtbls_parquet_to_sml(empty, files$toml), "has no rows")
    ## the MAF does not belong to the study description
    toml <- modified_toml(files$toml, function(x) {
        for (i in seq_along(x$run)) x$run[[i]]$maf_column <- NULL
        x
    })
    expect_error(suppressMessages(mtbls_parquet_to_sml(files$parquet, toml)),
                 "abundance column")
})

test_that("runs without a MAF column have missing abundances", {
    files <- study_files("MTBLS103")
    toml <- modified_toml(files$toml, function(x) {
        x$run[[3L]]$maf_column <- NULL
        x
    })
    sml <- suppressMessages(mtbls_parquet_to_sml(files$parquet, toml))
    expect_true(all(is.na(sml[["abundance_assay[3]"]])))
    expect_false(anyNA(sml[["abundance_assay[4]"]]))
    ## the mean and variation of the group of the run come from its other
    ## assays
    grp <- vapply(.toml_read(toml)$sample, function(s) s$factors$HIAE, "")
    ab <- as.matrix(sml[, paste0("abundance_assay[", 1:32, "]")])
    expect_false(anyNA(sml[["abundance_study_variable[1]"]]))
    expect_equal(sml[["abundance_study_variable[1]"]],
                 unname(rowMeans(ab[, grp == "CTR"], na.rm = TRUE)))
    cv <- apply(ab[, grp == "CTR"], 1L, function(v)
        sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE))
    expect_equal(sml[["abundance_variation_study_variable[1]"]], unname(cv))
    expect_false(anyNA(sml[["abundance_study_variable[2]"]]))
})

test_that("a study variable without any abundance is missing", {
    files <- study_files("MTBLS103")
    toml <- modified_toml(files$toml, function(x) {
        for (i in seq_along(x$run)) {
            s <- x$sample[[match(x$run[[i]]$sample,
                                 vapply(x$sample, `[[`, "", "name"))]]
            if (s$factors$HIAE == "CTR")
                x$run[[i]]$maf_column <- NULL
        }
        x
    })
    sml <- suppressMessages(mtbls_parquet_to_sml(files$parquet, toml))
    expect_type(sml[["abundance_study_variable[1]"]], "double")
    expect_true(all(is.na(sml[["abundance_study_variable[1]"]])))
    expect_false(any(is.nan(sml[["abundance_study_variable[1]"]])))
    expect_false(anyNA(sml[["abundance_study_variable[2]"]]))
})

test_that("only the average function mean ignores missing abundances", {
    files <- study_files("MTBLS103")
    toml <- modified_toml(files$toml, function(x) {
        x$run[[3L]]$maf_column <- NULL
        x
    })
    mtd <- mtbls_toml_to_mtd(toml)
    i <- mtd[, 1L] == "study_variable[1]-average_function"
    mtd[i, 2L] <- "[MS, MS:1002883, median, ]"
    sml <- suppressMessages(mtbls_parquet_to_sml(files$parquet, toml, mtd = mtd))
    expect_true(all(is.na(sml[["abundance_study_variable[1]"]])))
})
