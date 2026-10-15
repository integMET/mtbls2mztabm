## Intermediate files of the bundled studies, written once per test run.
.study_files_cache <- new.env()

study_files <- function(study = c("MTBLS103", "MTBLS40")) {
    study <- match.arg(study)
    if (is.null(.study_files_cache[[study]])) {
        dir <- tempfile(paste0(study, "-"))
        dir.create(dir)
        path <- system.file("extdata", study, package = "mtbls2mztabm")
        .study_files_cache[[study]] <- suppressMessages(list(
            toml = mtbls_factors_to_toml(
                study, file.path(dir, paste0(study, ".toml")), path = path),
            parquet = mtbls_maf_to_parquet(
                study, file.path(dir, paste0(study, ".parquet")),
                path = path)))
    }
    .study_files_cache[[study]]
}

## Value(s) of a metadata field of an MTD matrix.
mtd_value <- function(mtd, key) unname(mtd[mtd[, 1L] == key, 2L])

## Read a TOML study description, change it and write it to a new file.
modified_toml <- function(toml, fun) {
    x <- .toml_read(toml)
    f <- tempfile(fileext = ".toml")
    .toml_write(fun(x), f)
    f
}

## The official mzTab-M validator is downloaded by RmzTabM on first use: tests
## using it only run where it is available (or can be fetched).
skip_if_no_validator <- function() {
    skip_on_cran()
    path <- tryCatch(RmzTabM::mzTabMValidator_download(),
                     error = function(e) NULL)
    skip_if(is.null(path) || !file.exists(path),
            "the mzTab-M validator is not available")
}
