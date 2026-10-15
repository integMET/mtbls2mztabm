## Reading and writing of the TOML study description. The files are parsed
## with RcppTOML. RcppTOML cannot write TOML and the serializers available on
## CRAN either need a Rust toolchain to build or lose arrays of tables, so
## the few constructs of the study description are written here directly:
## strings, integers, logicals, arrays of strings, tables and arrays of
## tables (see `.toml_write()` for the accepted R structure).

## Basic string (double quotes) with all characters TOML requires escaped.
.toml_string <- function(x) {
    x <- enc2utf8(as.character(x))
    x <- gsub("\\", "\\\\", x, fixed = TRUE)
    x <- gsub("\"", "\\\"", x, fixed = TRUE)
    if (any(grepl("[[:cntrl:]]", x))) {
        x <- gsub("\b", "\\b", x, fixed = TRUE)
        x <- gsub("\t", "\\t", x, fixed = TRUE)
        x <- gsub("\n", "\\n", x, fixed = TRUE)
        x <- gsub("\f", "\\f", x, fixed = TRUE)
        x <- gsub("\r", "\\r", x, fixed = TRUE)
        ## remaining control characters
        for (cp in c(1:7, 11L, 14:31, 127L))
            x <- gsub(intToUtf8(cp), sprintf("\\u%04X", cp), x, fixed = TRUE)
    }
    paste0("\"", x, "\"")
}

## Bare key if possible, quoted key otherwise.
.toml_key <- function(x) {
    bare <- grepl("^[A-Za-z0-9_-]+$", x)
    x[!bare] <- .toml_string(x[!bare])
    x
}

## TOML representation of the elements of an atomic vector (no NA).
.toml_values <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) .toml_string(x)
    else if (is.logical(x)) ifelse(x, "true", "false")
    else if (is.integer(x)) as.character(x)
    else if (is.numeric(x))
        vapply(x, function(v) {
            if (is.nan(v)) "nan"
            else if (is.infinite(v)) if (v > 0) "inf" else "-inf"
            else {
                s <- format(v, digits = 15L, trim = TRUE)
                ## a float needs a fractional part or an exponent
                if (grepl("^-?[0-9]+$", s)) paste0(s, ".0") else s
            }
        }, character(1))
    else stop("Cannot write a value of class '", class(x)[1L], "' to TOML")
}

## `key = value` line(s) for an atomic vector. A vector of length one is
## written as a scalar unless it is marked with I() (an array of length
## one); longer vectors are arrays, split over several lines if long.
.toml_assignment <- function(key, x) {
    vals <- .toml_values(x)
    if (length(x) == 1L && !inherits(x, "AsIs"))
        return(paste(key, "=", vals))
    one <- paste0(key, " = [", paste(vals, collapse = ", "), "]")
    if (nchar(one, type = "width") <= 100L)
        return(one)
    c(paste(key, "= ["), paste0("  ", vals, ","), "]")
}

.toml_is_atomic <- function(x) is.atomic(x) && !is.null(x)

## An array of tables: an unnamed list whose elements are named lists.
.toml_is_aot <- function(x) {
    is.list(x) && is.null(names(x)) && length(x) &&
        all(vapply(x, function(e) is.list(e) && !is.null(names(e)), logical(1)))
}

## Lines of the table `x` (a named list) at the dotted key path `path`.
.toml_table_lines <- function(x, path = character()) {
    stopifnot(is.list(x), !is.null(names(x)) || !length(x))
    ## NULL, NA and zero-length entries carry no value and are left out
    keep <- vapply(x, function(e) {
        if (.toml_is_atomic(e)) length(e) > 0L &&
            !(length(e) == 1L && is.na(e))
        else length(e) > 0L
    }, logical(1))
    x <- x[keep]
    atomic <- vapply(x, .toml_is_atomic, logical(1))
    aot <- vapply(x, .toml_is_aot, logical(1))
    for (e in x[atomic])
        if (anyNA(e)) stop("Missing values cannot be written to TOML arrays")
    lns <- character()
    for (nm in names(x)[atomic])
        lns <- c(lns, .toml_assignment(.toml_key(nm), x[[nm]]))
    for (nm in names(x)[!atomic & !aot]) {
        p <- c(path, .toml_key(nm))
        lns <- c(lns, "", paste0("[", paste(p, collapse = "."), "]"),
                 .toml_table_lines(x[[nm]], p))
    }
    for (nm in names(x)[aot]) {
        p <- c(path, .toml_key(nm))
        header <- paste0("[[", paste(p, collapse = "."), "]]")
        lns <- c(lns, unlist(lapply(x[[nm]], function(el)
            c("", header, .toml_table_lines(el, p))), use.names = FALSE))
    }
    lns
}

## Write the named list `x` to `file` as TOML (UTF-8). Entries are written
## as follows: atomic vectors as scalars or arrays, named lists as tables and
## unnamed lists of named lists as arrays of tables. `comment` lines are put
## at the top of the file.
.toml_write <- function(x, file, comment = character()) {
    lns <- c(if (length(comment)) c(paste("#", comment), ""),
             .toml_table_lines(x))
    .ensure_dir(file)
    con <- file(file, open = "wb")
    on.exit(close(con))
    writeLines(enc2utf8(lns), con, useBytes = TRUE)
    invisible(file)
}

## Parse the TOML file `file` into a named list.
.toml_read <- function(file) {
    if (!is.character(file) || length(file) != 1L || !file.exists(file))
        stop("TOML file '", paste(file, collapse = ", "), "' does not exist")
    txt <- paste(readLines(file, encoding = "UTF-8", warn = FALSE),
                 collapse = "\n")
    res <- tryCatch(parseTOML(txt, fromFile = FALSE, escape = FALSE),
                    error = function(e)
                        stop("Could not parse the TOML file '", file, "': ",
                             conditionMessage(e), call. = FALSE))
    unclass(res)
}
