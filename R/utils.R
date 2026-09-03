## Internal helpers shared across the package. None of these touch the
## network; they operate on values already read from ISA-TAB files.

## Strip surrounding double quotes and whitespace.
.unquote <- function(x) {
    sub('^"(.*)"$', "\\1", trimws(x))
}

## Turn empty strings of the character columns of a data.frame into NA and
## strip quotes that `read.table(quote = "")` leaves in place.
.clean_table <- function(x) {
    x[] <- lapply(x, function(v) {
        if (is.character(v)) {
            v <- .unquote(v)
            v[!nzchar(v)] <- NA_character_
        }
        v
    })
    x
}

## Column key derived from an ISA-TAB name: "Factor Value[Age]" -> "age".
.key <- function(x) {
    gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", tolower(x)))
}

## A column of a data.frame as character, or NA if it does not exist.
.col_or_na <- function(x, nm) {
    if (nm %in% names(x)) as.character(x[[nm]])
    else rep(NA_character_, nrow(x))
}

## The qualifier column that follows column `nm` if it is named `qual`
## ("Unit", "Term Source REF", "Term Accession Number"), else NA.
.qualifier <- function(x, nm, qual, within = 3L) {
    i <- which(names(x) == nm)[1L]
    if (is.na(i)) return(rep(NA_character_, nrow(x)))
    idx <- seq.int(i + 1L, length.out = within)
    idx <- idx[idx <= ncol(x)]
    idx <- idx[names(x)[idx] == qual]
    ## stop at the next non-qualifier column
    for (j in seq.int(i + 1L, length.out = within)) {
        if (j > ncol(x) || !names(x)[j] %in%
            c("Unit", "Term Source REF", "Term Accession Number")) {
            idx <- idx[idx < j]
            break
        }
    }
    if (!length(idx)) return(rep(NA_character_, nrow(x)))
    as.character(x[[idx[1L]]])
}

## xcms-style feature identifiers: "FT" + zero-padded index.
.feature_ids <- function(n, prefix = "FT") {
    if (!n) return(character())
    digits <- max(ceiling(log10(n + 1L)), 1L)
    sprintf(paste0(prefix, "%0", digits, "d"), seq_len(n))
}

## mzTab-M controlled-vocabulary parameter "[LABEL, LABEL:ACC, name, ]" from
## an ontology term accession (usually an OBO PURL) and the term name. Falls
## back to "[,, name, ]" when no accession is available.
.cv_param <- function(accession, name, source = NA_character_) {
    n <- max(length(accession), length(name))
    accession <- rep_len(as.character(accession), n)
    name <- rep_len(as.character(name), n)
    source <- rep_len(as.character(source), n)
    res <- rep(NA_character_, n)
    has_name <- !is.na(name)
    has_acc <- !is.na(accession) & nzchar(accession)
    b <- basename(accession)
    label <- ifelse(grepl("_", b, fixed = TRUE), sub("_.*$", "", b),
                    ifelse(is.na(source), basename(dirname(accession)),
                           source))
    acc <- ifelse(grepl("_", b, fixed = TRUE), sub("_", ":", b, fixed = TRUE),
                  paste0(label, ":", b))
    ## accession given as "PREFIX:123" rather than a URL
    plain <- has_acc & grepl("^[A-Za-z]+:[A-Za-z0-9_]+$", accession)
    label[plain] <- sub(":.*$", "", accession[plain])
    acc[plain] <- accession[plain]
    both <- has_acc & has_name
    res[both] <- sprintf("[%s, %s, %s, ]", label[both], acc[both], name[both])
    res[!has_acc & has_name] <- sprintf("[,, %s, ]", name[!has_acc & has_name])
    res
}

## PSI-MS CV parameter for an instrument component: use the accession from
## the ISA-TAB file when there is one, otherwise the generic parent term with
## the free-text value in the value slot.
.ms_cv <- function(accession, value, parent_acc, parent_name) {
    accession <- accession[!is.na(accession)][1L]
    value <- value[!is.na(value)][1L]
    if (!is.na(accession) && grepl("MS_[0-9]+$", accession))
        return(sprintf("[MS, %s, %s, ]", sub("_", ":", basename(accession)),
                       if (is.na(value)) "" else value))
    sprintf("[MS, %s, %s, %s]", parent_acc, parent_name,
            if (is.na(value)) "" else value)
}

## Scan polarity normalised to the two values RmzTabM accepts ("positive",
## "negative"), which it turns into the PSI-MS CV parameters itself.
.polarity_plain <- function(x) {
    l <- tolower(trimws(as.character(x)))
    res <- rep(NA_character_, length(x))
    res[!is.na(l) & grepl("^pos", l)] <- "positive"
    res[!is.na(l) & grepl("^neg", l)] <- "negative"
    if (anyNA(res) && any(!is.na(x)))
        stop("Unsupported scan polarity value(s): ",
             paste(unique(x[is.na(res) & !is.na(x)]), collapse = ", "),
             ". mzTab-M export needs \"positive\" or \"negative\".")
    res
}

## Data file format as a PSI-MS CV parameter, from the file extension.
.file_format_cv <- function(f) {
    ext <- tolower(file_ext(as.character(f)))
    res <- rep(NA_character_, length(f))
    res[ext == "mzml"] <- "[MS, MS:1000584, mzML format, ]"
    res[ext == "mzxml"] <- "[MS, MS:1000566, ISB mzXML format, ]"
    res[ext == "mzdata"] <- "[MS, MS:1000564, PSI mzData format, ]"
    res[ext %in% c("cdf", "nc")] <- "[MS, MS:1000825, Andi-MS format, ]"
    res
}

## Native spectrum identifier format for a data file format, where it is fixed
## by the format itself (mzXML and Andi-MS/netCDF: scan numbers). NA where it
## depends on the vendor (mzML) or is unknown.
.id_format_cv <- function(f) {
    ext <- tolower(file_ext(as.character(f)))
    res <- rep(NA_character_, length(f))
    res[ext %in% c("mzxml", "cdf", "nc")] <-
        "[MS, MS:1000776, scan number only nativeID format, ]"
    res
}

## Numeric coercion without warnings: anything that is not a number is NA.
.as_numeric <- function(v) {
    v <- trimws(as.character(v))
    pat <- "^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$"
    ok <- !is.na(v) & grepl(pat, v)
    out <- rep(NA_real_, length(v))
    out[ok] <- as.numeric(v[ok])
    out
}
