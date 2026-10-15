roundtrip <- function(x, comment = character()) {
    f <- tempfile(fileext = ".toml")
    .toml_write(x, f, comment = comment)
    .toml_read(f)
}

test_that("TOML strings are escaped and survive a round trip", {
    s <- c("plain", "with \"quotes\" and \\ backslash", "tab\tnewline\ncr\r",
           "unicode: μ-β ≥ 5 °C", "bell\a null-ish \001 \177",
           "[brackets], commas, = signs # hashes", "'single'", "")
    for (v in s) {
        x <- roundtrip(list(a = v))
        expect_identical(x$a, v)
    }
    expect_identical(roundtrip(list(a = s))$a, s)
    expect_identical(.toml_string("a\"b"), "\"a\\\"b\"")
})

test_that("TOML scalars, arrays, tables and arrays of tables are written", {
    x <- list(
        version = 1L, flag = TRUE, ratio = 1.5, whole = 2, big = 1e10,
        name = "x", names = c("a", "b", "c"), single = I("only"),
        none = character(), missing = NA_character_, nothing = NULL,
        study = list(id = "S1", authors = c("A", "B"),
                     sub = list(deep = "yes")),
        run = list(list(assay = "a1", n = 1L),
                   list(assay = "a2", factors = list(g = "x", h = "y"))),
        empty = list())
    y <- roundtrip(x)
    expect_identical(y$version, 1L)
    expect_true(y$flag)
    expect_identical(y$ratio, 1.5)
    expect_identical(y$whole, 2)
    expect_identical(y$big, 1e10)
    expect_identical(y$names, c("a", "b", "c"))
    expect_identical(y$single, "only")
    ## zero-length, missing and empty entries are not written
    expect_false(any(c("none", "missing", "nothing", "empty") %in% names(y)))
    expect_identical(y$study$authors, c("A", "B"))
    expect_identical(y$study$sub$deep, "yes")
    expect_length(y$run, 2L)
    expect_identical(y$run[[1L]], list(assay = "a1", n = 1L))
    expect_identical(y$run[[2L]]$factors, list(g = "x", h = "y"))
})

test_that("an array with one element is written as an array if marked with I()", {
    f <- tempfile(fileext = ".toml")
    .toml_write(list(a = "x", b = I("x"), c = I(c("x", "y"))), f)
    lns <- readLines(f)
    expect_true("a = \"x\"" %in% lns)
    expect_true("b = [\"x\"]" %in% lns)
    expect_true("c = [\"x\", \"y\"]" %in% lns)
})

test_that("long arrays are written over several lines", {
    v <- paste0("level number ", 1:30)
    f <- tempfile(fileext = ".toml")
    .toml_write(list(levels = I(v)), f)
    lns <- readLines(f)
    expect_gt(length(lns), 30L)
    expect_identical(lns[1L], "levels = [")
    expect_identical(.toml_read(f)$levels, v)
})

test_that("keys that are not bare keys are quoted", {
    x <- list(study = list(`Age (years)` = "53", `a.b` = "1", plain_key = "2"))
    f <- tempfile(fileext = ".toml")
    .toml_write(x, f)
    lns <- readLines(f)
    expect_true("\"Age (years)\" = \"53\"" %in% lns)
    expect_true("plain_key = \"2\"" %in% lns)
    expect_identical(.toml_read(f)$study[["Age (years)"]], "53")
    expect_identical(.toml_read(f)$study[["a.b"]], "1")
    ## and in table headers
    y <- roundtrip(list(`my table` = list(a = "x")))
    expect_identical(y[["my table"]]$a, "x")
})

test_that("comments are written at the top and the output directory is made", {
    f <- file.path(tempfile(), "sub", "x.toml")
    .toml_write(list(a = "x"), f, comment = c("first line", "second line"))
    lns <- readLines(f)
    expect_identical(lns[1:2], c("# first line", "# second line"))
    expect_identical(.toml_read(f)$a, "x")
})

test_that("floats keep their type and special values are written", {
    expect_identical(.toml_values(c(1, 1.5, -2, 1e-7)),
                     c("1.0", "1.5", "-2.0", "1e-07"))
    expect_identical(.toml_values(c(NaN, Inf, -Inf)), c("nan", "inf", "-inf"))
    expect_identical(.toml_values(3L), "3")
    expect_identical(.toml_values(c(TRUE, FALSE)), c("true", "false"))
    expect_identical(.toml_values(factor("a")), "\"a\"")
    x <- roundtrip(list(a = 1, b = 1e-7, c = c(0.1, 0.25)))
    expect_identical(x$a, 1)
    expect_identical(x$b, 1e-7)
    expect_identical(x$c, c(0.1, 0.25))
})

test_that("values that TOML cannot hold are rejected", {
    f <- tempfile(fileext = ".toml")
    expect_error(.toml_write(list(a = c("x", NA)), f), "Missing values")
    expect_error(.toml_values(Sys.Date()), "Cannot write")
})

test_that(".toml_read reports unreadable files", {
    expect_error(.toml_read(tempfile()), "does not exist")
    f <- tempfile(fileext = ".toml")
    writeLines("this is = = not toml", f)
    expect_error(.toml_read(f), "Could not parse")
})
