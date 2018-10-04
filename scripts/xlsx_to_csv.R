#!/usr/bin/env Rscript
library(readr)
library(readxl)

# disable scientific notation
# see https://github.com/tidyverse/readr/issues/671

format_numeric <- function(x, ...) {
    numeric_cols <- vapply(x, is.numeric, logical(1))
  x[numeric_cols] <- lapply(x[numeric_cols], format, ...)
    x
}

args = commandArgs(trailingOnly=TRUE)
in_file = args[1]

s = read_excel(in_file)
cat(format_csv(format_numeric(s)))
