# Main Script to run all the sub Rmd files in order.
# the order is determined by the names of the files we give in order.txt


suppressPackageStartupMessages({
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required. Install with install.packages('rmarkdown').")
  }
})


args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) && nzchar(file_arg)) {
  script_dir <- dirname(normalizePath(file_arg))
  setwd(script_dir)
}

message("Working directory: ", getwd())

read_order <- function(path = "order.txt") {
  if (!file.exists(path)) stop("order.txt not found in current working directory.")
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x) & !startsWith(x, "#")]
  if (!length(x)) stop("order.txt is empty after removing blanks/comments.")
  missing <- x[!file.exists(x)]
  if (length(missing)) {
    msg <- paste0("These files from order.txt do not exist:\n- ", paste(missing, collapse = "\n- "),
                  "\n\nFiles present here:\n", paste(list.files(".", recursive = TRUE), collapse = "\n"))
    stop(msg)
  }
  normalizePath(x, winslash = "/", mustWork = TRUE)
}

render_one <- function(f) {
  message("\nRendering: ", f)
  rmarkdown::render(
    input = f,
    clean = TRUE,
    envir = new.env(parent = globalenv())
  )
}

main <- function() {
  rmds <- read_order("order.txt")
  message("Planned order:\n- ", paste(basename(rmds), collapse = "\n- "))
  for (f in rmds) render_one(f)
  message("\nAll done.")
}

# Only run main() if sourced as a script, not if imported
invisible(main())
