#' @export
print.mercRel <- function(x, ...) {
  cat("\nCall:\n",
      paste(deparse(attr(x, "call")),
            sep = "\n",
            collapse = "\n"),
      "\n",
      sep = "")
  if (length(x$Uncorrected)) {
    cat("\nCoefficients Uncorrected Model:\n")
    print(x$Uncorrected)
  }
  if (length(x$Corrected)) {
    cat("\nCoefficients Corrected Model:\n")
    print(x$Corrected)
  }

  cat("\n")
  invisible(x)
}
