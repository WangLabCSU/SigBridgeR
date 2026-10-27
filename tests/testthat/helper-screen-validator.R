# Shared fixtures for the `Validate*()` test files.

# Satisfies every documented interface requirement of a screen function
compliant_screen <- function(
  sc_data,
  matched_bulk,
  phenotype,
  label_type,
  phenotype_class,
  verbose = TRUE,
  ...
) {
  if (verbose) {
    message("screening")
  }
  list(scRNA_data = sc_data)
}
