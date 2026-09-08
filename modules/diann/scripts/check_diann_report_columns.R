# Verify that a DIA-NN report carries every column the downstream readers need.
#
# The CLI probe in modules/_shared/diann_env.py checks that the user's DIA-NN
# accepts the flags we pass; it cannot see what the search actually writes.
# DIA-NN has renamed and dropped report columns between releases, and both
# conduitR::diann_to_qfeatures() and the peptidotyping presence scripts select
# columns by name — a rename there surfaces as a confusing NULL/"object not
# found" error deep inside an R script, hours into a run. This rule turns that
# into one clear failure naming the missing columns.
#
# The required-column list is NOT duplicated here: it lives in
# modules/_shared/diann_env.py (DIANN_REQUIRED_REPORT_COLUMNS) and arrives via
# snakemake@params, so there is a single writer for the schema contract.

suppressPackageStartupMessages({
  library(arrow)
  library(jsonlite)
})

report_fp <- snakemake@input[["report"]]
required <- as.character(snakemake@params[["required_columns"]])
profile <- as.character(snakemake@params[["profile"]])
out_json <- snakemake@output[["report_schema_check"]]

# Read the schema only — never the rows. A monolithic InfiniDIA report can be
# tens of GB and we only need column names.
present <- names(arrow::open_dataset(report_fp, format = "parquet"))

missing <- setdiff(required, present)

result <- list(
  profile = profile,
  report = report_fp,
  ok = length(missing) == 0,
  n_columns_present = length(present),
  required_columns = required,
  missing_columns = missing,
  present_columns = present
)

dir.create(dirname(out_json), recursive = TRUE, showWarnings = FALSE)
write(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE), out_json)

if (length(missing) > 0) {
  stop(
    sprintf(
      paste0(
        "DIA-NN report %s is missing %d column(s) the '%s' consumers require:\n",
        "  %s\n\n",
        "Columns present in the report:\n  %s\n\n",
        "This almost always means the DIA-NN version in use writes a different\n",
        "report schema than this workflow expects. Check the DIA-NN version\n",
        "recorded in the run manifest against the tested versions listed in\n",
        "modules/_shared/diann_env.py (DIANN_TESTED_VERSIONS)."
      ),
      report_fp, length(missing), profile,
      paste(missing, collapse = ", "),
      paste(present, collapse = ", ")
    ),
    call. = FALSE
  )
}

cat(sprintf(
  "DIA-NN report schema OK: all %d required '%s' columns present in %s\n",
  length(required), profile, report_fp
))
