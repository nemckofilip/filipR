#' Run AlphaFold 3
#'
#' @description Builds a shell command that runs AlphaFold 3 (Apptainer container) on every
#' input JSON in a directory, in one process (the model is compiled once). The MSA search
#' (`run.data.pipeline`) is CPU only and can be run apart from the inference (`run.inference`,
#' needs a GPU). Submit with [fn_submit()], with `gpus = 1` for the inference.
#' @param input.dir Character. Directory with AlphaFold 3 input JSON files.
#' @param output.dir Output directory. Default 'db/af3/'.
#' @param base.name Name of this batch; the done-file is `<output.dir>/<base.name>.done`.
#' @param run.data.pipeline Logical. Run the MSA/template search. Default TRUE.
#' @param run.inference Logical. Run the structure prediction. Default TRUE.
#' @param jax.backend Character. "gpu" (default) or "cpu" for the inference; "cpu" is about 100x slower and not officially supported by AlphaFold 3.
#' @param num.seeds Integer. Seeds per input (only used if the JSON has a single `modelSeeds`). Default NULL (use the JSON).
#' @param sif Path to the AlphaFold 3 container.
#' @param model.dir Directory with the model parameters (af3.bin).
#' @param db.dir Directory with the AlphaFold 3 databases.
#' @param cores Threads for the MSA search. Default 8.
#' @return data.table with columns: input.dir, output.dir, cmd, path (the done-file).
#' @import data.table
#' @export
fn_alphafold3 <- function(input.dir,
                          output.dir = "db/af3/",
                          base.name,
                          run.data.pipeline = TRUE,
                          run.inference = TRUE,
                          jax.backend = "gpu",
                          num.seeds = NULL,
                          sif = "/large_storage/hornslab/userspace/filip/references/alphafold3/alphafold3_v3.0.4.sif",
                          model.dir = "/common_datasets/alphafold3/models/af3_weights",
                          db.dir = "/common_datasets/alphafold3/databases/af3_dbs",
                          cores = 8) {
  if (!dir.exists(input.dir)) stop("input.dir does not exist: ", input.dir)
  if (!file.exists(sif)) stop("Container not found: ", sif)
  if (!run.data.pipeline && !run.inference) stop("Nothing to run.")
  if (!jax.backend %in% c("gpu", "cpu")) stop("jax.backend must be 'gpu' or 'cpu'.")
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

  # ---- Command ----
  input.dir <- normalizePath(input.dir); output.dir <- normalizePath(output.dir)
  done <- file.path(output.dir, paste0(base.name, ".done"))
  cmd <- paste("apptainer exec", if (run.inference && jax.backend == "gpu") "--nv",
               "--bind", shQuote(input.dir), "--bind", shQuote(output.dir),
               "--bind", shQuote(paste0(model.dir, ":", model.dir, ":ro")),
               "--bind", shQuote(paste0(db.dir, ":", db.dir, ":ro")),
               shQuote(sif), "python /app/alphafold/run_alphafold.py",
               paste0("--input_dir=", shQuote(input.dir)),
               paste0("--output_dir=", shQuote(output.dir)),
               paste0("--model_dir=", shQuote(model.dir)),
               paste0("--db_dir=", shQuote(db.dir)),
               paste0("--run_data_pipeline=", tolower(run.data.pipeline)),
               paste0("--run_inference=", tolower(run.inference)),
               paste0("--jackhmmer_n_cpu=", cores), paste0("--nhmmer_n_cpu=", cores))
  if (jax.backend == "cpu") cmd <- paste(cmd, "--jax_backend=cpu --flash_attention_implementation=xla")
  if (!is.null(num.seeds)) cmd <- paste(cmd, paste0("--num_seeds=", num.seeds))
  cmd <- paste(cmd, "&&", "touch", shQuote(done))

  data.table::data.table(input.dir = input.dir, output.dir = output.dir, cmd = cmd, path = done)
}
