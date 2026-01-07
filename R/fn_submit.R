#' Submit Commands to Arc Slurm Cluster
#'
#' @description
#' Submits shell commands to the Arc Slurm cluster using sbatch.
#' Handles resource allocation, log directory creation, and script generation.
#'
#' @param cmd A data.table containing the commands. Intended columns:
#'   - `cmd`: Shell command to execute.
#'   - `path`: (Optional) Output file path (used for existence checks).
#'   - `file.type`: (Optional) Label.
#' @param job.name Name of the job. Default: "myJob".
#' @param partition Slurm partition (e.g., "cpu", "cpu_batch", "gpu"). Default: "cpu".
#' @param cores Number of CPUs per task. Default: 4.
#' @param mem Memory allocation (e.g., "16G"). Default: "16G".
#' @param time Time limit (HH:MM:SS). Default: "04:00:00".
#' @param logs Directory for log files and submission scripts. Default: "db/logs".
#' @param conda_env Name of the conda environment to activate. Default: "/home/filip.nemcko/miniconda/envs/r_env".
#' @param overwrite If TRUE, overwrites existing output files. Default: FALSE.
#' @param execute If TRUE, submits the job. If FALSE, only creates the script. Default: TRUE.
#' @param create.output.dirs Create output directories before running? Default: TRUE.
#'
#' @import data.table
#' @export
fn_submit <- function(cmd,
                      job.name= "myJob",
                      partition= "cpu_batch",
                      cores= 4,
                      mem= "16G",
                      time= "04:00:00",
                      logs= "db/logs",
                      conda_env= "/home/filip.nemcko/miniconda/envs/r_env",
                      overwrite= FALSE,
                      execute= TRUE,
                      create.output.dirs= TRUE)
{
  # Input checks
  if(!is.data.table(cmd) || !"cmd" %in% names(cmd))
    stop("cmd must be a data.table with a 'cmd' column.")
  
  # Filter existing files if not overwriting
  if(!overwrite && "path" %in% names(cmd)) {
    cmd <- cmd[is.na(path) | !file.exists(path)]
  }
  
  if(nrow(cmd) == 0) {
    message("No commands to run (files exist or empty input).")
    return(invisible(NULL))
  }

  # Create directories
  if(create.output.dirs && "path" %in% names(cmd)) {
    paths <- cmd$path[!is.na(cmd$path)]
    dirs <- unique(dirname(paths))
    lapply(dirs, function(d) if(!dir.exists(d)) dir.create(d, recursive=TRUE, showWarnings=FALSE))
  }
  
  if(!dir.exists(logs)) dir.create(logs, recursive=TRUE)

  # Script generation
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  script_path <- file.path(logs, paste0(job.name, "_", timestamp, ".sh"))
  
  # Slurm Header
  header <- c(
    "#!/bin/bash",
    paste0("#SBATCH --job-name=", job.name),
    paste0("#SBATCH --output=", file.path(logs, paste0(job.name, "_%j.out"))),
    paste0("#SBATCH --error=", file.path(logs, paste0(job.name, "_%j.err"))),
    paste0("#SBATCH --partition=", partition),
    paste0("#SBATCH --cpus-per-task=", cores),
    paste0("#SBATCH --mem=", mem),
    paste0("#SBATCH --time=", time)
  )
  
  # Environment setup
  env_setup <- c("source /opt/conda/etc/profile.d/conda.sh")
  
  if(!is.null(conda_env)) {
    env_setup <- c(env_setup, paste("conda activate", conda_env))
  }
  
  # Commands
  # We write all commands to the script. They will run sequentially.
  body <- cmd$cmd
  
  # Write the full script
  writeLines(c(header, "", env_setup, "", body), script_path)
  
  if(execute) {
    # Submit using sbatch
    # --parsable returns only the job ID
    submission <- system(paste("sbatch --parsable", script_path), intern=TRUE)
    message(paste("Submitted batch job", submission))
    return(invisible(submission))
  } else {
    message(paste("Script created at:", script_path))
    return(script_path)
  }
}

