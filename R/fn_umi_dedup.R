#' Deduplicate Reads using UMI-tools
#'
#' @name fn_umi_dedup
#' @description
#' Removes PCR duplicates using UMI-tools.
#' Assumes UMIs are in the read header (read_id).
#' Requires the input BAM to be sorted.
#' Automatically creates a temporary index (.bai) if one is missing, then deletes it.
#' Automatically detects paired-end vs single-end from the BAM file.
#'
#' @param bam Character string. Path to the sorted BAM file.
#' @param base.name Base name for output files.
#' @param output.dir Directory for deduplicated BAMs. Default: 'db/bam_dedup/'.
#' @param stats.dir Directory for UMI stats. Default: 'db/umi_stats/'.
#' @param umi_separator The character separating the Read ID from the UMI. Default '_'.
#' @param paired Logical. If TRUE, use paired-end mode. If FALSE, single-end. If NULL (default), auto-detect from BAM.
#' @param cores Number of threads. Used for the temporary indexing step.
#'
#' @return A data.table with columns: bam_in, bam_dedup, stats, cmd, job.name.
#' @export
fn_umi_dedup <- function(bam, 
                         base.name,
                         output.dir = 'db/bam_dedup/', 
                         stats.dir = 'db/umi_stats/',
                         umi_separator = '_',
                         paired = NULL,
                         cores = 4) {

  # 1. Checks
  if(length(bam) > 1) stop('fn_umi_dedup processes one file at a time.')
  if(!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  if(!dir.exists(stats.dir)) dir.create(stats.dir, recursive = TRUE)

  # 2. Define Paths
  bam_dedup <- file.path(output.dir, paste0(base.name, '_dedup.bam'))
  stats_prefix <- file.path(stats.dir, paste0(base.name, '_dedup'))
  
  # 3. Build Command
  
  # Step A: Check for index, create if missing (umi_tools needs it)
  cmd_index <- paste0(
    'if [ ! -f ', bam, '.bai ]; then ',
    'echo "Creating temporary index..."; ',
    'samtools index -@ ', cores, ' ', shQuote(bam), '; ',
    'INDEX_CREATED=1; ',
    'else INDEX_CREATED=0; fi'
  )

  # Step B: Detect paired-end if not specified
  if(is.null(paired)) {
    # Auto-detect: check if BAM has paired reads (flag 0x1)
    cmd_detect <- paste0(
      'if samtools view -c -f 1 ', shQuote(bam), ' | grep -q "^0$"; then ',
      'PAIRED_FLAG=""; ',
      'else PAIRED_FLAG="--paired"; fi'
    )
  } else {
    cmd_detect <- paste0('PAIRED_FLAG="', if(paired) '--paired' else '', '"')
  }

  # Step C: Run umi_tools dedup
  cmd_dedup <- paste(
    'umi_tools dedup',
    '-I', shQuote(bam),
    '-S', shQuote(bam_dedup),
    '$PAIRED_FLAG',
    '--extract-umi-method=read_id',
    paste0('--umi-separator=', shQuote(umi_separator)),
    '--method=directional',
    '--spliced-is-unique',
    '--output-stats', shQuote(stats_prefix),
    '--log', shQuote(paste0(stats_prefix, '.log'))
  )

  # Step D: Cleanup Index (only if we created it)
  cmd_cleanup <- paste0(
    'if [ $INDEX_CREATED -eq 1 ]; then ',
    'rm ', shQuote(paste0(bam, '.bai')), '; ',
    'fi'
  )

  # Combine
  full_cmd <- paste(cmd_index, cmd_detect, cmd_dedup, cmd_cleanup, sep = ' && ')

  # 4. Return metadata
  data.table::data.table(
    bam_in = bam,
    bam_dedup = bam_dedup,
    stats = paste0(stats_prefix, '_edit_distance.tsv'),
    cmd = full_cmd,
    job.name = paste0('dedup_', base.name),
    path = bam_dedup
  )
}