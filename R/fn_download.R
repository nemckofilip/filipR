#' Download a file from a URL if not already present
#'
#' @description
#' Downloads a file using `wget`, skipping if the destination file already
#' exists. Prints a message for each download attempt and on success.
#'
#' @param url Character. URL to download from.
#' @param dest Character. Destination file path.
#'
#' @return Invisibly returns `dest`.
#' @export
#' @examples
#' \dontrun{
#' fn_download("https://example.com/file.gz", "data/file.gz")
#' }
fn_download <- function(url, dest) {
  if (!file.exists(dest)) {
    message("Downloading: ", basename(dest))
    system(paste("wget -O", shQuote(dest), shQuote(url)))
    if (file.exists(dest)) message("  Done: ", dest)
  } else {
    message("Already exists, skipping: ", dest)
  }
  invisible(dest)
}
