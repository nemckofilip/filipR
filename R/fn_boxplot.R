#' Enhanced Boxplot
#'
#' @description
#' Wrapper around base boxplot with cleaner defaults, optional violin overlay,
#' jittered outlier points, and p-value brackets between groups.
#'
#' @param x Formula, list, matrix, or data passed to boxplot.
#' @param compute.pval List of length-2 vectors with group indices to compare.
#'   E.g., list(c(1,2), c(2,3)) compares groups 1v2 and 2v3.
#' @param pval.FUN Function to compute p-value. Takes two numeric vectors, returns p-value.
#'   Default: function(x, y) wilcox.test(x, y)$p.value.
#' @param pval.cex Cex for p-value text. Default: 0.7.
#' @param pval.stars Show significance stars? Default: TRUE.
#' @param pval.values Show numeric p-values? Default: FALSE.
#' @param tilt.names Tilt x-axis labels? Default: FALSE.
#' @param srt Rotation angle for tilted names. Default: 45.
#' @param violin Show violin density polygons behind boxes? Default: FALSE.
#' @param viocol Violin fill color. Default: transparent.
#' @param viowex Violin width expansion factor. Default: 0.4.
#' @param outline Show outliers? Default: FALSE.
#' @param pts.col Outlier point color. Default: adjustcolor("lightgrey", 0.4).
#' @param pts.cex Outlier point size. Default: 0.6.
#' @param ... Extra parameters passed to boxplot().
#'
#' @examples
#' # Basic usage with p-values
#' fn_boxplot(len ~ dose, data = ToothGrowth,
#'            col = c("#e76f51", "#2a9d8f", "#264653"),
#'            compute.pval = list(c(1,2), c(2,3), c(1,3)))
#'
#' # Violin mode
#' fn_boxplot(len ~ dose, data = ToothGrowth,
#'            violin = TRUE, viocol = adjustcolor("#75AAD8", 0.3),
#'            col = "#75AAD8")
#'
#' @export
fn_boxplot <- function(x, ...) UseMethod("fn_boxplot")

#' @describeIn fn_boxplot Default method for lists and vectors
#' @export
fn_boxplot.default <- function(x, ...,
                               compute.pval = NULL,
                               pval.FUN = function(x, y) wilcox.test(x, y)$p.value,
                               pval.cex = 0.7,
                               pval.stars = TRUE,
                               pval.values = FALSE,
                               tilt.names = FALSE,
                               srt = 45,
                               range = 1.5,
                               width = NULL,
                               varwidth = FALSE,
                               notch = FALSE,
                               outline = FALSE,
                               pts.col = adjustcolor("lightgrey", 0.4),
                               pts.cex = 0.6,
                               names,
                               plot = TRUE,
                               border = par("fg"),
                               col = NULL,
                               log = "",
                               pars = list(boxwex = ifelse(violin, 0.2, 0.4),
                                           staplewex = NA,
                                           outwex = NA,
                                           outpch = NA),
                               horizontal = FALSE,
                               add = FALSE,
                               at = NULL,
                               frame = FALSE,
                               whisklty = ifelse(violin, 2, 1),
                               lwd = par("lwd"),
                               ylim = NULL,
                               xaxt = "s",
                               violin = FALSE,
                               viocol = "#FFFFFF00",
                               viowex = 0.4) {
  
  # Boxplot stats
  if(!missing(names) && is.function(names)) {
    box <- boxplot(x, ..., plot = FALSE)
    box$names <- names(box$names)
  } else {
    box <- boxplot(x, ..., names = names, plot = FALSE)
  }
  
  # Extract groups
  args <- list(x, ...)
  namedargs <- if(!is.null(attributes(args)$names)) attributes(args)$names != "" else rep_len(FALSE, length(args))
  groups <- if(is.list(x)) x else args[!namedargs]
  if(0L == (n <- length(groups))) stop("invalid first argument")
  if(length(class(groups))) groups <- unclass(groups)
  attr(groups, "names") <- box$names
  
  if(plot) {
    # Initial plot (transparent if violin mode)
    .out <- boxplot(x, ..., range = range, width = width, varwidth = varwidth,
                    notch = notch, outline = outline,
                    names = if(tilt.names && !horizontal) NA else box$names,
                    plot = plot, border = if(violin) NA else border,
                    col = if(violin) NA else col, log = log,
                    pars = pars, horizontal = horizontal, add = add, at = at,
                    frame = frame, whisklty = if(violin) 0 else whisklty,
                    lwd = lwd, ylim = ylim, xaxt = xaxt)
    
    # Violin overlay
    if(violin) {
      xpos <- if(is.null(at)) seq_along(groups) else at
      viocols <- rep(viocol, length.out = length(groups))
      
      for(j in seq_along(groups)) {
        vals <- groups[[j]]
        if(length(vals) < 2) next
        .d <- density(vals, from = box$stats[1, j], to = box$stats[5, j], na.rm = TRUE)
        vx <- .d$y / max(.d$y) * viowex / 2
        vx <- xpos[j] - c(vx, rev(-vx))
        vy <- c(.d$x, rev(.d$x))
        if(horizontal) polygon(vy, vx, col = viocols[j], lwd = lwd)
        else polygon(vx, vy, col = viocols[j], lwd = lwd)
      }
      
      # Redraw boxes on top of violins
      boxplot(x, ..., range = range, width = width, varwidth = varwidth,
              notch = notch, outline = outline,
              names = if(tilt.names && !horizontal) NA else box$names,
              plot = plot, border = border, col = col, log = log,
              pars = pars, horizontal = horizontal, add = TRUE, at = at,
              frame = frame, whisklty = whisklty, lwd = lwd, ylim = ylim,
              xaxt = "n", yaxt = "n")
    }
    
    # Jittered outlier points
    if(outline && length(.out$out)) {
      x.out <- if(!is.null(at)) at[.out$group] else .out$group
      set.seed(1)
      x.out <- jitter(x.out, amount = pars$boxwex / 2)
      if(horizontal) points(.out$out, x.out, col = pts.col, pch = 19, cex = pts.cex)
      else points(x.out, .out$out, col = pts.col, pch = 19, cex = pts.cex)
    }
    
    # P-value brackets
    if(!is.null(compute.pval)) {
      if(!is.list(compute.pval) || !all(lengths(compute.pval) == 2))
        stop("compute.pval must be a list of length-2 vectors")
      if(any(unlist(compute.pval) > length(groups)))
        stop("compute.pval indices exceed number of groups")
      
      pval <- .fn_compute_pval_brackets(groups, box, compute.pval, pval.FUN,
                                         outline, at, horizontal, pval.values)
      if(nrow(pval))
        .fn_draw_pval_brackets(pval, horizontal, pval.cex, pval.stars, pval.values)
    }
    
    # Tilted names
    if(tilt.names && !horizontal && xaxt != "n") {
      at_pos <- if(is.null(at)) seq_along(box$names) else at
      line.width <- diff(grconvertY(c(0, 1), "line", "user"))
      adj <- line.width * par("mgp")[2]
      y <- par("usr")[3] - adj
      text(at_pos, y, box$names, srt = srt, offset = -0.1, pos = 2,
           xpd = NA, cex = par("cex.axis"))
    }
  }
  
  invisible(box)
}

#' @describeIn fn_boxplot Method for matrices
#' @export
fn_boxplot.matrix <- function(x, use.cols = TRUE, ...) {
  groups <- if(use.cols) {
    split(c(x), rep.int(1L:ncol(x), rep.int(nrow(x), ncol(x))))
  } else {
    split(c(x), seq(nrow(x)))
  }
  if(length(nam <- dimnames(x)[[1 + use.cols]])) names(groups) <- nam
  invisible(fn_boxplot(groups, ...))
}

#' @describeIn fn_boxplot Method for formulas
#' @export
fn_boxplot.formula <- function(formula, data = NULL, ..., subset, na.action = NULL) {
  if(missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$na.action <- na.action
  m[[1L]] <- quote(stats::model.frame)
  mf <- eval(m, parent.frame())
  response <- attr(attr(mf, "terms"), "response")
  fn_boxplot(split(mf[[response]], mf[-response]), ...)
}


# =========================================================================
# Internal helpers
# =========================================================================

#' Compute p-value bracket positions with overlap detection
#' @keywords internal
.fn_compute_pval_brackets <- function(groups, box, compute.pval, pval.FUN,
                                       outline, at, horizontal, pval.values) {
  
  dat <- data.table::data.table(
    dat = groups,
    x = if(is.null(at)) seq_along(groups) else at
  )
  dat[, dat := lapply(dat, na.omit)]
  if(outline) dat[, max := sapply(dat, max, na.rm = TRUE)]
  else dat[, max := box$stats[5, ]]
  
  # Order pairs so x0 <= x1
  compute.pval <- lapply(compute.pval, function(i) i[order(dat[i, x])])
  
  pval <- cbind(
    dat[sapply(compute.pval, `[`, 1), !"max"],
    dat[sapply(compute.pval, `[`, 2), !"max"]
  )
  data.table::setnames(pval, c("dat0", "x0", "dat1", "x1"))
  
  # Compute p-values
  pval[, pvalue := mapply(function(x, y) pval.FUN(unlist(x), unlist(y)), x = dat0, y = dat1)]
  
  # Text x position
  pval[, x := rowMeans(.SD), .SDcols = c("x0", "x1")]
  
  # Compute y positions in inches with overlap detection
  adj <- if(horizontal) strwidth("M", units = "inch") else strheight("M", units = "inch")
  if(pval.values) adj <- adj * 1.5
  
  # Base y: max value between the two bracket endpoints
  pval[, y := max(dat[.BY, max, on = c("x>=x0", "x<=x1")]), .(x0, x1)]
  if(horizontal) pval[, y := grconvertX(y, "user", "inch")]
  else pval[, y := grconvertY(y, "user", "inch")]
  
  # Stack brackets to avoid overlaps
  data.table::setorderv(pval, c("y", "x0", "x1"))
  pval[, idx := .I]
  for(i in seq_len(nrow(pval))) {
    .c <- sort(pval[pval[i], y, on = c("x0<=x1", "x1>=x0", "idx<=idx")])
    .c <- .c[diff(c(.c, Inf)) > 2 * adj & .c >= pval[i, y]]
    pval[i, y := data.table::first(.c) + adj]
  }
  
  return(pval)
}

#' Draw p-value brackets with formatted labels
#' @keywords internal
.fn_draw_pval_brackets <- function(pval, horizontal, pval.cex, pval.stars, pval.values) {
  
  # Convert to user coordinates
  if(horizontal) {
    pval[, y0 := grconvertX(y, "inch", "user")]
    pval[, y1 := y0]
    pval[, y := grconvertX(y + strheight("M", "inch") * 0.45, "inch", "user")]
    data.table::setnames(pval,
                         c("x", "y", "x0", "x1", "y0", "y1"),
                         c("y", "x", "y0", "y1", "x0", "x1"))
  } else {
    pval[, y0 := grconvertY(y, "inch", "user")]
    pval[, y1 := y0]
    pval[, y := grconvertY(y + strwidth("M", "inch") * 0.45, "inch", "user")]
  }
  
  pval[, {
    # Draw bracket lines
    segments(x0, y0, x1, y1, xpd = NA)
    
    # Format and draw labels (addPval logic)
    mapply(function(px, py, p) {
      star <- if(pval.stars) {
        if(p < 1e-5) "****"
        else if(p < 1e-3) "***"
        else if(p < 1e-2) "**"
        else if(p < 5e-2) "*"
        else "N.S"
      } else ""
      
      lab <- if(pval.values) formatC(p, digits = 1, format = "e") else ""
      
      label <- if(pval.values) {
        if(p < 2.2e-308) bquote(italic(P) < "2.2e-308" * .(star))
        else if(p > 0.05) bquote(italic(P) == .(lab)^"N.S")
        else bquote(italic(P) == .(lab) * .(star))
      } else {
        if(p > 0.05) bquote(.(lab)^"N.S")
        else bquote(.(lab) * .(star))
      }
      
      text(px, py, labels = label,
           offset = ifelse(pval.values, -0.2, -0.35),
           pos = 3, cex = pval.cex,
           srt = ifelse(horizontal, -90, 0),
           xpd = NA)
    }, px = x, py = y, p = pvalue)
  }]
}