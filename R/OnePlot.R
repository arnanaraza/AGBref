# ================================================================
# OnePlot.R
# Robust PVIR plot-to-map diagnostic plots
# ================================================================

OnePlot <- function(x, y, resultsFolder, caption = "", fname = "", save = TRUE) {
  
  ok <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2) {
    warning("OnePlot skipped: fewer than 2 valid plot-map pairs.")
    return(invisible(NULL))
  }
  
  x <- ifelse(x > 350, 400, x)
  
  ct <- findInterval(
    x,
    c(0:12 * 25, 7:12 * 50, Inf),
    left.open = TRUE
  )
  
  ux <- aggregate(x, list(ct), FUN = mean, na.rm = TRUE)[, 2]
  uy <- aggregate(y, list(ct), FUN = mean, na.rm = TRUE)[, 2]
  nu <- aggregate(y, list(ct), FUN = function(z) length(na.omit(z)))[, 2]
  q1 <- aggregate(y, list(ct), FUN = quantile, probs = 0.25, na.rm = TRUE)[, 2]
  q3 <- aggregate(y, list(ct), FUN = quantile, probs = 0.75, na.rm = TRUE)[, 2]
  
  if (length(ux) == 0 || length(uy) == 0) {
    warning("OnePlot skipped: no valid binned values.")
    return(invisible(NULL))
  }
  
  cx <- (2:7 * 0.25)[findInterval(nu, c(0, 10, 20, 50, 100, 200))]
  r <- c(0, 425)
  
  n <- length(x)
  rmse <- sqrt(mean((y - x)^2, na.rm = TRUE))
  bias <- mean(y - x, na.rm = TRUE)
  r2 <- suppressWarnings(cor(x, y, use = "complete.obs")^2)
  
  caption2 <- paste0(
    caption,
    "\nN = ", n,
    " | RMSE = ", round(rmse, 1),
    " | Bias = ", round(bias, 1),
    " | R² = ", round(r2, 2)
  )
  
  draw_plot <- function() {
    
    plot(
      ux, uy,
      las = 1,
      main = caption2,
      pch = 16,
      cex = cx,
      xlab = "Mean reference AGB (Mg/ha)",
      ylab = "Mapped AGB (Mg/ha)",
      xlim = r,
      ylim = r,
      xaxt = "n"
    )
    
    arrows(ux, q1, ux, q3, length = 0)
    abline(0, 1, lty = 2, lwd = 1)
    
    legend(
      "topleft",
      bty = "n",
      pch = 16,
      ncol = 2,
      pt.cex = 2:7 * 0.25,
      cex = 0.9,
      legend = c("< 10", "10-20", "20-50", "50-100", "100-200", "> 200"),
      title = "#/bin"
    )
    
    if (length(ux) > 0 && tail(ux, 1) == 400) {
      axis(1, at = 0:4 * 100, labels = c(0:3 * 100, ">300"))
      if (requireNamespace("plotrix", quietly = TRUE)) {
        plotrix::axis.break(breakpos = 350)
      }
    } else {
      axis(1, at = 0:4 * 100, labels = c(0:4 * 100))
    }
  }
  
  if (save) {
    if (!dir.exists(resultsFolder)) {
      dir.create(resultsFolder, recursive = TRUE, showWarnings = FALSE)
    }
    
    png(file.path(resultsFolder, fname), 1000, 1000, res = 150)
    draw_plot()
    dev.off()
  } else {
    draw_plot()
  }
  
  invisible(NULL)
}


OnePlotScatter <- function(x, y, resultsFolder, caption = "", fname = "") {
  
  ok <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2) {
    warning("OnePlotScatter skipped: fewer than 2 valid plot-map pairs.")
    return(invisible(NULL))
  }
  
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder, recursive = TRUE, showWarnings = FALSE)
  }
  
  n <- length(x)
  rmse <- sqrt(mean((y - x)^2, na.rm = TRUE))
  bias <- mean(y - x, na.rm = TRUE)
  r2 <- suppressWarnings(cor(x, y, use = "complete.obs")^2)
  
  caption2 <- paste0(
    caption,
    "\nN = ", n,
    " | RMSE = ", round(rmse, 1),
    " | Bias = ", round(bias, 1),
    " | R² = ", round(r2, 2)
  )
  
  r <- c(0, 800)
  
  png(file.path(resultsFolder, fname), 1000, 1000, res = 150)
  
  plot(
    x, y,
    las = 1,
    main = caption2,
    pch = 16,
    cex = 0.5,
    xlab = "Reference AGB (Mg/ha)",
    ylab = "Mapped AGB (Mg/ha)",
    xlim = r,
    ylim = r,
    xaxt = "n"
  )
  
  abline(0, 1, lty = 2, lwd = 1)
  axis(1, at = 0:6 * 200, labels = c(0:6 * 200))
  
  dev.off()
  
  invisible(NULL)
}