
#' add GIXRD style to ggplot of experimental data
#'
#' @author K. Juraic
#' @param p ggplot2 with experimental data
#' @return
#' @import ggplot2
#' @export
#' @examples \dontrun{xrd_ggplot_style(p)}
xrd_ggplot_style <- function(p) {
  pp <- p + theme_bw(base_size = 20) +
            theme(legend.position = "none") +
            xlab(expression(2*theta~'['*degree*']')) +
            ylab('Intensity [a.u.]')
  pp
}


#' add marks for XRD peak positions (tth) to gglot XRD graph
#'
#' @author K. Juraic
#' @param p ggplot XRD graph
#' @param peaks_df data.frame with column tth (xrd peaks possitions)
#' @param stick_y (y_min, y(max)) y coordinates of stick stard and end
#' @return
#' @export
#' @iport ggplot2
#' @examples \dontrun{xrd_ggplot_add_peak_pos(p, df, c(0,100))}
xrd_ggplot_add_peak_pos <- function(p, peaks_df, stick_y = c(-100, 0)) {
  pp <- p + geom_linerange(data = peaks_df, mapping = aes_(x = ~ tth, ymin = stick_y[1], ymax = stick_y[2]))
  pp
}
