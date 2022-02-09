
#' Shift intensity in steps for presenting data ain shifted plot
#' @author K. Juraić
#' @description Shift of intensity data for shifted plot. If shift factor is
#'              single number it is constructed intensity shift array with
#'              equidistant shift factor for each diffractogram in data.frame.
#'              Array lengt is equial to number of diffractograms in data frames.
#'              If shift array is provided as aaray of same length as number of
#'              diffractograms in data.frames is directly uses to shift data in
#'              data.frames.
#' @param dat_df data.frame(tth, intensity, name) with XRD data
#' @param shift_factor single number or array with values to shift intensity data
#'
#' @return data.frame(tth, intensity, name)
#' @export
#'
#' @examples \dontrun{xrd_intensity_shift(dat_df= xrd_data, shift_factor = 100)}
xrd_intensity_shift <- function(dat_df, shift_factor = 0) {
  df_lst <- split(dat_df, dat_df$name)
  if (length(shift_factor) == 1) {
    intensity_shift <- seq(from = 0, to = length(df_lst) - 1) * shift_factor
  } else {
    intensity_shift = shift_factor
  }
  df <- map2_df(.x = df_lst, .y = intensity_shift, .f = ~mutate(.x, intensity = intensity + .y))
  return(df)
}
