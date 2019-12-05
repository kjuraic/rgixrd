
#' d_2_tth(d, lambda)
#'
#' @description calculate position of XRD peak in 2*theta scale for d-spacing and
#'              wavelength lambda. d-spacing and wavelength shoud be given in same
#'              units: angstroms or nanometers
#' @author K. Juraić
#' @param d d-spacing in angstroms
#' @param lambda X-ray wavelength in angstroms
#' @return tth 2*theta position of XRD peak
#' @export
#' @examples d_2_tth(1.59, 1.54056)
d_2_tth <- function(d, lambda = 1.54056) {
  tth <- 360 / pi * asin(.5 * lambda / d)
  tth
}


#' d_tth_2_d(tth, lambda)
#'
#' @description calculate d-spacing for tth position of XRD peak in 2*theta scale
#'              tth should be given in degrees
#' @author K. Juraić
#' @param tth 2*theta position of XRD peak in degrees
#' @param lambda X-ray wavelength in angstroms or nanometers
#' @return d-spacing in same units as wavelength
#' @export
#' @examples tth_2_d(28.1, 1.54056)
tth_2_d <- function(tth, lambda = 1.54056) {
  d <- .5 * lambda / sin(tth * pi / 360.)
  d
}


#' tetragonal_hkl_2_d(h, k, l, a, c)
#'
#' @description For tetragona crystall structure calculate d-spacing by using
#'              known lattice paramteres (a, c) and Muller indexes
#' @author K. Juraic
#' @param h millerov indeks
#' @param k millerov indeks
#' @param l millerov index
#' @param a lattice parameter
#' @param c lattice parameter
#' @return d-spacing : razmak mrežnih ravnina
#' @export
#' @examples \dontrun{tetragonal_hkl_2_d(h = 1, k = 1, l = 0
#'                                       a = 4.73727, c = 3.186383)}
tetragonal_hkl_2_d <- function(h, k, l, a, c) {
  d <- 1 / sqrt((h^2 + l^2) / a^2 + l^2 / c^2)
  d
}


#' tetragonal_d_2_c(h, k, l, d, a)
#'
#' @description calculate lattice parameter c for tetragona crystall structure
#'              from d-spacing from known Muller indexes (h,k,l) and lattice
#'              paramter a. For specioal case when Muller index h = 0, c can be
#'              calculated independently on Muller index a
#' @author K. Juraić
#' @param h Muller index
#' @param k Muller index
#' @param l Muller index
#' @param d d-spacing
#' @param a lattice paramerr
#'
#' @return c lattice parameter of hexagonal crystall structure
#' @export
#'
#' @examples tetragonal_d_2_c(0, 0, 1, 1.59, 1)
tetragonal_d_2_c <- function(h, k, l, d, a) {
  c <- l / sqrt(1/d^2 - (h^2 + k^2) / a^2)
  c
}


#' tetragonal_d_2_a(h, k, l, d, c)
#'
#' @description calculate lattice parameter a for tetragona crystall structure
#'              from d-spacing from known Muller indexes (h,k,l) and lattice
#'              paramter c. For specioal case when Muller index l = 0, a can be
#'              calculated independently on Muller index c
#' @author K. Juraić
#' @param h Muller index
#' @param k Muller index
#' @param l Muller index
#' @param d d-spacing
#' @param c lattice paramerr
#'
#' @return a lattice parameter of hexagonal crystall structure
#' @export
#'
#' @examples tetragonal_d_2_c(1, 0, 0, 1.59, 1)
tetragonal_d_2_a <- function(h, k, l, d, c) {
  a <- sqrt((h^2 + k^2) / (1/d^2 - l^2 / c^2))
  a
}
