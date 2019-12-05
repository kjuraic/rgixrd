

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
