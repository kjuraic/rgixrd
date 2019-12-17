
#' read GIXRD data from MCX file
#'
#' @param xrd_file name of filme with MCX GIXRD data
#' @return list(xrd_file, alpha_i, data.frame(tth, intensity, background))
#' @export
#' @examples \dontrun{read_mcx_xrd()}
#' @import tcltk
read_mcx_xrd <- function(xrd_file = tcltk::tk_choose.files(multi = FALSE)){
  if (file.exists(xrd_file)) {
    # angle of incidence
    tmp <- strsplit(readLines(con = xrd_file, n = 1), split = " ")[[1]]
    if (tmp[1] == "#") {
      alpha_i <- as.numeric(tmp[4])
    } else {
      alpha_i <- NA
    }
    # xrd data
    xrd_dat <- utils::read.table(xrd_file, skip = 1)
    names(xrd_dat) <- c("tth","intensity","background")
    cat("[", xrd_file, "] Reeading ... OK!\n")
    cat("\t alpha.i =", alpha_i, "\n")
    cat("\t 2*theta = [", min(xrd_dat$tth), " - ", max(xrd_dat$tth), "]   N-points = ", length(xrd_dat$tth), "\n" )
  } else {
    cat("[", xrd_file, "] does not exist\n")
  }
  list(file = xrd_file, alpha_i = alpha_i, data = xrd_dat)
}



#' read multiple MCX XRD files
#'
#' @param xrd_files array with full file namse
#' @return list of xrd data (for each file name single list element)
#' @export
#' @examples \dontrun{read_mcx_xrd_multi}
#' @import plyr
read_mcx_xrd_multi <- function(xrd_files = tcltk::tk_choose.files()) {
  xrd_lst <- plyr::alply(.data = xrd_files, .margins = 1, .fun = read_mcx_xrd)
  xrd_lst
}


# read_spec(file_name) ----------------------------------------------
#' Read GIXRD data from file (Buljan instrument)
#' @description Read GIXRD data from file (Buljan instrument) to list.
#'              In one file are writen all scans for one single sample.
#'              Each list element coresponds to one scan in file.
#'              Columns description for file format: MUST BE ADDED
#' @author K. Juraic
#' @param file_name GIXRD data file name
#' @return list with GIXRD data
#' @importFrom stringr str_replace
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @export
#' @examples
#'             \dontrun{read_spec("gixrd_name")}
read_spec <- function(file_name){
  tmp <- readLines(file_name)
  blank_line <- which(tmp == "")
  blank_line <- c(blank_line, length(tmp))
  #diff(blank_line)
  xrd_dat <- list()
  for(i in 1:(length(blank_line)-1)){
    scan_lines <- tmp[(blank_line[i]+1):blank_line[i+1]]
    scan_header <- scan_lines[startsWith(x = scan_lines, "#")]
    scan_points <- scan_lines[!startsWith(x = scan_lines, "#")]
    xrd_command <- scan_lines[startsWith(x = scan_lines, "#S")]
    xrd_date <- scan_lines[startsWith(x = scan_lines, "#D")]
    xrd_acqTime <- scan_lines[startsWith(x = scan_lines, "#T")]
    n_column <- scan_lines[startsWith(x = scan_lines, "#N")] %>%
      str_replace(pattern = "#N ", replacement = "") %>%
      as.numeric()
    n_points <- scan_lines[startsWith(x = scan_lines, "#P1")] %>%
      str_replace(pattern = "#P1 ", replacement = "") %>%
      as.numeric()
    init_values <-scan_lines[startsWith(x = scan_lines, "#P0")] %>%
      str_replace(pattern = "#P0 ", replacement = "") %>%
      str_split(pattern = " +") %>% unlist() %>% as.numeric()
    init_values_name <-tmp[startsWith(x = tmp, "#O0")] %>% str_split(pattern = " +") %>% unlist()
    names(init_values) <- init_values_name[-1]
    col_nms <- strsplit(x = scan_lines[startsWith(x = scan_lines, "#L")], split = " +") %>% unlist()
    dat <- read.table(textConnection(scan_points))
    colnames(dat) <- col_nms[-1]
    xrd_dat[[i]] <- list(command = xrd_command, date = xrd_date, acqTime = xrd_acqTime, init_values, dat = dat)
  }
  return(xrd_dat)
}

# gixrd_average(xrd_dat) -------------------------------------------------
#' gixrd_average(xrd_dat)
#' @description calculate average difractogram from list of GIXRD scans read by
#'              gixrd_read_file function (Maja difractomerer file format).
#' @param xrd_dat list with GIDRD scans data
#' @param scans_to_average array of scans ids to average
#' @return data frame with calculated average
#' @examples \dontrun{gixrd_average(gixrd_data)}
#' @importFrom stats sd median
#' @importFrom utils read.table
#' @export
gixrd_average <- function(xrd_dat, scans_to_average = 1:length(xrd_dat)){
  xrd_dat <- xrd_dat[scans_to_average]
  tth <- xrd_dat[[1]]$dat$TwoTheta
  intensity_mat <- matrix(NA, nrow = length(xrd_dat), ncol = length(tth))
  for (i in 1:length(xrd_dat)) {
    intensity_mat[i,1:length(xrd_dat[[i]]$dat$Detector)] <- xrd_dat[[i]]$dat$Detector
  }
  intensity_mean <- colMeans(intensity_mat, na.rm = TRUE)
  intensity_sd <- apply(X = intensity_mat, MARGIN = 2, FUN = sd, na.rm = TRUE)
  intensity_median <- apply(X = intensity_mat, MARGIN = 2, FUN = median, na.rm = TRUE)
  xrd_av <- data.frame(tth, Imean = intensity_mean, Isd = intensity_sd, Imedian = intensity_median)
  return(xrd_av)
}


#' read_mercury_hkl(file_name)
#' @description read mercury hkl file with d spacing [Angstrom],
#'              structure factor and multiplicit for (hkl) difraction peaks
#' @author K. Juraic
#' @param file_name mercurt *.hkl file name
#' @return data.frame(h, k, l, d, f2, mult)
#' @export
#' @examples
#' \dontrun{read_mercury_hkl("ZnO.hkl")}
#'
read_mercury_hkl <- function(file_name = file.choose()) {
  if (file.exists(file_name)) {
    dat = read.table(file = file_name, header = TRUE)
    if (dim(dat)[2] != 6) {
      cat("File format not valid!\n\t File should have 6 columns!\n")
    } else {
      names(dat) <- c("h", "k", "l", "d", "f2", "mult")
      return(dat)
    }
  } else {
    cat("File does not exist!\n")
  }
}



