
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
    xrd_dat <- data.frame()
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

#' file_name_2_sample_name(file_name)
#' @author K. Juraić
#' @description extract sample name from file full path, removing directory
#'              structure and file extension
#' @param file_name file full path
#'
#' @return sample_name extracted sample name
#' @export
#'
#' @examples file_name_2_sample_name("../ime_uzorka.dat")
file_name_2_sample_name <- function(file_name) {
  sample_name <- tools::file_path_sans_ext(base::basename(file_name))
  return(sample_name)
}

#' read_xrd_data(data_dir =  tcltk::tk_choose.dir(), pattern = "")
#' @author K. Juraić
#' @description read xrd experimental data from folder selected by file extension
#'              At the moment applicable for .xy and Bruker .raw format
#' @param data_dir folder path with xrd files
#' @param pattern  file extension (.xy, .raw)
#'
#' @return data.frame(tth, intesity)
#' @export
#' @importFrom purrr map
#' @importFrom purrr map2_df
#'
#' @examples \dontrun{read_xrd_data("../data)}
read_xrd_data <- function(data_dir =  tcltk::tk_choose.dir(), pattern = ""){
  fnms <- list.files(path = data_dir, pattern = pattern, full.names = TRUE)
  sample_names <- file_name_2_sample_name(fnms)
  if (pattern == ".xy") {
    dat_lst <- purrr::map(.x = fnms, .f = read_xy)
    dat_df <- purrr::map2_df(.x = dat_lst, .y = sample_names, .f = ~ mutate(.x, name = .y))
  } else if (pattern == ".raw") {
    dat_lst <- purrr::map(.x = fnms, .f = read_bruker_raw4)
    dat_df <- purrr::map2_df(.x = dat_lst, .y = sample_names, .f = ~ mutate(.x, name = .y))
  }
  cat(paste("Found", length(fnms), "files:\n"))
  print(fnms)
  return(dat_df)
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
    xrd_dat[[i]] <- list(command = xrd_command,
                         date = xrd_date,
                         acqTime = xrd_acqTime,
                         init_values = init_values,
                         dat = dat)
  }
  return(xrd_dat)
}


#' read XRD data from xy file
#' @author K. Juraić
#' @description read XRD data from xy file. Can be generated from raw instrument
#'              data with powDLL software. Data are stored as two
#'              column (tth, intensity).
#' @param file_name XRD data filename (full path)
#'
#' @return data.frame(tth, intensity)
#' @importFrom readr read_table2
#' @export
#'
#' @examples \dontrun{read_xy(""XRD_pattern.xy)}
#'
read_xy <- function(file_name) {
  xy_data <- read_table2(file_name, col_names = c("tth", "intensity"))
  return(xy_data)
}


#' read Bruker raw file (version 4)
#' @author K. Juraić
#' @description read Bruker XRD raw file (version 4.0)
#' @param file_name full path of file with XRD data
#'
#' @return data.frame(tth, intensity)
#' @export
#'
#' @examples \dontrun{read_bruker_raw4("XRD_pattern.raw")}
#'
read_bruker_raw4 <- function(file_name) {
  file_size <- file.size(file_name)
  con = file(file_name, "rb")
  # Bruker raw file version
  bruker_version <- readBin(con = con, character(), size = 5, n = 1)
  bruker_version
  # Date
  seek(con, 12)
  bruker_date <- readBin(con = con, character(), size = 10, n = 1)
  bruker_date
  #Time
  seek(con, 24)
  bruker_time <- readBin(con = con, character(), size = 8, n = 1)
  bruker_time
  # N points
  seek(con, 471)
  n_points <- readBin(con = con, integer(), n = 1)
  n_points
  # data tth start and step
  seek(con, 539)
  dat <- readBin(con = con, double(), n = 2)
  dat[1:2]
  # tth scale
  tth <- seq(from = dat[1], by = dat[2], length.out = n_points)
  # byte where intensity data starts (32bit double)
  data_start_byte <- file_size - n_points*4
  seek(con, data_start_byte)
  # read intensity data
  intensity <- readBin(con = con, double(), size = 4, n = n_points)
  close(con)
  bruker_data <- data.frame(tth = tth, intensity = intensity)
  return(bruker_data)
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


# gixrd_read_file_lst(file_names, save_to_file = TRUE) ----------------------
#' gixrd_read_file_lst(file_names, save_to_file = TRUE)
#' @description Automaticaly read multiple files with GIXRD data and
#'              average them and write averaged data to files
#' @param file_names array of file names with GIXRD data (full paths)
#' @param save_to_file Should be the averaged data saved to file
#' @return list of data frames with averaged data
#' @examples \dontrun{gixrd_aveeage_lst(fnms)}
#' @export
gixrd_average_lst <- function(file_names, save_to_file = FALSE) {
  xrd_lst <- list()
  for (i in 1:length(file_names)) {
    cat("[",i,"] ", file_names[i],"\n")
    xrd_dat <- read_spec(file_names[i])
    xrd_av <- gixrd_average(xrd_dat)
    xrd_lst[[i]] <- xrd_av
    if (save_to_file == TRUE) {
      fnm <- paste0(file_names[i], ".dat")
      utils::write.table(x = xrd_av,
                  file = fnm,
                  sep = "\t",
                  row.names = FALSE, col.names = TRUE)
      cat("\tWriting to file:", fnm, "\n")
    }
  }
  return(xrd_lst)
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





