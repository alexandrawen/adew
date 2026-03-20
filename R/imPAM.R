#' @importFrom magrittr "%>%"
#' @export
imPAM <- function(dir, key, ipam.pattern="*.csv") {
  ipam <- list.files(path=dir, pattern=ipam.pattern, full.names=TRUE)
  names <- list.files(path=dir, pattern=ipam.pattern, full.names=FALSE)

  data <- NULL
  # Read in each file
  for (i in 1:length(ipam)) {
    # Read in IPAM data file (should be a single row, ;-delimited)
    df <- suppressMessages(readr::read_delim(ipam[i], ";",
                                             escape_double = FALSE, trim_ws = TRUE,
                                             show_col_types = FALSE, name_repair = 'unique')) %>%
      purrr::discard(~all(is.na(.))) %>%
      dplyr::arrange(desc(`No.`)) %>%
      dplyr::slice(1)
    date <- lubridate::dmy(as.character(df[1,1]))
    PAR <- as.numeric(df[1,4])
    file <- substr(names[i],0,(nchar(names[i])-nchar(ipam.pattern)+1))

    # Make data long form and recast with rows for each AOI (based on column names)
    df2 <- df[ , -c(1:4)] %>%
      tidyr::pivot_longer(cols = everything(),
                          names_to = "variable",
                          values_to = "value") %>%
      dplyr::mutate(AOI = as.numeric(stringr::str_extract(variable, pattern = "[0-9]+")),
                    var = stringr::str_extract(variable, pattern = "[A-Z,a-z]+")) %>%
      tidyr::pivot_wider(id_cols = AOI,
                         names_from = var,
                         values_from = value) %>%
      dplyr::mutate(Date = date,
                    PAR = PAR,
                    file = file)

    data <- plyr::rbind.fill(data,df2)
  }

  if (is.null(key)) {
    data %>%
      return(.)
  } else {
    key %<>%
      dplyr::group_by(file) %>%
      dplyr::mutate(N = dplyr::n())

    data %>%
      dplyr::group_by(file) %>%
      dplyr::mutate(N = dplyr::n()) %>%
      dplyr::full_join(key) %>%
      dplyr::select(-c(N)) %>%
      return(.)
  }
}
