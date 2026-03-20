#' @importFrom magrittr "%>%"
#' @export
RC_importIPAM <- function(dir, key, ipam.pattern="*.csv") {
  ipam <- list.files(path=dir, pattern=ipam.pattern, full.names=TRUE)
  basenames <- gsub("\\..*$", "", basename(ipam))

  # Read in each file
  for (i in 1:length(ipam)) {

    # Read in IPAM data file (should be a single row, ;-delimited)
    df <- suppressMessages(readr::read_delim(ipam[i], ";",
                            escape_double = FALSE, trim_ws = TRUE,
                            show_col_types = FALSE, name_repair = 'unique')) %>%
      dplyr::arrange(desc(`No.`)) %>%
      purrr::discard(~all(is.na(.)))
    df <- df[1, ]
    date = dmy(as.character(df[1,1]))

    # Make data long form and recast with rows for each AOI (based on column names)
    df2 <- df %>%
      dplyr::select(-c(1:4)) %>%
      tidyr::pivot_longer(everything(),
                   names_to = "variable",
                   values_to = "value") %>%
      dplyr::mutate(AOI = as.numeric(stringr::str_extract(variable, pattern = "[0-9]+")),
             var = stringr::str_extract(variable, pattern = "[A-Z,a-z]+"))
    df3 <- reshape2::dcast(na.omit(df2), AOI ~ var)

    df4 <- df3 %>%
      dplyr::mutate(Date = date)

    # Merge data with IDs and assign to global environment
    df4$file <- basenames[i]
    assign(basenames[i], df4)

  }

  # Merge all data frames together
  df5 <- do.call(rbind, mget(basenames))
  rownames(df5) <- NULL  # gets rid of rownames

  # Return final data frame
  if (is.null(key)) {
    result <- df5
    return(result)
  } else {
    result <- df5 %>%
      merge(key)
    return(result)
  }
}
