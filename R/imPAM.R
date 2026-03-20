#' Import Imaging-PAM AOI summary files
#'
#' Import one or more Imaging-PAM `.csv` files from a directory, extract the
#' most recent summary row from each file, reshape AOI metrics into tidy format,
#' and optionally join to a metadata key.
#'
#' @param dir Character. Path to directory containing Imaging-PAM files.
#' @param key Data frame with metadata to join by `file`. Default = NULL.
#' @param ipam.pattern Character regex for file matching. Default = "\\.csv$".
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom rlang .data
#' @export
imPAM <- function(dir, key = NULL, ipam.pattern = "\\.csv$") {

  # ===== Argument checks =====
  if (!is.character(dir) || length(dir) != 1) {
    rlang::abort("`dir` must be a single character string specifying a directory path.")
  }

  if (!dir.exists(dir)) {
    rlang::abort(paste0("Directory does not exist: ", dir))
  }

  if (!is.character(ipam.pattern) || length(ipam.pattern) != 1) {
    rlang::abort("`ipam.pattern` must be a single character string (regular expression).")
  }

  if (!is.null(key) && !is.data.frame(key)) {
    rlang::abort("`key` must be a data.frame or NULL.")
  }

  if (!is.null(key) && !"file" %in% names(key)) {
    rlang::abort("`key` must contain a column named `file` for joining.")
  }

  # ===== Identify files =====
  ipam_files <- list.files(path = dir, pattern = ipam.pattern, full.names = TRUE)
  file_names <- list.files(path = dir, pattern = ipam.pattern, full.names = FALSE)

  if (length(ipam_files) == 0) {
    rlang::abort(
      paste0(
        "No files found in directory '", dir,
        "' matching pattern '", ipam.pattern, "'."
      )
    )
  }

  message("Found ", length(ipam_files), " Imaging-PAM files.")

  data <- NULL

  # ===== Read and process each file =====
  for (i in seq_along(ipam_files)) {

    message("Processing file ", i, "/", length(ipam_files), ": ", file_names[i])

    df <- tryCatch(
      suppressMessages(
        readr::read_delim(
          file = ipam_files[i],
          delim = ";",
          escape_double = FALSE,
          trim_ws = TRUE,
          show_col_types = FALSE,
          name_repair = "unique"
        )
      ),
      error = function(e) {
        rlang::abort(
          paste0(
            "Failed to read file: ", file_names[i],
            "\nUnderlying error: ", e$message
          )
        )
      }
    )

    # Check minimum expected structure
    if (ncol(df) < 5) {
      rlang::abort(
        paste0("File appears malformed (fewer than 5 columns): ", file_names[i])
      )
    }

    # Remove rows that are entirely NA, then keep the most recent summary row
    df <- df %>%
      dplyr::filter(!dplyr::if_all(dplyr::everything(), is.na)) %>%
      dplyr::arrange(dplyr::desc(.data$`No.`)) %>%
      dplyr::slice(1)

    if (nrow(df) == 0) {
      rlang::abort(
        paste0("No valid data rows found after filtering in file: ", file_names[i])
      )
    }

    # ===== Extract metadata safely =====
    date <- tryCatch(
      lubridate::dmy(as.character(df[[1, 1]])),
      error = function(e) {
        rlang::abort(
          paste0("Failed to parse date in file: ", file_names[i])
        )
      }
    )

    if (is.na(date)) {
      rlang::abort(
        paste0("Date could not be parsed in file: ", file_names[i])
      )
    }

    par_value <- suppressWarnings(as.numeric(df[[1, 4]]))

    if (is.na(par_value)) {
      rlang::abort(
        paste0("PAR value could not be converted to numeric in file: ", file_names[i])
      )
    }

    file <- tools::file_path_sans_ext(file_names[i])

    # ===== Reshape AOI data =====
    df2 <- df[, -c(1:4), drop = FALSE] %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "variable",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        AOI = as.numeric(stringr::str_extract(.data$variable, pattern = "[0-9]+")),
        var = stringr::str_extract(.data$variable, pattern = "[A-Za-z]+")
      )

    if (all(is.na(df2$AOI))) {
      rlang::abort(
        paste0(
          "Failed to extract AOI identifiers from column names in file: ",
          file_names[i]
        )
      )
    }

    df2 <- df2 %>%
      tidyr::pivot_wider(
        id_cols = "AOI",
        names_from = "var",
        values_from = "value"
      ) %>%
      dplyr::mutate(
        Date = date,
        PAR = par_value,
        file = file
      )

    data <- plyr::rbind.fill(data, df2)
  }

  message("Finished processing all files.")

  # ===== Optional join with key =====
  if (is.null(key)) {
    return(data)
  } else {

    message("Joining with key...")

    key %<>%
      dplyr::group_by(.data$file) %>%
      dplyr::mutate(N = dplyr::n()) %>%
      dplyr::ungroup()

    result <- data %>%
      dplyr::group_by(.data$file) %>%
      dplyr::mutate(N = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::full_join(key, by = c("file", "N")) %>%
      dplyr::select(-.data$N)

    if (nrow(result) == 0) {
      rlang::abort("Join with `key` resulted in an empty dataset. Check matching `file` values.")
    }

    message("Join complete.")

    return(result)
  }
}
