#' Import HOBO logger files and calculate daily heat stress metrics
#'
#' `hoboDHWs()` imports HOBO logger exports from a data frame, a single file,
#' multiple files, or a directory of files; standardizes common HOBO column
#' formats across export styles; optionally calculates daily heat stress metrics
#' and Degree Heating Weeks (DHWs); and returns either the full time series or a
#' daily summary table.
#'
#' The function is designed to accommodate variation among HOBO exports,
#' including:
#' \itemize{
#'   \item combined or separate date and time columns
#'   \item temperature recorded in either degrees Celsius or Fahrenheit
#'   \item light columns recorded as Light, Lux, Intensity, or lum/ft²
#'   \item CSV and Excel export formats
#'   \item directories that contain multiple files
#' }
#'
#' When `calc = TRUE`, daily hotspot values are calculated relative to a
#' bleaching threshold of `MMM + anomaly`, where hotspot magnitude is
#' calculated as `Temperature \u00B0C - MMM` for observations at or above the
#' threshold. Daily DHWs are then calculated across an explicit 84-day rolling
#' window using actual calendar dates rather than row position. Missing days are
#' completed as `NA` and propagate uncertainty through DHW calculations rather
#' than being assumed to be zero.
#'
#' If `groupingVariable` is supplied, the named column is used for summaries,
#' DHWs, and plotting without requiring that column to be literally named
#' `"Group"`. The original column name is retained in the output.
#'
#' @param path A data frame, a single file path, a character vector of file
#'   paths, or a directory containing supported HOBO files.
#' @param MMM Numeric. Maximum Monthly Mean temperature used as the bleaching
#'   baseline. Required when `calc = TRUE`. May be left as `NA` when
#'   `calc = FALSE`.
#' @param anomaly Numeric. Hotspot threshold above the MMM. Default is `1`.
#' @param calc Logical. If `TRUE`, calculate daily hotspot values and DHWs. If
#'   `FALSE`, only standardize the data and generate daily temperature summaries.
#' @param groupingVariable Optional character string naming a grouping column.
#'   This column is used directly if present in the input data. If not supplied
#'   and multiple files are imported, `Source_File` is used as the grouping
#'   variable.
#' @param summary Logical. If `TRUE`, return the daily summary table. If
#'   `FALSE`, return the full standardized time series with appended daily
#'   metrics when available.
#' @param plot Logical. If `TRUE`, generate a temperature plot and, when
#'   `calc = TRUE`, a companion DHW plot.
#' @param plotFile Optional file path for saving the plot. If `NULL`, the plot
#'   is printed but not saved.
#'
#' @return A data frame. Returns either the full standardized time series or a
#'   daily summary table, depending on `summary`.
#'
#' @examples
#' \dontrun{
#' # ===== Import one file without calculating DHWs =====
#' tmp_hobo <- hoboDHWs(
#'   path = "logger1.xlsx",
#'   calc = FALSE
#' )
#'
#' # ===== Import all supported files in a directory and calculate daily DHWs =====
#' tmp_summary <- hoboDHWs(
#'   path = "hobo_exports/",
#'   MMM = 29.5,
#'   anomaly = 1,
#'   calc = TRUE,
#'   summary = TRUE
#' )
#'
#' # ===== Import multiple files, add a grouping column, and summarize by that column =====
#' tmp_hobo <- hoboDHWs(
#'   path = c("logger1.xlsx", "logger2.xlsx"),
#'   calc = FALSE
#' )
#'
#' tmp_hobo <- tmp_hobo %>%
#'   dplyr::mutate(
#'     site = dplyr::case_when(
#'       grepl("Pier", Source_File) ~ "RSMAS Pier",
#'       grepl("DiLido", Source_File) ~ "DiLido",
#'       TRUE ~ "Other"
#'     )
#'   )
#'
#' tmp_site_summary <- hoboDHWs(
#'   path = tmp_hobo,
#'   MMM = 29.5,
#'   anomaly = 1,
#'   calc = TRUE,
#'   groupingVariable = "site",
#'   summary = TRUE,
#'   plot = TRUE
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import patchwork
#' @export
hoboDHWs <- function(path,
                     MMM = NA,
                     anomaly = 1,
                     calc = TRUE,
                     groupingVariable = NA,
                     summary = FALSE,
                     plot = FALSE,
                     plotFile = NULL) {

  temp_col_name <- "Temperature \u00B0C"

  # ===== Validate inputs =====
  if (!is.logical(calc) || length(calc) != 1 || is.na(calc)) {
    stop("`calc` must be TRUE or FALSE.")
  }

  if (calc) {
    if (!is.numeric(MMM) || length(MMM) != 1 || is.na(MMM)) {
      stop("`MMM` must be a single non-missing numeric value when `calc = TRUE`.")
    }
  } else {
    if (!(is.numeric(MMM) && length(MMM) == 1) && !is.na(MMM)) {
      stop("`MMM` must be NA or a single numeric value when `calc = FALSE`.")
    }
  }

  if (!is.numeric(anomaly) || length(anomaly) != 1 || is.na(anomaly)) {
    stop("`anomaly` must be a single non-missing numeric value.")
  }

  if (!is.logical(summary) || length(summary) != 1 || is.na(summary)) {
    stop("`summary` must be TRUE or FALSE.")
  }

  if (!is.logical(plot) || length(plot) != 1 || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.")
  }

  if (!is.na(groupingVariable) &&
      (!is.character(groupingVariable) || length(groupingVariable) != 1)) {
    stop("`groupingVariable` must be NA or a single character string.")
  }

  if (!is.null(plotFile) &&
      (!is.character(plotFile) || length(plotFile) != 1 || !nzchar(plotFile))) {
    stop("`plotFile` must be NULL or a single non-empty character string.")
  }

  # ===== Return the first non-missing value in a vector =====
  safe_first_non_na <- function(x) {
    x <- x[!is.na(x)]

    if (length(x) == 0) {
      return(NA)
    }

    x[[1]]
  }

  # ===== Extract logger serial number from header-like text =====
  extract_serial_number <- function(x) {
    if (length(x) == 0 || all(is.na(x))) {
      return(NA_character_)
    }

    tmp_match <- stringr::str_match(
      paste(x, collapse = " | "),
      "(?:LGR\\s*S/N:|SEN\\s*S/N:|Serial\\s*Number:?|S/N:?)[[:space:]]*([0-9]+)"
    )

    tmp_match[, 2]
  }

  # ===== Parse HOBO datetimes across common export styles =====
  parse_hobo_datetime <- function(x) {
    tmp_raw_examples <- utils::head(unique(x[!is.na(x)]), 5)

    try_return <- function(value) {
      if (sum(!is.na(value)) > 0) {
        return(value)
      }

      NULL
    }

    if (inherits(x, c("POSIXct", "POSIXt"))) {
      return(as.POSIXct(x, tz = "UTC"))
    }

    if (inherits(x, "Date")) {
      return(as.POSIXct(x, tz = "UTC"))
    }

    if (is.factor(x)) {
      x <- as.character(x)
    }

    if (is.numeric(x)) {
      tmp_out <- try_return(
        as.POSIXct(x * 86400, origin = "1899-12-30", tz = "UTC")
      )

      if (!is.null(tmp_out)) {
        return(tmp_out)
      }

      message(
        "Failed to parse numeric DateTime values. Example raw values: ",
        paste(tmp_raw_examples, collapse = ", ")
      )

      stop(
        "DateTime values could not be parsed from numeric input. ",
        "Check whether the source column contains valid Excel date-time values."
      )
    }

    x <- as.character(x)
    x <- stringr::str_squish(x)
    x[x == ""] <- NA_character_

    tmp_numeric <- suppressWarnings(as.numeric(x))
    tmp_out <- try_return(
      as.POSIXct(tmp_numeric * 86400, origin = "1899-12-30", tz = "UTC")
    )

    if (!is.null(tmp_out)) {
      return(tmp_out)
    }

    tmp_formats <- c(
      "%Y-%m-%d %H:%M:%S",
      "%Y-%m-%d %H:%M",
      "%m/%d/%y %I:%M:%S %p",
      "%m/%d/%Y %I:%M:%S %p",
      "%m/%d/%y %I:%M %p",
      "%m/%d/%Y %I:%M %p",
      "%m/%d/%y %H:%M:%S",
      "%m/%d/%Y %H:%M:%S",
      "%m/%d/%y %H:%M",
      "%m/%d/%Y %H:%M"
    )

    for (tmp_format in tmp_formats) {
      tmp_out <- try_return(as.POSIXct(x, format = tmp_format, tz = "UTC"))

      if (!is.null(tmp_out)) {
        return(tmp_out)
      }
    }

    message(
      "Failed to parse DateTime values. Example raw values: ",
      paste(tmp_raw_examples, collapse = ", ")
    )

    stop(
      "DateTime values could not be parsed. ",
      "Check whether the source column contains a supported datetime format."
    )
  }

  # ===== Read one HOBO file and recover header metadata =====
  read_hobo_file <- function(file_path) {
    tmp_ext <- tolower(tools::file_ext(file_path))

    if (!tmp_ext %in% c("csv", "xlsx", "xls")) {
      stop("Unsupported file type: ", basename(file_path))
    }

    if (tmp_ext == "csv") {
      tmp_first_line <- readr::read_lines(file_path, n_max = 1)

      if (length(tmp_first_line) == 0) {
        stop("File appears to be empty: ", basename(file_path))
      }

      tmp_plot_title <- ifelse(
        grepl("Plot Title", tmp_first_line, ignore.case = TRUE),
        stringr::str_remove(tmp_first_line, '^"?Plot Title:\\s*'),
        NA_character_
      )

      tmp_skip <- ifelse(
        grepl("Plot Title", tmp_first_line, ignore.case = TRUE),
        1,
        0
      )

      tmp_df <- suppressMessages(
        readr::read_csv(
          file = file_path,
          skip = tmp_skip,
          show_col_types = FALSE,
          progress = FALSE
        )
      ) %>%
        as.data.frame(stringsAsFactors = FALSE)

      return(list(
        data = tmp_df,
        logger_sn = extract_serial_number(c(names(tmp_df), tmp_first_line)),
        plot_title = tmp_plot_title
      ))
    }

    tmp_sheet <- readxl::excel_sheets(file_path)[1]

    tmp_df <- suppressMessages(
      readxl::read_excel(
        path = file_path,
        sheet = tmp_sheet,
        col_names = FALSE
      )
    ) %>%
      as.data.frame(stringsAsFactors = FALSE)

    if (nrow(tmp_df) < 2) {
      stop("Excel file does not contain enough rows to identify headers: ", basename(file_path))
    }

    tmp_first_cell <- as.character(tmp_df[1, 1])

    if (!is.na(tmp_first_cell) && grepl("Plot Title", tmp_first_cell, ignore.case = TRUE)) {
      tmp_header_row <- 2
      tmp_plot_title <- tmp_first_cell
      tmp_logger_sn <- extract_serial_number(unlist(tmp_df[tmp_header_row, ]))
      names(tmp_df) <- as.character(unlist(tmp_df[tmp_header_row, ]))
      tmp_df <- tmp_df[-c(1, 2), , drop = FALSE]
    } else {
      tmp_header_row <- 1
      tmp_plot_title <- NA_character_
      tmp_logger_sn <- extract_serial_number(unlist(tmp_df[tmp_header_row, ]))
      names(tmp_df) <- as.character(unlist(tmp_df[tmp_header_row, ]))
      tmp_df <- tmp_df[-1, , drop = FALSE]
    }

    tmp_df <- janitor::remove_empty(tmp_df, which = c("rows", "cols"))

    list(
      data = tmp_df,
      logger_sn = tmp_logger_sn,
      plot_title = tmp_plot_title
    )
  }

  # ===== Standardize input columns to a common HOBO schema =====
  standardize_hobo <- function(hoboFile, groupingVariable = NA) {
    tmp_names <- names(hoboFile)

    if (is.null(tmp_names) || length(tmp_names) == 0) {
      stop("Imported data do not contain column names.")
    }

    tmp_map <- tibble::tibble(
      original_name = tmp_names
    ) %>%
      dplyr::mutate(
        function_name = dplyr::case_when(
          grepl("date", .data$original_name, ignore.case = TRUE) &
            grepl("time", .data$original_name, ignore.case = TRUE) ~ "DateTime",
          grepl("date", .data$original_name, ignore.case = TRUE) &
            !grepl("time", .data$original_name, ignore.case = TRUE) ~ "Date",
          !grepl("date", .data$original_name, ignore.case = TRUE) &
            grepl("time", .data$original_name, ignore.case = TRUE) ~ "Time",
          grepl("temp", .data$original_name, ignore.case = TRUE) ~ "Temperature",
          grepl("lux", .data$original_name, ignore.case = TRUE) |
            grepl("light", .data$original_name, ignore.case = TRUE) |
            grepl("lum", .data$original_name, ignore.case = TRUE) ~ "Light",
          !is.na(groupingVariable) & .data$original_name == groupingVariable ~ "GroupingVariable",
          TRUE ~ NA_character_
        ),
        units = dplyr::case_when(
          .data$function_name == "Temperature" &
            grepl("f", .data$original_name, ignore.case = TRUE) ~ "F",
          .data$function_name == "Temperature" &
            grepl("c", .data$original_name, ignore.case = TRUE) ~ "C",
          .data$function_name == "Light" &
            grepl("lux", .data$original_name, ignore.case = TRUE) ~ "lux",
          .data$function_name == "Light" &
            grepl("lum", .data$original_name, ignore.case = TRUE) ~ "lum",
          TRUE ~ NA_character_
        )
      )

    tmp_pick <- function(label) {
      tmp_match <- tmp_map %>%
        dplyr::filter(.data$function_name == label)

      if (nrow(tmp_match) == 0) {
        return(NA_character_)
      }

      tmp_match$original_name[[1]]
    }

    tmp_datetime_col <- tmp_pick("DateTime")
    tmp_date_col <- tmp_pick("Date")
    tmp_time_col <- tmp_pick("Time")
    tmp_temp_col <- tmp_pick("Temperature")
    tmp_light_col <- tmp_pick("Light")
    tmp_group_col <- tmp_pick("GroupingVariable")

    if (is.na(tmp_temp_col)) {
      stop("No temperature column could be identified.")
    }

    if (is.na(tmp_datetime_col) && is.na(tmp_date_col)) {
      stop("No date or datetime column could be identified.")
    }

    if (!is.na(tmp_datetime_col)) {
      tmp_datetime <- parse_hobo_datetime(hoboFile[[tmp_datetime_col]])
    } else if (!is.na(tmp_date_col) && !is.na(tmp_time_col)) {
      tmp_datetime <- parse_hobo_datetime(
        paste(hoboFile[[tmp_date_col]], hoboFile[[tmp_time_col]])
      )
    } else {
      tmp_datetime <- parse_hobo_datetime(hoboFile[[tmp_date_col]])
    }

    tmp_units <- tmp_map %>%
      dplyr::filter(.data$original_name == tmp_temp_col) %>%
      dplyr::pull(.data$units)

    tmp_units <- tmp_units[!is.na(tmp_units)][1]

    tmp_out <- tibble::tibble(
      DateTime = tmp_datetime
    )

    tmp_out[[temp_col_name]] <- suppressWarnings(as.numeric(hoboFile[[tmp_temp_col]]))

    if (!is.na(tmp_light_col)) {
      tmp_out$Light <- suppressWarnings(as.numeric(hoboFile[[tmp_light_col]]))
    }

    if (!is.na(tmp_group_col)) {
      tmp_out[[groupingVariable]] <- as.character(hoboFile[[tmp_group_col]])
    }

    if (all(is.na(tmp_out[[temp_col_name]]))) {
      stop("Temperature column could not be converted to numeric.")
    }

    if (is.na(tmp_units)) {
      tmp_units <- ifelse(mean(tmp_out[[temp_col_name]], na.rm = TRUE) > 45, "F", "C")
    }

    if (!is.na(tmp_units) && tmp_units == "F") {
      tmp_out <- tmp_out %>%
        dplyr::mutate(
          !!temp_col_name := weathermetrics::fahrenheit.to.celsius(.data[[temp_col_name]])
        )
    }

    tmp_out %>%
      dplyr::filter(!is.na(.data$DateTime), !is.na(.data[[temp_col_name]]))
  }

  # ===== Decide which column should act as the grouping variable =====
  determine_group_column <- function(df, groupingVariable = NA) {
    if (!is.na(groupingVariable)) {
      if (!groupingVariable %in% names(df)) {
        stop(
          "`groupingVariable = '",
          groupingVariable,
          "'` was supplied, but that column was not found in the data."
        )
      }

      return(groupingVariable)
    }

    if ("Source_File" %in% names(df) &&
        dplyr::n_distinct(df$Source_File, na.rm = TRUE) > 1) {
      return("Source_File")
    }

    NA_character_
  }

  # ===== Calculate rolling DHWs using actual dates and an 84-day window =====
  calc_dhw_by_date <- function(dates, dhDay) {
    vapply(
      seq_along(dates),
      FUN.VALUE = numeric(1),
      FUN = function(i) {
        tmp_window_start <- dates[i] - 83
        tmp_window_idx <- dates >= tmp_window_start & dates <= dates[i]
        tmp_window_values <- dhDay[tmp_window_idx]

        if (any(is.na(tmp_window_values))) {
          return(NA_real_)
        }

        sum(tmp_window_values) / 7
      }
    )
  }

  # ===== Summarize daily rows without assuming Logger_SN exists =====
  summarise_daily_core <- function(df_daily) {
    tmp_has_logger <- "Logger_SN" %in% names(df_daily)
    tmp_has_hotspot_weighted <- "hotspot_weighted" %in% names(df_daily)

    tmp_result <- df_daily %>%
      dplyr::group_by(.data$.group_internal, .data$Date) %>%
      dplyr::summarise(
        Source_File = safe_first_non_na(.data$Source_File),
        Temperature_Average = mean(.data[[temp_col_name]], na.rm = TRUE),
        Temperature_StDev = stats::sd(.data[[temp_col_name]], na.rm = TRUE),
        Temperature_Min = min(.data[[temp_col_name]], na.rm = TRUE),
        Temperature_Max = max(.data[[temp_col_name]], na.rm = TRUE),
        N_Records = dplyr::n(),
        dhDay = if (tmp_has_hotspot_weighted) {
          if (all(is.na(.data$hotspot_weighted))) {
            NA_real_
          } else {
            sum(.data$hotspot_weighted, na.rm = TRUE)
          }
        } else {
          NA_real_
        },
        .groups = "drop"
      )

    if (tmp_has_logger) {
      tmp_logger <- df_daily %>%
        dplyr::group_by(.data$.group_internal, .data$Date) %>%
        dplyr::summarise(
          Logger_SN = safe_first_non_na(.data$Logger_SN),
          .groups = "drop"
        )

      tmp_result <- tmp_result %>%
        dplyr::left_join(tmp_logger, by = c(".group_internal", "Date")) %>%
        dplyr::relocate(dplyr::any_of("Logger_SN"), .after = .data$Date)
    }

    tmp_result
  }

  # ===== Calculate daily temperature summaries, hotspots, and DHWs =====
  calculate_daily_metrics <- function(df, MMM, anomaly, group_col = NA_character_) {
    if (is.na(group_col)) {
      df <- df %>%
        dplyr::mutate(.group_internal = "All_Data")
    } else {
      df <- df %>%
        dplyr::mutate(.group_internal = as.character(.data[[group_col]]))
    }

    tmp_daily_raw <- df %>%
      dplyr::arrange(.data$.group_internal, .data$DateTime) %>%
      dplyr::group_by(.data$.group_internal) %>%
      dplyr::mutate(
        Date = as.Date(.data$DateTime),
        next_time = dplyr::lead(.data$DateTime),
        interval_minutes = as.numeric(
          difftime(.data$next_time, .data$DateTime, units = "mins")
        )
      ) %>%
      dplyr::group_by(.data$.group_internal, .data$Date, .add = FALSE) %>%
      dplyr::mutate(
        interval_minutes = dplyr::if_else(
          is.na(.data$interval_minutes) | .data$interval_minutes <= 0,
          stats::median(
            .data$interval_minutes[.data$interval_minutes > 0],
            na.rm = TRUE
          ),
          .data$interval_minutes
        ),
        interval_minutes = dplyr::if_else(
          is.na(.data$interval_minutes) | !is.finite(.data$interval_minutes),
          NA_real_,
          .data$interval_minutes
        ),
        hotspot = dplyr::if_else(
          .data[[temp_col_name]] >= (MMM + anomaly),
          .data[[temp_col_name]] - MMM,
          0
        ),
        hotspot_weighted = .data$hotspot * (.data$interval_minutes / (24 * 60))
      ) %>%
      dplyr::ungroup()

    tmp_daily_summary <- summarise_daily_core(tmp_daily_raw) %>%
      dplyr::group_by(.data$.group_internal) %>%
      tidyr::complete(
        Date = seq(min(.data$Date), max(.data$Date), by = "day")
      ) %>%
      dplyr::arrange(.data$Date, .by_group = TRUE) %>%
      tidyr::fill(dplyr::any_of(c("Logger_SN", "Source_File")), .direction = "downup") %>%
      dplyr::mutate(
        DHWs = calc_dhw_by_date(.data$Date, .data$dhDay)
      ) %>%
      dplyr::ungroup()

    if (is.na(group_col)) {
      tmp_daily_summary_out <- tmp_daily_summary %>%
        dplyr::select(-.data$.group_internal)
    } else {
      tmp_daily_summary_out <- tmp_daily_summary
      tmp_daily_summary_out[[group_col]] <- tmp_daily_summary_out$.group_internal
      tmp_daily_summary_out <- tmp_daily_summary_out %>%
        dplyr::select(-.data$.group_internal)
    }

    if (is.na(group_col)) {
      tmp_complete <- tmp_daily_raw %>%
        dplyr::left_join(
          tmp_daily_summary_out %>%
            dplyr::select(.data$Date, .data$dhDay, .data$DHWs),
          by = "Date"
        ) %>%
        dplyr::select(
          -dplyr::any_of(c(
            "next_time",
            "Date",
            "interval_minutes",
            "hotspot",
            "hotspot_weighted",
            ".group_internal"
          ))
        )
    } else {
      tmp_complete <- tmp_daily_raw %>%
        dplyr::left_join(
          tmp_daily_summary_out %>%
            dplyr::select(dplyr::all_of(group_col), .data$Date, .data$dhDay, .data$DHWs),
          by = c(group_col, "Date")
        ) %>%
        dplyr::select(
          -dplyr::any_of(c(
            "next_time",
            "Date",
            "interval_minutes",
            "hotspot",
            "hotspot_weighted",
            ".group_internal"
          ))
        )
    }

    tmp_complete <- tmp_complete %>%
      dplyr::select(-dplyr::any_of("Logger_SN"))

    list(
      full = tmp_complete,
      summary = tmp_daily_summary_out
    )
  }

  # ===== Summarize daily temperatures without DHW calculations =====
  summarize_daily_only <- function(df, group_col = NA_character_) {
    if (is.na(group_col)) {
      df <- df %>%
        dplyr::mutate(.group_internal = "All_Data")
    } else {
      df <- df %>%
        dplyr::mutate(.group_internal = as.character(.data[[group_col]]))
    }

    tmp_summary <- df %>%
      dplyr::mutate(Date = as.Date(.data$DateTime)) %>%
      summarise_daily_core()

    if (is.na(group_col)) {
      tmp_summary_out <- tmp_summary %>%
        dplyr::select(-.data$.group_internal)
    } else {
      tmp_summary_out <- tmp_summary
      tmp_summary_out[[group_col]] <- tmp_summary_out$.group_internal
      tmp_summary_out <- tmp_summary_out %>%
        dplyr::select(-.data$.group_internal)
    }

    tmp_complete <- df %>%
      dplyr::select(-dplyr::any_of(c("Logger_SN", ".group_internal")))

    list(
      full = tmp_complete,
      summary = tmp_summary_out
    )
  }

  # ===== Import and standardize data =====
  if (is.data.frame(path)) {
    hoboFile <- standardize_hobo(
      hoboFile = path,
      groupingVariable = groupingVariable
    )

    tmp_extra_cols <- setdiff(names(path), names(hoboFile))

    if (length(tmp_extra_cols) > 0) {
      hoboFile <- dplyr::bind_cols(
        hoboFile,
        path %>% dplyr::select(dplyr::all_of(tmp_extra_cols))
      )
    }

  } else if (length(path) == 1 && dir.exists(path)) {
    tmp_files <- list.files(
      path,
      pattern = "\\.(csv|xlsx|xls)$",
      full.names = TRUE
    )

    if (length(tmp_files) == 0) {
      stop("No supported HOBO files were found in the supplied directory.")
    }

    tmp_list <- lapply(tmp_files, function(tmp_file) {
      tmp_raw <- read_hobo_file(tmp_file)

      tmp_std <- standardize_hobo(
        hoboFile = tmp_raw$data,
        groupingVariable = groupingVariable
      )

      tmp_std$Source_File <- basename(tmp_file)

      if (!is.na(tmp_raw$logger_sn)) {
        tmp_std$Logger_SN <- tmp_raw$logger_sn
      }

      tmp_std
    })

    hoboFile <- dplyr::bind_rows(tmp_list) %>%
      dplyr::distinct()

  } else if (length(path) == 1 && file.exists(path)) {
    tmp_raw <- read_hobo_file(path)

    hoboFile <- standardize_hobo(
      hoboFile = tmp_raw$data,
      groupingVariable = groupingVariable
    )

    hoboFile$Source_File <- basename(path)

    if (!is.na(tmp_raw$logger_sn)) {
      hoboFile$Logger_SN <- tmp_raw$logger_sn
    }

  } else if (is.character(path) && length(path) > 1) {
    tmp_files <- path[file.exists(path)]

    if (length(tmp_files) == 0) {
      stop("None of the supplied file paths could be found.")
    }

    tmp_list <- lapply(tmp_files, function(tmp_file) {
      tmp_raw <- read_hobo_file(tmp_file)

      tmp_std <- standardize_hobo(
        hoboFile = tmp_raw$data,
        groupingVariable = groupingVariable
      )

      tmp_std$Source_File <- basename(tmp_file)

      if (!is.na(tmp_raw$logger_sn)) {
        tmp_std$Logger_SN <- tmp_raw$logger_sn
      }

      tmp_std
    })

    hoboFile <- dplyr::bind_rows(tmp_list) %>%
      dplyr::distinct()

  } else {
    stop("Import failed.")
  }

  if (nrow(hoboFile) == 0) {
    stop("No valid HOBO observations remained after import and cleaning.")
  }

  # ===== Determine the grouping column to use downstream =====
  tmp_group_col <- determine_group_column(
    df = hoboFile,
    groupingVariable = groupingVariable
  )

  # ===== Calculate daily metrics =====
  if (calc) {
    tmp_metrics <- calculate_daily_metrics(
      df = hoboFile,
      MMM = MMM,
      anomaly = anomaly,
      group_col = tmp_group_col
    )

    hoboComplete <- tmp_metrics$full
    hoboSummary <- tmp_metrics$summary
  } else {
    warning("`calc = FALSE`, so DHW metrics were not calculated.")

    tmp_metrics <- summarize_daily_only(
      df = hoboFile,
      group_col = tmp_group_col
    )

    hoboComplete <- tmp_metrics$full
    hoboSummary <- tmp_metrics$summary
  }

  # ===== Plot output =====
  if (plot) {
    tmp_facet_formula <- if (!is.na(tmp_group_col)) {
      stats::as.formula(paste("~", tmp_group_col))
    } else {
      NULL
    }

    p1 <- ggplot2::ggplot(
      hoboComplete,
      ggplot2::aes(
        x = .data$DateTime,
        y = .data[[temp_col_name]]
      )
    ) +
      ggplot2::geom_line(color = "#2C77B8", alpha = 0.8) +
      ggplot2::labs(
        y = temp_col_name,
        x = NULL
      ) +
      ggthemes::theme_few(base_size = 13)

    if (calc) {
      p1 <- p1 +
        ggplot2::geom_hline(
          yintercept = MMM + anomaly,
          color = "darkred",
          linetype = "dashed",
          linewidth = 0.6
        )
    }

    if (!is.null(tmp_facet_formula)) {
      p1 <- p1 + ggplot2::facet_wrap(tmp_facet_formula, scales = "free_x")
    }

    if (calc) {
      p2 <- ggplot2::ggplot(
        hoboSummary,
        ggplot2::aes(
          x = .data$Date,
          y = .data$DHWs
        )
      ) +
        ggplot2::geom_line(color = "#F05D5E", linewidth = 0.8, na.rm = TRUE) +
        ggplot2::labs(
          title = "Degree Heating Weeks (DHWs)",
          y = "DHWs",
          x = "Date"
        ) +
        ggthemes::theme_few(base_size = 13)

      if (!is.null(tmp_facet_formula)) {
        p2 <- p2 + ggplot2::facet_wrap(tmp_facet_formula, scales = "free_x")
      }

      full_plot <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(2, 1))
    } else {
      full_plot <- p1
    }

    if (is.null(plotFile)) {
      print(full_plot)
    } else {
      ggplot2::ggsave(
        filename = plotFile,
        plot = full_plot,
        width = 10,
        height = if (calc) 15 else 10,
        dpi = 300
      )

      message("Plot saved to: ", normalizePath(plotFile))
    }
  }

  # ===== Return output =====
  if (summary) {
    return(hoboSummary)
  } else {
    return(hoboComplete)
  }
}
