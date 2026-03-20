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
#'   \item directories that contain multiple files or nested folders of files
#' }
#'
#' Provenance metadata are retained in standardized columns, including
#' `Source_File`, `Source_Series`, `Logger_SN`, and `Plot_Title`, so imported
#' data can be cleaned and rerun through the function later without depending on
#' a literal column named `Group`.
#'
#' When `calc = TRUE`, daily hotspot values are calculated relative to a
#' bleaching threshold of `MMM + anomaly`, where the hotspot magnitude is
#' calculated as `Temperature - MMM` for observations at or above the threshold.
#' Daily DHWs are then calculated across an explicit 84-day rolling window using
#' actual calendar dates rather than row position. Missing days are completed as
#' `NA` and propagate uncertainty through DHW calculations rather than being
#' assumed to be zero.
#'
#' @param path A data frame, a single file path, a character vector of file
#'   paths, or a directory containing supported HOBO files.
#' @param MMM Numeric. Maximum Monthly Mean temperature used as the bleaching
#'   baseline.
#' @param anomaly Numeric. Hotspot threshold above the MMM. Default is `1`,
#'   corresponding to the common MMM + 1 °C threshold.
#' @param calc Logical. If `TRUE`, calculate daily hotspot values and DHWs. If
#'   `FALSE`, only standardize the data and generate daily temperature summaries.
#' @param groupingVariable Optional character string naming a grouping column to
#'   use for daily summaries, DHWs, and plotting. If not supplied, the function
#'   uses `Source_Series` when available; otherwise all data are treated as one
#'   series.
#' @param recursive Logical. If `TRUE` and `path` is a directory, search
#'   subdirectories recursively for supported files.
#' @param summary Logical. If `TRUE`, return the daily summary table. If
#'   `FALSE`, return the full standardized time series with appended daily
#'   metrics when available.
#' @param plot Logical. If `TRUE`, generate a simple temperature plot and, when
#'   `calc = TRUE`, a companion DHW plot.
#' @param plotFile Optional file path for saving the plot. If `NULL`, the plot
#'   is printed but not saved.
#'
#' @return A data frame. Returns either the full standardized time series or a
#'   daily summary table, depending on `summary`.
#'
#' @examples
#' \dontrun{
#' tmp_hobo <- hoboDHWs(
#'   path = c("logger1.csv", "logger2.csv"),
#'   MMM = 29.5,
#'   calc = FALSE
#' )
#'
#' tmp_hobo <- tmp_hobo %>%
#'   dplyr::mutate(site = dplyr::case_when(
#'     grepl("Pier", Source_File) ~ "RSMAS Pier",
#'     grepl("DiLido", Source_File) ~ "DiLido",
#'     TRUE ~ "Other"
#'   ))
#'
#' tmp_summary <- hoboDHWs(
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
                     MMM,
                     anomaly = 1,
                     calc = TRUE,
                     groupingVariable = NA,
                     recursive = TRUE,
                     summary = FALSE,
                     plot = FALSE,
                     plotFile = NULL) {

  # ===== Validate Inputs =====
  if (!is.numeric(MMM) || length(MMM) != 1 || is.na(MMM)) {
    stop("`MMM` must be a single non-missing numeric value.")
  }

  if (!is.numeric(anomaly) || length(anomaly) != 1 || is.na(anomaly)) {
    stop("`anomaly` must be a single non-missing numeric value.")
  }

  if (!is.logical(calc) || length(calc) != 1 || is.na(calc)) {
    stop("`calc` must be TRUE or FALSE.")
  }

  if (!is.logical(recursive) || length(recursive) != 1 || is.na(recursive)) {
    stop("`recursive` must be TRUE or FALSE.")
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

  # ===== Helper Functions =====
  safe_first_non_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA)
    }
    x[[1]]
  }

  extract_serial_number <- function(x) {
    if (length(x) == 0 || all(is.na(x))) {
      return(NA_character_)
    }

    tmp_match <- stringr::str_match(
      string = paste(x, collapse = " | "),
      pattern = "(?:LGR\\s*S/N:|SEN\\s*S/N:|Serial\\s*Number:?|S/N:?)[[:space:]]*([0-9]+)"
    )

    tmp_match[, 2]
  }

  clean_hobo_excel <- function(df_raw) {
    df_raw <- as.data.frame(df_raw, stringsAsFactors = FALSE)

    if (nrow(df_raw) < 2) {
      stop("Excel file does not contain enough rows to identify headers.")
    }

    tmp_first_cell <- as.character(df_raw[1, 1])

    if (!is.na(tmp_first_cell) && grepl("Plot Title", tmp_first_cell, ignore.case = TRUE)) {
      tmp_header_row <- 2
      tmp_serial <- extract_serial_number(unlist(df_raw[tmp_header_row, ]))
      tmp_plot_title <- tmp_first_cell

      names(df_raw) <- as.character(unlist(df_raw[tmp_header_row, ]))
      df_raw <- df_raw[-c(1, 2), , drop = FALSE]
    } else {
      tmp_header_row <- 1
      tmp_serial <- extract_serial_number(unlist(df_raw[tmp_header_row, ]))
      tmp_plot_title <- NA_character_

      names(df_raw) <- as.character(unlist(df_raw[tmp_header_row, ]))
      df_raw <- df_raw[-1, , drop = FALSE]
    }

    df_raw <- janitor::remove_empty(df_raw, which = c("rows", "cols"))

    list(
      data = df_raw,
      logger_sn = tmp_serial,
      plot_title = tmp_plot_title
    )
  }

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

      tmp_skip <- ifelse(grepl("Plot Title", tmp_first_line, ignore.case = TRUE), 1, 0)

      tmp_df <- suppressMessages(
        readr::read_csv(
          file = file_path,
          skip = tmp_skip,
          show_col_types = FALSE,
          progress = FALSE
        )
      ) %>%
        as.data.frame(stringsAsFactors = FALSE)

      tmp_serial <- extract_serial_number(c(names(tmp_df), tmp_first_line))

      return(list(
        data = tmp_df,
        logger_sn = tmp_serial,
        plot_title = tmp_plot_title
      ))
    }

    if (tmp_ext %in% c("xlsx", "xls")) {
      tmp_sheet <- readxl::excel_sheets(file_path)[1]

      tmp_df <- suppressMessages(
        readxl::read_excel(
          path = file_path,
          sheet = tmp_sheet,
          col_names = FALSE
        )
      )

      return(clean_hobo_excel(tmp_df))
    }
  }

  parse_datetime_vector <- function(x) {

    if (inherits(x, c("POSIXct", "POSIXlt", "POSIXt"))) {
      return(as.POSIXct(x, tz = "UTC"))
    }

    if (inherits(x, "Date")) {
      return(as.POSIXct(x, tz = "UTC"))
    }

    if (is.numeric(x)) {
      tmp_dt <- suppressWarnings(as.POSIXct(x, origin = "1899-12-30", tz = "UTC"))
      if (!all(is.na(tmp_dt))) {
        return(tmp_dt)
      }
    }

    x <- as.character(x)

    tmp_parsed <- suppressWarnings(
      lubridate::parse_date_time(
        x,
        orders = c(
          "mdy HMS p",
          "mdy HM p",
          "mdy IMS p",
          "mdy IM p",
          "mdy HMS",
          "mdy HM",
          "ymd HMS",
          "ymd HM",
          "Ymd HMS",
          "Ymd HM"
        ),
        tz = "UTC"
      )
    )

    if (all(is.na(tmp_parsed))) {
      stop(
        "DateTime values could not be parsed. ",
        "Check whether the source column contains a supported datetime format."
      )
    }

    as.POSIXct(tmp_parsed, tz = "UTC")
  }

  standardize_hobo_data <- function(hoboFile,
                                    logger_sn = NA_character_,
                                    plot_title = NA_character_) {

    tmp_original_names <- names(hoboFile)

    if (is.null(tmp_original_names) || length(tmp_original_names) == 0) {
      stop("Imported data do not contain column names.")
    }

    tmp_name_map <- tibble::tibble(
      original_name = tmp_original_names
    ) %>%
      dplyr::mutate(
        clean_name = stringr::str_squish(original_name),
        function_name = dplyr::case_when(
          grepl("date[ -]?time|datetime", clean_name, ignore.case = TRUE) ~ "DateTime",
          grepl("^date$|date", clean_name, ignore.case = TRUE) &
            !grepl("time", clean_name, ignore.case = TRUE) ~ "Date",
          grepl("^time$|time", clean_name, ignore.case = TRUE) &
            !grepl("date", clean_name, ignore.case = TRUE) ~ "Time",
          grepl("temp|temperature", clean_name, ignore.case = TRUE) ~ "Temperature",
          grepl("lux|light|intensity|lum", clean_name, ignore.case = TRUE) ~ "Light",
          TRUE ~ NA_character_
        ),
        units = dplyr::case_when(
          function_name == "Temperature" &
            grepl("°f|\\(f\\)|fahrenheit|temp,\\s*f", clean_name, ignore.case = TRUE) ~ "°F",
          function_name == "Temperature" &
            grepl("°c|\\(c\\)|celsius|temp,\\s*c", clean_name, ignore.case = TRUE) ~ "°C",
          function_name == "Light" &
            grepl("lux", clean_name, ignore.case = TRUE) ~ "lux",
          function_name == "Light" &
            grepl("lum", clean_name, ignore.case = TRUE) ~ "lum/ft²",
          TRUE ~ NA_character_
        )
      )

    if (!any(tmp_name_map$function_name == "Temperature")) {
      stop("No temperature column could be identified.")
    }

    if (!any(tmp_name_map$function_name %in% c("DateTime", "Date"))) {
      stop("No date or datetime column could be identified.")
    }

    tmp_keep <- tmp_name_map %>%
      dplyr::filter(!is.na(function_name))

    hoboFile <- hoboFile[, tmp_keep$original_name, drop = FALSE]
    names(hoboFile) <- tmp_keep$function_name

    if ("DateTime" %in% names(hoboFile)) {
      hoboFile <- hoboFile %>%
        dplyr::mutate(DateTime = parse_datetime_vector(DateTime))
    } else if (all(c("Date", "Time") %in% names(hoboFile))) {
      hoboFile <- hoboFile %>%
        dplyr::mutate(DateTime = parse_datetime_vector(paste(Date, Time)))
    } else if ("Date" %in% names(hoboFile)) {
      hoboFile <- hoboFile %>%
        dplyr::mutate(DateTime = parse_datetime_vector(Date))
    }

    if (all(is.na(hoboFile$DateTime))) {
      stop("No DateTime values could be parsed.")
    }

    hoboFile <- hoboFile %>%
      dplyr::filter(!is.na(DateTime))

    tmp_temp_units <- tmp_name_map %>%
      dplyr::filter(function_name == "Temperature") %>%
      dplyr::pull(units)

    tmp_temp_units <- tmp_temp_units[!is.na(tmp_temp_units)][1]

    if (is.na(tmp_temp_units)) {
      tmp_temp_values <- suppressWarnings(as.numeric(hoboFile$Temperature))

      if (all(is.na(tmp_temp_values))) {
        stop("Temperature column could not be converted to numeric.")
      }

      tmp_temp_units <- ifelse(mean(tmp_temp_values, na.rm = TRUE) > 45, "°F", "°C")
      message("Temperature units were not explicit in the source header. Inferred units as ", tmp_temp_units, ".")
    }

    hoboFile <- hoboFile %>%
      dplyr::mutate(
        Temperature = suppressWarnings(as.numeric(Temperature)),
        Temperature = dplyr::if_else(
          tmp_temp_units == "°F",
          weathermetrics::fahrenheit.to.celsius(Temperature),
          Temperature
        )
      ) %>%
      dplyr::filter(!is.na(Temperature)) %>%
      dplyr::rename(`Temperature °C` = Temperature)

    if ("Light" %in% names(hoboFile)) {
      hoboFile <- hoboFile %>%
        dplyr::mutate(Light = suppressWarnings(as.numeric(Light)))
    }

    hoboFile %>%
      dplyr::mutate(
        Logger_SN = logger_sn,
        Plot_Title = plot_title
      ) %>%
      dplyr::select(dplyr::any_of(c(
        "DateTime",
        "Temperature °C",
        "Light",
        "Logger_SN",
        "Plot_Title"
      )))
  }

  infer_file_groups <- function(path, recursive = TRUE) {
    if (length(path) > 1) {
      tmp_files <- path[file.exists(path)]

      if (length(tmp_files) == 0) {
        stop("None of the supplied file paths could be found.")
      }

      tmp_index <- tibble::tibble(
        file_path = tmp_files,
        file_name = basename(tmp_files),
        folder_name = basename(dirname(tmp_files)),
        file_stem = tools::file_path_sans_ext(file_name),
        Source_Series = file_stem,
        grouping_rule = "vector_of_files_each_file_is_series"
      )

      return(tmp_index)
    }

    if (length(path) == 1 && file.exists(path)) {
      tmp_index <- tibble::tibble(
        file_path = path,
        file_name = basename(path),
        folder_name = basename(dirname(path)),
        file_stem = tools::file_path_sans_ext(basename(path)),
        Source_Series = tools::file_path_sans_ext(basename(path)),
        grouping_rule = "single_file_single_series"
      )

      return(tmp_index)
    }

    if (length(path) == 1 && dir.exists(path)) {
      tmp_files <- list.files(
        path = path,
        pattern = "\\.(csv|xlsx|xls)$",
        full.names = TRUE,
        recursive = recursive
      )

      if (length(tmp_files) == 0) {
        stop("No supported HOBO files were found in the supplied directory.")
      }

      tmp_index <- tibble::tibble(
        file_path = tmp_files,
        file_name = basename(tmp_files),
        rel_dir = dirname(fs::path_rel(tmp_files, start = path)),
        folder_name = basename(dirname(tmp_files)),
        file_stem = tools::file_path_sans_ext(file_name)
      ) %>%
        dplyr::mutate(
          Source_Series = dplyr::case_when(
            rel_dir == "." ~ file_stem,
            TRUE ~ folder_name
          ),
          grouping_rule = dplyr::case_when(
            rel_dir == "." ~ "top_level_files_each_file_is_series",
            TRUE ~ "subfolder_files_combined_within_folder"
          )
        )

      return(tmp_index)
    }

    stop("`path` must be a data frame, an existing file path, a vector of existing file paths, or a directory.")
  }

  report_parsed_groups <- function(file_index) {
    tmp_report <- file_index %>%
      dplyr::count(Source_Series, name = "n_files") %>%
      dplyr::arrange(Source_Series)

    message("Parsed ", nrow(file_index), " HOBO file(s) into ", nrow(tmp_report), " inferred series.")

    purrr::pwalk(
      list(tmp_report$Source_Series, tmp_report$n_files),
      function(Source_Series, n_files) {
        if (n_files == 1) {
          message("  - `", Source_Series, "`: 1 file treated as a standalone series.")
        } else {
          message("  - `", Source_Series, "`: ", n_files, " files combined as one continuous series.")
        }
      }
    )
  }

  determine_grouping_column <- function(df, groupingVariable = NA) {
    if (!is.na(groupingVariable)) {
      if (!groupingVariable %in% names(df)) {
        stop("`groupingVariable = '", groupingVariable, "'` was supplied, but that column was not found in the data.")
      }
      return(groupingVariable)
    }

    if ("Source_Series" %in% names(df)) {
      return("Source_Series")
    }

    return(NA_character_)
  }

  prepare_grouped_data <- function(df, groupingVariable = NA) {
    tmp_group_var <- determine_grouping_column(df, groupingVariable = groupingVariable)

    if (!is.na(tmp_group_var)) {
      df <- df %>%
        dplyr::mutate(
          .group_internal = as.character(.data[[tmp_group_var]]),
          Plot_Group = as.character(.data[[tmp_group_var]])
        )

      message("Using `", tmp_group_var, "` as the grouping variable for summaries, DHWs, and plotting.")
    } else {
      df <- df %>%
        dplyr::mutate(
          .group_internal = "All_Data",
          Plot_Group = "All_Data"
        )

      message("No grouping variable detected. Treating all data as one continuous series.")
    }

    df
  }

  calc_dhw_by_date <- function(dates, dhDay) {
    purrr::map_dbl(seq_along(dates), function(i) {
      tmp_window_start <- dates[i] - 83
      tmp_window_idx <- dates >= tmp_window_start & dates <= dates[i]
      tmp_window_values <- dhDay[tmp_window_idx]

      if (any(is.na(tmp_window_values))) {
        return(NA_real_)
      }

      sum(tmp_window_values) / 7
    })
  }

  calculate_daily_metrics <- function(df, MMM, anomaly) {
    tmp_daily_raw <- df %>%
      dplyr::arrange(.group_internal, DateTime) %>%
      dplyr::group_by(.group_internal, Plot_Group) %>%
      dplyr::mutate(
        Date = as.Date(DateTime),
        next_time = dplyr::lead(DateTime),
        interval_minutes = as.numeric(difftime(next_time, DateTime, units = "mins"))
      ) %>%
      dplyr::group_by(.group_internal, Plot_Group, Date) %>%
      dplyr::mutate(
        interval_minutes = dplyr::if_else(
          is.na(interval_minutes) | interval_minutes <= 0,
          stats::median(interval_minutes[interval_minutes > 0], na.rm = TRUE),
          interval_minutes
        ),
        interval_minutes = dplyr::if_else(
          is.na(interval_minutes) | !is.finite(interval_minutes),
          NA_real_,
          interval_minutes
        ),
        hotspot = dplyr::if_else(
          `Temperature °C` >= (MMM + anomaly),
          `Temperature °C` - MMM,
          0
        ),
        hotspot_weighted = hotspot * (interval_minutes / (24 * 60))
      ) %>%
      dplyr::ungroup()

    tmp_daily_summary <- tmp_daily_raw %>%
      dplyr::group_by(.group_internal, Plot_Group, Date) %>%
      dplyr::summarise(
        Logger_SN = safe_first_non_na(Logger_SN),
        Plot_Title = safe_first_non_na(Plot_Title),
        Source_File = safe_first_non_na(Source_File),
        Source_Series = safe_first_non_na(Source_Series),
        Temperature_Average = mean(`Temperature °C`, na.rm = TRUE),
        Temperature_StDev = stats::sd(`Temperature °C`, na.rm = TRUE),
        Temperature_Min = min(`Temperature °C`, na.rm = TRUE),
        Temperature_Max = max(`Temperature °C`, na.rm = TRUE),
        N_Records = dplyr::n(),
        dhDay = if (all(is.na(hotspot_weighted))) NA_real_ else sum(hotspot_weighted, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(.group_internal, Plot_Group) %>%
      tidyr::complete(Date = seq(min(Date), max(Date), by = "day")) %>%
      dplyr::arrange(Date, .by_group = TRUE) %>%
      tidyr::fill(
        Logger_SN,
        Plot_Title,
        Source_File,
        Source_Series,
        .direction = "downup"
      ) %>%
      dplyr::mutate(
        DHWs = calc_dhw_by_date(Date, dhDay)
      ) %>%
      dplyr::ungroup()

    tmp_complete <- tmp_daily_raw %>%
      dplyr::select(-next_time) %>%
      dplyr::mutate(Date = as.Date(DateTime)) %>%
      dplyr::left_join(
        tmp_daily_summary %>%
          dplyr::select(.group_internal, Plot_Group, Date, dhDay, DHWs),
        by = c(".group_internal", "Plot_Group", "Date")
      ) %>%
      dplyr::select(-Date, -interval_minutes, -hotspot, -hotspot_weighted)

    list(
      full = tmp_complete,
      summary = tmp_daily_summary
    )
  }

  summarize_daily_only <- function(df) {
    df %>%
      dplyr::mutate(Date = as.Date(DateTime)) %>%
      dplyr::group_by(.group_internal, Plot_Group, Date) %>%
      dplyr::summarise(
        Logger_SN = safe_first_non_na(Logger_SN),
        Plot_Title = safe_first_non_na(Plot_Title),
        Source_File = safe_first_non_na(Source_File),
        Source_Series = safe_first_non_na(Source_Series),
        Temperature_Average = mean(`Temperature °C`, na.rm = TRUE),
        Temperature_StDev = stats::sd(`Temperature °C`, na.rm = TRUE),
        Temperature_Min = min(`Temperature °C`, na.rm = TRUE),
        Temperature_Max = max(`Temperature °C`, na.rm = TRUE),
        N_Records = dplyr::n(),
        .groups = "drop"
      )
  }

  # ===== Import and Standardize Data =====
  if (is.data.frame(path)) {
    message("Importing HOBO data from an in-memory data frame.")

    hoboFile <- standardize_hobo_data(
      hoboFile = path,
      logger_sn = if ("Logger_SN" %in% names(path)) dplyr::first(path$Logger_SN) else NA_character_,
      plot_title = if ("Plot_Title" %in% names(path)) dplyr::first(path$Plot_Title) else NA_character_
    )

    for (tmp_col in c("Source_File", "Source_Series")) {
      if (tmp_col %in% names(path) && !tmp_col %in% names(hoboFile)) {
        hoboFile[[tmp_col]] <- path[[tmp_col]]
      }
    }

    tmp_extra_cols <- setdiff(names(path), names(hoboFile))

    if (length(tmp_extra_cols) > 0) {
      hoboFile <- dplyr::bind_cols(
        hoboFile,
        path %>% dplyr::select(dplyr::all_of(tmp_extra_cols))
      )
    }

  } else {
    tmp_index <- infer_file_groups(path = path, recursive = recursive)
    report_parsed_groups(tmp_index)

    tmp_imported <- purrr::pmap(
      list(tmp_index$file_path, tmp_index$Source_Series),
      function(file_path, Source_Series) {
        message("Reading: ", basename(file_path))

        tmp_raw <- read_hobo_file(file_path)

        tmp_std <- standardize_hobo_data(
          hoboFile = tmp_raw$data,
          logger_sn = tmp_raw$logger_sn,
          plot_title = tmp_raw$plot_title
        )

        tmp_std %>%
          dplyr::mutate(
            Source_File = basename(file_path),
            Source_Series = Source_Series
          )
      }
    )

    hoboFile <- dplyr::bind_rows(tmp_imported) %>%
      dplyr::arrange(Source_Series, DateTime) %>%
      dplyr::distinct()
  }

  if (nrow(hoboFile) == 0) {
    stop("No valid HOBO observations remained after import and cleaning.")
  }

  # ===== Prepare Grouping Structure =====
  tmp_group_var <- determine_grouping_column(
    df = hoboFile,
    groupingVariable = groupingVariable
  )

  hoboFile <- prepare_grouped_data(
    df = hoboFile,
    groupingVariable = groupingVariable
  )

  plot_group_label <- dplyr::case_when(
    is.na(tmp_group_var) ~ "Series",
    tmp_group_var == "Source_Series" ~ "Series",
    TRUE ~ tmp_group_var
  )

  # ===== Calculate Daily Metrics =====
  if (calc) {
    tmp_metrics <- calculate_daily_metrics(
      df = hoboFile,
      MMM = MMM,
      anomaly = anomaly
    )

    hoboComplete <- tmp_metrics$full
    hoboSummary <- tmp_metrics$summary
  } else {
    warning("`calc = FALSE`, so DHW metrics were not calculated.")
    hoboComplete <- hoboFile
    hoboSummary <- summarize_daily_only(hoboFile)
  }

  # ===== Plot Outputs =====
  if (plot) {
    if ("Plot_Group" %in% names(hoboComplete) &&
        dplyr::n_distinct(hoboComplete$Plot_Group, na.rm = TRUE) > 1) {

      p_temp <- ggplot2::ggplot(
        hoboComplete,
        ggplot2::aes(
          x = DateTime,
          y = `Temperature °C`,
          color = Plot_Group,
          fill = Plot_Group,
          group = Plot_Group
        )
      ) +
        ggplot2::geom_line(linewidth = 0.6, alpha = 0.8) +
        ggplot2::geom_point(size = 1.5, alpha = 0.8) +
        ggplot2::geom_hline(
          yintercept = MMM + anomaly,
          linetype = "dashed",
          linewidth = 0.6
        ) +
        ggplot2::labs(
          title = "Temperature time series",
          x = NULL,
          y = "Temperature (°C)",
          color = plot_group_label,
          fill = plot_group_label
        ) +
        ggthemes::theme_few(base_size = 13)

    } else {
      p_temp <- ggplot2::ggplot(
        hoboComplete,
        ggplot2::aes(
          x = DateTime,
          y = `Temperature °C`
        )
      ) +
        ggplot2::geom_line(linewidth = 0.6, alpha = 0.8) +
        ggplot2::geom_point(size = 1.5, alpha = 0.8) +
        ggplot2::geom_hline(
          yintercept = MMM + anomaly,
          linetype = "dashed",
          linewidth = 0.6
        ) +
        ggplot2::labs(
          title = "Temperature time series",
          x = NULL,
          y = "Temperature (°C)"
        ) +
        ggthemes::theme_few(base_size = 13)
    }

    if (calc) {
      if ("Plot_Group" %in% names(hoboSummary) &&
          dplyr::n_distinct(hoboSummary$Plot_Group, na.rm = TRUE) > 1) {

        p_dhw <- ggplot2::ggplot(
          hoboSummary,
          ggplot2::aes(
            x = Date,
            y = DHWs,
            color = Plot_Group,
            fill = Plot_Group,
            group = Plot_Group
          )
        ) +
          ggplot2::geom_line(linewidth = 0.7, alpha = 0.8, na.rm = TRUE) +
          ggplot2::geom_point(size = 1.6, alpha = 0.8, na.rm = TRUE) +
          ggplot2::labs(
            title = "Degree Heating Weeks (DHWs)",
            x = "Date",
            y = "DHWs",
            color = plot_group_label,
            fill = plot_group_label
          ) +
          ggthemes::theme_few(base_size = 13)

      } else {
        p_dhw <- ggplot2::ggplot(
          hoboSummary,
          ggplot2::aes(
            x = Date,
            y = DHWs
          )
        ) +
          ggplot2::geom_line(linewidth = 0.7, alpha = 0.8, na.rm = TRUE) +
          ggplot2::geom_point(size = 1.6, alpha = 0.8, na.rm = TRUE) +
          ggplot2::labs(
            title = "Degree Heating Weeks (DHWs)",
            x = "Date",
            y = "DHWs"
          ) +
          ggthemes::theme_few(base_size = 13)
      }

      full_plot <- p_temp + p_dhw + patchwork::plot_layout(ncol = 2)
    } else {
      full_plot <- p_temp
    }

    if (!is.null(plotFile)) {
      ggplot2::ggsave(
        filename = plotFile,
        plot = full_plot,
        width = ifelse(calc, 16, 10),
        height = 7,
        dpi = 300
      )

      message("Plot saved to: ", normalizePath(plotFile))
    } else {
      print(full_plot)
    }
  }

  # ===== Clean Output =====
  hoboComplete <- hoboComplete %>%
    dplyr::select(-.group_internal)

  hoboSummary <- hoboSummary %>%
    dplyr::select(-.group_internal)

  # ===== Return Requested Output =====
  if (summary) {
    return(hoboSummary)
  } else {
    return(hoboComplete)
  }
}
