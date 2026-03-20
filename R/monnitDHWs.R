#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
monnitDHWs <- function(path, MMM, anomaly = 1, calc = TRUE, groupingVariable = NA, summary = FALSE) {

  readIn <- function(f){
    file <- f %>%
      dplyr::rename(`Temperature °C` = `Raw Data`,
                    `Temperature °F` = Value) %>%
      dplyr::mutate(DateTime = lubridate::mdy_hm(Date)) %>%
      dplyr::select(SensorID,`Sensor Name`,DateTime,`Temperature °C`,`Temperature °F`)
  }

  reReadIn <- function(f){
    file <- f %>%
      dplyr::select(SensorID,`Sensor Name`,Group,DateTime,`Temperature °C`,`Temperature °F`)
  }

  if (is.data.frame(path) == TRUE){
    monnitFile <- reReadIn(path)
  }else if (dir.exists(path) == TRUE) {
    list <- list.files(path, full.names = TRUE) %>%
      lapply(rio::import) %>%
      lapply(readIn)
    names(list) <- basename(list.files(path,full.names = TRUE))
    monnitFile <- list %>%
      purrr:::map_dfr(plyr::rbind.fill) %>%
      unique() %>%
      dplyr::mutate(Group = `Sensor Name`)
  } else {print("Import failed.")}

  ####Calculation Math####
  if (calc == TRUE) {
    if (is.na(groupingVariable) == FALSE) {
      calcTemperatures <- monnitFile %>%
        dplyr::group_by(Group) %>%
        dplyr::arrange(DateTime) %>%
        dplyr::mutate(lagTime = dplyr::lag(DateTime,1),
                      diffTime = as.integer(difftime(DateTime,lagTime, units = "mins")),
                      Date = lubridate::date(DateTime)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(Group, Date) %>%
        dplyr::mutate(sumMinutes = sum(diffTime,na.rm=TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(percentMinutes = diffTime/sumMinutes,
                      diffTemp = ifelse(`Temperature °C` >= (MMM+anomaly), `Temperature °C` - MMM, 0),
                      percentTemp = diffTemp * percentMinutes) %>%
        dplyr::group_by(Group, Date) %>%
        dplyr::arrange(Date) %>%
        dplyr::mutate(dhDay = sum(percentTemp, na.rm=TRUE)) %>%
        dplyr::ungroup()

      monnitSummary <- monnitFile %>%
        dplyr::mutate(Date = lubridate::date(DateTime)) %>%
        dplyr::group_by(Group, Date) %>%
        dplyr::summarise(Temperature_Average = ave(`Temperature °C`,na.rm=TRUE),
                       Temperature_StDev = sd(`Temperature °C`,na.rm=TRUE)) %>%
        unique() %>%
        dplyr::ungroup() %>%
        dplyr::left_join(unique(calcTemperatures[c("Date", "dhDay", "Group")]), by = c("Date", "Group")) %>%
        dplyr::mutate(dhWeek=dhDay*(1/7)) %>%
        dplyr::group_by(Group) %>%
        dplyr::arrange(Date) %>%
        dplyr::mutate(DHWs = zoo::rollapplyr(dhWeek, (7*12), sum, partial = TRUE))

      monnitComplete <- monnitFile %>%
        dplyr::mutate(Date = lubridate::date(DateTime)) %>%
        dplyr::left_join(monnitSummary[c("Date","dhDay","DHWs", "Group")], by = c("Date", "Group")) %>%
        dplyr::select(-c(Date))
    } else {
      calcTemperatures <- monnitFile %>%
        dplyr::arrange(DateTime) %>%
        dplyr::mutate(lagTime = dplyr::lag(DateTime,1),
                      diffTime = as.integer(difftime(DateTime,lagTime, units = "mins")),
                      Date = lubridate::date(DateTime)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(Date) %>%
        dplyr::mutate(sumMinutes = sum(diffTime,na.rm=TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(percentMinutes = diffTime/sumMinutes,
                      diffTemp = ifelse(`Temperature °C` >= (MMM+anomaly), `Temperature °C` - MMM, 0),
                      percentTemp = diffTemp * percentMinutes) %>%
        dplyr::group_by(Date) %>%
        dplyr::arrange(Date) %>%
        dplyr::mutate(dhDay = sum(percentTemp, na.rm=TRUE)) %>%
        dplyr::ungroup()

      monnitSummary <- monnitFile %>%
        dplyr::mutate(Date = lubridate::date(DateTime)) %>%
        dplyr::group_by(Date) %>%
        dplyr::summarise(Temperature_Average = ave(`Temperature °C`,na.rm=TRUE),
                       Temperature_StDev = sd(`Temperature °C`,na.rm=TRUE)) %>%
        unique() %>%
        dplyr::ungroup() %>%
        dplyr::left_join(unique(calcTemperatures[c("Date", "dhDay")]), by = c("Date")) %>%
        dplyr::mutate(dhWeek=dhDay*(1/7)) %>%
        dplyr::arrange(Date) %>%
        dplyr::mutate(DHWs = zoo::rollapplyr(dhWeek, (7*12), sum, partial = TRUE))

      monnitComplete <- monnitFile %>%
        dplyr::mutate(Date = lubridate::date(DateTime)) %>%
        dplyr::left_join(monnitSummary[c("Date","dhDay","DHWs")], by = "Date") %>%
        dplyr::select(-c(Date))
    }
  } else if (calc == FALSE){
    if (is.na(groupingVariable) == FALSE) {
      monnitSummary <- monnitFile %>%
        dplyr::mutate(Date = lubridate::date(DateTime)) %>%
        dplyr::group_by(Group, Date) %>%
        dplyr::summarise(Temperature_Average = ave(`Temperature °C`,na.rm=TRUE),
                       Temperature_StDev = sd(`Temperature °C`,na.rm=TRUE)) %>%
        unique()

      monnitComplete <- monnitFile
    } else {
      monnitSummary <- monnitFile %>%
        dplyr::mutate(Date = lubridate::date(DateTime)) %>%
        dplyr::group_by(Date) %>%
        dplyr::summarise(Temperature_Average = ave(`Temperature °C`,na.rm=TRUE),
                       Temperature_StDev = sd(`Temperature °C`,na.rm=TRUE)) %>%
        unique()

      monnitComplete <- monnitFile
    }
  }

  if (summary == TRUE){
    monnitSummary %>%
      return(.)
  }else if (summary == FALSE){
    monnitComplete %>%
      return(.)
  }
}
