#' Calculate symbiont community proportions and symbiont:host signal from QuantStudio output
#'
#' This function processes qPCR summary output exported from QuantStudio, applies
#' replicate-based quality control, calculates relative symbiont proportions, and,
#' when host assays are present, estimates total symbiont:host signal.
#'
#' @param df A data frame containing QuantStudio summary output. Expected columns
#'   include `Sample.Name`, `File.Name`, and target-specific summary columns such as
#'   `.reps`, `.CT.mean`, and `.CT.sd`.
#' @param sdThreshold Numeric threshold used to determine whether duplicate Ct values
#'   pass replicate consistency QC.
#' @param QC Logical; if `TRUE`, return only the QC/rerun table.
#' @param ctCutoff Numeric Ct threshold for pass/fail filtering. Use `NA` to ignore
#'   Ct cutoff and use only replicate count and SD criteria.
#' @param singles Logical; if `TRUE`, allow single replicates to pass under the
#'   specified Ct rules.
#' @param returnAll Logical; if `TRUE`, return the processed output joined back to
#'   the original input data.
#' @param rerunDetail Logical; if `TRUE` and `QC = TRUE`, return the detailed
#'   target-level rerun table rather than the collapsed sample-level rerun table.
#' @param ctWindow Numeric Ct window above the minimum symbiont Ct within a sample
#'   used to classify failed symbiont amplifications as relatively important for
#'   rerun review.
#'
#' @return By default, a processed data frame containing symbiont proportions,
#'   dominant/background symbiont identity, summed symbiont signal, and host-normalized
#'   signal when host assays are present. If `QC = TRUE`, returns the rerun/QC table.
#'   If `QC = TRUE` and `rerunDetail = TRUE`, returns the detailed target-level rerun
#'   table. If `returnAll = TRUE`, returns the processed output joined back to the
#'   original input data.
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @export
qPCRcalc <- function(
    df,
    sdThreshold,
    QC = FALSE,
    ctCutoff = 40,
    singles = FALSE,
    returnAll = FALSE,
    rerunDetail = FALSE,
    ctWindow = 1
) {

  # ===== Identify assay targets =====
  targets <- data.frame(
    target = strex::str_before_last_dot(
      colnames(
        df %>%
          dplyr::select(dplyr::contains(".reps"))
      )
    )
  ) %>%
    dplyr::mutate(
      type = dplyr::case_when(
        .data$target %in% c(
          "A", "B", "C", "D", "E", "F", "G", "H", "M",
          "a", "b", "c", "d", "e", "f", "g", "h", "m",
          "Symbiodinium", "Breviolum", "Cladocopium", "Durusdinium", "Effrenium", "Fugacium",
          "Sym", "sym", "symbiont"
        ) ~ "Symbiont",
        TRUE ~ "Coral"
      )
    )

  if (length(unique(targets$type)) == 1) {
    plateType <- "symComm"
  } else {
    plateType <- "symHost"
  }

  # ===== Helper to parse QuantStudio summary columns =====
  prep_qpcr_table <- function(df, targets, keep_zero = FALSE, summarize_values = TRUE) {

    out <- df %>%
      tidyr::pivot_longer(
        cols = -dplyr::all_of(c("Sample.Name", "File.Name"))
      ) %>%
      tidyr::separate(.data$name, into = c("target", "calc"), extra = "merge") %>%
      dplyr::full_join(targets, by = "target")

    if (summarize_values == TRUE) {
      out <- out %>%
        dplyr::mutate(
          value = dplyr::case_when(
            is.nan(.data$value) ~ NA_real_,
            TRUE ~ .data$value
          )
        ) %>%
        dplyr::filter(!is.na(.data$target) & .data$target != "NA") %>%
        dplyr::group_by(.data$Sample.Name, .data$File.Name, .data$target, .data$calc) %>%
        dplyr::mutate(value = mean(.data$value, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        unique()
    }

    if (keep_zero == TRUE) {
      out <- out %>%
        dplyr::filter(.data$value >= 0)
    } else {
      out <- out %>%
        dplyr::filter(.data$value > 0)
    }

    out %>%
      tidyr::pivot_wider(
        names_from = "calc",
        values_from = "value"
      )
  }

  # ===== Apply standard QC rules =====
  if (is.na(ctCutoff)) {

    dff2 <- prep_qpcr_table(
      df = df,
      targets = targets,
      keep_zero = TRUE,
      summarize_values = TRUE
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        qc = dplyr::case_when(
          .data$reps == 2 & .data$CT.sd < sdThreshold ~ "Pass",
          .data$reps == 0 ~ "NotPresent",
          TRUE ~ "Fail"
        )
      ) %>%
      dplyr::filter(.data$type == "Symbiont" | (.data$type == "Coral" & .data$qc != "Fail")) %>%
      dplyr::relocate(.data$qc)

    dff <- dff2 %>%
      dplyr::filter(.data$reps > 0)

  } else {

    dff2 <- prep_qpcr_table(
      df = df,
      targets = targets,
      keep_zero = FALSE,
      summarize_values = TRUE
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        qc = dplyr::case_when(
          .data$reps == 2 & .data$CT.sd < sdThreshold & .data$CT.mean <= ctCutoff ~ "Pass",
          .data$reps == 0 ~ "NotPresent",
          TRUE ~ "Fail"
        )
      ) %>%
      dplyr::filter(.data$type == "Symbiont" | (.data$type == "Coral" & .data$qc != "Fail")) %>%
      dplyr::relocate(.data$qc)

    dff <- dff2 %>%
      dplyr::filter(.data$reps > 0)
  }

  # ===== Optionally allow single-replicate targets =====
  if (singles == TRUE) {

    dff <- prep_qpcr_table(
      df = df,
      targets = targets,
      keep_zero = FALSE,
      summarize_values = FALSE
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        qc = dplyr::case_when(
          .data$reps > 0 ~ "Pass",
          .data$reps == 0 ~ "NotPresent",
          TRUE ~ "Fail"
        )
      ) %>%
      dplyr::filter(.data$type == "Symbiont" | (.data$type == "Coral" & .data$qc == "Pass")) %>%
      dplyr::relocate(.data$qc)

    dff2 <- prep_qpcr_table(
      df = df,
      targets = targets,
      keep_zero = FALSE,
      summarize_values = TRUE
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        qc = dplyr::case_when(
          .data$reps == 1 & .data$CT.mean <= ctCutoff ~ "Pass",
          .data$reps == 2 & .data$CT.sd < sdThreshold & .data$CT.mean <= ctCutoff ~ "Pass",
          .data$reps == 0 ~ "NotPresent",
          TRUE ~ "Fail"
        )
      ) %>%
      dplyr::filter(.data$type == "Symbiont" | (.data$type == "Coral" & .data$qc != "Fail")) %>%
      dplyr::relocate(.data$qc)
  }

  # ===== Build denominator table for symbiont composition =====
  denoms <- dff %>%
    dplyr::filter(.data$qc != "Fail") %>%
    dplyr::filter(.data$type == "Symbiont") %>%
    dplyr::select(-c(.data$qc, .data$CT.mean, .data$CT.sd, .data$reps, .data$type)) %>%
    dplyr::rename(numerator = .data$target)

  reps <- dff2 %>%
    dplyr::select(.data$Sample.Name, .data$File.Name, .data$target, .data$reps, .data$qc)

  # ===== Sum relative values across symbiont targets =====
  denoms2 <- denoms %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(c("Sample.Name", "File.Name", "numerator")),
      values_to = "value",
      names_to = "target"
    ) %>%
    dplyr::left_join(targets, by = "target") %>%
    dplyr::filter(.data$type == "Symbiont") %>%
    dplyr::group_by(.data$Sample.Name, .data$File.Name, .data$target) %>%
    dplyr::mutate(
      N = dplyr::n(),
      valSum = sum(.data$value, na.rm = TRUE)
    ) %>%
    dplyr::left_join(reps, by = c("Sample.Name", "File.Name", "target")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      value2 = dplyr::case_when(
        .data$qc == "Pass" ~ .data$valSum
      )
    ) %>%
    dplyr::select(-c(.data$reps, .data$qc, .data$value, .data$valSum, .data$numerator, .data$N)) %>%
    unique()

  # ===== Calculate symbiont proportions =====
  symPropsLong <- denoms2 %>%
    dplyr::arrange(.data$target) %>%
    dplyr::filter(.data$type == "Symbiont") %>%
    dplyr::mutate(
      prop = 1 / (.data$value2 + 1)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$Sample.Name, .data$File.Name, .data$target, .data$prop)

  maxProps <- symPropsLong %>%
    dplyr::group_by(.data$Sample.Name, .data$File.Name) %>%
    dplyr::arrange(dplyr::desc(.data$prop), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(prop = round(.data$prop, 3)) %>%
    dplyr::rename(
      maxProp = .data$prop,
      maxSym = .data$target
    ) %>%
    dplyr::filter(!is.na(.data$maxProp))

  backgroundProps <- symPropsLong %>%
    dplyr::group_by(.data$Sample.Name, .data$File.Name) %>%
    dplyr::arrange(dplyr::desc(.data$prop), .by_group = TRUE) %>%
    dplyr::slice(2:2) %>%
    dplyr::mutate(prop = round(.data$prop, 3)) %>%
    dplyr::rename(
      backProp = .data$prop,
      backSym = .data$target
    ) %>%
    dplyr::filter(!is.na(.data$backProp))

  symProps <- symPropsLong %>%
    tidyr::pivot_wider(
      names_from = "target",
      values_from = "prop",
      names_prefix = "prop",
      names_sep = ""
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      sumSyms = round(sum(dplyr::c_across(dplyr::starts_with("prop")), na.rm = TRUE), 3)
    ) %>%
    dplyr::ungroup()

  # ===== Calculate symbiont:host signal when host assays are present =====
  if (plateType == "symHost") {

    denoms <- dff %>%
      dplyr::filter(.data$qc != "Fail") %>%
      dplyr::select(-c(.data$qc, .data$CT.mean, .data$CT.sd, .data$reps, .data$type)) %>%
      dplyr::rename(numerator = .data$target)

    SHtotals <- denoms %>%
      tidyr::pivot_longer(
        cols = -dplyr::all_of(c("Sample.Name", "File.Name", "numerator")),
        values_to = "value",
        names_to = "target"
      ) %>%
      dplyr::left_join(targets, by = "target") %>%
      dplyr::filter(.data$type == "Coral", !is.na(.data$value)) %>%
      tidyr::pivot_wider(
        names_from = "numerator",
        names_glue = "{numerator}.Host",
        values_from = "value"
      ) %>%
      dplyr::select(-c(.data$type, .data$target)) %>%
      dplyr::group_by(.data$Sample.Name, .data$File.Name) %>%
      dplyr::mutate(
        symHost = sum(dplyr::c_across(dplyr::ends_with(".Host")), na.rm = TRUE)
      ) %>%
      dplyr::relocate(.data$Sample.Name, .data$File.Name, .data$symHost) %>%
      dplyr::ungroup()

    df2 <- symProps %>%
      dplyr::full_join(SHtotals, by = c("Sample.Name", "File.Name")) %>%
      dplyr::full_join(maxProps, by = c("Sample.Name", "File.Name")) %>%
      dplyr::full_join(backgroundProps, by = c("Sample.Name", "File.Name"))

  } else {

    df2 <- dplyr::full_join(symProps, maxProps, by = c("Sample.Name", "File.Name")) %>%
      dplyr::full_join(backgroundProps, by = c("Sample.Name", "File.Name"))
  }

  # ===== Build rerun recommendation table =====
  tmp_ct_context <- dff2 %>%
    dplyr::filter(.data$reps > 0, !is.na(.data$CT.mean)) %>%
    dplyr::group_by(.data$File.Name, .data$Sample.Name) %>%
    dplyr::mutate(
      sampleMinCtAll = min(.data$CT.mean, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$File.Name, .data$Sample.Name, .data$type) %>%
    dplyr::mutate(
      sampleMinCtType = min(.data$CT.mean, na.rm = TRUE),
      sampleMeanCtType = mean(.data$CT.mean, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      File.Name = .data$File.Name,
      Sample.Name = .data$Sample.Name,
      target = .data$target,
      type = .data$type,
      sampleMinCtAll = .data$sampleMinCtAll,
      sampleMinCtType = .data$sampleMinCtType,
      sampleMeanCtType = .data$sampleMeanCtType,
      isLowestCtAll = .data$CT.mean == .data$sampleMinCtAll,
      isLowestCtType = .data$CT.mean == .data$sampleMinCtType
    )

  tmp_failed_detail <- dff2 %>%
    dplyr::filter(.data$qc == "Fail") %>%
    dplyr::select(
      .data$Sample.Name,
      .data$File.Name,
      .data$target,
      .data$type,
      .data$reps,
      .data$CT.mean,
      .data$CT.sd
    ) %>%
    dplyr::left_join(
      tmp_ct_context,
      by = c("File.Name", "Sample.Name", "target", "type")
    ) %>%
    dplyr::mutate(
      isLowestCtAll = dplyr::coalesce(.data$isLowestCtAll, FALSE),
      isLowestCtType = dplyr::coalesce(.data$isLowestCtType, FALSE),

      isImportantSymFail = dplyr::case_when(
        .data$type != "Symbiont" ~ FALSE,
        is.na(.data$CT.mean) ~ FALSE,
        .data$reps == 0 ~ FALSE,
        .data$isLowestCtAll ~ TRUE,
        .data$isLowestCtType ~ TRUE,
        !is.na(.data$sampleMinCtType) & .data$CT.mean <= .data$sampleMinCtType + ctWindow ~ TRUE,
        !is.na(.data$sampleMeanCtType) & .data$CT.mean <= .data$sampleMeanCtType ~ TRUE,
        TRUE ~ FALSE
      ),

      failedReason = dplyr::case_when(
        .data$type == "Coral" ~ "failed host amplification",
        .data$isImportantSymFail & .data$isLowestCtAll ~ "lowest Ct amplification failed QC",
        .data$isImportantSymFail & .data$isLowestCtType ~ "lowest symbiont Ct failed QC",
        .data$isImportantSymFail ~ "important symbiont amplification failed QC",
        .data$reps == 2 & !is.na(.data$CT.sd) & .data$CT.sd >= sdThreshold ~ "duplicate amplification failed SD",
        .data$reps == 1 & singles == FALSE ~ "single amplification not counted",
        TRUE ~ "failed amplification"
      ),

      rerunPriorityNum = dplyr::case_when(
        .data$type == "Coral" ~ 1,
        .data$isImportantSymFail & .data$isLowestCtAll ~ 2,
        .data$isImportantSymFail & .data$isLowestCtType ~ 3,
        .data$isImportantSymFail ~ 4,
        .data$reps == 2 & !is.na(.data$CT.sd) & .data$CT.sd >= sdThreshold ~ 5,
        .data$reps == 1 & singles == FALSE ~ 6,
        TRUE ~ 7
      ),

      rerunPriority = dplyr::case_when(
        .data$rerunPriorityNum %in% c(1, 2, 3, 4) ~ "1_High",
        .data$rerunPriorityNum %in% c(5, 6) ~ "2_Moderate",
        TRUE ~ "3_Low"
      ),

      rerunTargets = .data$target,
      rerunSummary = paste0(.data$failedReason, ": ", .data$target)
    ) %>%
    dplyr::arrange(.data$File.Name, .data$Sample.Name, .data$rerunPriorityNum, .data$CT.mean, .data$target)

  tmp_failed_summary <- tmp_failed_detail %>%
    dplyr::group_by(
      .data$File.Name,
      .data$Sample.Name,
      .data$rerunPriorityNum,
      .data$rerunPriority,
      .data$failedReason
    ) %>%
    dplyr::summarise(
      rerunTargets = paste(.data$target, collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$File.Name, .data$Sample.Name, .data$rerunPriorityNum, .data$failedReason) %>%
    dplyr::group_by(.data$File.Name, .data$Sample.Name) %>%
    dplyr::summarise(
      rerunPriorityNum = min(.data$rerunPriorityNum),
      rerunPriority = paste(.data$rerunPriority, collapse = " | "),
      rerunReason = paste(.data$failedReason, collapse = " | "),
      rerunTargets = paste(.data$rerunTargets, collapse = " | "),
      rerunSummary = paste0(.data$failedReason, ": ", .data$rerunTargets, collapse = " || "),
      .groups = "drop"
    )

  noReps <- reps %>%
    dplyr::left_join(targets, by = "target") %>%
    dplyr::group_by(.data$Sample.Name, .data$File.Name, .data$type) %>%
    dplyr::summarise(sum = sum(.data$reps, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      rerunPriorityNum = dplyr::case_when(
        .data$type == "Coral" & .data$sum != 2 ~ 1,
        .data$type == "Symbiont" & .data$sum == 0 ~ 2,
        TRUE ~ NA_real_
      ),
      rerunPriority = dplyr::case_when(
        .data$type == "Coral" & .data$sum != 2 ~ "1_High",
        .data$type == "Symbiont" & .data$sum == 0 ~ "1_High",
        TRUE ~ NA_character_
      ),
      rerunReason = dplyr::case_when(
        .data$type == "Symbiont" & .data$sum == 0 ~ "no symbiont DNA",
        .data$type == "Coral" & .data$sum != 2 ~ "failed host amplification",
        TRUE ~ NA_character_
      ),
      rerunTargets = dplyr::case_when(
        .data$type == "Symbiont" & .data$sum == 0 ~ "sumSyms",
        .data$type == "Coral" & .data$sum != 2 ~ "symHost",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(.data$rerunReason)) %>%
    dplyr::group_by(.data$Sample.Name, .data$File.Name) %>%
    dplyr::mutate(
      nProblems = dplyr::n(),
      omit = dplyr::case_when(
        .data$nProblems == 1 & .data$rerunTargets == "symHost" ~ "Y"
      )
    ) %>%
    dplyr::filter(is.na(.data$omit)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      rerunSummary = paste0(.data$rerunReason, ": ", .data$rerunTargets)
    ) %>%
    dplyr::select(
      .data$File.Name,
      .data$Sample.Name,
      .data$rerunPriorityNum,
      .data$rerunPriority,
      .data$rerunReason,
      .data$rerunTargets,
      .data$rerunSummary
    )

  reruns_detail <- dplyr::bind_rows(
    tmp_failed_detail %>%
      dplyr::select(
        .data$File.Name,
        .data$Sample.Name,
        .data$rerunPriorityNum,
        .data$rerunPriority,
        rerunReason = .data$failedReason,
        .data$rerunTargets,
        .data$rerunSummary,
        .data$type,
        .data$reps,
        .data$CT.mean,
        .data$CT.sd,
        .data$sampleMinCtAll,
        .data$sampleMinCtType,
        .data$sampleMeanCtType,
        .data$isLowestCtAll,
        .data$isLowestCtType,
        .data$isImportantSymFail
      ),
    noReps %>%
      dplyr::mutate(
        type = NA_character_,
        reps = NA_real_,
        CT.mean = NA_real_,
        CT.sd = NA_real_,
        sampleMinCtAll = NA_real_,
        sampleMinCtType = NA_real_,
        sampleMeanCtType = NA_real_,
        isLowestCtAll = NA,
        isLowestCtType = NA,
        isImportantSymFail = NA
      ) %>%
      dplyr::select(
        .data$File.Name,
        .data$Sample.Name,
        .data$rerunPriorityNum,
        .data$rerunPriority,
        .data$rerunReason,
        .data$rerunTargets,
        .data$rerunSummary,
        .data$type,
        .data$reps,
        .data$CT.mean,
        .data$CT.sd,
        .data$sampleMinCtAll,
        .data$sampleMinCtType,
        .data$sampleMeanCtType,
        .data$isLowestCtAll,
        .data$isLowestCtType,
        .data$isImportantSymFail
      )
  ) %>%
    dplyr::arrange(.data$File.Name, .data$Sample.Name, .data$rerunPriorityNum, .data$CT.mean, .data$rerunTargets) %>%
    dplyr::left_join(df, by = c("Sample.Name", "File.Name")) %>%
    dplyr::filter(!.data$Sample.Name %in% c("NTC", "-")) %>%
    dplyr::arrange(.data$rerunPriority, .data$File.Name)

  reruns <- dplyr::bind_rows(
    tmp_failed_summary,
    noReps
  ) %>%
    dplyr::arrange(.data$File.Name, .data$Sample.Name, .data$rerunPriorityNum) %>%
    dplyr::group_by(.data$File.Name, .data$Sample.Name) %>%
    dplyr::summarise(
      rerunPriority = paste(.data$rerunPriority, collapse = " | "),
      rerunReason = paste(.data$rerunReason, collapse = " | "),
      rerunTargets = paste(.data$rerunTargets, collapse = " | "),
      rerunSummary = paste(.data$rerunSummary, collapse = " || "),
      .groups = "drop"
    ) %>%
    dplyr::left_join(df, by = c("Sample.Name", "File.Name")) %>%
    dplyr::filter(!.data$Sample.Name %in% c("NTC", "-")) %>%
    dplyr::arrange(.data$rerunPriority, .data$File.Name)

  df3 <- df2 %>%
    dplyr::left_join(
      reruns %>%
        dplyr::mutate(checkQC = "checkQC") %>%
        dplyr::select(.data$File.Name, .data$Sample.Name, .data$checkQC),
      by = c("File.Name", "Sample.Name")
    )

  # ===== Return requested output =====
  if (QC != TRUE & returnAll == FALSE) {
    return(df3)
  } else if (QC != TRUE & returnAll == TRUE) {
    df4 <- df3 %>%
      dplyr::left_join(df, by = c("Sample.Name", "File.Name"))
    return(df4)
  } else if (QC == TRUE & rerunDetail == TRUE) {
    return(reruns_detail)
  } else {
    return(reruns)
  }
}
