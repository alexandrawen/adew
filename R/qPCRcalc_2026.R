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
#' @export
qPCRcalc <- function(df, sdThreshold, QC = FALSE, ctCutoff = 40, singles = FALSE, returnAll = FALSE, rerunDetail = FALSE, ctWindow = 1) {

  # ===== Identify assay targets =====
  # Reconstruct target names from QuantStudio summary columns ending in ".reps"
  # and classify each assay as either symbiont or coral host.
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
        target %in% c(
          "A", "B", "C", "D", "E", "F", "G", "H", "M",
          "a", "b", "c", "d", "e", "f", "g", "h", "m",
          "Symbiodinium", "Breviolum", "Cladocopium", "Durusdinium", "Effrenium", "Fugacium",
          "Sym", "sym", "symbiont"
        ) ~ "Symbiont",
        TRUE ~ "Coral"
      )
    )

  # Determine whether the plate contains only symbiont assays or both symbiont
  # and host assays, since host-containing plates require an additional
  # symbiont:host calculation later.
  if (length(unique(targets$type)) == 1) {
    plateType <- "symComm"
  } else {
    plateType <- "symHost"
  }

  # ===== Helper to parse QuantStudio summary columns =====
  # This helper keeps the repeated reshaping steps in one place:
  # 1) pivot all target summary columns long,
  # 2) split target from summary statistic,
  # 3) join target metadata,
  # 4) optionally average duplicated entries,
  # 5) return to wide format for QC evaluation.
  prep_qpcr_table <- function(df, targets, keep_zero = FALSE, summarize_values = TRUE) {

    out <- df %>%
      tidyr::pivot_longer(cols = -c(Sample.Name, File.Name)) %>%
      tidyr::separate(name, into = c("target", "calc"), extra = "merge") %>%
      dplyr::full_join(targets, by = "target")

    if (summarize_values == TRUE) {
      out <- out %>%
        dplyr::mutate(
          value = dplyr::case_when(
            is.nan(value) ~ NA,
            TRUE ~ value
          )
        ) %>%
        dplyr::filter(!is.na(target) & target != "NA") %>%
        dplyr::group_by(Sample.Name, File.Name, target, calc) %>%
        dplyr::mutate(value = mean(value, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        unique()
    }

    if (keep_zero == TRUE) {
      out <- out %>%
        dplyr::filter(value >= 0)
    } else {
      out <- out %>%
        dplyr::filter(value > 0)
    }

    out %>%
      tidyr::pivot_wider(
        names_from = "calc",
        values_from = "value"
      )
  }

  # ===== Apply standard QC rules =====
  # dff2 stores the full target-level QC table.
  # dff stores only amplified targets for downstream calculations.
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
          reps == 2 & CT.sd < sdThreshold ~ "Pass",
          reps == 0 ~ "NotPresent",
          TRUE ~ "Fail"
        )
      ) %>%
      dplyr::filter(type == "Symbiont" | (type == "Coral" & qc != "Fail")) %>%
      dplyr::relocate(qc)

    dff <- dff2 %>%
      dplyr::filter(reps > 0)

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
          reps == 2 & CT.sd < sdThreshold & CT.mean <= ctCutoff ~ "Pass",
          reps == 0 ~ "NotPresent",
          TRUE ~ "Fail"
        )
      ) %>%
      dplyr::filter(type == "Symbiont" | (type == "Coral" & qc != "Fail")) %>%
      dplyr::relocate(qc)

    dff <- dff2 %>%
      dplyr::filter(reps > 0)
  }

  # ===== Optionally allow single-replicate targets =====
  # When singles = TRUE, overwrite dff and dff2 using more permissive rules
  # that allow one successful replicate to pass under the Ct cutoff.
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
          reps > 0 ~ "Pass",
          reps == 0 ~ "NotPresent",
          TRUE ~ "Fail"
        )
      ) %>%
      dplyr::filter(type == "Symbiont" | (type == "Coral" & qc == "Pass")) %>%
      dplyr::relocate(qc)

    dff2 <- prep_qpcr_table(
      df = df,
      targets = targets,
      keep_zero = FALSE,
      summarize_values = TRUE
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        qc = dplyr::case_when(
          reps == 1 & CT.mean <= ctCutoff ~ "Pass",
          reps == 2 & CT.sd < sdThreshold & CT.mean <= ctCutoff ~ "Pass",
          reps == 0 ~ "NotPresent",
          TRUE ~ "Fail"
        )
      ) %>%
      dplyr::filter(type == "Symbiont" | (type == "Coral" & qc != "Fail")) %>%
      dplyr::relocate(qc)
  }

  # ===== Build denominator table for symbiont composition =====
  # Keep only passing symbiont assays, then rename the current target as the
  # numerator so the pairwise ratio table can be reassembled downstream.
  denoms <- dff %>%
    dplyr::filter(qc != "Fail") %>%
    dplyr::filter(type == "Symbiont") %>%
    dplyr::select(-c(qc, CT.mean, CT.sd, reps, type)) %>%
    dplyr::rename(numerator = target)

  # Retain target-level replicate and QC metadata so failed or absent assays can
  # be flagged later during the rerun summary step.
  reps <- dff2 %>%
    dplyr::select(Sample.Name, File.Name, target, reps, qc)

  # ===== Sum relative values across symbiont targets =====
  # This step rebuilds target totals from the pairwise ratio structure and keeps
  # only values associated with targets that passed QC.
  denoms2 <- denoms %>%
    tidyr::pivot_longer(
      cols = -c(Sample.Name, File.Name, numerator),
      values_to = "value",
      names_to = "target"
    ) %>%
    dplyr::left_join(targets, by = "target") %>%
    dplyr::filter(type == "Symbiont") %>%
    dplyr::group_by(Sample.Name, File.Name, target) %>%
    dplyr::mutate(
      N = dplyr::n(),
      valSum = sum(value, na.rm = TRUE)
    ) %>%
    dplyr::left_join(reps, by = c("Sample.Name", "File.Name", "target")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      value2 = dplyr::case_when(
        qc == "Pass" ~ valSum
      )
    ) %>%
    dplyr::select(-c(reps, qc, value, valSum, numerator, N)) %>%
    unique()

  # ===== Calculate symbiont proportions =====
  # Convert summed relative values into within-sample symbiont proportions and
  # retain both long and wide forms for downstream summaries.
  symPropsLong <- denoms2 %>%
    dplyr::arrange(target) %>%
    dplyr::filter(type == "Symbiont") %>%
    dplyr::mutate(
      prop = 1 / (value2 + 1)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(Sample.Name, File.Name, target, prop)

  # Extract the dominant symbiont identity and proportion for each sample.
  maxProps <- symPropsLong %>%
    dplyr::group_by(Sample.Name, File.Name) %>%
    dplyr::arrange(-prop, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(prop = round(prop, 3)) %>%
    dplyr::rename(
      maxProp = prop,
      maxSym = target
    ) %>%
    dplyr::filter(!is.na(maxProp))

  # Extract the second-highest symbiont proportion to describe background
  # symbiont signal within each sample.
  backgroundProps <- symPropsLong %>%
    dplyr::group_by(Sample.Name, File.Name) %>%
    dplyr::arrange(-prop, .by_group = TRUE) %>%
    dplyr::slice(2:2) %>%
    dplyr::mutate(prop = round(prop, 3)) %>%
    dplyr::rename(
      backProp = prop,
      backSym = target
    ) %>%
    dplyr::filter(!is.na(backProp))

  # Convert the long-format symbiont proportions to wide format and compute the
  # summed symbiont proportion as a quick internal check.
  symProps <- symPropsLong %>%
    tidyr::pivot_wider(
      names_from = "target",
      values_from = "prop",
      names_prefix = "prop",
      names_sep = ""
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      sumSyms = round(sum(dplyr::across(dplyr::starts_with("prop")), na.rm = TRUE), 3)
    )

  # ===== Calculate symbiont:host signal when host assays are present =====
  # On mixed symbiont/host plates, calculate an aggregate host-normalized signal
  # and append it to the symbiont composition output.
  if (plateType == "symHost") {

    denoms <- dff %>%
      dplyr::filter(qc != "Fail") %>%
      dplyr::select(-c(qc, CT.mean, CT.sd, reps, type)) %>%
      dplyr::rename(numerator = target)

    SHtotals <- denoms %>%
      tidyr::pivot_longer(
        cols = -c(Sample.Name, File.Name, numerator),
        values_to = "value",
        names_to = "target"
      ) %>%
      dplyr::left_join(targets, by = "target") %>%
      dplyr::filter(type == "Coral", !is.na(value)) %>%
      tidyr::pivot_wider(
        names_from = "numerator",
        names_glue = "{numerator}.Host",
        values_from = "value"
      ) %>%
      dplyr::select(-c(type, target)) %>%
      dplyr::group_by(Sample.Name, File.Name) %>%
      dplyr::mutate(
        symHost = sum(dplyr::across(dplyr::ends_with(".Host")), na.rm = TRUE)
      ) %>%
      dplyr::relocate(Sample.Name, File.Name, symHost) %>%
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
  # Prioritize failed targets that likely represent biologically meaningful
  # amplification within each sample, while always flagging failed host targets.

  # Build within-sample Ct context from all amplifying targets. For symbionts,
  # use both the minimum and mean Ct within each sample to identify failed
  # targets that still look relatively strong compared with the rest of the sample.
  tmp_ct_context <- dff2 %>%
    dplyr::filter(reps > 0, !is.na(CT.mean)) %>%
    dplyr::group_by(File.Name, Sample.Name) %>%
    dplyr::mutate(
      sampleMinCtAll = min(CT.mean, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(File.Name, Sample.Name, type) %>%
    dplyr::mutate(
      sampleMinCtType = min(CT.mean, na.rm = TRUE),
      sampleMeanCtType = mean(CT.mean, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      File.Name,
      Sample.Name,
      target,
      type,
      sampleMinCtAll,
      sampleMinCtType,
      sampleMeanCtType,
      isLowestCtAll = CT.mean == sampleMinCtAll,
      isLowestCtType = CT.mean == sampleMinCtType
    )

  # Keep failed targets and classify rerun importance. Failed host assays are
  # always prioritized. Failed symbiont assays are prioritized when they are
  # relatively low-Ct within the sample, particularly when they are the lowest
  # amplifier overall or among symbionts.
  tmp_failed_detail <- dff2 %>%
    dplyr::filter(qc == "Fail") %>%
    dplyr::select(Sample.Name, File.Name, target, type, reps, CT.mean, CT.sd) %>%
    dplyr::left_join(
      tmp_ct_context,
      by = c("File.Name", "Sample.Name", "target", "type")
    ) %>%
    dplyr::mutate(
      isLowestCtAll = dplyr::coalesce(isLowestCtAll, FALSE),
      isLowestCtType = dplyr::coalesce(isLowestCtType, FALSE),

      isImportantSymFail = dplyr::case_when(
        type != "Symbiont" ~ FALSE,
        is.na(CT.mean) ~ FALSE,
        reps == 0 ~ FALSE,
        isLowestCtAll ~ TRUE,
        isLowestCtType ~ TRUE,
        !is.na(sampleMinCtType) & CT.mean <= sampleMinCtType + ctWindow ~ TRUE,
        !is.na(sampleMeanCtType) & CT.mean <= sampleMeanCtType ~ TRUE,
        TRUE ~ FALSE
      ),

      failedReason = dplyr::case_when(
        type == "Coral" ~ "failed host amplification",
        isImportantSymFail & isLowestCtAll ~ "lowest Ct amplification failed QC",
        isImportantSymFail & isLowestCtType ~ "lowest symbiont Ct failed QC",
        isImportantSymFail ~ "important symbiont amplification failed QC",
        reps == 2 & !is.na(CT.sd) & CT.sd >= sdThreshold ~ "duplicate amplification failed SD",
        reps == 1 & singles == FALSE ~ "single amplification not counted",
        TRUE ~ "failed amplification"
      ),

      rerunPriorityNum = dplyr::case_when(
        type == "Coral" ~ 1,
        isImportantSymFail & isLowestCtAll ~ 2,
        isImportantSymFail & isLowestCtType ~ 3,
        isImportantSymFail ~ 4,
        reps == 2 & !is.na(CT.sd) & CT.sd >= sdThreshold ~ 5,
        reps == 1 & singles == FALSE ~ 6,
        TRUE ~ 7
      ),

      rerunPriority = dplyr::case_when(
        rerunPriorityNum %in% c(1, 2, 3, 4) ~ "1_High",
        rerunPriorityNum %in% c(5, 6) ~ "2_Moderate",
        TRUE ~ "3_Low"
      ),

      rerunTargets = target,
      rerunSummary = paste0(failedReason, ": ", target)
    ) %>%
    dplyr::arrange(File.Name, Sample.Name, rerunPriorityNum, CT.mean, target)

  # Summarize failed targets at the sample level while preserving priority order
  # and explicit target identities.
  tmp_failed_summary <- tmp_failed_detail %>%
    dplyr::group_by(File.Name, Sample.Name, rerunPriorityNum, rerunPriority, failedReason) %>%
    dplyr::summarise(
      rerunTargets = paste(target, collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(File.Name, Sample.Name, rerunPriorityNum, failedReason) %>%
    dplyr::group_by(File.Name, Sample.Name) %>%
    dplyr::summarise(
      rerunPriorityNum = min(rerunPriorityNum),
      rerunPriority = paste(rerunPriority, collapse = " | "),
      rerunReason = paste(failedReason, collapse = " | "),
      rerunTargets = paste(rerunTargets, collapse = " | "),
      rerunSummary = paste0(failedReason, ": ", rerunTargets, collapse = " || "),
      .groups = "drop"
    )

  # Detect samples lacking symbiont or host amplification at the sample level.
  noReps <- reps %>%
    dplyr::left_join(targets, by = "target") %>%
    dplyr::group_by(Sample.Name, File.Name, type) %>%
    dplyr::summarise(sum = sum(reps, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      rerunPriorityNum = dplyr::case_when(
        type == "Coral" & sum != 2 ~ 1,
        type == "Symbiont" & sum == 0 ~ 2,
        TRUE ~ NA_real_
      ),
      rerunPriority = dplyr::case_when(
        type == "Coral" & sum != 2 ~ "1_High",
        type == "Symbiont" & sum == 0 ~ "1_High",
        TRUE ~ NA_character_
      ),
      rerunReason = dplyr::case_when(
        type == "Symbiont" & sum == 0 ~ "no symbiont DNA",
        type == "Coral" & sum != 2 ~ "failed host amplification",
        TRUE ~ NA_character_
      ),
      rerunTargets = dplyr::case_when(
        type == "Symbiont" & sum == 0 ~ "sumSyms",
        type == "Coral" & sum != 2 ~ "symHost",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(rerunReason)) %>%
    dplyr::group_by(Sample.Name, File.Name) %>%
    dplyr::mutate(
      nProblems = dplyr::n(),
      omit = dplyr::case_when(
        nProblems == 1 & rerunTargets == "symHost" ~ "Y"
      )
    ) %>%
    dplyr::filter(is.na(omit)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      rerunSummary = paste0(rerunReason, ": ", rerunTargets)
    ) %>%
    dplyr::select(
      File.Name,
      Sample.Name,
      rerunPriorityNum,
      rerunPriority,
      rerunReason,
      rerunTargets,
      rerunSummary
    )

  # Combine detailed target-level failed amplifications with sample-level
  # missing-signal flags for optional QC review.
  reruns_detail <- dplyr::bind_rows(
    tmp_failed_detail %>%
      dplyr::select(
        File.Name,
        Sample.Name,
        rerunPriorityNum,
        rerunPriority,
        rerunReason = failedReason,
        rerunTargets,
        rerunSummary,
        type,
        reps,
        CT.mean,
        CT.sd,
        sampleMinCtAll,
        sampleMinCtType,
        sampleMeanCtType,
        isLowestCtAll,
        isLowestCtType,
        isImportantSymFail
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
        File.Name,
        Sample.Name,
        rerunPriorityNum,
        rerunPriority,
        rerunReason,
        rerunTargets,
        rerunSummary,
        type,
        reps,
        CT.mean,
        CT.sd,
        sampleMinCtAll,
        sampleMinCtType,
        sampleMeanCtType,
        isLowestCtAll,
        isLowestCtType,
        isImportantSymFail
      )
  ) %>%
    dplyr::arrange(File.Name, Sample.Name, rerunPriorityNum, CT.mean, rerunTargets) %>%
    dplyr::left_join(df, by = c("Sample.Name", "File.Name")) %>%
    dplyr::filter(!Sample.Name %in% c("NTC", "-")) %>%
    dplyr::arrange(rerunPriority, File.Name)

  # Collapse rerun recommendations to one row per sample for the standard QC output.
  reruns <- dplyr::bind_rows(
    tmp_failed_summary,
    noReps
  ) %>%
    dplyr::arrange(File.Name, Sample.Name, rerunPriorityNum) %>%
    dplyr::group_by(File.Name, Sample.Name) %>%
    dplyr::summarise(
      rerunPriority = paste(rerunPriority, collapse = " | "),
      rerunReason = paste(rerunReason, collapse = " | "),
      rerunTargets = paste(rerunTargets, collapse = " | "),
      rerunSummary = paste(rerunSummary, collapse = " || "),
      .groups = "drop"
    ) %>%
    dplyr::left_join(df, by = c("Sample.Name", "File.Name")) %>%
    dplyr::filter(!Sample.Name %in% c("NTC", "-")) %>%
    dplyr::arrange(rerunPriority, File.Name)

  # Append a simple QC flag back to the main processed output.
  df3 <- df2 %>%
    dplyr::left_join(
      reruns %>%
        dplyr::mutate(checkQC = "checkQC") %>%
        dplyr::select(File.Name, Sample.Name, checkQC),
      by = c("File.Name", "Sample.Name")
    )

  # ===== Return requested output =====
  if (QC != TRUE & returnAll == FALSE) {
    return(df3)
  } else if (QC != TRUE & returnAll == TRUE) {
    df4 <- df3 %>%
      dplyr::left_join(df, by = dplyr::join_by(Sample.Name, File.Name))
    return(df4)
  } else if (QC == TRUE & rerunDetail == TRUE) {
    return(reruns_detail)
  } else {
    return(reruns)
  }
}
