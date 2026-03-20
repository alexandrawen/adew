#' @importFrom magrittr "%>%"
#' @export
qPCRcalc <- function(df, sdThreshold, QC = FALSE, ctCutoff = 40, singles = FALSE, returnAll = FALSE) {
  #reconstruct list of targets from dataExtract and classify as host or symbiont
  targets <- data.frame(target = strex::str_before_last_dot(colnames(df %>%
                                                                       dplyr::select(c(contains(".reps")))))) %>%
    dplyr::mutate(type = dplyr::case_when(target %in% c("A","B","C","D","E","F","G","a","b","c","d","e","f","g","Sym","sym","symbiont") ~ "Symbiont",
                                          TRUE ~ "Coral"))

  if (length(unique(targets$type)) == 1) {
    plateType = "symComm"
  } else{
    plateType = "symHost"}

  #determine pass or fail for each target
  if (is.na(ctCutoff)) {
    dff2 <- df %>%
      tidyr::pivot_longer(cols = -c(Sample.Name,File.Name)) %>%
      tidyr::separate(name, into = c("target","calc"),extra="merge") %>%
      dplyr::full_join(targets, by = "target") %>%
      dplyr::mutate(value = dplyr::case_when(is.nan(value) ~ NA,
                                             TRUE ~ value)) %>%
      dplyr::filter(!is.na(target) & target != "NA") %>%
      dplyr::group_by(Sample.Name,File.Name,target,calc) %>%
      dplyr::mutate(value = mean(value, na.rm=TRUE)) %>%
      dplyr::ungroup() %>%
      unique() %>%
      dplyr::filter(value >= 0) %>%
      tidyr::pivot_wider(names_from = "calc",
                         values_from = "value") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(qc = dplyr::case_when((reps == 2 & CT.sd < sdThreshold) ~ "Pass",
                                          reps == 0 ~ "NotPresent",
                                          TRUE ~ "Fail")) %>%
      dplyr::filter(type == "Symbiont" |
                      type == "Coral" & qc != "Fail") %>%
      dplyr::relocate(qc)

    dff <- dff2 %>%
      dplyr::filter(reps > 0)
  } else {
    dff2 <- df %>%
      tidyr::pivot_longer(cols = -c(Sample.Name,File.Name)) %>%
      tidyr::separate(name, into = c("target","calc"),extra="merge") %>%
      dplyr::full_join(targets, by = "target") %>%
      dplyr::mutate(value = case_when(is.nan(value) ~ NA,
                                      TRUE ~ value)) %>%
      dplyr::filter(!is.na(target) & target != "NA") %>%
      dplyr::group_by(Sample.Name,File.Name,target,calc) %>%
      dplyr::mutate(value = mean(value, na.rm=TRUE)) %>%
      dplyr::ungroup() %>%
      unique() %>%
      dplyr::filter(value > 0) %>%
      tidyr::pivot_wider(names_from = "calc",
                         values_from = "value") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(qc = dplyr::case_when((reps == 2 & CT.sd < sdThreshold & CT.mean <= ctCutoff) ~ "Pass",
                                          reps == 0 ~ "NotPresent",
                                          TRUE ~ "Fail")) %>%
      dplyr::filter(type == "Symbiont" |
                      type == "Coral" & qc != "Fail") %>%
      dplyr::relocate(qc)

    dff <- dff2 %>%
      dplyr::filter(reps > 0)
  }


  singles = singles
  if (singles == TRUE) {
    dff <- df %>%
      tidyr::pivot_longer(cols = -c(Sample.Name,File.Name)) %>%
      tidyr::separate(name, into = c("target","calc"),extra="merge") %>%
      dplyr::full_join(targets, by = "target") %>%
      tidyr::pivot_wider(names_from = "calc",
                         values_from = "value") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(qc = dplyr::case_when((reps > 0) ~ "Pass",
                                          reps == 0 ~ "NotPresent",
                                          TRUE ~ "Fail")) %>%
      dplyr::filter(type == "Symbiont" |
                      type == "Coral" & qc == "Pass") %>%
      dplyr::relocate(qc)

    dff2 <- df %>%
      tidyr::pivot_longer(cols = -c(Sample.Name,File.Name)) %>%
      tidyr::separate(name, into = c("target","calc"),extra="merge") %>%
      dplyr::full_join(targets, by = "target") %>%
      dplyr::mutate(value = case_when(is.nan(value) ~ NA,
                                      TRUE ~ value)) %>%
      dplyr::filter(!is.na(target) & target != "NA") %>%
      dplyr::group_by(Sample.Name,File.Name,target,calc) %>%
      dplyr::mutate(value = mean(value, na.rm=TRUE)) %>%
      dplyr::ungroup() %>%
      unique() %>%
      dplyr::filter(value > 0) %>%
      tidyr::pivot_wider(names_from = "calc",
                         values_from = "value") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(qc = dplyr::case_when((reps == 1 & CT.mean <= ctCutoff) ~ "Pass",
                                          (reps == 2 & CT.sd < sdThreshold & CT.mean <= ctCutoff) ~ "Pass",
                                          reps == 0 ~ "NotPresent",
                                          TRUE ~ "Fail")) %>%
      dplyr::filter(type == "Symbiont" |
                      type == "Coral" & qc != "Fail") %>%
      dplyr::relocate(qc)
  }

  #aggregate all ratio denominators
  denoms <- dff %>%
    dplyr::filter(qc != "Fail") %>%
    dplyr::filter(type == "Symbiont") %>%
    dplyr::select(-c(qc,CT.mean,CT.sd,reps,type)) %>%
    dplyr::rename(numerator = target)

  #check for single dominance
  reps <- dff2 %>%
    dplyr::select(c(Sample.Name,File.Name,target,reps,qc))

  #obtain totals for symbiont relative values
  denoms2 <- denoms %>%
    tidyr::pivot_longer(cols = -c(Sample.Name,File.Name,numerator),
                        values_to = "value",
                        names_to = "target") %>%
    dplyr::left_join(targets, by = "target") %>%
    dplyr::filter(type == "Symbiont") %>%
    dplyr::group_by(Sample.Name,File.Name,target) %>%
    dplyr::mutate(N = dplyr::n(),
                  valSum = sum(value, na.rm=TRUE)) %>%
    dplyr::left_join(reps, by = c("Sample.Name","File.Name","target")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(value2 = dplyr::case_when(qc == "Pass" ~ valSum)) %>%
    dplyr::select(-c(reps,qc,value,valSum,numerator,N)) %>%
    unique()

  #calculate symbiont proportions relative to one another
  symPropsLong <- denoms2 %>%
    dplyr::arrange(target) %>%
    dplyr::filter(type == "Symbiont") %>%
    dplyr::mutate(#value = na_if(value, "0"),
      prop = (1/(value2+1))) %>%
    dplyr::ungroup() %>%
    dplyr::select(c(Sample.Name,File.Name,target,prop))

  maxProps <- symPropsLong %>%
    dplyr::group_by(Sample.Name,File.Name) %>%
    dplyr::arrange(-prop, .by_group=TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(prop = round(prop,3)) %>%
    dplyr::rename(maxProp = prop,
                  maxSym = target) %>%
    dplyr::filter(!is.na(maxProp))

  backgroundProps <- symPropsLong %>%
    dplyr::group_by(Sample.Name,File.Name) %>%
    dplyr::arrange(-prop, .by_group=TRUE) %>%
    dplyr::slice(2:2) %>%
    dplyr::mutate(prop = round(prop,3)) %>%
    dplyr::rename(backProp = prop,
                  backSym = target) %>%
    dplyr::filter(!is.na(backProp))

  symProps <- symPropsLong %>%
    tidyr::pivot_wider(names_from = "target",
                       values_from = "prop",
                       names_prefix = "prop",
                       names_sep = "") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sumSyms = round(sum(dplyr::across(starts_with("prop")), na.rm=TRUE),3))

  if (plateType == "symHost") {
    #calculate symbiont:host proportions for each symbiont
    denoms <- dff %>%
      dplyr::filter(qc != "Fail") %>%
      # dplyr::filter(type == "Symbiont") %>%
      dplyr::select(-c(qc,CT.mean,CT.sd,reps,type)) %>%
      dplyr::rename(numerator = target)
    #obtain totals for symbiont:host values
    SHtotals <- denoms %>%
      tidyr::pivot_longer(cols = -c(Sample.Name,File.Name,numerator),
                          values_to = "value",
                          names_to = "target") %>%
      dplyr::left_join(targets, by = "target") %>%
      dplyr::filter(type == "Coral",
                    !is.na(value)) %>%
      tidyr::pivot_wider(names_from = "numerator",
                         names_glue = "{numerator}.Host",
                         values_from = "value") %>%
      dplyr::select(-c(type,target)) %>%
      dplyr::group_by(Sample.Name,File.Name) %>%
      dplyr::mutate(symHost = sum(dplyr::across(dplyr::ends_with(".Host")), na.rm=TRUE)) %>%
      dplyr::relocate(Sample.Name,File.Name,symHost) %>%
      dplyr::ungroup()

    df2 <- symProps %>%
      dplyr::full_join(SHtotals, by = c("Sample.Name", "File.Name")) %>%
      dplyr::full_join(maxProps, by = c("Sample.Name", "File.Name")) %>%
      dplyr::full_join(backgroundProps, by = c("Sample.Name", "File.Name"))

  } else {
    df2 <- dplyr::full_join(symProps,maxProps, by = c("Sample.Name", "File.Name")) %>%
      dplyr::full_join(backgroundProps, by = c("Sample.Name", "File.Name"))
  }

  #return to this later
  failedqc <- dff2 %>%
    dplyr::filter(qc == "Fail") %>%
    dplyr::select(c(Sample.Name:reps,CT.sd)) %>%
    dplyr::select(-c(type)) %>%
    dplyr::mutate(Reason = dplyr::case_when(reps == "1" ~ "singleRep",
                                            !is.na(CT.sd) ~ "SD",
                                            TRUE ~ NA)) %>%
    dplyr::select(c(Sample.Name,File.Name,target,Reason))

  noReps <- reps %>%
    dplyr::left_join(targets) %>%
    dplyr::group_by(Sample.Name,File.Name,type) %>%
    dplyr::summarise(sum = sum(reps, na.rm=T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Reason = dplyr::case_when((type == "Symbiont" & sum == 0) ~ "noSymDNA",
                                            (type == "Coral" & sum != 2) ~ "noHostDNA"),
                  target = dplyr::case_when((type == "Symbiont" & sum == 0) ~ "sumSyms",
                                            (type == "Coral" & sum != 2) ~ "symHost")) %>%
    dplyr::select(-c(sum,type)) %>%
    dplyr::filter(!is.na(Reason)) %>%
    dplyr::group_by(Sample.Name,File.Name) %>%
    dplyr::mutate(nProblems = dplyr::n(),
                  omit = dplyr::case_when(nProblems == 1 & target == "symHost" ~ "Y")) %>%
    dplyr::filter(is.na(omit)) %>%
    dplyr::select(-c(nProblems,omit))

  if(nrow(failedqc)!=0){
    rerunSummary <- dff %>%
      dplyr::filter(Sample.Name %in% c(failedqc$Sample.Name)) %>%
      dplyr::group_by(File.Name,Sample.Name,type) %>%
      dplyr::mutate(minCT = min(CT.mean, na.rm=TRUE),
                    minCT = dplyr::case_when(minCT == Inf ~ NA,
                                             minCT == CT.mean ~ target,
                                             TRUE ~ NA),
                    Summary = "minCT") %>%
      # suppressWarnings(dplyr::mutate(minCT = min(CT.mean, na.rm=TRUE)),
      #                                minCT = dplyr::case_when(minCT == Inf ~ NA,
      #                                                         minCT == CT.mean ~ target,
      #                                                         TRUE ~ NA),
      #                                Summary = "minCT")) %>%
      dplyr::ungroup() %>%
      dplyr::select(c(File.Name,Sample.Name,target = minCT,Summary)) %>%
      dplyr::filter(!is.na(target)) %>%
      unique() %>%
      dplyr::right_join(failedqc) %>%
      dplyr::mutate(Summary = dplyr::case_when(Summary == "minCT" ~ "failedTarget = firstAmplified",
                                               TRUE ~ NA)) %>%
      dplyr::select(c(File.Name,Sample.Name,Summary)) %>%
      dplyr::filter(!is.na(Summary))

    reruns <- plyr::rbind.fill(failedqc,noReps) %>%
      dplyr::arrange(target) %>%
      tidyr::pivot_wider(names_from = target,
                         values_from = Reason) %>%
      dplyr::left_join(rerunSummary) %>%
      dplyr::left_join(df, by = c("Sample.Name","File.Name")) %>%
      dplyr::filter(!Sample.Name %in% c("NTC","-")) %>%
      dplyr::arrange(File.Name)

    df3 <- df2 %>%
      dplyr::mutate(checkQC = case_when(File.Name %in% reruns$File.Name
                                        & Sample.Name %in% reruns$Sample.Name ~ "checkQC",
                                        TRUE ~ NA))
  }else{
    reruns <- plyr::rbind.fill(failedqc,noReps) %>%
      dplyr::arrange(target) %>%
      tidyr::pivot_wider(names_from = target,
                         values_from = Reason) %>%
      dplyr::left_join(rerunSummary) %>%
      dplyr::left_join(df, by = c("Sample.Name","File.Name")) %>%
      dplyr::filter(!Sample.Name %in% c("NTC","-")) %>%
      dplyr::arrange(File.Name)

    df3 <- df2 %>%
      dplyr::left_join(reruns %>%
                         dplyr::mutate(checkQC = "checkQC") %>%
                         dplyr::select(c(File.Name,Sample.Name,checkQC)))
  }



  QC = QC
  returnAll = returnAll
  # df2
  # reruns
  if (QC != TRUE & returnAll == FALSE) {
    return(df3)
  }else if (QC != TRUE & returnAll == TRUE){
    df4 <- df3 %>%
      dplyr::left_join(df, by = join_by(Sample.Name, File.Name))
    return(df4)
  }else {
    return(reruns)
  }
}
