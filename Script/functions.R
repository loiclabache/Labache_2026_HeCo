################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Shared helpers for the manuscript analyses....................................
################################################################################
require_columns <- function(data, columns, object_name = deparse(substitute(data))) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(object_name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

drop_index_columns <- function(data) {
  data[, !grepl("^(X|Unnamed:.*)$", names(data)), drop = FALSE]
}

robust_scale <- function(data) {
  data <- as.data.frame(data)
  
  out <- lapply(data, function(x) {
    med <- stats::median(x)
    mad_value <- stats::mad(x, center = med)
    (x - med) / mad_value
  })
  as.data.frame(out, check.names = FALSE)
}

standardize_demographics <- function(data) {
  aliases <- list(Participant = c("Participant", "Subject", "subject", "id"),
                  Age_years = c("Age_years", "Age_in_Yrs"),
                  Educational_years = c("Educational_years", "SSAGA_Educ"),
                  Total_Intracranial_Volume_eTIV_mm3 = c("Total_Intracranial_Volume_eTIV_mm3",
                                                         "FS_IntraCranial_Vol",
                                                         "FS_InterCranial_Vol"),
                  Gender = c("Gender", "sex"),
                  Manual_Preference = c("Manual_Preference", "PM_strict_30ES"),
                  phenotype_SENSAAS = c( "phenotype_SENSAAS",
                                         "Cluster...MMHC..category.",
                                         "Cluster - MMHC (category)",
                                         "Cluster MMHC category"))
  
  for (target in names(aliases)) {
    if (target %in% names(data)) next
    source <- aliases[[target]][aliases[[target]] %in% names(data)][1]
    if (!is.na(source)) names(data)[names(data) == source] <- target
  }
  
  data
}

attention_labels <- function(x) {
  x <- as.character(x)
  out <- dplyr::case_when(x %in% c("1", "Typical", "Strong Typical",
                                   "Typical (A)") ~ "Typical",
                          x %in% c("2", "Bilateral", "Mild Typical", 
                                   "Bilateral (A)") ~ "Bilateral",
                          x %in% c("3", "Atypical", "Weak Typical",
                                   "Atypical (A)") ~ "Atypical",
                          TRUE ~ NA_character_)
  factor(out, levels = c("Typical", "Bilateral", "Atypical"))
}

language_labels <- function(x) {
  x <- as.character(x)
  out <- dplyr::case_when(x %in% c("strong_typical", "Strong Typical", 
                                   "Strong Typical (L)") ~ "Strong Typical",
                          x %in% c("mild_typical", "Mild Typical",
                                   "Mild Typical (L)") ~ "Mild Typical",
                          x %in% c("Atypical", "Atypical (L)") ~ "Atypical",
                          TRUE ~ NA_character_)
  factor(out, levels = c("Strong Typical", "Mild Typical", "Atypical"))
}

build_crowding_phenotype <- function(attention, language) {
  attention <- as.character(attention)
  language <- as.character(language)
  
  out <- dplyr::case_when(attention == "Bilateral" ~ "Mixed",
                          attention == "Typical" & language %in% c("Strong Typical", 
                                                                   "Mild Typical") ~ "Split",
                          attention == "Atypical" & language == "Atypical" ~ "Split",
                          attention == "Typical" & language == "Atypical" ~ "Crowded",
                          attention == "Atypical" & language %in% c("Strong Typical",
                                                                    "Mild Typical") ~ "Crowded",
                          TRUE ~ NA_character_)
  factor(out, levels = c("Split", "Mixed", "Crowded"))
}

load_functional_markers <- function(path_folder, include_demographics = TRUE) {
  rest <- read.csv(file.path(path_folder, "Data",
                             "HCP_987_participants_restMarkers_ALANs.csv"),
                   check.names = FALSE)[, -1] |>
    drop_index_columns()
  task <- read.csv(file.path(path_folder, "Data",
                             "HCP_987_participants_taskMarkers_ALANs.csv"),
                   check.names = FALSE)[, -1] |>
    drop_index_columns()
  
  require_columns(rest, "subject", "resting-state marker file")
  require_columns(task, "Participant", "task marker file")
  task$Participant <- as.character(task$Participant)
  rest$subject <- as.character(rest$subject)
  
  rest_vars <- c("ParietoFrontal_strength_Asym",
                 "TemporoFrontal_strength_Asym",
                 "ParietoFrontal_strength_Sum",
                 "TemporoFrontal_strength_Sum",
                 "ParietoFrontal_ihhc",
                 "TemporoFrontal_ihhc")
  task_vars <- c("PF", "PF_hub", "TF", "TF_hub")
  require_columns(rest, rest_vars, "resting-state marker file")
  require_columns(task, task_vars, "task marker file")
  
  data <- dplyr::inner_join(rest[, c("subject", rest_vars)],
                            task[, c("Participant", task_vars)],
                            by = c("subject" = "Participant"))
  
  left_positive <- c("ParietoFrontal_strength_Asym",
                     "TemporoFrontal_strength_Asym",
                     "PF", "PF_hub", "TF", "TF_hub")
  data[left_positive] <- lapply(data[left_positive], function(x) -x)
  
  if (!include_demographics) return(data)
  
  homotopic_vars <- c("ParietoFrontal_ihhc", "TemporoFrontal_ihhc")

  phenotype <- read.csv(file.path(path_folder, "Data",
                                  "HCP_987_participants.csv"),
                        check.names = FALSE)[, -1] |>
    drop_index_columns()
  require_columns(phenotype, c("subject", "phenotype_ALANs"),
                  "attention phenotype file")
  phenotype$subject <- as.character(phenotype$subject)
  
  demographics <- read.csv(file.path(path_folder, "Data", 
                                     "HCP_994_participants_Demographics_Cognition.csv"),
                           check.names = FALSE)[, -1] |>
    drop_index_columns() |>
    standardize_demographics()
  require_columns(demographics,
                  c("Participant", "Age_years", "Educational_years",
                    "Total_Intracranial_Volume_eTIV_mm3", "Gender", 
                    "Manual_Preference",
                    "phenotype_SENSAAS"),
                  "demographic file")
  demographics$Participant <- as.character(demographics$Participant)
  
  data <- data |>
    dplyr::inner_join(phenotype[, c("subject", "phenotype_ALANs")],
                      by = "subject") |>
    dplyr::inner_join(demographics, by = c("subject" = "Participant")) |>
    dplyr::mutate(attention_phenotype = attention_labels(phenotype_ALANs),
                  language_phenotype = language_labels(phenotype_SENSAAS),
                  Gender = factor(Gender, levels = c("F", "M")),
                  Manual_Preference = factor(Manual_Preference, 
                                             levels = c("LH", "RH")))
  
  data
}

fit_functional_models <- function(data) {
  outcomes <- c("ParietoFrontal_strength_Asym",
                "TemporoFrontal_strength_Asym",
                "ParietoFrontal_strength_Sum",
                "TemporoFrontal_strength_Sum",
                "ParietoFrontal_ihhc",
                "TemporoFrontal_ihhc",
                "PF", "PF_hub", "TF", "TF_hub")
  
  model_data <- data |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(outcomes,
                                                "attention_phenotype", 
                                                "Manual_Preference",
                                                "Age_years",
                                                "Educational_years",
                                                "Total_Intracranial_Volume_eTIV_mm3", 
                                                "Gender")),
                                ~ !is.na(.x))) |>
    droplevels()
  
  cell_counts <- table(model_data$attention_phenotype,
                       model_data$Manual_Preference)
  
  fits <- setNames(vector("list", length(outcomes)), outcomes)
  for (outcome in outcomes) {
    formula <- stats::as.formula(paste0(outcome,
                                        " ~ attention_phenotype * Manual_Preference + ",
                                        "Age_years + Educational_years + ",
                                        "Total_Intracranial_Volume_eTIV_mm3 + Gender"))
    model <- stats::lm(formula, data = model_data)
    anova <- car::Anova(model, type = "III")
    emm_grid <- suppressMessages(emmeans::emmeans(model, ~ attention_phenotype,
                                                  weights = "equal"))
    fits[[outcome]] <- list(model = model,
                            anova = anova,
                            emmeans = emm_grid,
                            pairwise = emmeans::contrast(emm_grid,
                                                         method = "pairwise",
                                                         adjust = "tukey"))
  }
  
  mainEffect <- purrr::imap_dfr(fits, function(fit, outcome) {
    row <- fit$anova["attention_phenotype", , drop = FALSE]
    tibble::tibble(outcome = outcome,
                   df_num = as.numeric(row$Df),
                   df_den = stats::df.residual(fit$model),
                   F = as.numeric(row$`F value`),
                   p = as.numeric(row$`Pr(>F)`))}) |>
    dplyr::mutate(p_FDR = stats::p.adjust(p, method = "BH"))
  
  pairwise <- purrr::imap_dfr(fits, function(fit, outcome) {
    as.data.frame(summary(fit$pairwise, infer = c(TRUE, TRUE))) |>
      dplyr::mutate(outcome = outcome, .before = 1)})
  
  emmeans <- purrr::imap_dfr(fits, function(fit, outcome) {
    as.data.frame(summary(fit$emmeans, infer = c(TRUE, TRUE))) |>
      dplyr::mutate(outcome = outcome, .before = 1)})
  
  interaction <- purrr::imap_dfr(fits, function(fit, outcome) {
    row <- fit$anova[ "attention_phenotype:Manual_Preference", , drop = FALSE]
    tibble::tibble(outcome = outcome,
                   df_num = as.numeric(row$Df),
                   df_den = stats::df.residual(fit$model),
                   F = as.numeric(row$`F value`),
                   p = as.numeric(row$`Pr(>F)`))})
  
  list(data = model_data,
       fits = fits,
       mainEffect = mainEffect,
       interaction = interaction,
       emmeans = emmeans,
       pairwise = pairwise)
}

extract_type3_effect <- function(anova_table, model, term) {
  row <- anova_table[term, , drop = FALSE]
  df_num <- as.numeric(row$Df)
  df_den <- stats::df.residual(model)
  F_value <- as.numeric(row$`F value`)
  tibble::tibble(df_num = df_num,
                 df_den = df_den,
                 F = F_value,
                 p = as.numeric(row$`Pr(>F)`),
                 eta2_partial = (F_value * df_num) / (F_value * df_num + df_den))
}
