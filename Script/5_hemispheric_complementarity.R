################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Population Variation in Hemispheric Complementarity...........................
################################################################################
library(dplyr)
library(tidyr)


# Data..........................................................................
################################################################################
path_folder <- "/Labache_2026_HeCo"
code_root <- file.path(path_folder, 
                       "Script")
source(file.path(code_root, "functions.R"))

phenotype <- read.csv(file.path(path_folder, "Data", 
                                "HCP_987_participants.csv"),
                      check.names = FALSE)[, -1] |>
  drop_index_columns()
demographics <- read.csv(file.path(path_folder, "Data", 
                                   "HCP_994_participants_Demographics_Cognition.csv"),
                         check.names = FALSE)[, -1] |>
  drop_index_columns() |>
  standardize_demographics()

require_columns(phenotype, c("subject", "phenotype_ALANs"), "attention phenotype file")
require_columns(demographics,
                c("Participant", "phenotype_SENSAAS", "Age_years", "Educational_years",
                  "Total_Intracranial_Volume_eTIV_mm3", "Gender", "Manual_Preference"),
                "demographic file")

analysis_data <- inner_join(transform(phenotype, subject = as.character(subject)),
                            transform(demographics, Participant = as.character(Participant)),
                            by = c("subject" = "Participant")) |>
  mutate(attention_phenotype = attention_labels(phenotype_ALANs),
         language_phenotype = language_labels(phenotype_SENSAAS)) |>
  filter(!is.na(attention_phenotype), !is.na(language_phenotype))


# Analysis: 3 by 3..............................................................
################################################################################
tab_3x3 <- table(analysis_data$attention_phenotype, 
                 analysis_data$language_phenotype)
chi_3x3 <- chisq.test(tab_3x3)
cramers_v <- sqrt(unname(chi_3x3$statistic) / (sum(tab_3x3) * min(nrow(tab_3x3) - 1, ncol(tab_3x3) - 1)))

expected <- chi_3x3$expected
row_prop <- rowSums(tab_3x3) / sum(tab_3x3)
col_prop <- colSums(tab_3x3) / sum(tab_3x3)
adj_resid <- (tab_3x3 - expected) / sqrt(expected * (1 - row_prop) %*% t(1 - col_prop))
cell_p <- 2 * pnorm(-abs(adj_resid))
cell_fdr <- matrix(p.adjust(as.vector(cell_p), method = "BH"),
                   nrow = nrow(tab_3x3),
                   dimnames = dimnames(tab_3x3))
cell_results <- as.data.frame(as.table(adj_resid)) |>
  rename(attention_phenotype = Var1, language_phenotype = Var2, adjusted_residual = Freq) |>
  mutate(p = as.vector(cell_p),
         p_FDR = as.vector(cell_fdr),
         direction = case_when(
           p_FDR < 0.05 & adjusted_residual > 0 ~ "enriched",
           p_FDR < 0.05 & adjusted_residual < 0 ~ "depleted",
           TRUE ~ "not significant"))

print(tab_3x3)
print(chi_3x3)
cat("Cramer's V:", cramers_v, "\n")
print(cell_results)


# Analysis: 2 by 2 - no bilateral individuals...................................
################################################################################
direct_data <- analysis_data |>
  filter(attention_phenotype != "Bilateral") |>
  mutate(attention_dominance = if_else(attention_phenotype == "Typical", "Right", "Left"),
         language_dominance = if_else(language_phenotype == "Atypical", "Right", "Left"))
direct_table <- table(direct_data$attention_dominance, direct_data$language_dominance)
direct_chi <- chisq.test(direct_table, correct = FALSE)
phi <- sqrt(unname(direct_chi$statistic) / sum(direct_table))
# Odds ratio for complementary vs same-hemisphere organization..................
a <- direct_table["Right", "Left"]
b <- direct_table["Right", "Right"]
c <- direct_table["Left", "Left"]
d <- direct_table["Left", "Right"]
odds_ratio <- (a * d) / (b * c)
se_log_or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
odds_ratio_ci <- exp(log(odds_ratio) + c(-1, 1) * qnorm(0.975) * se_log_or)

print(direct_table)
print(direct_chi)
cat("Phi:", phi, "\nOdds ratio:", odds_ratio, "\n95% CI:", odds_ratio_ci, "\n")

analysis_data <- analysis_data |>
  mutate(crowding_phenotype = build_crowding_phenotype(attention_phenotype, language_phenotype),
         Gender = factor(Gender, levels = c("F", "M")),
         Manual_Preference = factor(Manual_Preference, levels = c("LH", "RH")))
print(table(analysis_data$crowding_phenotype, useNA = "ifany"))

crowding_export <- analysis_data |>
  transmute(subject,
            crowding_phenotype,
            attention_phenotype,
            language_phenotype,
            Age_years,
            Total_Intracranial_Volume_eTIV_mm3,
            Gender,
            Manual_Preference,
            Educational_years)
head(crowding_export)

association_summary <- data.frame(test = c("3x3 Pearson chi-square",
                                           "Directional 2x2 Pearson chi-square"),
                                  statistic = c(unname(chi_3x3$statistic), 
                                                unname(direct_chi$statistic)),
                                  df = c(unname(chi_3x3$parameter),
                                         unname(direct_chi$parameter)),
                                  p = c(chi_3x3$p.value,
                                        direct_chi$p.value),
                                  effect_size = c(cramers_v, phi),
                                  effect_size_name = c("Cramer's V", "Phi"))
print(association_summary)
odds_ratio_summary <- data.frame(odds_ratio = odds_ratio,
                                 CI_low = odds_ratio_ci[1],
                                 CI_high = odds_ratio_ci[2])
print(odds_ratio_summary)