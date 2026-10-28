################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Comparable Cognitive Performance Regardless of Crowding.......................
################################################################################
library(dplyr)
library(purrr)
library(tibble)
library(car)
library(emmeans)

options(contrasts = c("contr.sum", "contr.poly"))

path_folder <- "/Labache_2026_HeCo"
code_root <- file.path(path_folder, 
                       "Script")
source(file.path(code_root, "functions.R"))

crowding <- read.csv(file.path(path_folder, "Data", 
                               "HCP_987_participants_phenotype_crowding.csv"))[, -1] |>
  drop_index_columns()
cognition <- read.csv(file.path(path_folder, "Data", 
                                "HCP_1206_participants_unrestricted.csv"))[, -1]

cog_map <- c(Total_Composite = "CogTotalComp_Unadj",
             Crystal_Composite = "CogCrystalComp_Unadj",
             Fluid_Composite = "CogFluidComp_Unadj",
             Visual_Episodic_Memory = "PicSeq_Unadj",
             Cognitive_Flexibility_CardSort = "CardSort_Unadj",
             Inhibition_Flanker = "Flanker_Unadj",
             Reading_Decoding = "ReadEng_Unadj",
             Vocabulary_Comprehension = "PicVocab_Unadj",
             Processing_Speed = "ProcSpeed_Unadj",
             Working_Memory_ListSorting = "ListSort_Unadj")

cognition <- cognition[, c("Subject", unname(cog_map))]
names(cognition) <- c("subject", names(cog_map))
cognition$subject <- as.character(cognition$subject)
crowding$subject <- as.character(crowding$subject)

outcomes <- names(cog_map)
data <- inner_join(crowding, cognition, by = "subject") |>
  mutate(crowding_phenotype = factor(crowding_phenotype, 
                                     levels = c("Split", "Mixed", "Crowded")),
         Manual_Preference = factor(Manual_Preference,
                                    levels = c("LH", "RH")),
         Gender = factor(Gender, levels = c("F", "M"))) |>
  filter(if_all(all_of(c(outcomes, "crowding_phenotype", "Manual_Preference", "Age_years",
                         "Educational_years", "Total_Intracranial_Volume_eTIV_mm3", "Gender")),
                ~ !is.na(.x))) |>
  droplevels()
cat("Complete-case N:", nrow(data), "\n")

fits <- setNames(vector("list", length(outcomes)), outcomes)
for (outcome in outcomes) {
  main_formula <- reformulate(c("crowding_phenotype", "Manual_Preference", 
                                "Age_years", "Educational_years",
                                "Total_Intracranial_Volume_eTIV_mm3", 
                                "Gender"),
                              response = outcome)
  interaction_formula <- as.formula(paste0(outcome,
                                           " ~ crowding_phenotype * Manual_Preference + Age_years + Educational_years + ",
                                           "Total_Intracranial_Volume_eTIV_mm3 + Gender"))
  
  main_model <- lm(main_formula, data = data)
  interaction_model <- lm(interaction_formula, data = data)
  marginal_means <- emmeans(main_model, ~ crowding_phenotype)
  
  fits[[outcome]] <- list(main_model = main_model,
                          main_anova = Anova(main_model, type = "III"),
                          interaction_model = interaction_model,
                          interaction_anova = Anova(interaction_model, type = "III"),
                          emmeans = marginal_means,
                          pairwise = pairs(marginal_means, adjust = "tukey"))
}

main_results <- imap_dfr(fits, function(fit, outcome) {
  extract_type3_effect(fit$main_anova, fit$main_model, "crowding_phenotype") |>
    mutate(outcome = outcome, .before = 1)
}) |>
  mutate(p_FDR = p.adjust(p, method = "BH"))

interaction_results <- imap_dfr(fits, function(fit, outcome) {
  extract_type3_effect(fit$interaction_anova,
                       fit$interaction_model,
                       "crowding_phenotype:Manual_Preference") |>
    mutate(outcome = outcome, .before = 1)
}) |>
  mutate(p_FDR = p.adjust(p, method = "BH"))

emmeans_results <- imap_dfr(fits, function(fit, outcome) {
  as.data.frame(summary(fit$emmeans, infer = c(TRUE, TRUE))) |>
    mutate(outcome = outcome, .before = 1)})
pairwise_results <- imap_dfr(fits, function(fit, outcome) {
  as.data.frame(summary(fit$pairwise, infer = c(TRUE, TRUE))) |>
    mutate(outcome = outcome, .before = 1)})

print(main_results)
print(interaction_results)
