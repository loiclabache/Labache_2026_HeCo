################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Attention Phenotype Associations with Self-Reported Sex and Handedness........
################################################################################
library(dplyr)
library(nnet)
library(car)
library(marginaleffects)

options(contrasts = c("contr.sum", "contr.poly"))

path_folder <- "/Labache_2026_HeCo"
code_root <- file.path(path_folder, 
                       "Script")
source(file.path(code_root, "functions.R"))


# Load Data.....................................................................
################################################################################
source_data <- load_functional_markers(file.path(path_folder),
                                       include_demographics = TRUE)
complete_fields <- c("ParietoFrontal_strength_Asym",
                     "TemporoFrontal_strength_Asym",
                     "ParietoFrontal_strength_Sum",
                     "TemporoFrontal_strength_Sum",
                     "ParietoFrontal_ihhc",
                     "TemporoFrontal_ihhc",
                     "PF", "PF_hub", "TF", "TF_hub", 
                     "attention_phenotype", 
                     "language_phenotype", "Gender", "Manual_Preference", 
                     "Age_years", "Educational_years", 
                     "Total_Intracranial_Volume_eTIV_mm3")
analysis_data <- source_data |>
  filter(if_all(all_of(complete_fields), ~ !is.na(.x))) |>
  transmute(subject,
            attention_phenotype = droplevels(attention_phenotype),
            Gender = relevel(factor(Gender), ref = "F"),
            Manual_Preference = relevel(factor(Manual_Preference), ref = "RH")) |>
  droplevels()
cat("Multinomial-model N:", nrow(analysis_data), "\n")


# Model with Interaction........................................................
################################################################################
model_interaction <- multinom(attention_phenotype ~ Gender * Manual_Preference,
                              data = analysis_data,
                              trace = FALSE)
interaction_test <- Anova(model_interaction, type = "III")
print(interaction_test)


# Model without Interaction.....................................................
################################################################################
contrasts(analysis_data$Gender) <- contr.treatment(levels(analysis_data$Gender), base = 1)
contrasts(analysis_data$Manual_Preference) <- contr.treatment(levels(analysis_data$Manual_Preference), base = 1)

model_additive <- multinom(attention_phenotype ~ Gender + Manual_Preference,
                           data = analysis_data,
                           trace = FALSE)
additive_test <- Anova(model_additive, type = "II")
print(additive_test)

beta <- coef(model_additive)
se <- summary(model_additive)$standard.errors
coef_table <- expand.grid(term = colnames(beta),
                          outcome = rownames(beta),
                          stringsAsFactors = FALSE) |>
  mutate(estimate = as.vector(t(beta)),
         SE = as.vector(t(se)),
         z = estimate / SE,
         p = 2 * pnorm(-abs(z)),
         RRR = exp(estimate),
         CI_low = exp(estimate - qnorm(0.975) * SE),
         CI_high = exp(estimate + qnorm(0.975) * SE))
print(coef_table)

pred_gender <- avg_predictions(model_additive,
                               newdata = analysis_data,
                               variables = list(Gender = c("F", "M")),
                               type = "probs",
                               conf_level = 0.95) |>
  as.data.frame()

pred_handedness <- avg_predictions(model_additive,
  newdata = analysis_data,
  variables = list(Manual_Preference = c("LH", "RH")),
  type = "probs",
  conf_level = 0.95) |>
  as.data.frame()

print(pred_gender)
print(pred_handedness)