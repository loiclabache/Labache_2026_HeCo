################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Intrinsic connectivity distinguishes bilateral attention organization.........
################################################################################
library(dplyr)
library(car)
library(emmeans)
library(purrr)
library(tibble)

options(contrasts = c("contr.sum", "contr.poly"))

path_folder <- "/Labache_2026_HeCo"
code_root <- file.path(path_folder, 
                       "Script")
source(file.path(code_root, "functions.R"))


# Analysis of Covariance........................................................
################################################################################
analysis <- fit_functional_models(load_functional_markers(file.path(path_folder)))
rest_vars <- c("ParietoFrontal_strength_Asym",
               "TemporoFrontal_strength_Asym",
               "ParietoFrontal_strength_Sum",
               "TemporoFrontal_strength_Sum",
               "ParietoFrontal_ihhc",
               "TemporoFrontal_ihhc")

rest_mainEffect <- analysis$mainEffect |> filter(outcome %in% rest_vars)
rest_emmeans <- analysis$emmeans |> filter(outcome %in% rest_vars)
rest_pairwise <- analysis$pairwise |> filter(outcome %in% rest_vars)
rest_interaction <- analysis$interaction |> filter(outcome %in% rest_vars)

print(rest_mainEffect)
print(rest_emmeans)
print(rest_pairwise)
print(rest_interaction)