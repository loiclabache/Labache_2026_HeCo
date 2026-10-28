################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Task-evoked asymmetry differentiates the attention phenotypes.................
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
task_vars <- c("PF", "PF_hub", "TF", "TF_hub")

task_mainEffect <- analysis$mainEffect |> filter(outcome %in% task_vars)
task_emmeans <- analysis$emmeans |> filter(outcome %in% task_vars)
task_pairwise <- analysis$pairwise |> filter(outcome %in% task_vars)
task_interaction <- analysis$interaction |> filter(outcome %in% task_vars)

print(task_mainEffect)
print(task_emmeans)
print(task_pairwise)
print(task_interaction)