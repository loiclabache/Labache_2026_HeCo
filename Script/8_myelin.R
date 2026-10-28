################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Greater Cortical Myelination in Mixed Hemispheric Dominance...................
################################################################################
library(dplyr)
library(car)
library(emmeans)
library(progress)
library(VertexWiseR)

options(contrasts = c("contr.sum", "contr.poly"))

path_folder <- "/Labache_2026_HeCo"
code_root <- file.path(path_folder, 
                       "Script")
source(file.path(code_root, "functions.R"))

crowding <- read.csv(file.path(path_folder, "Data", 
                               "HCP_987_participants_phenotype_crowding.csv")) |>
  drop_index_columns()

rest_ids <- read.csv(file.path(path_folder, "Data", 
                               "HCP_987_participants_restMarkers_ALANs.csv"))$subject
task_ids <- read.csv(file.path(path_folder, "Data", 
                               "HCP_987_participants_taskMarkers_ALANs.csv"))$Participant
source_ids <- as.character(intersect(rest_ids, task_ids))

crowding <- crowding |>
  mutate(subject = as.character(subject),
         crowding_phenotype = factor(crowding_phenotype, levels = c("Split", "Mixed", "Crowded")),
         Gender = factor(Gender, levels = c("F", "M")),
         Manual_Preference = factor(Manual_Preference, levels = c("LH", "RH"))) |>
  filter(subject %in% source_ids,
         subject %in% source_ids,
         complete.cases(crowding_phenotype, Age_years, Educational_years,
                        Total_Intracranial_Volume_eTIV_mm3, Gender, Manual_Preference)) |>
  droplevels()
cat("Myelin-analysis N:", nrow(crowding), "\n")

read_surface_data <- function(file, fallback_ids) {
  object <- readRDS(file)
  if (is.list(object) && length(object) >= 2L) {
    ids <- as.character(unlist(object[[1]], use.names = FALSE))
    data <- as.matrix(object[[2]])
    if (length(ids) != nrow(data)) stop("Embedded surface IDs do not match matrix rows: ", basename(file))
    return(list(ids = ids, data = data))
  }
  data <- as.matrix(object)
  if (nrow(data) != length(fallback_ids)) {
    stop("Legacy surface-matrix rows do not match the participant order: ", basename(file))
  }
  list(ids = fallback_ids, data = data)
}

myelin_object <- read_surface_data(file.path(path_folder, "Data", 
                                             "HCP_987_participants_SmoothedMyelinMap.rds"),
                                   source_ids)
curvature_object <- read_surface_data(file.path(path_folder, "Data", 
                                                "HCP_987_participants_Curvature.rds"),
                                      source_ids)
thickness_object <- read_surface_data(file.path(path_folder, "Data", 
                                                "HCP_987_participants_Thickness.rds"),
                                      source_ids)

myelin <- myelin_object$data
curvature <- curvature_object$data
thickness <- thickness_object$data
surface_ids <- myelin_object$ids

row_index <- match(crowding$subject, surface_ids)
myelin <- myelin[row_index, , drop = FALSE]
curvature <- curvature[row_index, , drop = FALSE]
thickness <- thickness[row_index, , drop = FALSE]


# Residualize myelin for local cortical thickness and curvature at each vertex..
################################################################################
myelin_resid <- matrix(NA_real_, nrow(myelin), ncol(myelin))
pb <- progress_bar$new(format = "Residualizing [:bar] :percent in :elapsed",
                       total = ncol(myelin), clear = FALSE, width = 60)
for (vertex in seq_len(ncol(myelin))) {
  fit <- lm(myelin[, vertex] ~ thickness[, vertex] + curvature[, vertex])
  myelin_resid[, vertex] <- residuals(fit)
  pb$tick()
}

crowding$crowding_trend <- dplyr::recode(as.character(crowding$crowding_phenotype),
                                         "Crowded" = -1,
                                         "Split" = 0,
                                         "Mixed" = 1)
crowding$Total_Intracranial_Volume_eTIV_mm3 <- crowding$Total_Intracranial_Volume_eTIV_mm3 / 1e6

covariates <- crowding[, c("Age_years", "Gender", "Total_Intracranial_Volume_eTIV_mm3",
                           "Educational_years", "Manual_Preference", "crowding_trend")]

rft_model <- RFT_vertex_analysis(model = covariates,
                                 contrast = crowding$crowding_trend,
                                 surf_data = myelin_resid,
                                 atlas = 3,
                                 p = 0.05)
print(rft_model$cluster_level_results)


# "Post-Hoc" Analysis...........................................................
################################################################################
crowding$cluster_myelin <- as.numeric(myelin_resid %*% rft_model$pos_mask) /
  sum(rft_model$pos_mask > 0)
followup_model <- lm(cluster_myelin ~ crowding_phenotype + Age_years + Gender + 
                       Educational_years + Total_Intracranial_Volume_eTIV_mm3 + 
                       Manual_Preference,
                     data = crowding)
followup_anova <- Anova(followup_model, type = "III")
followup_emmeans <- emmeans(followup_model, ~ crowding_phenotype)
followup_pairs <- pairs(followup_emmeans, adjust = "tukey")

print(followup_anova)
print(followup_emmeans)
print(followup_pairs)
