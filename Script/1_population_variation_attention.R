################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Population Variation in Sustained-Attention Lateralization....................
################################################################################
library(dplyr)
library(NbClust)

options(contrasts = c("contr.sum", "contr.poly"))

path_folder <- "/Labache_2026_HeCo"
code_root <- file.path(path_folder, 
                       "Script")
source(file.path(code_root, "functions.R"))


# Read Data.....................................................................
################################################################################
sent_core <- read.csv(file.path(path_folder, 
                                "Atlas", "SENSAAS_description.csv"))
alans <- read.csv(file.path(path_folder, 
                            "Atlas", "ALANs_description.txt"))


# Overlap Language & Attention..................................................
################################################################################
dice <- function(a, b) {
  a <- unique(a[!is.na(a)])
  b <- unique(b[!is.na(b)])
  if (!length(a) && !length(b)) return(NA_real_)
  2 * length(intersect(a, b)) / (length(a) + length(b))
}

language_right <- sent_core |>
  filter(Network == "SENT_CORE", Hemisphere == "Right") |>
  pull(Index)
pf_right <- alans |>
  filter(Network == "ParietoFrontal", Atlas_Homotopic == "Atlas", 
         Hemisphere == "Right") |>
  pull(Index)
tf_right <- alans |>
  filter(Network == "TemporoFrontal", Atlas_Homotopic == "Atlas", 
         Hemisphere == "Right") |>
  pull(Index)

overlap <- data.frame(comparison = c("Language vs full attention",
                                     "Language vs parieto-frontal",
                                     "Language vs temporo-frontal"),
                      dice = c(dice(language_right, c(pf_right, tf_right)), 
                               dice(language_right, pf_right), 
                               dice(language_right, tf_right)))
print(overlap)


# Hierarchical Classification...................................................
################################################################################
markers <- load_functional_markers(file.path(path_folder),
                                   include_demographics = TRUE)
cluster_vars <- c("ParietoFrontal_strength_Asym",
                  "TemporoFrontal_strength_Asym",
                  "ParietoFrontal_strength_Sum",
                  "TemporoFrontal_strength_Sum",
                  "ParietoFrontal_ihhc",
                  "TemporoFrontal_ihhc",
                  "PF", "PF_hub", "TF", "TF_hub")

markers_scaled <- robust_scale(markers[, cluster_vars])
nb <- NbClust(markers_scaled, distance = "euclidean",
              min.nc = 2, max.nc = ceiling((dim(markers_scaled)[1])^0.3),
              method = "ward.D2", index = "all")
recommended_k <- as.integer(names(which.max(table(nb$Best.nc["Number_clusters", ]))))
cat("NbClust majority recommendation:", recommended_k, "clusters\n")

hc <- hclust(dist(markers_scaled, method = "euclidean"), method = "ward.D2")
cluster_id <- cutree(hc, k = 3L)

task_vars <- c("PF", "PF_hub", "TF", "TF_hub")
profile <- data.frame(cluster_id = cluster_id,
                      task_score = rowMeans(markers_scaled[, task_vars])) |>
  group_by(cluster_id) |>
  summarise(task_score = mean(task_score), n = n(), .groups = "drop")

typical_id <- profile$cluster_id[which.min(profile$task_score)]
atypical_id <- profile$cluster_id[which.max(profile$task_score)]
bilateral_id <- setdiff(profile$cluster_id, c(typical_id, atypical_id))

label_map <- c(setNames("Typical", typical_id),
               setNames("Bilateral", bilateral_id),
               setNames("Atypical", atypical_id))
attention_phenotype <- factor(label_map[as.character(cluster_id)], 
                              levels = c("Typical", "Bilateral", "Atypical"))

phenotype <- data.frame(subject = markers$subject,
                        phenotype_ALANs = match(attention_phenotype,
                                                c("Typical", 
                                                  "Bilateral",
                                                  "Atypical")),
                        attention_phenotype = attention_phenotype)

phenotype_counts <- as.data.frame(table(phenotype$attention_phenotype),
                                  stringsAsFactors = FALSE)
names(phenotype_counts) <- c("attention_phenotype", "n")
print(phenotype_counts)
print(profile)

# Gender & Handedness by attention phenotype....................................
################################################################################
gender_var <- grep("sex|gender", names(markers), value = TRUE, ignore.case = TRUE)[1]
hand_var <- grep("hand|manual", names(markers), value = TRUE, ignore.case = TRUE)[1]
phenotype_demographics <- data.frame(subject = markers$subject,
                                     attention_phenotype = attention_phenotype,
                                     gender = markers[[gender_var]],
                                     handedness = markers[[hand_var]])
gender_table <- phenotype_demographics |>
  filter(!is.na(gender)) |>
  count(attention_phenotype, gender, name = "n") |>
  group_by(attention_phenotype) |>
  mutate(value = sprintf("%d (%.1f%%)", n, 100 * n / sum(n)),
         Characteristic = as.character(gender)) |>
  ungroup() |>
  select(attention_phenotype, Characteristic, value)
handedness_table <- phenotype_demographics |>
  filter(!is.na(handedness)) |>
  mutate(left_handed = grepl("LH", as.character(handedness), ignore.case = TRUE)) |>
  group_by(attention_phenotype) |>
  summarise(n = sum(left_handed), N = n(), .groups = "drop") |>
  mutate(value = sprintf("%d (%.1f%%)", n, 100 * n / N),
         Characteristic = "Left-handed") |>
  select(attention_phenotype, Characteristic, value)
demographic_table <- bind_rows(gender_table, handedness_table) |>
  mutate(attention_phenotype = factor(attention_phenotype,
                                      levels = c("Typical", "Bilateral", "Atypical"))) |>
  tidyr::pivot_wider(names_from = attention_phenotype, values_from = value)
print(demographic_table)