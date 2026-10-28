################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Heritability of Attention and Crowding Phenotype..............................
################################################################################
library(dplyr)
library(solarius)

path_folder <- "/Labache_2026_HeCo"
code_root <- file.path(path_folder, 
                       "Script")
source(file.path(code_root, "h2_multi.R"))

hcp <- read.csv(file.path(path_folder, "Data", 
                          "HCP_984_participants_FamilyData.csv"),
                check.names = FALSE)

hcp <- hcp |>
  mutate(age2 = Age_in_Yrs^2,
         age2_scale = as.numeric(scale(age2)),
         Age_in_Yrs_scale = as.numeric(scale(Age_in_Yrs)),
         FS_IntraCranial_Vol_scale = as.numeric(scale(FS_InterCranial_Vol)),
         PM_bin = if_else(PM_strict_30ES == "RH", 1, 0)) |>
  arrange(id)

K_big <- solarKinship2(hcp)
row_idx <- match(hcp$id, rownames(K_big))
col_idx <- match(hcp$id, colnames(K_big))
K <- K_big[row_idx, col_idx, drop = FALSE]

X <- model.matrix(~ PM_bin + Age_in_Yrs_scale + sex + age2_scale +
                    FS_IntraCranial_Vol_scale + Age_in_Yrs_scale:sex + age2_scale:sex,
                  data = hcp)

phenotypes <- c(Attention = "ALANs_num", Crowding = "Crowding_num")
results <- lapply(names(phenotypes), function(label) {
  y <- factor(hcp[[phenotypes[[label]]]])
  if (anyNA(y)) stop(label, " phenotype contains missing values.")
  
  P <- outer(y, y, FUN = "==") * 1
  diag(P) <- 1
  
  estimate <- h2_multi_R(P = P, K = K, X = X, n_perm = 1000L, seed = 13L)
  jackknife <- jackknife_family_se(P = P, K = K, X = X, famid = hcp$famid)
  
  data.frame(phenotype = label,
             h2 = estimate$h2,
             family_jackknife_SE = jackknife$se,
             permutation_p = estimate$p_perm)
}) |>
  bind_rows()

print(results)
