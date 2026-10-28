################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Crowded Organization Is Associated with Distributed Increases in Fiber Density
################################################################################
library(dplyr)
library(tidyr)
library(R.matlab)
library(car)
library(emmeans)
library(progress)

options(contrasts = c("contr.sum", "contr.poly"))

path_folder <- "/Labache_2026_HeCo"
code_root <- file.path(path_folder, "Script")
source(file.path(code_root, "functions.R"))

crowding <- read.csv(file.path(path_folder, "Data", 
                               "HCP_987_participants_phenotype_crowding.csv"))[, -1] |>
  drop_index_columns()
require_columns(crowding,
                c("subject", "crowding_phenotype", "Age_years", "Educational_years",
                  "Total_Intracranial_Volume_eTIV_mm3", "Gender", "Manual_Preference"),
                "crowding file")

crowding <- crowding |>
  mutate(subject = as.character(subject),
         crowding_phenotype = factor(crowding_phenotype, levels = c("Split", "Mixed", "Crowded")),
         Gender = factor(Gender, levels = c("F", "M")),
         Manual_Preference = factor(Manual_Preference, levels = c("LH", "RH"))) |>
  filter(complete.cases(crowding_phenotype, Age_years, Educational_years,
                        Total_Intracranial_Volume_eTIV_mm3, Gender, Manual_Preference)) |>
  droplevels()

sc <- readRDS(file.path(path_folder, "Data", 
                        "HCP_981_participants_400YanROIs_FiberDensity.Rds"))

crowding <- crowding[match(dimnames(sc)[[3]], crowding$subject), , drop = FALSE]

n_roi <- dim(sc)[1]
upper <- which(upper.tri(matrix(FALSE, n_roi, n_roi)), arr.ind = TRUE)
all_zero <- vapply(seq_len(nrow(upper)), function(i) {
  all(sc[upper[i, 1], upper[i, 2], ] == 0)}, 
  logical(1))
upper <- upper[!all_zero, , drop = FALSE]
cat("Nonzero edges tested:", nrow(upper), "\n")

base_formula <- paste("crowding_phenotype * Manual_Preference +",
                      "Age_years + Educational_years +",
                      "Total_Intracranial_Volume_eTIV_mm3 + Gender")
p_values <- numeric(nrow(upper))
pb <- progress_bar$new(format = "Edgewise models [:bar] :percent in :elapsed",
                       total = nrow(upper), clear = FALSE, width = 60)

for (i in seq_len(nrow(upper))) {
  y <- sc[upper[i, 1], upper[i, 2], ]
  edge_data <- crowding
  edge_data$HMOA <- y
  model <- lm(as.formula(paste("HMOA ~", base_formula)), data = edge_data)
  type3 <- Anova(model, type = "III")
  p_values[i] <- type3["crowding_phenotype", "Pr(>F)"]
  pb$tick()
}

edge_results <- data.frame(Index1 = upper[, 1],
                           Index2 = upper[, 2],
                           p = p_values,
                           p_FDR = p.adjust(p_values, method = "BH"))
significant_edges <- edge_results |> filter(p_FDR < 0.05)
cat("FDR-significant edges:", nrow(significant_edges), "\n")
if (!nrow(significant_edges)) stop("No FDR-significant structural-connectivity edges were found.")

atlas <- read.csv(file.path(path_folder, "Atlas", 
                            "400Parcels_Yeo2011_7Networks.txt"))
require_columns(atlas, c("Index", "Network"), "Yan atlas description")

significant_edges <- significant_edges |>
  mutate(Network1 = atlas$Network[match(Index1, atlas$Index)],
         Network2 = atlas$Network[match(Index2, atlas$Index)],
         Network_Pair = mapply(function(a, b) paste(sort(c(a, b)), collapse = "_"),
                               Network1, Network2,
                               USE.NAMES = FALSE))

# Network-pair summaries characterize the selected significant edges............
################################################################################
network_pairs <- unique(significant_edges$Network_Pair)
emmeans_list <- vector("list", length(network_pairs))
pairwise_list <- vector("list", length(network_pairs))
names(emmeans_list) <- names(pairwise_list) <- network_pairs

for (pair in network_pairs) {
  selected <- significant_edges |> filter(Network_Pair == pair)
  edge_values <- vapply(seq_len(nrow(selected)),
                        function(i) sc[selected$Index1[i], selected$Index2[i], ],
                        numeric(dim(sc)[3]))
  mean_hmoa <- if (is.null(dim(edge_values))) edge_values else rowMeans(edge_values)
  
  pair_data <- crowding
  pair_data$mean_HMOA <- mean_hmoa
  model <- lm(as.formula(paste("mean_HMOA ~", base_formula)), data = pair_data)
  emm <- emmeans(model, ~ crowding_phenotype)
  emmeans_list[[pair]] <- emm
  pairwise_list[[pair]] <- pairs(emm, adjust = "tukey")
}

network_emmeans <- bind_rows(lapply(network_pairs, function(pair) {
  as.data.frame(emmeans_list[[pair]]) |> mutate(Network_Pair = pair)
}))
network_pairwise <- bind_rows(lapply(network_pairs, function(pair) {
  as.data.frame(pairwise_list[[pair]]) |> mutate(Network_Pair = pair)
}))
print(network_emmeans)
print(network_pairwise)

# Normalize significant-edge counts by the number of anatomically possible edges
################################################################################
network_sizes <- atlas |> count(Network, name = "n_regions")
network_density <- significant_edges |>
  transmute(Network_A = mapply(function(a, b) sort(c(a, b))[1], Network1, Network2),
            Network_B = mapply(function(a, b) sort(c(a, b))[2], Network1, Network2)) |>
  count(Network_A, Network_B, name = "n_significant") |>
  mutate(Network_Pair = paste(Network_A, Network_B, sep = "_")) |>
  left_join(network_sizes, by = c("Network_A" = "Network")) |>
  rename(nA = n_regions) |>
  left_join(network_sizes, by = c("Network_B" = "Network")) |>
  rename(nB = n_regions) |>
  mutate(possible_edges = if_else(Network_A == Network_B, nA * (nA - 1) / 2, nA * nB),
         normalized_density = n_significant / possible_edges)
print(network_density)

network_involvement <- bind_rows(network_density |> transmute(Network = Network_A, normalized_density),
                                 network_density |> filter(Network_B != Network_A) |> transmute(Network = Network_B, normalized_density)) |>
  group_by(Network) |>
  summarise(total_normalized_density = sum(normalized_density), .groups = "drop") |>
  arrange(desc(total_normalized_density))
print(network_involvement)