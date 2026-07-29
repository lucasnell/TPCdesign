



library(MASS)  # must come before dplyr
library(ape)
library(tidyverse)
library(taxize)  # for finding alternative species names

# Cleans up species names
spp_cleaner <- function(species) {
    species |>
        # remove everything after first space:
        str_remove_all("\\ .*") |>
        # remove everything before and including first _:
        str_remove("^.*?_") |>
        # remove everything after and including first _:
        str_remove("_.*")
}


tpc_df <- "~/Library/CloudStorage/Box-Box/TPCdesign/phylogeny_input_briere_params.csv" |>
    read_csv(col_types ="ccccdddddc") |>
    rename(ctmin = Tmin, ctmax = Tmax) |>
    # remove genus name from species, then clean it:
    mutate(species = species |> str_remove(".*\\ ") |> spp_cleaner()) |>
    filter(!is.na(genus)) |>
    group_by(genus) |>
    mutate(n_genus_spp = n()) |>
    ungroup() |>
    select(-source)


loo <- function(x) {
    map_dbl(1:length(x), \(i) x[i] - mean(x[-i]))
}


tpc_diff_df <- tpc_df |>
    select(-b) |>
    filter(n_genus_spp >= 3) |>
    select(-n_genus_spp) |>
    group_by(genus) |>
    mutate(across(ctmin:log_b, loo)) |>
    ungroup()

tpc_diff_df |>
    (\(x) {
        par(mfrow = c(1, 3))
        hist(x$ctmin)
        hist(x$ctmax)
        hist(x$log_b)
        invisible(NULL)
    })()




distr_df <- crossing(varbl = c("ctmin", "ctmax", "log_b"),
                     distr = c("cauchy", "normal")) |>
    mutate(fit = map2(distr, varbl, \(d, v) suppressWarnings(fitdistr(tpc_diff_df[[v]], d))),
           LL = map_dbl(fit, logLik)) |>
    group_by(varbl) |>
    filter(LL == max(LL)) |>
    ungroup()

distr_df$fit |> set_names(distr_df$varbl)
# $ctmax
#     location      scale
#   -0.1048662    1.2966241
#  ( 0.2081358) ( 0.2099444)
#
# $ctmin
#        mean             sd
#   -3.460435e-17    3.325073e+00
#  ( 3.789273e-01) ( 2.679420e-01)
#
# $log_b
#        mean            sd
#   1.225571e-17   3.602843e-01
#  (4.105820e-02) (2.903253e-02)


tpc_diff_df |>
    (\(x) {
        par(mfrow = c(1, 3))
        hist(x$ctmin, breaks = 10)
        curve(19 / 0.12 * dnorm(x, distr_df$fit[[2]][["estimate"]][["mean"]],
                      distr_df$fit[[2]][["estimate"]][["sd"]]),
              col = "red", lwd = 1.5, add = TRUE)
        hist(x$ctmax, breaks = 10)
        curve(26 / 0.25 * dcauchy(x, distr_df$fit[[1]][["estimate"]][["location"]],
                      distr_df$fit[[1]][["estimate"]][["scale"]]),
              col = "red", lwd = 1.5, add = TRUE)
        hist(x$log_b, breaks = 10)
        curve(16 / 1.1 * dnorm(x, distr_df$fit[[3]][["estimate"]][["mean"]],
                      distr_df$fit[[3]][["estimate"]][["sd"]]),
              col = "red", lwd = 1.5, add = TRUE)
        invisible(NULL)
    })()


# For use in testing sims:
tpc_diff_df |>
    select(-species_for_tree) |>
    write_csv(file = "_testing/interm-data/tpc-genus-diffs.csv")




# =============================================================================*
# =============================================================================*
# Phylogeny-based estimates ----
# =============================================================================*
# =============================================================================*


# Instructions for using `ape::ace`:
# # 2. Create a named vector of trait values for the known tips
# trait_data <- c(2.1, 1.8, 3.5, 3.2, 0.9)
# names(trait_data) <- c("A", "B", "C", "D", "E")
#
# # 3. Predict/estimate ancestral states using ace()
# # For continuous traits (default model is Brownian Motion "BM")
# recon_cont <- ace(trait_data, tree, type = "continuous", method = "REML")

# Errors output from timetree when making tree from list of species:
error_df <- "~/Library/CloudStorage/Box-Box/TPCdesign/tpc-timetree-errors.txt" |>
    read_lines() |>
    str_remove_all("\\)") |>
    str_split(" \\(") |>
    map(\(x) c(str_split(x[[1]], " ")[[1]], x[[2]])) |>
    (\(x) tibble(genus = map_chr(x, \(z) z[[1]]),
                 species = map_chr(x, \(z) z[[2]]),
                 reason = map_chr(x, \(z) z[[3]])))()

repl_df <- error_df |>
    filter(grepl("replaced with", reason)) |>
    mutate(reason = str_remove_all(reason, "replaced with ")) |>
    mutate(repl_genus = str_remove_all(reason, "\\ .*"),
           repl_species = str_remove(reason, "^.*?\\ ")) |>
    select(-reason)


phy <- "~/Library/CloudStorage/Box-Box/TPCdesign/tpc-species.nwk" |>
    read.tree()


# LEFT OFF ----
# Match species names between sources

phy$tip.label[!phy$tip.label %in% tpc_df$species_for_tree]

tpc_df |> filter(species == "flos-aquae")

# Produce a new "species_for_tree" that matches `phy`
# spp_matcher <- function(genus, species, phy) {}

genus <- tpc_df$genus
species <- tpc_df$species


tip_labels <- phy$tip.label
species_for_tree <- paste(genus, species, sep = "_")

# unmatched tip labels:
un_tip_labs <- tip_labels[!tip_labels %in% species_for_tree]

un_tip_genus <- str_split(un_tip_labs, "_") |> map_chr(\(x) x[1])
un_tip_species <- str_split(un_tip_labs, "_") |> map_chr(\(x) x[2])

un_df_idx <- which(!species_for_tree %in% tip_labels)
un_df_genus <- genus[un_df_idx]
un_df_species <- species[un_df_idx]

sum(!un_tip_species %in% un_df_species)
sum(!un_tip_genus %in% un_df_genus)

un_tip_labs[!un_tip_species %in% un_df_species & !un_tip_genus %in% un_df_genus]




# un_tip_labs[!un_tip_species %in% un_df_species & !un_tip_genus %in% un_df_genus] |>
#     str_replace_all("_", " ") |>
#     synonyms(db = "itis")



taxadb::filter_name("Zoramia leptacantha", provider = "ncbi")
taxadb::filter_name("Podarcis atrata", provider = "ncbi")
taxadb::filter_name("Cosmocomoidea triguttata", provider = "ncbi")



read_csv("~/Library/CloudStorage/Box-Box/TPCdesign/tpcs.csv",
         col_types = "cccdddcddcc") |>
    # filter(!is.na(species), !is.na(genus)) |>
    mutate(species = spp_cleaner(species)) |>
    filter(grepl("flos", species))




# =============================================================================*
# =============================================================================*
# Phylogeny ----
# =============================================================================*
# =============================================================================*



# https://github.com/HuckleyLab/ThermalStress/blob/master/data/tpcs.csv
read_csv("~/Library/CloudStorage/Box-Box/TPCdesign/tpcs.csv",
         col_types = "cccdddcddcc") |>
    filter(!is.na(species), !is.na(genus)) |>
    select(genus, species) |>
    mutate(species = spp_cleaner(species)) |>
    distinct(genus, species) |>
    arrange(genus, species) |>
    pmap_chr(\(genus, species) paste(genus, species)) |>
    write_lines("~/Library/CloudStorage/Box-Box/TPCdesign/tpc-species.txt")



