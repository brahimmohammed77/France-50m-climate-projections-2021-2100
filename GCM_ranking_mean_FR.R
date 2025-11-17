# ============================================================
# Climate Model Rankings — France Average (R script)
# Output: CSV with rankings by Variable × Scenario × Season
#
# Variable codes (observational side):
#   TM = mean temperature (°C)
#   TN = minimum temperature (°C)
#   TX = maximum temperature (°C)
#   RR = precipitation (mm)
#
# Model variable names:
#   tas     = near-surface air temperature (K)             → corresponds to TM
#   tasmin  = minimum near-surface air temperature (K) → corresponds to TN
#   tasmax  = maximum near-surface air temperature (K) → corresponds to TX
#   pr      = precipitation (mm)       → corresponds to RR
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(stringr)
  library(readr)
})

# ---------- CONFIG ----------
mod_dir <- "extractions_points"
obs_dir <- "variables_separees"

# Mapping observational variable -> model variable name + observation file
# TM/TN/TX = temperatures in °C; RR = precipitation (mm)
# tas/tasmin/tasmax = temperatures in Kelvin; pr = precipitation flux
var_map <- tibble(
  variable = c("TM",            "TN",            "TX",            "RR"),
  var_mod  = c("tas",           "tasmin",        "tasmax",        "pr"),  # model variable names
  obs_file = c("donnees_TM.csv","donnees_TN.csv","donnees_TX.csv","donnees_RR.csv")
)

# Global filters
date_start <- as.Date("2015-01-01")
date_end   <- as.Date("2024-01-01")
alt_range  <- c(0, 4000)                # station altitude range in meters
seasons_keep <- c("All")
sort_by <- "RMSE"                       # "RMSE" or "|BIAS|"

# Output file
output_csv <- file.path(getwd(), sprintf("model_rankings_france_average_%s.csv", Sys.Date()))

# ---------- UTILITY FUNCTIONS ----------
label_season <- function(d) {
  m <- month(d)
  case_when(
    m %in% c(12,1,2) ~ "Winter",
    m %in% 3:5       ~ "Spring",
    m %in% 6:8       ~ "Summer",
    m %in% 9:11      ~ "Autumn",
    TRUE             ~ NA_character_
  )
}

# Read observations + compute national average per date (after altitude filter)
make_obs_nat <- function(obs_path, alt_rng) {
  obs_full <- read_csv(obs_path, show_col_types = FALSE) %>%
    mutate(DATE = ymd(DATE), OBS = VALEUR)
  
  valid_stations <- obs_full %>%
    select(POSTE, ALTI) %>% distinct() %>%
    filter(ALTI >= alt_rng[1], ALTI <= alt_rng[2]) %>%
    pull(POSTE)
  
  list(
    obs_nat = obs_full %>%
      filter(POSTE %in% valid_stations) %>%
      group_by(DATE) %>%
      summarise(OBS = mean(OBS, na.rm = TRUE), .groups = "drop"),
    valid_stations = valid_stations
  )
}

# Build (for one variable) the long data.frame obs vs model for France average
# For temperatures (tas, tasmin, tasmax): convert Kelvin → Celsius
# For precipitation (pr): units depend on model output (commonly mm/day)
build_var_df <- function(vm, mod_dir, obs_dir, date_rng, alt_rng) {
  obs_path <- file.path(obs_dir, vm$obs_file)
  if (!file.exists(obs_path)) {
    message("⚠️ Observation file not found: ", obs_path)
    return(tibble())
  }
  
  o <- make_obs_nat(obs_path, alt_rng)
  obs_nat <- o$obs_nat
  valid_stations <- o$valid_stations
  
  files_var <- list.files(
    mod_dir,
    pattern = paste0("^", vm$var_mod, ".*\\.csv$"),
    full.names = TRUE, recursive = TRUE
  )
  if (length(files_var) == 0) {
    message("⚠️ No model file found for ", vm$var_mod)
    return(tibble())
  }
  
  tibble(file = files_var) %>%
    mutate(
      SCENARIO = str_extract(file, "ssp[0-9]+"),
      MODEL_ID = basename(file) %>%
        str_extract("_(.*?)_") %>% str_remove_all("_")
    ) %>%
    distinct() %>%
    mutate(Variable = vm$variable, VarMod = vm$var_mod) %>%
    pmap_dfr(function(file, SCENARIO, MODEL_ID, Variable, VarMod) {
      
      dfm <- read_csv(file, show_col_types = FALSE) %>%
        mutate(DATE = ymd(DATE), MODEL = VALEUR) %>%
        filter(POSTE %in% valid_stations) %>%
        group_by(DATE) %>%
        summarise(MODEL = mean(MODEL, na.rm = TRUE), .groups = "drop")
      
      # Convert temperatures: tas/tasmin/tasmax (K) → Celsius
      if (VarMod %in% c("tas","tasmin","tasmax")) {
        dfm <- dfm %>% mutate(MODEL = MODEL - 273.15)
      }
      
      left_join(obs_nat, dfm, by = "DATE") %>%
        filter(DATE >= date_rng[1], DATE <= date_rng[2]) %>%
        mutate(
          ERROR = MODEL - OBS,
          Variable = Variable,
          SCENARIO = SCENARIO,
          MODEL_ID = MODEL_ID
        )
    })
}

# ---------- MAIN PIPELINE ----------
message("⏳ Computing national-average errors (all variables / scenarios / models)...")

df_all <- var_map %>%
  split(.$variable) %>%
  map_dfr(~ build_var_df(.x, mod_dir, obs_dir, c(date_start, date_end), alt_range))

if (nrow(df_all) == 0) stop("No combined data produced. Check paths and files.")

# Stats by season
df_seas <- df_all %>% mutate(SEASON = label_season(DATE))

stats_season <- df_seas %>%
  group_by(Variable, SCENARIO, SEASON, MODEL_ID) %>%
  summarise(
    BIAS = mean(ERROR, na.rm = TRUE),
    RMSE = sqrt(mean(ERROR^2, na.rm = TRUE)),
    N    = n(),
    .groups = "drop"
  )

# Stats "All seasons"
stats_all <- df_all %>%
  group_by(Variable, SCENARIO, MODEL_ID) %>%
  summarise(
    BIAS = mean(ERROR, na.rm = TRUE),
    RMSE = sqrt(mean(ERROR^2, na.rm = TRUE)),
    N    = n(),
    .groups = "drop"
  ) %>% mutate(SEASON = "All")

stats <- bind_rows(stats_season, stats_all)

# Filter selected seasons
stats <- stats %>% filter(SEASON %in% seasons_keep)

# Ranking (by group Variable × Scenario × Season)
ranked <- stats %>%
  group_by(Variable, SCENARIO, SEASON) %>%
  arrange(if (identical(sort_by, "RMSE")) RMSE else abs(BIAS), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup() %>%
  transmute(
    Variable, SCENARIO, Season = SEASON, MODEL_ID,
    BIAS = round(BIAS, 2), RMSE = round(RMSE, 2), N, Rank
  ) %>%
  arrange(Variable, SCENARIO, Season, Rank)

# ---------- EXPORT ----------
write_csv(ranked, output_csv)
message("✅ Rankings exported: ", output_csv)
