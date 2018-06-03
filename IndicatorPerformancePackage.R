#From Saskia Otto, indicator analysis tool
devtools::install_github("saskiaotto/INDperform")

library(INDperform)

# load the CC data
load("CC.CW.redv2.0.RData")

# Using the demo data
head(ind_ex)
head(press_ex)
head(press_type_ex)
# Scoring template:
crit_scores_tmpl


# Trend modelling -------------

m_trend <- model_trend(ind_tbl = ind_ex[ ,-1],
                       time = ind_ex$Year)
# Model diagnostics
pd <- plot_diagnostics(model_list = m_trend$model)
pd$all_plots[[1]] # first indicator
# Inspect trends
pt <- plot_trend(m_trend)
pt$TZA # shows trend of TZA indicator


# Indicator response modelling ------------

### Initialize data (combining IND with pressures)
dat_init <- ind_init(ind_tbl = ind_ex[ ,-1],
                     press_tbl = press_ex[ ,-1], time = ind_ex$Year)

### Model responses
m_gam <- model_gam(init_tbl = dat_init)

# Model diagnostics (e.g. first model)
plot_diagnostics(model_list = m_gam$model[[1]])$all_plots[[1]]
# Any outlier? 
m_gam$pres_outlier %>% purrr::compact(.)
# - get number of models with outliers detected
purrr::map_lgl(m_gam$pres_outlier, ~!is.null(.)) %>% sum() 
# - which models and what observations?
m_gam %>%
  dplyr::select(id, ind, press, pres_outlier) %>%
  dplyr::filter(!purrr::map_lgl(m_gam$pres_outlier, .f = is.null)) %>%
  tidyr::unnest(pres_outlier)
# Exclude outlier in models
m_gam <- model_gam(init_tbl = dat_init, excl_outlier = m_gam$pres_outlier)
# Any temporal autocorrelation
sum(m_gam$tac)
# - which models 
m_gam %>%
  dplyr::select(id, ind, press, tac) %>%
  dplyr::filter(tac)

# If temporal autocorrelation present
m_gamm <- model_gamm(init_tbl = dat_init,
                     filter = m_gam$tac)
# Again, any outlier?
purrr::map_lgl(m_gamm$pres_outlier, ~!is.null(.)) %>% sum() 

# Select best GAMM from different correlation structures
# (based on AIC)
best_gamm <- select_model(gam_tbl = m_gam,
                          gamm_tbl = m_gamm)
plot_diagnostics(model_list = best_gamm$model[[1]])$all_plots[[1]]
# Merge GAM and GAMMs
m_merged <- merge_models(m_gam[m_gam$tac == FALSE, ], best_gamm)

# Calculate derivatives
m_calc <- calc_deriv(init_tbl = dat_init,
                     mod_tbl = m_merged)

# Test for pressure interactions
it <- select_interaction(mod_tbl = m_calc)
# (creates combinations to test for)
m_all <- test_interaction(init_tbl = dat_init, mod_tbl = m_calc,
                          interactions = it)


# Scoring based on model output ------------
scores <- scoring(trend_tbl = m_trend, mod_tbl = m_all, press_type = press_type_ex)
# Runs a shiny app to modify the score for the subcriterion 10.1:
#scores <- expect_resp(mod_tbl = m_all, scores_tbl = scores)
sum_sc <- summary_sc(scores)
spie <- plot_spiechart(sum_sc)
spie$TZA # shows the spiechart of the indicator TZA