#From Saskia Otto, indicator analysis tool

# Saskia: Install the latest development version
# as I fixed some things after analysing your data
# devtools::install_github("saskiaotto/INDperform")

library(INDperform)
library(tidyverse)
library(gridExtra)

# Saskia: To save output
dir.create("./INDperform_output")

# load the CC data
load("CC.CW.redv2.0.RData")

# Saskia: colnomes should have no hyphens! -----
# This is why you had error messages when modelling trends or pressure relationships
# with the new INDperform version >=0.2.2 this shouldn't be a problem anymore
# check_names <- names(dat.red1)
# check_names <- str_replace(check_names, pattern = "-", replacement = "_")
# names(dat.red1)<- str_replace(check_names, pattern = "__", replacement = "_")
# --------------

data.full<-dat.red1
str(data.full)
#View(data.full)

#remove NAs
dat.red1<-na.omit(dat.red1)
dim(dat.red1)

envind<-c(12:15,17,18)
humind<-c(3,4,5,9,16)
ecosysind<-c(2,6,7,8,10,11,19)
cantcontrolind<-c(3,13,17,18)

envind<-c(grep("NOI",names(dat.red1)),grep("NPGO",names(dat.red1)),grep("PDO",names(dat.red1)))
ecosysind<-c(grep("lion",names(dat.red1)),grep("GF",names(dat.red1)),grep("Scav",names(dat.red1)),grep("Cop",names(dat.red1)))
humind<-c(grep("Commercial",names(dat.red1)),grep("Coastal",names(dat.red1)),grep("pollution",names(dat.red1)),grep("Nutrient",names(dat.red1)),grep("landings",names(dat.red1)),grep("Dredging",names(dat.red1)),grep("Habitat",names(dat.red1)))

ind_ex <- dat.red1[,ecosysind]
press_ex <- dat.red1[,c(envind,humind)]
names(press_ex)

press_type_ex <- data.frame(
  press = names(press_ex),
  press_type = NA
  )
press_type_ex$press_type[press_type_ex$press %in% names(dat.red1)[envind]] <- "Climate"
press_type_ex$press_type[press_type_ex$press %in% names(dat.red1)[humind]] <- "Human"
press_type_ex

# Using the demo data
head(ind_ex)
head(press_ex)
head(press_type_ex)
# Scoring template:
crit_scores_tmpl




# Trend modelling -------------

m_trend <- model_trend(ind_tbl = ind_ex,
                       time = dat.red1$year)
# Model diagnostics
pd <- plot_diagnostics(model_list = m_trend$model)
pd$all_plots[[1]] # first indicator

# check for outliers in all models
grid.arrange(grobs = pd$cooks_dist, ncol = 3) # 4 models have outliers
# check normality in all models
grid.arrange(grobs = pd$qq_plot, ncol = 3)
# check homogeneity in all models
grid.arrange(grobs = pd$resid_plot, ncol = 3)
# check for autocorrelation in all models
grid.arrange(grobs = pd$acf_plot, ncol = 3)
# check for partial autocorrelation in all models
grid.arrange(grobs = pd$pacf_plot, ncol = 3)

# Save diagnostic plots per indicator in pdf
ml <- marrangeGrob(grobs = pd$all_plots, ncol = 1, nrow = 1)
ggsave("INDperform_output/Trend_diagnostics.pdf", ml, height = 8, width = 12)


# Inspect trends
pt <- plot_trend(m_trend)
pt[[1]] # shows first indicator
# show all together:
ml <- grid.arrange(grobs = pt, ncol = 3)
ml
# save as pdf
ggsave("INDperform_output/Trend_plots.pdf", ml, height = 12, width = 15)

# show only significant trends
grid.arrange(
  grobs = pt[which(m_trend$p_val <= 0.05)],
  ncol = 2)


# Indicator response modelling ------------

### 1. Initialize data (combining IND with pressures)
# Saskia: check if you want to combine all pressures with all indicators?
dat_init <- ind_init(ind_tbl = ind_ex,
                     press_tbl = press_ex,
                     time = dat.red1$year,
                     train = 0.8) # Saskia: since we have only 8 observations! --> 6 training obs, 2 test obs.
                      # or choose train = 1 and then remove crit in scoring

### 2a. Model responses with GAMs
m_gam <- model_gam(init_tbl = dat_init)

# Model diagnostics  --> Saskia: many models to look at!!!!
pd <- plot_diagnostics(model_list = m_gam$model)
pd$all_plots[[1]] # all diagnostics of first indicator
# Save diagnostic plots per indicator in pdf
ml <- marrangeGrob(grobs = pd$all_plots, ncol = 1, nrow = 1)
ggsave("INDperform_output/GAM_diagnostics.pdf", ml, height = 8, width = 12)

# Inspect diagnostics already saved in output tibble:
# Any outlier?
m_gam$pres_outlier %>% purrr::compact(.) # Saskia: yes but then you have only few data points!!
# - get number of models with outliers detected
purrr::map_lgl(m_gam$pres_outlier, ~!is.null(.)) %>% sum()
# - which models and what observations?
m_gam %>%
  dplyr::select(id, ind, press, pres_outlier) %>%
  dplyr::filter(!purrr::map_lgl(m_gam$pres_outlier, .f = is.null)) %>%
  tidyr::unnest(pres_outlier) %>%
  print(n = 70)
# Exclude outlier in models  --> Saskia: I wouldn't advise to do this
#m_gam <- model_gam(init_tbl = dat_init, excl_outlier = m_gam$pres_outlier)
# Any temporal autocorrelation
sum(m_gam$tac)
# - which models
m_gam %>%
  dplyr::select(id, ind, press, tac) %>%
  dplyr::filter(tac) %>%
  print(n = 70)

# 2b. Aply GAMMs to account for TAC
# If temporal autocorrelation present --> Saskia: not all cor structures work due to limited sample size
m_gamm <- model_gamm(init_tbl = dat_init,
                     filter = m_gam$tac)
# Again, any outlier?  --> again, I wouldn't like to remove data
purrr::map_lgl(m_gamm$pres_outlier, ~!is.null(.)) %>% sum()

# Select best GAMM from different correlation structures
# (based on AIC) --> Saskia: or choose manually one corr structure per indicator
best_gamm <- select_model(gam_tbl = m_gam,
                          gamm_tbl = m_gamm)
View(best_gamm)

# Still any temporal autocorrelation?
sum(best_gamm$tac) # --> Saskia: still some models show TAC
# - which models
best_gamm %>%
  dplyr::select(id, ind, press, tac, corrstruc) %>%
  dplyr::filter(tac) %>%
  print(n = 70)

# GAM(M) diagnostics  --> Saskia: many models to look at!!!!
pd <- plot_diagnostics(model_list = best_gamm$model)
# Save diagnostic plots per indicator in pdf
ml <- marrangeGrob(grobs = pd$all_plots, ncol = 1, nrow = 1)
ggsave("INDperform_output/GAM(M)_diagnostics.pdf", ml, height = 8, width = 12)


# 2c. Merge GAM and GAMMs
m_merged <- merge_models(m_gam[m_gam$tac == FALSE, ], best_gamm)



# 3. Calculate derivatives of GAM/GAMMs tibble ---
m_deriv <- calc_deriv(init_tbl = dat_init,
                     mod_tbl = m_merged)
# View sign. models
filter(m_deriv, p_val <= 0.05) %>%
  View() # --> 30 sign. models
# Show final significant GAM/GAMM model plots:
sel <- which(m_deriv$p_val <= 0.05)
pm <- plot_model(init_tbl = dat_init[sel, ], mod_tbl = m_deriv[sel, ])
# Save all sign. indicator plots
ml <- gridExtra::marrangeGrob(grobs = pm$all_plots, ncol = 1, nrow = 1)
ggplot2::ggsave("INDperform_output/Sign_GAM(M)_results.pdf", ml, height = 10, width = 12)


# 3. Calculate derivatives of GAMs only ---
m_gam_deriv <- calc_deriv(init_tbl = dat_init,
                     mod_tbl = m_gam)
# View sign. GAMs
filter(m_gam_deriv, p_val <= 0.05) %>%
  View() # --> 7 sign. models
# Show final significant GAM-only model plots:
sel <- which(m_gam_deriv$p_val <= 0.05)
pm <- plot_model(init_tbl = dat_init[sel, ], mod_tbl = m_gam_deriv[sel, ])
# save all sign. indicator plots
ml <- gridExtra::marrangeGrob(grobs = pm$all_plots, ncol = 1, nrow = 1)
ggplot2::ggsave("INDperform_output/Sign_GAM_results.pdf", ml, height = 10, width = 12)



# 4. Test for pressure interactions with threshold GAMs
# Saskia: Here, check what you really want to combine...

# the following creates all combinations to test for:
# it <- select_interaction(mod_tbl = m_deriv)
# --> with this you would test for 1890 combinations!!!
# maybe filter some out from the 'it' tibble?

# m_all <- test_interaction(init_tbl = dat_init, mod_tbl = m_deriv,
#                           interactions = it)


# Saskia: I will continue with the scoring without the interactions for now
# (hence we will not score crit. 10_4 pressure interactions)
my_crit_scores <- filter(crit_scores_tmpl, subcrit != "C10_4")


# Scoring based on model output ------------
m_fin <- m_gam_deriv
scores <- scoring(trend_tbl = m_trend, mod_tbl = m_fin, press_type = press_type_ex,
  crit_scores = my_crit_scores)
# Runs a shiny app to modify the score for the subcriterion 10.1:
# Saskia: You have to go through every response and see if it as expected:
scores <- expect_resp(mod_tbl = m_fin, scores_tbl = scores)
sum_sc <- summary_sc(scores)
spie <- plot_spiechart(sum_sc, title_size = 4, lab_size = 3)
spie[[1]] # shows the spiechart of the first indicator

# Show all spiecharts:
ml <- grid.arrange(grobs = spie, ncol = 3)
ggsave("INDperform_output/Spiecharts.pdf", ml, height = 10, width = 15)


# Score-based cluster analysis -----

# Unweighted version:
scores_mat <- summary_sc(scores)$scores_matrix
dist_matrix <- dist_sc(scores_mat)
clust_scores <- clust_sc(dist_matrix)
plot_clust_sc(clust_scores)


# Weighted by the different criteria and pressure types: (see also ?dist_sc_group)
# Split the scores by pressure-independent criteria and pressure types
mat_nam <- names(scores_mat)
env_nam <- names(dat.red1)[envind]
hum_nam <- names(dat.red1)[humind]

col_clim <- purrr::map_lgl(1:length(mat_nam),
  ~ any(str_detect(mat_nam[.], env_nam))) %>%
  which()
col_hum <- purrr::map_lgl(1:length(mat_nam),
  ~ any(str_detect(mat_nam[.], hum_nam))) %>%
  which()

dist_w_matrix <- dist_sc_group(
  x = list(
    scores_mat[,1:2], # crit 8 and 11
    scores_mat[ ,col_clim],  # crit 9 and 10 - climate
    scores_mat[ ,col_hum] # crit 9 and 10 - humans
  )
)

clust_w_scores <- clust_sc(dist_w_matrix)
plot_clust_sc(clust_w_scores)

ggsave("INDperform_output/Weighted_cluster_analysis.pdf")

save(m_trend, dat_init, m_gam, m_gamm, best_gamm, m_merged, m_deriv, m_gam_deriv,
  file = "INDperform_output/Saskia_all_model4checking.Rdata")

# load("INDperform_output/Saskia_all_model4checking.Rdata")
