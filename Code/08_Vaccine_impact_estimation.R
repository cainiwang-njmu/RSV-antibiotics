# 08 Estimating the effect of vaccine intervention
## 08-1 Estimate mono-infection proportion and the prevalence of antibiotics use within mono-infection by WHO region
file_path <- "Data/data_analysis.xlsx"
data_5 <- read_excel(file_path, sheet = "mono-infection")

# Filter valid RSV cases
data5 <- data_5 %>%
  filter(RSV != 0 & !is.na(RSV))

filtered_data5 <- data5 %>%
  mutate(ID = as.character(ID)) %>%
  left_join(data_1, by = "ID") #data_1: 1. Data Cataloguing

# Recode LMICs to UMICs (limited LMIC data; combined for stable estimation)
filtered_data5$WHO_region <- gsub("LMICs", "UMICs", filtered_data5$WHO_region)

# Meta-analysis: mono-infection proportion by region
study_labels5 <- paste(filtered_data5$"Author & year")
meta_region5 <- rma.glmm(
  measure = "PLO",
  data = filtered_data5,
  xi = `RSV-mono`,
  ni = RSV,
  mods = ~ WHO_region - 1
)
meta_region5$slab5 <- study_labels5
summary(meta_region5)

res.df5 <- data.frame(
  est = meta_region5$yi,
  se = sqrt(meta_region5$vi),
  n.all = nrow(filtered_data5)
)
res.df5$prop.est <- transf.ilogit(res.df5$est)
res.df5$prop.lci <- transf.ilogit(res.df5$est - 1.96 * res.df5$se)
res.df5$prop.uci <- transf.ilogit(res.df5$est + 1.96 * res.df5$se)
res.df5$slab <- meta_region5$slab5

res.df_combine5 <- data.frame(
  est = as.numeric(meta_region5$b),
  se = meta_region5$se,
  n.all = nrow(filtered_data5)
)
res.df_combine5$prop.est <- transf.ilogit(res.df_combine5$est)
res.df_combine5$prop.lci <- transf.ilogit(res.df_combine5$est - 1.96 * res.df_combine5$se)
res.df_combine5$prop.uci <- transf.ilogit(res.df_combine5$est + 1.96 * res.df_combine5$se)

 -----------------------------------------------------------
# Meta-analysis: antibiotics use in mono-infection (HICs and UMICs) 
hic_umic_dat <- subset(filtered_mono, WHO_region %in% c("HICs", "UMICs"))
study_labels6 <- paste(hic_umic_dat$"Author & year")

meta_region6 <- rma.glmm(
  measure = "PLO",
  data = hic_umic_dat,
  xi = Use,
  ni = Case,
  mods = ~ WHO_region - 1 # Stratified by region; remove for overall estimate
)
meta_region6$slab6 <- study_labels6
summary(meta_region6)
table(hic_umic_dat$WHO_region)

res.df6 <- data.frame(
  est = meta_region6$yi,
  se = sqrt(meta_region6$vi),
  n.all = nrow(hic_umic_dat)
)
res.df6$prop.est <- transf.ilogit(res.df6$est)
res.df6$prop.lci <- transf.ilogit(res.df6$est - 1.96 * res.df6$se)
res.df6$prop.uci <- transf.ilogit(res.df6$est + 1.96 * res.df6$se)
res.df6$slab <- meta_region6$slab6

res.df_combine6 <- data.frame(
  est = as.numeric(meta_region6$b),
  se = meta_region6$se,
  n.all = nrow(hic_umic_dat)
)
res.df_combine6$prop.est <- transf.ilogit(res.df_combine6$est)
res.df_combine6$prop.lci <- transf.ilogit(res.df_combine6$est - 1.96 * res.df_combine6$se)
res.df_combine6$prop.uci <- transf.ilogit(res.df_combine6$est + 1.96 * res.df_combine6$se)

 -----------------------------------------------------------
# Monte Carlo simulation from logit-scale estimates
#' @param logit_est Logit-scale point estimate
#' @param logit_se Logit-scale standard error
#' @param param_name Parameter name for labeling
#' @param n_sim Number of simulations (default: 1000)
#' @param seed Random seed for reproducibility
#' @return List with summary and full simulation draws
logit_prop_mc <- function(logit_est, logit_se, param_name, n_sim = 1000, seed = 123) {
  set.seed(seed)
  
  logit_sim <- rnorm(n_sim, mean = logit_est, sd = logit_se)
  
  prop_sim <- plogis(logit_sim)
  
  summary_df <- tibble(
    parameter = param_name,
    point_estimate = plogis(logit_est),
    sim_median = median(prop_sim),
    sim_mean = mean(prop_sim),
    sim_lower = quantile(prop_sim, 0.025),
    sim_upper = quantile(prop_sim, 0.975)
  )
  
  draws_df <- tibble(
    parameter = param_name,
    draw_value = prop_sim,
    draw_id = 1:n_sim,
    logit_value = logit_sim
  )
  
  return(list(
    summary = summary_df,
    full_simulations = setNames(list(prop_sim), param_name),
    model_info = list(n_simulations = n_sim, distribution = "logit-normal")
  ))
}

 -----------------------------------------------------------
# Define parameters for simulation
# Estimates from meta-analysis models:
# mono_HIC: mono-infection proportion in HICs
# mono_LMICs: mono-infection proportion in LMICs/UMICs
# HICs: antibiotics use in HICs (mono-infection)
# UMICs: antibiotics use in UMICs (mono-infection)
# Global: global antibiotics use in mono-infection
params_in <- list(
  logit_est = c(2.0945,      
                0.1006,     
                -0.5033,     
                1.6808,     
                1.4166),     
  logit_se  = c(0.4572,      
                0.1897,      
                2.6788,      
                1.0867,     
                1.0481),     
  param_name = c("mono_HIC", "mono_LMICs", "HICs", "LMICs", "Global")
)

# Run Monte Carlo simulations
mc_prop <- pmap(params_in, logit_prop_mc, n_sim = 1000, seed = 123)

logit_draws_df <- map_dfr(mc_prop, function(x) {
  tibble(
    parameter = x$summary$parameter,
    draw_value = x$full_simulations[[1]],
    draw_id = 1:length(x$full_simulations[[1]])
  )
})

summary_prop <- map_dfr(mc_prop, ~ .x$summary)
print(summary_prop)

split_prop <- split(logit_draws_df, logit_draws_df$parameter)

HICs_use    <- as.data.frame(split_prop$HICs)
LMICs_use   <- as.data.frame(split_prop$LMICs)
Global_use  <- as.data.frame(split_prop$Global)
mono_HIC    <- as.data.frame(split_prop$mono_HIC)
mono_LMICs  <- as.data.frame(split_prop$mono_LMICs)

combined_mono <- cbind(
  HICs_use,
  LMICs_use,
  Global_use,
  mono_HIC,
  mono_LMICs
)

names(combined_mono) <- c(
  "HICs", "HICs_use", "draw_id1",
  "LMICs", "LMICs_use", "draw_id2",
  "Global", "Global_use", "draw_id3",
  "mono_HIC", "mono_HIC_rate", "draw_id4",
  "mono_LMICs", "mono_LMICs_rate", "draw_id5"
)

  -----------------------------------------------------------
# Vaccine effectiveness (VE) of RSV vaccine against infection and antibiotics use avoidance
#' Calculate IRR and 95% CI from vaccine and placebo data
#'
#' @param vaccine_cases Number of cases in vaccine group
#' @param vaccine_py Person-years in vaccine group
#' @param placebo_cases Number of cases in placebo group
#' @param placebo_py Person-years in placebo group
#' @param alpha Significance level (default: 0.05 for 95% CI)
#' @return List containing IRR, 95% CI, log IRR, and standard error
irr_ci_simple <- function(vaccine_cases, vaccine_py, 
                          placebo_cases, placebo_py, 
                          alpha = 0.05) {
  
  IRR <- (vaccine_cases / vaccine_py) / (placebo_cases / placebo_py)
  
  logIRR <- log(IRR)
  se_log <- sqrt(1/vaccine_cases + 1/placebo_cases)
  
  z <- qnorm(1 - alpha/2)
  logCI <- logIRR + c(-z, z) * se_log
  irrCI <- exp(logCI)
  
  return(list(
    IRR = IRR,
    IRR_95CI = setNames(irrCI, c("lower", "upper")),
    log_IRR = logIRR,
    se_log_IRR = se_log,
    parameters = list(
      vaccine_cases = vaccine_cases,
      vaccine_py = vaccine_py,
      placebo_cases = placebo_cases,
      placebo_py = placebo_py,
      alpha = alpha
    )
  ))
}

 -----------------------------------------------------------
# Define scenarios for VE calculation 
# Scenario definitions:
# - Global_antiVE: Global VE against antibiotics use
# - HICs_antiVE: HICs VE against antibiotics use
# - LMICs_antiVE: LMICs VE against antibiotics use
# - VE/VE_fnano: VE against infection

scenarios <- list(
  list(name = "Global_antiVE", 
       vac_cases = 518, vac_py = 730, 
       pla_cases = 311, pla_py = 379),
  
  list(name = "HICs_antiVE", 
       vac_cases = 25, vac_py = 242, 
       pla_cases = 27, pla_py = 132),
  
  list(name = "LMICs_antiVE", 
       vac_cases = 493, vac_py = 488, 
       pla_cases = 284, pla_py = 247),
  
  list(name = "VE", 
       vac_cases = 57, vac_py = 3495, 
       pla_cases = 117, pla_py = 3480),
  
  list(name = "VE_fnano", 
       vac_cases = 57, vac_py = 2765, 
       pla_cases = 53, pla_py = 1430)
  
)

# Calculate IRR for each scenario 
results <- list()
for (scen in scenarios) {
  res <- irr_ci_simple(
    vaccine_cases = scen$vac_cases, 
    vaccine_py = scen$vac_py,
    placebo_cases = scen$pla_cases, 
    placebo_py = scen$pla_py
  )
  results[[scen$name]] <- res
  cat(sprintf("%s: IRR=%.3f [%.3f, %.3f]\n", 
              scen$name, res$IRR, res$IRR_95CI[1], res$IRR_95CI[2]))
}

 -----------------------------------------------------------
# Monte Carlo simulation for VE
irr_log_normal_mc <- function(irr_point, irr_lower, irr_upper, 
                              param_name, n_sim = 1000, seed = 123) {
  set.seed(seed)

  mu <- log(irr_point)
  sigma <- (log(irr_upper) - log(irr_lower)) / (2 * 1.96)
  
  # Sample from normal distribution on log scale
  log_irr_sim <- rnorm(n_sim, mean = mu, sd = sigma)
  
  irr_sim <- exp(log_irr_sim)
  
  # Calculate VE (proportion scale)
  ve_sim <- (1 - irr_sim)

  summary_df <- tibble(
    parameter = param_name,
    IRR_point = irr_point,
    IRR_sim_median = median(irr_sim),
    IRR_sim_mean = mean(irr_sim),
    IRR_sim_lower = quantile(irr_sim, 0.025),
    IRR_sim_upper = quantile(irr_sim, 0.975),
    VE_point = 1 - irr_point,
    VE_sim_median = median(ve_sim),
    VE_sim_mean = mean(ve_sim),
    VE_sim_lower = quantile(ve_sim, 0.025),
    VE_sim_upper = quantile(ve_sim, 0.975)
  )

  draws_df <- tibble(
    parameter = param_name,
    IRR_draw = irr_sim,
    VE_draw = ve_sim,
    draw_id = 1:n_sim
  )
  
  return(list(
    summary = summary_df,
    full_simulations = setNames(list(irr_sim), param_name),
    ve_simulations = setNames(list(ve_sim), paste0(param_name, "_VE")),
    draws_df = draws_df,
    model_info = list(
      n_simulations = n_sim, 
      distribution = "log-normal",
      log_mu = mu,
      log_sigma = sigma
    )
  ))
}

# Define parameters for Monte Carlo simulation
params_in <- list(
  irr_point = c(0.865, 0.505, 0.879, 0.485, 0.556),
  irr_lower = c(0.751, 0.293, 0.759, 0.353, 0.383),
  irr_upper = c(0.995, 0.870, 1.017, 0.666, 0.808),
  param_name = c("Global_antiVE", "HICs_antiVE", "LMICs_antiVE", "VE", "VE_fnano")
)

# Run Monte Carlo simulation

mc_VE <- pmap(params_in, irr_log_normal_mc, seed = 123)

irr_draws_df <- map_dfr(mc_VE, ~ .x$draws_df)
summary_VE <- map_dfr(mc_VE, ~ .x$summary)
print(summary_VE)

# Split draws by parameter 

split_VE <- split(irr_draws_df, irr_draws_df$parameter)
table(irr_draws_df$parameter)

# Extract individual parameter draws
Global_antiVE <- as.data.frame(split_VE$Global_antiVE)
HICs_antiVE <- as.data.frame(split_VE$HICs_antiVE)
LMICs_antiVE <- as.data.frame(split_VE$LMICs_antiVE)
VE <- as.data.frame(split_VE$VE)

# Combine into a single data frame 
combined_VE <- cbind(HICs_antiVE, LMICs_antiVE, Global_antiVE, VE)
names(combined_VE) <- c(
  "HICs", "HICs_IRR", "HICs_antiVE", "draw_id1",
  "LMICs", "LMICs_IRR", "LMICs_antiVE", "draw_id2",
  "Global", "Global_IRR", "Global_antiVE", "draw_id3",
  "VE", "VE_IRR", "VE_inf", "draw_id4"
)


----------------------------------------------------------------------------------------------------------------
#08-2 Hospital admission estimates for infants aged 0-6 months
# Monte Carlo simulation for hospital admissions (0-6 months) 

rsv_hospital_mc <- function(scenarios, R = 1000, seed = 123) {
  set.seed(seed)
  
  results_list <- list()
  full_simulations <- list()
  
  for (i in seq_len(nrow(scenarios))) {
    cur <- scenarios[i, ]
  
    mu_log <- log(cur$hospital_admission)
    sigma_log <- (log(cur$UCL) - log(cur$LCL)) / 3.92
    
    hospital_sim <- rlnorm(R, meanlog = mu_log, sdlog = sigma_log)
 
    results_list[[i]] <- data.frame(
      income_group = cur$income_level,
      hosp_median = round(median(hospital_sim), 0),
      hosp_lower = round(quantile(hospital_sim, 0.025), 0),
      hosp_upper = round(quantile(hospital_sim, 0.975), 0),
      scenario_id = i,
      stringsAsFactors = FALSE
    )
    
    full_simulations[[i]] <- hospital_sim
  }
  
  summary_table <- do.call(rbind, results_list) %>%
    `row.names<-`(NULL)
  
  names(full_simulations) <- scenarios$income_level
  
  return(list(
    summary = summary_table,
    full_simulations = full_simulations,
    model_info = list(n_simulations = R, distribution = "lognormal")
  ))
}

hospital_scenarios <- data.frame(
  income_level = c("HICs", "LMICs", "Global"),
  hospital_admission = c(194000, 1188000, 1376000),
  LCL = c(133000, 802000, 1017000),
  UCL = c(283000, 1759000, 1982000)
) #data source: Li et al., 2022

mc_results_hosp <- rsv_hospital_mc(hospital_scenarios, R = 1000)

print(mc_results_hosp$summary)

hospital_draws_df <- map2_dfr(
  mc_results_hosp$full_simulations,
  names(mc_results_hosp$full_simulations),
  ~ tibble(
    income_group = .y,
    hosp_draw = .x,
    draw_id = 1:length(.x)
  )
)

# Split by income group and combine 
split_df <- split(hospital_draws_df, hospital_draws_df$income_group)

HICs_draw  <- as.data.frame(split_df$HICs)
LMICs_draw <- as.data.frame(split_df$LMICs)
Global_draw <- as.data.frame(split_df$Global)

# Combine into single data frame 
combined_hosp <- cbind(HICs_draw, LMICs_draw, Global_draw)
names(combined_hosp) <- c(
  "HICs", "HICs_pop_draw", "draw_id_HICs",
  "LMICs", "LMICs_pop_draw", "draw_id_LMICs",
  "Global", "Global_pop_draw", "draw_id_Global"
)


--------------------------------------------------------------------------------------------------------------------------------------------
## 08-3 Estimate the impact of maternal RSV vaccination on reducing antibiotics use in infants aged 0–6 months using an indirect approach
# Inputs:
# combined_hosp: Hospital admission draws (HICs_pop_draw, LMICs_pop_draw, Global_pop_draw)
# combined_mono: Mono-infection and antibiotic use draws (mono_HIC_rate, mono_LMICs_rate, HICs_use, LMICs_use)
# combined_VE: Vaccine effectiveness draws (VE_inf)
# vaccine_coverage: Vaccine coverage scenarios (default: 20%, 40%, 60%)

calculate_avoidable_antibiotics_indirect <- function(combined_hosp, combined_mono, combined_VE,
                                                     vaccine_coverage = c(0.2, 0.4, 0.6)) {
  
  hosp_cols <- c("HICs_pop_draw", "LMICs_pop_draw", "Global_pop_draw")
  missing_hosp <- hosp_cols[!hosp_cols %in% names(combined_hosp)]
  if (length(missing_hosp) > 0) {
    stop(sprintf("combined_hosp missing columns: %s", 
                 paste(missing_hosp, collapse = ", ")))
  }
  
  mono_cols <- c("mono_HIC_rate", "mono_LMICs_rate", "HICs_use", "LMICs_use")
  missing_mono <- mono_cols[!mono_cols %in% names(combined_mono)]
  if (length(missing_mono) > 0) {
    stop(sprintf("combined_mono missing columns: %s", 
                 paste(missing_mono, collapse = ", ")))
  }
  
  ve_cols <- c("VE_inf")
  missing_ve <- ve_cols[!ve_cols %in% names(combined_VE)]
  if (length(missing_ve) > 0) {
    stop(sprintf("combined_VE missing columns: %s", 
                 paste(missing_ve, collapse = ", ")))
  }
  
  n_rows <- c(nrow(combined_hosp), nrow(combined_mono), nrow(combined_VE))
  if (length(unique(n_rows)) != 1) {
    stop(sprintf("Row counts differ: combined_hosp=%d, combined_mono=%d, combined_VE=%d", 
                 n_rows[1], n_rows[2], n_rows[3]))
  }
  
  hosp_hics <- combined_hosp[["HICs_pop_draw"]]
  hosp_lmics <- combined_hosp[["LMICs_pop_draw"]]
  hosp_global <- combined_hosp[["Global_pop_draw"]]
  
  mono_hic <- combined_mono[["mono_HIC_rate"]]
  mono_lmic <- combined_mono[["mono_LMICs_rate"]]
  
  hic_use <- combined_mono[["HICs_use"]]
  lmic_use <- combined_mono[["LMICs_use"]]
  
  ve_inf <- combined_VE[["VE_inf"]]
  
  results_list <- list()
  scenario_id <- 1
  
  for (coverage in vaccine_coverage) {
    
    # Absolute avoidable cases
    avoidable_hics <- hosp_hics * mono_hic * coverage * ve_inf * hic_use
    avoidable_lmics <- hosp_lmics * mono_lmic * coverage * ve_inf * lmic_use
    avoidable_global <- avoidable_hics + avoidable_lmics
    
    total_hospitalizations <- hosp_global
    
    # Global VE (weighted average)
    global_ve <- avoidable_global / (coverage * hosp_global)
    
    # Percentage of hospitalizations avoided
    percent_global <- coverage * global_ve
    percent_hics <- mono_hic * coverage * ve_inf * hic_use      
    percent_lmics <- mono_lmic * coverage * ve_inf * lmic_use    
    
    # Anti-VE (vaccine effectiveness against antibiotics use)
    antiVE_hics <- mono_hic * hic_use * ve_inf
    antiVE_lmics <- mono_lmic * lmic_use * ve_inf
    antiVE_global <- global_ve  
    
    results_list[[scenario_id]] <- data.frame(
      coverage = paste0(coverage * 100, "%"),
      scenario_id = scenario_id,
      
      # Global total
      global_total_median = round(median(avoidable_global) / 100) * 100,
      global_total_lower = round(quantile(avoidable_global, 0.025) / 100) * 100,
      global_total_upper = round(quantile(avoidable_global, 0.975) / 100) * 100,
      
      # Global percentage
      global_percent_median = round(median(percent_global) * 100, 2),
      global_percent_lower = round(quantile(percent_global, 0.025) * 100, 2),
      global_percent_upper = round(quantile(percent_global, 0.975) * 100, 2),
      
      # HICs total
      hics_total_median = round(median(avoidable_hics) / 100) * 100,
      hics_total_lower = round(quantile(avoidable_hics, 0.025) / 100) * 100,
      hics_total_upper = round(quantile(avoidable_hics, 0.975) / 100) * 100,
      
      # HICs percentage
      hics_percent_median = round(median(percent_hics) * 100, 2),
      hics_percent_lower = round(quantile(percent_hics, 0.025) * 100, 2),
      hics_percent_upper = round(quantile(percent_hics, 0.975) * 100, 2),
      
      # LMICs total
      lmics_total_median = round(median(avoidable_lmics) / 100) * 100,
      lmics_total_lower = round(quantile(avoidable_lmics, 0.025) / 100) * 100,
      lmics_total_upper = round(quantile(avoidable_lmics, 0.975) / 100) * 100,
      
      # LMICs percentage
      lmics_percent_median = round(median(percent_lmics) * 100, 2),
      lmics_percent_lower = round(quantile(percent_lmics, 0.025) * 100, 2),
      lmics_percent_upper = round(quantile(percent_lmics, 0.975) * 100, 2),
      
      # Anti-VE (HICs)
      antiVE_hics_median = round(median(antiVE_hics) * 100, 2),
      antiVE_hics_lower = round(quantile(antiVE_hics, 0.025) * 100, 2),
      antiVE_hics_upper = round(quantile(antiVE_hics, 0.975) * 100, 2),
      
      # Anti-VE (LMICs)
      antiVE_lmics_median = round(median(antiVE_lmics) * 100, 2),
      antiVE_lmics_lower = round(quantile(antiVE_lmics, 0.025) * 100, 2),
      antiVE_lmics_upper = round(quantile(antiVE_lmics, 0.975) * 100, 2),
      
      # Anti-VE (Global)
      antiVE_global_median = round(median(antiVE_global) * 100, 2),
      antiVE_global_lower = round(quantile(antiVE_global, 0.025) * 100, 2),
      antiVE_global_upper = round(quantile(antiVE_global, 0.975) * 100, 2),
      
      stringsAsFactors = FALSE
    )
    
    scenario_id <- scenario_id + 1
  }
  
 
  summary_results <- do.call(rbind, results_list)
  row.names(summary_results) <- NULL
  
  return(list(
    summary = summary_results, 
    assumption = "Using separate mono-infection proportions: HICs use mono_HIC, LMICs use mono_LMICs. global_ve = avoidable_global / (coverage × hosp_global) (weighted average antiVE).",
    formulas = list(
      hics_absolute = "HICs avoidable = HICs hospitalizations × mono_HIC × coverage × VE_inf × HICs_use",
      lmics_absolute = "LMICs avoidable = LMICs hospitalizations × mono_LMICs × coverage × VE_inf × LMICs_use",
      global_absolute = "Global avoidable = HICs avoidable + LMICs avoidable",
      global_percent = "Global percentage = coverage × global_ve",
      global_ve = "global_ve = (HICs_hosp×mono_HIC×VE_inf×HICs_use + LMICs_hosp×mono_LMICs×VE_inf×LMICs_use) / Global_total_hosp",
      antiVE_hics = "antiVE_hics = mono_HIC × HICs_use × VE_inf",
      antiVE_lmics = "antiVE_lmics = mono_LMICs × LMICs_use × VE_inf",
      antiVE_global = "antiVE_global = global_ve (weighted average antiVE)"
    )
  ))
}

results_indirect <- calculate_avoidable_antibiotics_indirect(
  combined_hosp = combined_hosp,
  combined_mono = combined_mono,
  combined_VE = combined_VE
)
print(results_indirect$summary)


-----------------------------------------------------------------------------------------------------------------------------------------
## 08-4 Estimate the impact of maternal RSV vaccination on reducing antibiotics use in infants aged 0–6 months using a direct approach

# Inputs:
# combined_hosp: Hospital admission draws (HICs_pop_draw, LMICs_pop_draw, Global_pop_draw)
# combined_VE: Vaccine effectiveness draws (HICs_antiVE, LMICs_antiVE)
# vaccine_coverage: Vaccine coverage scenarios (default: 20%, 40%, 60%)
# adj_factor: Adjustment factor (default: 1.35)

calculate_avoidable_antibiotics_direct <- function(combined_hosp, combined_VE,
                                                   vaccine_coverage = c(0.2, 0.4, 0.6),
                                                   adj_factor = 1.35) {
  
  hosp_cols <- c("HICs_pop_draw", "LMICs_pop_draw", "Global_pop_draw")
  missing_hosp <- hosp_cols[!hosp_cols %in% names(combined_hosp)]
  if (length(missing_hosp) > 0) {
    stop(sprintf("combined_hosp missing columns: %s", 
                 paste(missing_hosp, collapse = ", ")))
  }
  
  ve_cols <- c("HICs_antiVE", "LMICs_antiVE")
  missing_ve <- ve_cols[!ve_cols %in% names(combined_VE)]
  if (length(missing_ve) > 0) {
    stop(sprintf("combined_VE missing columns: %s", 
                 paste(missing_ve, collapse = ", ")))
  }
  
  n_rows <- c(nrow(combined_hosp), nrow(combined_VE))
  if (length(unique(n_rows)) != 1) {
    stop(sprintf("Row counts differ: combined_hosp=%d, combined_VE=%d", 
                 n_rows[1], n_rows[2]))
  }
  
  hosp_hics <- combined_hosp[["HICs_pop_draw"]]
  hosp_lmics <- combined_hosp[["LMICs_pop_draw"]]
  hosp_global <- combined_hosp[["Global_pop_draw"]]
  
  hics_ve <- combined_VE[["HICs_antiVE"]]
  lmics_ve <- combined_VE[["LMICs_antiVE"]]
  
  # Calculate avoidable antibiotic use by coverage scenario
  
  results_list <- list()
  scenario_id <- 1
  
  for (coverage in vaccine_coverage) {
    
    # Absolute avoidable cases
    avoidable_hics <- hosp_hics * coverage * hics_ve * adj_factor
    avoidable_lmics <- hosp_lmics * coverage * lmics_ve * adj_factor
    avoidable_global <- avoidable_hics + avoidable_lmics
    
    total_hospitalizations <- hosp_global  
    
    # Global VE (back-calculated from avoidable cases)
    global_ve <- avoidable_global / (coverage * total_hospitalizations * adj_factor)
    
    # Percentage of hospitalizations avoided
    percent_global <- coverage * adj_factor * global_ve
    percent_hics <- coverage * hics_ve * adj_factor
    percent_lmics <- coverage * lmics_ve * adj_factor
    
    # Anti-VE (vaccine effectiveness against antibiotics use, adjusted)
    antiVE_hics <- hics_ve * adj_factor      
    antiVE_lmics <- lmics_ve * adj_factor    
    antiVE_global <- global_ve * adj_factor  
    
    results_list[[scenario_id]] <- data.frame(
      coverage = paste0(coverage * 100, "%"),
      adj_factor = adj_factor,
      scenario_id = scenario_id,
      
      # Global total
      global_total_median = round(median(avoidable_global) / 100) * 100,
      global_total_lower = round(quantile(avoidable_global, 0.025) / 100) * 100,
      global_total_upper = round(quantile(avoidable_global, 0.975) / 100) * 100,
      
      # Global percentage
      global_percent_median = round(median(percent_global) * 100, 2),
      global_percent_lower = round(quantile(percent_global, 0.025) * 100, 2),
      global_percent_upper = round(quantile(percent_global, 0.975) * 100, 2),
      
      # Global anti-VE
      antiVE_global_median = round(median(antiVE_global) * 100, 2),
      antiVE_global_lower = round(quantile(antiVE_global, 0.025) * 100, 2),
      antiVE_global_upper = round(quantile(antiVE_global, 0.975) * 100, 2),
      
      # HICs total
      hics_total_median = round(median(avoidable_hics) / 100) * 100,
      hics_total_lower = round(quantile(avoidable_hics, 0.025) / 100) * 100,
      hics_total_upper = round(quantile(avoidable_hics, 0.975) / 100) * 100,
      
      # HICs percentage
      hics_percent_median = round(median(percent_hics) * 100, 2),
      hics_percent_lower = round(quantile(percent_hics, 0.025) * 100, 2),
      hics_percent_upper = round(quantile(percent_hics, 0.975) * 100, 2),
      
      # HICs anti-VE
      antiVE_hics_median = round(median(antiVE_hics) * 100, 2),
      antiVE_hics_lower = round(quantile(antiVE_hics, 0.025) * 100, 2),
      antiVE_hics_upper = round(quantile(antiVE_hics, 0.975) * 100, 2),
      
      # LMICs total
      lmics_total_median = round(median(avoidable_lmics) / 100) * 100,
      lmics_total_lower = round(quantile(avoidable_lmics, 0.025) / 100) * 100,
      lmics_total_upper = round(quantile(avoidable_lmics, 0.975) / 100) * 100,
      
      # LMICs percentage
      lmics_percent_median = round(median(percent_lmics) * 100, 2),
      lmics_percent_lower = round(quantile(percent_lmics, 0.025) * 100, 2),
      lmics_percent_upper = round(quantile(percent_lmics, 0.975) * 100, 2),
      
      # LMICs anti-VE
      antiVE_lmics_median = round(median(antiVE_lmics) * 100, 2),
      antiVE_lmics_lower = round(quantile(antiVE_lmics, 0.025) * 100, 2),
      antiVE_lmics_upper = round(quantile(antiVE_lmics, 0.975) * 100, 2),
      
      stringsAsFactors = FALSE
    )
    
    scenario_id <- scenario_id + 1
  }
  
  summary_results <- do.call(rbind, results_list)
  row.names(summary_results) <- NULL
  
  return(list(
    summary = summary_results,
    assumption = "Global hospitalizations are read from the 'Global_pop_draw' column in combined_hosp. Global VE is back-calculated: global_ve = avoidable_global / (coverage × total_hospitalizations × adj_factor).",
    formulas = list(
      hics = "HICs avoidable = HICs hospitalizations × coverage × HICs_antiVE × adj_factor",
      lmics = "LMICs avoidable = LMICs hospitalizations × coverage × LMICs_antiVE × adj_factor",
      global_absolute = "Global avoidable = HICs avoidable + LMICs avoidable",
      global_percent = "Global percentage reduction = coverage × global_ve × adj_factor",
      global_ve = "Global VE = Global avoidable / (coverage × Global hospitalizations × adj_factor)",
      hics_percent = "HICs percentage reduction = coverage × HICs_antiVE × adj_factor",
      lmics_percent = "LMICs percentage reduction = coverage × LMICs_antiVE × adj_factor"
    )
  ))
}

results_direct <- calculate_avoidable_antibiotics_direct(
  combined_hosp = combined_hosp,
  combined_VE = combined_VE
)
print(results_direct$summary)
