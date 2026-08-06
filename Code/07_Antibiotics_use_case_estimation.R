# 07 Estimates of antibiotics use in children younger than 5 years hospitalised with RSV-associated acute lower respiratory infection by World Bank income regions
## 07-1 Predict the prevalence of antibiotics use based on meta-regression model
### Predict logit-scale effects for all scenarios
predict_effects_logit <- function(model, data) {
  
  study_levels <- levels(data$Study_duration)
  age_levels <- levels(data$Age_group)
  who_levels <- levels(data$WHO_region)
  diag_levels <- levels(data$Diagnosis)
  wards_levels <- levels(data$Wards_stand)
  
  scenarios <- expand.grid(
    WHO_region = c("HICs", "UMICs", "LMICs", "LICs"),
    Age_group = c("＜1", "＜5"),
    Study_duration = "2010-2019",
    Diagnosis = "ALRI",
    Wards_stand = "IP",
    stringsAsFactors = FALSE
  )
  
  scenarios$Study_duration <- factor(scenarios$Study_duration, levels = study_levels)
  scenarios$Age_group <- factor(scenarios$Age_group, levels = age_levels)
  scenarios$WHO_region <- factor(scenarios$WHO_region, levels = who_levels)
  scenarios$Diagnosis <- factor(scenarios$Diagnosis, levels = diag_levels)
  scenarios$Wards_stand <- factor(scenarios$Wards_stand, levels = wards_levels)
  
  ### Create model matrix for prediction
  mod_matrix <- model.matrix(
    ~ Study_duration + Age_group + WHO_region + Diagnosis + Wards_stand,
    data = scenarios
  )
  
  mod_matrix <- mod_matrix[, -1]
  
  ### Predict on logit scale 
  predictions <- predict(
    model,
    newmods = mod_matrix,
    transf = NULL
  )
  
  return(predictions)
}

 -----------------------------------------------------------
### Monte Carlo simulation to quantify prediction uncertainty
monte_carlo_from_predict <- function(pred_obj, n_sim = 1000) {

  mu_logit <- pred_obj$pred
  se_logit <- pred_obj$se
  n_scenarios <- length(mu_logit)
  
  logit_sim_matrix <- matrix(NA, nrow = n_sim, ncol = n_scenarios)
  prob_sim_matrix <- matrix(NA, nrow = n_sim, ncol = n_scenarios)
  
  for (j in 1:n_scenarios) {
    logit_sim <- rnorm(n_sim, mean = mu_logit[j], sd = se_logit[j])
    logit_sim_matrix[, j] <- logit_sim
    prob_sim_matrix[, j] <- plogis(logit_sim)  # Transform to probability scale
  }
  
  # Scenario labels
  scenario_names <- data.frame(
    income_group = rep(c("HICs", "UMICs", "LMICs", "LICs"), 2),
    age_group = rep(c("＜1", "＜5"), each = 4)
  )
  
  # Summary statistics for each scenario
  summary_stats <- data.frame(
    Scenario = 1:n_scenarios,
    Income_Group = scenario_names$income_group,
    Age_Group = scenario_names$age_group,
    
    Logit_Mean = mu_logit,
    Logit_SE = se_logit,
    Logit_CI_Lower = pred_obj$ci.lb,
    Logit_CI_Upper = pred_obj$ci.ub,
    
    # Probability scale (from simulation)
    Pred_Prob_Median = apply(prob_sim_matrix, 2, median),
    Pred_Prob_Mean = apply(prob_sim_matrix, 2, mean),
    Pred_Prob_SD = apply(prob_sim_matrix, 2, sd),
    Pred_Prob_2.5 = apply(prob_sim_matrix, 2, quantile, 0.025),
    Pred_Prob_25 = apply(prob_sim_matrix, 2, quantile, 0.25),
    Pred_Prob_75 = apply(prob_sim_matrix, 2, quantile, 0.75),
    Pred_Prob_97.5 = apply(prob_sim_matrix, 2, quantile, 0.975),
    
    # Probability scale (direct from logit prediction)
    Direct_Prob = plogis(mu_logit),
    Direct_CI_Lower = plogis(pred_obj$ci.lb),
    Direct_CI_Upper = plogis(pred_obj$ci.ub),
    
    stringsAsFactors = FALSE
  )
  
  numeric_cols <- sapply(summary_stats, is.numeric)
  summary_stats[numeric_cols] <- round(summary_stats[numeric_cols], 4)
  
  scenario_labels <- paste0(scenario_names$income_group, "_", scenario_names$age_group)
  
  raw_draws <- data.frame(
    scenario = rep(scenario_labels, each = n_sim),
    draw_id = rep(1:n_sim, times = n_scenarios),
    predicted_probability = as.vector(prob_sim_matrix),
    stringsAsFactors = FALSE
  )
  
  return(list(
    raw_draws = raw_draws,
    summary_stats = summary_stats,
    logit_sim_matrix = logit_sim_matrix,
    prob_sim_matrix = prob_sim_matrix,
    raw_pred = pred_obj,
    model_info = list(
      n_simulations = n_sim,
      n_scenarios = n_scenarios
    )
  ))
}

### Run prediction and simulation
set.seed(123)
pred_logit <- predict_effects_logit(glmm_model, reg_data)
mc_results <- monte_carlo_from_predict(pred_logit, n_sim = 1000)
rate_predict <- mc_results$raw_draws
lvl <- c("LICs_＜1", "LMICs_＜1", "UMICs_＜1", "HICs_＜1",
         "LICs_＜5", "LMICs_＜5", "UMICs_＜5", "HICs_＜5")
rate_predict$scenario <- factor(rate_predict$scenario, levels = lvl)
rate_predict <- arrange(rate_predict, scenario)

----------------------------------------------------------------------------------------------------------------
## 07-2 Monte Carlo simulation for hospital admission population estimates (data source: Li et al., 2022)
hospital_data <- read_excel(
  "Data/RSV_hospital_data.xlsx",
  col_names = c("income_level", "age_group", "hospital_admission", "LCL", "UCL"),
  skip = 1
)

### Monte Carlo simulation for population effects
pop_effects_mc <- function(pop_data, R = 1000, seed = 123) {
  set.seed(seed)
  pop_data <- pop_data %>%
    mutate(
      income_group = income_level,
      age_group = age_group
    )
  
  results_list <- list()
  full_simulations <- list()
  
  for (i in seq_len(nrow(pop_data))) {
    cur <- pop_data[i, ]
    
    mu_log <- log(cur$hospital_admission)
    sigma_log <- (log(cur$UCL) - log(cur$LCL)) / 3.92 
    
    pop_sim <- rlnorm(R, meanlog = mu_log, sdlog = sigma_log)
    
    # Summary table
    results_list[[i]] <- data.frame(
      income_group = as.character(cur$income_group),
      age_group = as.character(cur$age_group),
      pop_median = round(median(pop_sim), 0),
      pop_lower = round(quantile(pop_sim, 0.025), 0),
      pop_upper = round(quantile(pop_sim, 0.975), 0),
      scenario_id = i,
      stringsAsFactors = FALSE
    )
    
    full_simulations[[i]] <- pop_sim
  }
  
  # Combine results
  summary_table <- do.call(rbind, results_list) %>%
    `row.names<-`(NULL)
  
  names(full_simulations) <- paste0(
    pop_data$income_group, "_", pop_data$age_group
  )
  
  return(list(
    summary = summary_table[, c("income_group", "age_group", "pop_median", "pop_lower", "pop_upper")],
    full_simulations = full_simulations,
    scenarios = pop_data,
    model_info = list(n_simulations = R)
  ))
}

### Run simulation
pop_mc_results <- pop_effects_mc(hospital_data, R = 1000)
pop_draws_df <- pop_mc_results$full_simulations %>%
  imap_dfr(~ tibble(scenario = .y, population_draw = .x))

### View results
print(pop_mc_results$summary)
unique(pop_draws_df$scenario)

-------------------------------------------------------------------------------------------------------------------------------------------
## 07-3 Population estimates of antibiotics use
global_pop <- pop_draws_df[pop_draws_df$scenario %in% c("Global_＜1", "Global_＜5"), ] # Separate global population estimates for later use

pop_draws_df <- pop_draws_df[!grepl("Global", pop_draws_df$scenario, ignore.case = TRUE), ]#Remove global scenarios from main data, keeping income-age groups for aggregation

### Combine population draws with rate draws
df_final <- data.frame(
  pop_draw  = pop_draws_df$population_draw,
  rate_draw = rate_predict$predicted_probability,
  inf_draw  = pop_draws_df$population_draw * rate_predict$predicted_probability,
  group1 = rate_predict$scenario,
  group2 = pop_draws_df$scenario
)

### Summary by income-age group 
summary_tbl <- df_final %>%
  group_by(group1) %>%
  summarise(
    inf_median = round(median(inf_draw) / 1000) * 1000,
    inf_lower  = round(quantile(inf_draw, 0.025) / 1000) * 1000,
    inf_upper  = round(quantile(inf_draw, 0.975) / 1000) * 1000,
    .groups = "drop"
  )

print(summary_tbl)

### Check matching between population and rate scenarios
df_check <- data.frame(
  pop_row   = seq_along(pop_draws_df$population_draw),
  rate_row  = seq_along(rate_predict$predicted_probability),
  scenario_pop = pop_draws_df$scenario,
  scenario_rate = rate_predict$scenario,
  pop_draw  = pop_draws_df$population_draw,
  rate_draw = rate_predict$predicted_probability
) %>%
  mutate(matched = scenario_pop == scenario_rate)

table(df_check$matched)

### Calculate global antibiotics use 
global_pop_1_vec <- global_pop[global_pop$scenario == "Global_＜1", "population_draw", drop = TRUE]
global_pop_5_vec <- global_pop[global_pop$scenario == "Global_＜5", "population_draw", drop = TRUE]

split_sum <- split(df_final, df_final$group1)

### Children under 1 year
LICs_1 <- split_sum$"LICs_＜1"
LMICs_1 <- split_sum$"LMICs_＜1"
UMICs_1 <- split_sum$"UMICs_＜1"
HICs_1 <- split_sum$"HICs_＜1"

global_inf_1 <- LICs_1$inf_draw + LMICs_1$inf_draw + UMICs_1$inf_draw + HICs_1$inf_draw
global_rate_1 <- global_inf_1 / global_pop_1_vec

summary_inf_1 <- quantile(global_inf_1, c(0.025, 0.5, 0.975))
summary_rate_1 <- quantile(global_rate_1, c(0.025, 0.5, 0.975))

### Children under 5 years
LICs_5 <- split_sum$"LICs_＜5"
LMICs_5 <- split_sum$"LMICs_＜5"
UMICs_5 <- split_sum$"UMICs_＜5"
HICs_5 <- split_sum$"HICs_＜5"

global_inf_5 <- LICs_5$inf_draw + LMICs_5$inf_draw + UMICs_5$inf_draw + HICs_5$inf_draw
global_rate_5 <- global_inf_5 / global_pop_5_vec

summary_inf_5 <- quantile(global_inf_5, c(0.025, 0.5, 0.975))
summary_rate_5 <- quantile(global_rate_5, c(0.025, 0.5, 0.975))
