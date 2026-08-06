#04 Meta-regression: RSV any pathogens
## Prepare data (filter missing values)
reg_data <- filtered_anypathogens %>%
  filter(
    !is.na(Study_duration),
    !is.na(Diagnosis_stan)
  )

reg_data$WHO_region <- factor(
  reg_data$WHO_region,
  levels = c("HICs", "UMICs", "LMICs", "LICs")
)

reg_data$Study_duration <- factor(
  reg_data$Study_duration,
  levels = c("before1999", "2000-2009", "2010-2019", "after2020")
)

reg_data$Age_group <- factor(
  reg_data$Age_group,
  levels = c("＜1", "＜2", "＜5")
)

reg_data$Diagnosis_stan <- factor(
  reg_data$Diagnosis_stan,
  levels = c("ARI", "ALRI", "Bronchiolitis", "Pneumonia")
)

reg_data$Wards_stand <- factor(
  reg_data$Wards_stand,
  levels = c("OP", "IP", "ER", "ICU", "others")
)

reg_data <- data.frame(
  Study_duration = reg_data$Study_duration,
  Age_group = reg_data$Age_group,
  WHO_region = reg_data$WHO_region,
  Types = reg_data$Group,
  Diagnosis = reg_data$Diagnosis_stan,
  use = reg_data$Use,
  case = reg_data$Case,
  Wards_stand = reg_data$Wards_stand
)

## Meta-regression model
glmm_model <- rma.glmm(
  measure = "PLO",
  xi = use,
  ni = case,
  mods = ~ Study_duration + Age_group + WHO_region + Diagnosis + Wards_stand,
  data = reg_data
)
summary(glmm_model)
