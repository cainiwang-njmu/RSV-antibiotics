#05 Sensitivity analyses
file_path <- "Data/data_analysis.xlsx"
data_4 <- read_excel(file_path, sheet = "Quality_evaluation")
filtered_data4 <- data_4 %>%
  filter(ID != 0 & !is.na(ID) & Bias != 0)
filtered_data4$ID <- as.character(filtered_data4$ID)

-----------------------------------------------------------
#Sensitivity 1: excluding studies with high risk of biases
## Any pathogens
combined_pathogens <- left_join(filtered_anypathogens, filtered_data4, by = "ID") %>%
  filter(Bias %in% c("low", "medium"))
quality_pathogens <- rma.glmm(
  measure = "PLO",
  data = combined_pathogens,
  xi = Use,
  ni = Case,
  mods = ~ WHO_region - 1
)
summary(quality_pathogens)
res.qp <- data.frame(
  est = as.numeric(quality_pathogens$b),
  se = quality_pathogens$se,
  n.all = nrow(combined_pathogens)
)
res.qp$prop.est <- transf.ilogit(res.qp$est)
res.qp$prop.lci <- transf.ilogit(res.qp$est - 1.96 * res.qp$se)
res.qp$prop.uci <- transf.ilogit(res.qp$est + 1.96 * res.qp$se)

## Any bacteria
combined_bacteria <- left_join(filtered_anybacteria, filtered_data4, by = "ID") %>%
  filter(Bias %in% c("low", "medium"))
quality_bacteria <- rma.glmm(
  measure = "PLO",
  data = combined_bacteria,
  xi = Use,
  ni = Case,
  method = "FE",
  mods = ~ WHO_region - 1
)
summary(quality_bacteria)
res.qb <- data.frame(
  est = as.numeric(quality_bacteria$b),
  se = quality_bacteria$se,
  n.all = nrow(combined_bacteria)
)
res.qb$prop.est <- transf.ilogit(res.qb$est)
res.qb$prop.lci <- transf.ilogit(res.qb$est - 1.96 * res.qb$se)
res.qb$prop.uci <- transf.ilogit(res.qb$est + 1.96 * res.qb$se)

## Any viruses
combined_viruses <- left_join(filtered_anyviruses, filtered_data4, by = "ID") %>%
  filter(Bias %in% c("low", "medium"))
quality_viruses <- rma.glmm(
  measure = "PLO",
  data = combined_viruses,
  xi = Use,
  ni = Case,
  mods = ~ WHO_region - 1
)
summary(quality_viruses)
res.qv <- data.frame(
  est = as.numeric(quality_viruses$b),
  se = quality_viruses$se,
  n.all = nrow(combined_viruses)
)
res.qv$prop.est <- transf.ilogit(res.qv$est)
res.qv$prop.lci <- transf.ilogit(res.qv$est - 1.96 * res.qv$se)
res.qv$prop.uci <- transf.ilogit(res.qv$est + 1.96 * res.qv$se)

## Mono-infection
combined_mono <- left_join(filtered_mono, filtered_data4, by = "ID") %>%
  filter(Bias %in% c("low", "medium"))
quality_mono <- rma.glmm(
  measure = "PLO",
  data = combined_mono,
  xi = Use,
  ni = Case,
  mods = ~ WHO_region - 1
)
summary(quality_mono)
res.qm <- data.frame(
  est = as.numeric(quality_mono$b),
  se = quality_mono$se,
  n.all = nrow(combined_mono)
)
res.qm$prop.est <- transf.ilogit(res.qm$est)
res.qm$prop.lci <- transf.ilogit(res.qm$est - 1.96 * res.qm$se)
res.qm$prop.uci <- transf.ilogit(res.qm$est + 1.96 * res.qm$se)

--------------------------------------------------------------------------------
#Sensitivity 2: using pre-intervention data only (any pathogens, mono-infection)
## Any pathogens
combined_pathogens <- left_join(filtered_anypathogens, filtered_data4, by = "ID")
combined_pathogens2 <- combined_pathogens %>%
  filter(!(ID %in% c(25938, 28895, 28896, 211)))
pathogens_before <- data_2 %>%
  filter(Group == "anypathogens_before")
pathogens_before2 <- left_join(pathogens_before, data_1, by = "ID")
combine_pathogensbefore <- bind_rows(combined_pathogens2, pathogens_before2)
table(combine_pathogensbefore$WHO_region)

before_pathogens <- rma.glmm(
  measure = "PLO",
  data = combine_pathogensbefore,
  xi = Use,
  ni = Case,
  mods = ~ WHO_region - 1
)
summary(before_pathogens)
res.bp <- data.frame(
  est = as.numeric(before_pathogens$b),
  se = before_pathogens$se,
  n.all = nrow(combine_pathogensbefore)
)
res.bp$prop.est <- transf.ilogit(res.bp$est)
res.bp$prop.lci <- transf.ilogit(res.bp$est - 1.96 * res.bp$se)
res.bp$prop.uci <- transf.ilogit(res.bp$est + 1.96 * res.bp$se)

## Mono-infection
combined_mono <- left_join(filtered_mono, filtered_data4, by = "ID")
combined_mono2 <- combined_mono %>%
  filter(!(ID %in% c(30179, 30430)))
mono_before <- data_2 %>%
  filter(Group == "mono_before")
mono_before2 <- left_join(mono_before, data_1, by = "ID")
combine_monobefore <- bind_rows(combined_mono2, mono_before2)
table(combine_monobefore$WHO_region)
before_mono <- rma.glmm(
  measure = "PLO",
  data = combine_monobefore,
  xi = Use,
  ni = Case,
  mods = ~ WHO_region - 1
)
summary(before_mono)
res.mb <- data.frame(
  est = as.numeric(before_mono$b),
  se = before_mono$se,
  n.all = nrow(combine_monobefore)
)
res.mb$prop.est <- transf.ilogit(res.mb$est)
res.mb$prop.lci <- transf.ilogit(res.mb$est - 1.96 * res.mb$se)
res.mb$prop.uci <- transf.ilogit(res.mb$est + 1.96 * res.mb$se)

-----------------------------------------------------------
#Sensitivity 3: including only PCR-based studies
## Any pathogens
pathogens_PCR <-filtered_anypathogens %>%
  filter(V_method_Stan == "PCR")
table(pathogens_PCR$WHO_region)
PCR_pathogens<- rma.glmm(measure = "PLO", data = pathogens_PCR, 
                         xi = pathogens_PCR$Use,
                         ni = pathogens_PCR$Case,
                         mods = ~ WHO_region - 1
)
summary(PCR_pathogens)
res.pcrp <- data.frame(est = as.numeric(PCR_pathogens$b),
                       se = PCR_pathogens$se,
                       n.all = nrow(pathogens_PCR))
res.pcrp$prop.est <- transf.ilogit(res.pcrp$est)
res.pcrp$prop.lci <- transf.ilogit(res.pcrp$est - 1.96 * res.pcrp$se)
res.pcrp$prop.uci <- transf.ilogit(res.pcrp$est + 1.96 * res.pcrp$se)

bacteria_PCR <-filtered_anybacteria %>%
  filter(V_method_Stan == "PCR")
table(bacteria_PCR $WHO_region)

## Any bacteria
PCR_bacteria<- rma.glmm(measure = "PLO", data = bacteria_PCR, 
                        xi = bacteria_PCR$Use,
                        ni = bacteria_PCR$Case,
                        mods = ~ WHO_region - 1,
                        method ="FE"
                        
)
summary(PCR_bacteria)
res.pcrb <- data.frame(est = as.numeric(PCR_bacteria$b),
                       se = PCR_bacteria$se,
                       n.all = nrow(bacteria_PCR))
res.pcrb$prop.est <- transf.ilogit(res.pcrb$est)
res.pcrb$prop.lci <- transf.ilogit(res.pcrb$est - 1.96 * res.pcrb$se)
res.pcrb$prop.uci <- transf.ilogit(res.pcrb$est + 1.96 * res.pcrb$se)

## Any viruses
viruses_PCR <-filtered_anyviruses %>%
  filter(V_method_Stan == "PCR")
table(viruses_PCR $WHO_region)
PCR_viruses<- rma.glmm(measure = "PLO", data = viruses_PCR, 
                       xi = viruses_PCR$Use,
                       ni = viruses_PCR$Case,
                       mods = ~ WHO_region - 1
                       
)
summary(PCR_viruses)
res.pcrv <- data.frame(est = as.numeric(PCR_viruses$b),
                       se = PCR_viruses$se,
                       n.all = nrow(viruses_PCR))
res.pcrv$prop.est <- transf.ilogit(res.pcrv$est)
res.pcrv$prop.lci <- transf.ilogit(res.pcrv$est - 1.96 * res.pcrv$se)
res.pcrv$prop.uci <- transf.ilogit(res.pcrv$est + 1.96 * res.pcrv$se)

## Mono-infection
mono_PCR <-filtered_mono %>%
  filter(V_method_Stan == "PCR")
table(mono_PCR $WHO_region)
mono_PCR <- mono_PCR %>%
  filter(!(ID %in% c( 30430)))
PCR_mono<- rma.glmm(measure = "PLO", data = mono_PCR, 
                    xi = mono_PCR$Use,
                    ni = mono_PCR$Case
                    
)
summary(PCR_mono)
res.pcrm <- data.frame(est = as.numeric(PCR_mono$b),
                       se = PCR_mono$se,
                       n.all = nrow(mono_PCR))
res.pcrm$prop.est <- transf.ilogit(res.pcrm$est)
res.pcrm$prop.lci <- transf.ilogit(res.pcrm$est - 1.96 * res.pcrm$se)
res.pcrm$prop.uci <- transf.ilogit(res.pcrm$est + 1.96 * res.pcrm$se)

-----------------------------------------------------------
#05 Sensitivity analyses: Figure S3
dat_ses <- read_excel("Data/sensitivity_glmm.xlsx")

dat_ses$`Sensitivity analysis` <- factor(
  dat_ses$`Sensitivity analysis`,
  levels = c(
    "Main analysis",
    "Excluding studies with high risk of biases",
    "Including only PCR-based studies",
    "Using pre-intervention data only"
  )
)

dat_ses$Country <- factor(
  dat_ses$Country,
  levels = c("HICs", "UMICs", "LMICs", "LICs"),
  ordered = TRUE
)

## Create sensitivity plot
p_type2 <- ggplot(dat_ses, aes(x = Country, color = `Sensitivity analysis`)) +
  
  geom_point(aes(y = PP), position = position_dodge2(width = 0.5)) +
  geom_point(aes(y = PP, color = `Sensitivity analysis`),
             position = position_dodge2(width = 0.5),
             size = 3, alpha = 0.35) +

  geom_errorbar(aes(ymin = L95CI, ymax = H95CI, color = `Sensitivity analysis`),
                width = 0.5, linewidth = 0.6, position = position_dodge2(width = 0.5)) +
  
  scale_y_continuous(
    name = "Proportion (95% CI)",
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
    limits = c(0, 1),
    expand = c(0, 0.08)
  ) +
  
  labs(x = "Country Income Group") +
  
  scale_color_lancet(
    name = "Sensitivity analysis",
    breaks = c(
      "Main analysis",
      "Excluding studies with high risk of biases",
      "Including only PCR-based studies",
      "Using pre-intervention data only"
    )
  ) +
  
  facet_wrap(
    ~ factor(Type,
             levels = c("any病原体", "any细菌", "any病毒", "单独感染"),
             labels = c(
               "Overall RSV-infection group",
               "RSV-bacterial co-infection group",
               "RSV-viral co-infection group",
               "RSV mono-infection group"
             )
    ),
    ncol = 2
  ) +
  
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.margin = margin(t = 5, b = 5),
    plot.margin = margin(t = 15, r = 15, b = 5, l = 15, unit = "mm"),
    text = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3))) +
  
  geom_text(aes(x = Country, y = 0, label = Number, color = `Sensitivity analysis`),
            position = position_dodge(0.7),
            size = 3, fontface = "bold", vjust = 1.2, show.legend = FALSE)

print(p_type2)

ggsave(
  filename = "Output/Figures/sensitivity_glmm.pdf",
  plot = p_type2,
  width = 13,
  height = 8,
  dpi = 300
)
