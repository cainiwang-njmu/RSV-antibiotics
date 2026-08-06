#03 Subgroup analysis: RSV any pathogens
pathogens_sub<-filtered_anypathogens%>%
  filter(Age_stand != 0 ) #Select studies for subgroup analysis

## HICs
sub_HICs<-pathogens_sub%>%
  filter(WHO_region=="HICs")
labels_HICs <- paste(sub_HICs$"Author & year")
meta_HICs <- rma.glmm(measure = "PLO", data = sub_HICs, 
                      xi = sub_HICs$Use,
                      ni = sub_HICs$Case,
                      mods = ~ Age_stand - 1
)
meta_HICs$slab <- labels_HICs
summary(meta_HICs)

res.HICs <- data.frame(est = as.numeric(meta_HICs$b),
                       se = meta_HICs$se,
                       n.all = nrow(sub_HICs))
res.HICs$prop.est <- transf.ilogit(res.HICs$est)
res.HICs$prop.lci <- transf.ilogit(res.HICs$est - 1.96 * res.HICs$se)
res.HICs$prop.uci <- transf.ilogit(res.HICs$est + 1.96 * res.HICs$se)

## UMICs
sub_UMICs<-pathogens_sub%>%
  filter(WHO_region=="UMICs")
labels_UMICs <- paste(sub_UMICs$"Author & year")
meta_UMICs <- rma.glmm(measure = "PLO", data = sub_UMICs, 
                       xi = sub_UMICs$Use,
                       ni = sub_UMICs$Case,
                       mods = ~ Age_stand - 1
)
meta_UMICs$slab <- labels_UMICs
summary(meta_UMICs)

res.UMICs <- data.frame(est = as.numeric(meta_UMICs$b),
                        se = meta_UMICs$se,
                        n.all = nrow(sub_UMICs))
res.UMICs$prop.est <- transf.ilogit(res.UMICs$est)
res.UMICs$prop.lci <- transf.ilogit(res.UMICs$est - 1.96 * res.UMICs$se)
res.UMICs$prop.uci <- transf.ilogit(res.UMICs$est + 1.96 * res.UMICs$se)

## LMICs
sub_LMICs<-pathogens_sub%>%
  filter(WHO_region=="LMICs")
labels_LMICs <- paste(sub_LMICs$"Author & year")
meta_LMICs <- rma.glmm(measure = "PLO", data = sub_LMICs, 
                       xi = sub_LMICs$Use,
                       ni = sub_LMICs$Case,
                       mods = ~ Age_stand - 1
)
meta_LMICs$slab <- labels_LMICs
summary(meta_LMICs)

res.LMICs <- data.frame(est = as.numeric(meta_LMICs$b),
                        se = meta_LMICs$se,
                        n.all = nrow(sub_LMICs))
res.LMICs$prop.est <- transf.ilogit(res.LMICs$est)
res.LMICs$prop.lci <- transf.ilogit(res.LMICs$est - 1.96 * res.LMICs$se)
res.LMICs$prop.uci <- transf.ilogit(res.LMICs$est + 1.96 * res.LMICs$se)

## LICs
sub_LICs<-pathogens_sub%>%
  filter(WHO_region=="LICs")
labels_LICs <- paste(sub_LICs$"Author & year")
meta_LICs <- rma.glmm(measure = "PLO", data = sub_LICs, 
                      xi = sub_LICs$Use,
                      ni = sub_LICs$Case,
                      mods = ~ Age_stand - 1
)
meta_LICs$slab <- labels_LICs
summary(meta_LICs)

res.LICs <- data.frame(est = as.numeric(meta_LICs$b),
                       se = meta_LICs$se,
                       n.all = nrow(sub_LICs))
res.LICs$prop.est <- transf.ilogit(res.LICs$est)
res.LICs$prop.lci <- transf.ilogit(res.LICs$est - 1.96 * res.LICs$se)
res.LICs$prop.uci <- transf.ilogit(res.LICs$est + 1.96 * res.LICs$se)

-----------------------------------------------------------
#04 Subgroup analysis: RSV any pathogens (Figure S2)
dat <- read_excel("Data/subgroup_forest_glmm.xlsx")

dat$group <- factor(
  dat$group,
  levels = c("HICs", "UMICs", "LMICs", "LICs"),
  ordered = TRUE
)

# Create subgroup plot
p_type1 <- ggplot(
  dat,
  aes(x = group, y = PP, color = `Age group (years)`)
) +
  
  geom_errorbar(
    aes(ymin = L95CI, ymax = H95CI),
    position = position_dodge(.7),
    width = 0.5,
    linewidth = 1
  ) +
  
  geom_point(position = position_dodge(0.7), size = 1) +
  geom_point(position = position_dodge(0.7), size = 3, alpha = 0.6) +
  
  scale_color_lancet(
    name = "Age group (years)"
  ) +
  
  labs(
    x = "Country Income Group",
    y = "Proportion 95% (CI)"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7),
    legend.margin = margin(t = 1, b = 1),
    plot.margin = margin(t = 15, r = 15, unit = "mm")
  ) +
  
  guides(colour = guide_legend(nrow = 1)) +
  
  geom_text(
    aes(x = group, y = 0, label = number,
        color = `Age group (years)`),
    position = position_dodge(0.7),
    size = 3,
    fontface = "bold",
    vjust = 1.2,
    show.legend = FALSE
  )

ggsave(
  filename = "Output/Figures/Figure S2.pdf",
  plot = p_type1,
  width = 7,
  height = 5,
  dpi = 300
)
