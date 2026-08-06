# 02 Main analysis: filter data by Age_filter and Group
## Load data from Excel
file_path <- "Data/data_analysis.xlsx"

data_1 <- read_excel(
  file_path,
  sheet = "1. Data Cataloguing",
  skip = 1,
  col_names = TRUE
)

data_1 <- data_1 %>%
  filter(!duplicated(ID))

data_2 <- read_excel(
  file_path,
  sheet = "2.Data Extraction"
)

filtered_data2 <- data_2 %>%
  filter(Age != 0 & !is.na(Age))

filtered_data2$ID <- as.character(filtered_data2$ID)
data_combine <- left_join(filtered_data2, data_1, by = "ID")

## Filter for main analysis 
data_type <- data_combine %>%
  filter(Age_filter == "Overall")

## Define infection types of interest
group_list <- c("anypathogens", "anybacteria", "anyviruses", "mono-infection")

## Filter and split by Group
filtered_list <- data_type %>%
  filter(Group %in% group_list) %>%
  split(.$Group)

## Extract individual data frames
filtered_anypathogens <- filtered_list[["anypathogens"]]
filtered_anybacteria <- filtered_list[["anybacteria"]]
filtered_anyviruses <- filtered_list[["anyviruses"]]
filtered_mono <- filtered_list[["mono-infection"]]

-----------------------------------------------------------
#02 Main analysis: RSV-bacterial co-infection group
## Note: This analysis includes only high-income countries (HICs), Low- and middle-income countries (LMICs) were excluded due to insufficient data
hic_dat1  <- subset(filtered_anybacteria, WHO_region == "HICs")
HIC_labels <- paste(hic_dat1$"Author & year")
meta_hic1  <- rma.glmm(measure = "PLO", data = hic_dat1,
                       xi = Use, ni = Case, mods = ~ 1)
meta_hic1$slab <- HIC_labels 

## Extract individual study estimates (logit scale)
res.hic1 <- data.frame(est1 =meta_hic1$yi,
                       se1 = sqrt(meta_hic1$vi),
                       n.all = nrow(hic_dat1))
res.hic1$prop.est1 <- transf.ilogit(res.hic1$est1)
res.hic1$prop.lci1 <- transf.ilogit(res.hic1$est1 - 1.96 * res.hic1$se1)
res.hic1$prop.uci1 <- transf.ilogit(res.hic1$est1 + 1.96 * res.hic1$se1)
res.hic1$slab <-meta_hic1$slab

## Extract regional pooled estimates
res.df_hic1 <- data.frame(est = as.numeric(meta_hic1$b),
                          se = meta_hic1$se,
                          n.all = nrow(hic_dat1))
res.df_hic1$prop.est <- transf.ilogit(res.df_hic1$est)
res.df_hic1$prop.lci <- transf.ilogit(res.df_hic1$est - 1.96 * res.df_hic1$se)
res.df_hic1$prop.uci <- transf.ilogit(res.df_hic1$est + 1.96 * res.df_hic1$se)

-----------------------------------------------------------
# 02 Main analysis: RSV-viral co-infection group
study_labels2 <- paste(filtered_anyviruses$"Author & year")
meta_region2 <- rma.glmm(measure = "PLO", data = filtered_anyviruses, 
                         xi = filtered_anyviruses$Use,
                         ni = filtered_anyviruses$Case,
                         mods = ~ WHO_region - 1 
                         )
meta_region2$slab2 <- study_labels2

res.df2 <- data.frame(est =meta_region2$yi,
                      se = sqrt(meta_region2$vi),
                      n.all = nrow(filtered_anyviruses))
res.df2$prop.est <- transf.ilogit(res.df2$est)
res.df2$prop.lci <- transf.ilogit(res.df2$est - 1.96 * res.df2$se)
res.df2$prop.uci <- transf.ilogit(res.df2$est + 1.96 * res.df2$se)
res.df2$slab <-meta_region2$slab2

res.df_combine2 <- data.frame(est = as.numeric(meta_region2$b),
                              se = meta_region2$se,
                              n.all = nrow(filtered_anyviruses))
res.df_combine2$prop.est <- transf.ilogit(res.df_combine2$est)
res.df_combine2$prop.lci <- transf.ilogit(res.df_combine2$est - 1.96 * res.df_combine2$se)
res.df_combine2$prop.uci <- transf.ilogit(res.df_combine2$est + 1.96 * res.df_combine2$se)

## Run summary() for each model 
summary()

## Estimate heterogeneity indices for different income groups and overall
hic_dat2  <- subset(filtered_anyviruses, WHO_region == "HICs")
umic_dat2 <- subset(filtered_anyviruses, WHO_region == "UMICs")
meta_hic2  <- rma.glmm(measure = "PLO", data = hic_dat2,
                       xi = Use, ni = Case, mods = ~ 1 )
meta_umic2<- rma.glmm(measure = "PLO", data = umic_dat2,
                      xi = Use, ni = Case, mods = ~ 1)
meta_all2 <- rma.glmm(measure = "PLO", data = filtered_anyviruses, 
                      xi = filtered_anyviruses$Use,
                      ni = filtered_anyviruses$Case
)
summary()

-----------------------------------------------------------
# 02 Main analysis: RSV mono-infection group
hic_umic_dat <- subset(filtered_mono, WHO_region %in% c("HICs", "UMICs"))
study_labels3 <- paste(hic_umic_dat$"Author & year")
meta_region3 <- rma.glmm(measure = "PLO", data =hic_umic_dat, 
                         xi =hic_umic_dat$Use,
                         ni =hic_umic_dat$Case,
                         mods = ~ WHO_region - 1
)

meta_region3$slab3 <- study_labels3
summary(meta_region3)
print(meta_region3$slab3)
table(hic_umic_dat$WHO_region)

res.df3 <- data.frame(est =meta_region3$yi,
                      se = sqrt(meta_region3$vi),
                      n.all = nrow(hic_umic_dat))
res.df3$prop.est <- transf.ilogit(res.df3$est)
res.df3$prop.lci <- transf.ilogit(res.df3$est - 1.96 * res.df3$se)
res.df3$prop.uci <- transf.ilogit(res.df3$est + 1.96 * res.df3$se)
res.df3$slab <-meta_region3$slab3

res.df_combine3 <- data.frame(est = as.numeric(meta_region3$b),
                              se = meta_region3$se,
                              n.all = nrow(hic_umic_dat))
res.df_combine3$prop.est <- transf.ilogit(res.df_combine3$est)
res.df_combine3$prop.lci <- transf.ilogit(res.df_combine3$est - 1.96 * res.df_combine3$se)
res.df_combine3$prop.uci <- transf.ilogit(res.df_combine3$est + 1.96 * res.df_combine3$se)

hic_dat3  <- subset(filtered_mono, WHO_region == "HICs")
umic_dat3 <- subset(filtered_mono, WHO_region == "UMICs")
meta_hic3<- rma.glmm(measure = "PLO", data = hic_dat3,
                      xi = Use, ni = Case, mods = ~ 1)
meta_umic3<- rma.glmm(measure = "PLO", data = umic_dat3,
                      xi = Use, ni = Case, mods = ~ 1)
meta_all3 <- rma.glmm(measure = "PLO", data = hic_umic_dat, 
                      xi = hic_umic_dat$Use,
                      ni = hic_umic_dat$Case)
summary()

-----------------------------------------------------------
# 02 Main analysis: RSV with any pathogen group
study_labels4 <- paste(filtered_anypathogens$"Author & year")
meta_region4 <- rma.glmm(measure = "PLO", data = filtered_anypathogens, 
                         xi = filtered_anypathogens$Use,
                         ni = filtered_anypathogens$Case,
                         mods = ~ WHO_region - 1
)
meta_region4$slab4 <- study_labels4
summary(meta_region4)
print(meta_region4$slab4)
table(filtered_anypathogens$WHO_region)

res.df4 <- data.frame(est =meta_region4$yi,
                      se = sqrt(meta_region4$vi),
                      n.all = nrow(filtered_anypathogens))
res.df4$prop.est <- transf.ilogit(res.df4$est)
res.df4$prop.lci <- transf.ilogit(res.df4$est - 1.96 * res.df4$se)
res.df4$prop.uci <- transf.ilogit(res.df4$est + 1.96 * res.df4$se)
res.df4$slab <-meta_region4$slab4

res.df_combine4 <- data.frame(est = as.numeric(meta_region4$b),
                              se = meta_region4$se,
                              n.all = nrow(filtered_anypathogens))
res.df_combine4$prop.est <- transf.ilogit(res.df_combine4$est)
res.df_combine4$prop.lci <- transf.ilogit(res.df_combine4$est - 1.96 * res.df_combine4$se)
res.df_combine4$prop.uci <- transf.ilogit(res.df_combine4$est + 1.96 * res.df_combine4$se)

hic_dat4  <- subset(filtered_anypathogens, WHO_region == "HICs")
umic_dat4 <- subset(filtered_anypathogens, WHO_region == "UMICs")
lmic_dat4<- subset(filtered_anypathogens, WHO_region == "LMICs")
lic_dat4<- subset(filtered_anypathogens, WHO_region == "LICs")

meta_hic4  <- rma.glmm(measure = "PLO", data = hic_dat4,
                       xi = Use, ni = Case, mods = ~ 1 )
meta_umic4<- rma.glmm(measure = "PLO", data = umic_dat4,
                      xi = Use, ni = Case, mods = ~ 1)
meta_lmic4<- rma.glmm(measure = "PLO", data = lmic_dat4,
                      xi = Use, ni = Case, mods = ~ 1)
meta_lic4<- rma.glmm(measure = "PLO", data = lic_dat4,
                     xi = Use, ni = Case, mods = ~ 1)
meta_all4 <- rma.glmm(measure = "PLO", data = filtered_anypathogens, 
                      xi = filtered_anypathogens$Use,
                      ni = filtered_anypathogens$Case
)
summary()

-----------------------------------------------------------
#02 Main analysis: forest plot (Figure 3)
## RSV mono-infection group
file_path <- "Data/forest_mono.xlsx"
forest <- read_excel(file_path)

# Format text columns for display
forest$V1 <- ifelse(is.na(forest$V3),
                    forest$V1,
                    paste0("   ", forest$V1))
forest$V3 <- ifelse(is.na(forest$V3), "", forest$V3)
forest$V4 <- ifelse(is.na(forest$V4), "", forest$V4)

n_row   <- nrow(forest)
txt_col <- rep("#353338", n_row)
txt_col[c(2, 4)] <- "#1B325F"

# Ensure Output/Figures directory exists
if (!dir.exists("Output/Figures")) {
  dir.create("Output/Figures", recursive = TRUE)
}

# save as PDF 
pdf_file <- "Output/Figures/forest_mono_infection_glmm.png"
png(pdf_file,width = 3000, height = 2100, res = 300)

# Forest plot
forest1 <- forestplot(
  labeltext = as.matrix(forest[, 1:5]),
  mean = forest$V6,
  lower = forest$V7,
  upper = forest$V8,
  is.summary = c(T, T, F, T, F, F, F, F, F, F, F, T, F, T, F, F, F),
  hrzl_lines = list(
    '1' = gpar(lty = 1, lwd = 2),
    '2' = gpar(lty = 1, lwd = 2),
    '14' = gpar(lty = 1, lwd = 2)
  ),
  zero = NA,
  boxsize = 0.4,
  lineheight = unit(8, 'mm'),
  colgap = unit(2, 'mm'),
  lwd.zero = 2,
  lwd.ci = 2,
  col = fpColors(box = '#77a8c5', summary = "#285b90",
                 lines = 'black', text = txt_col, zero = '#4D4D4D'),
  lwd.xaxis = 2,
  lty.ci = "solid",
  graph.pos = 5,
  xticks = c(0.25, 0.5, 0.75, 1),
  txt_gp = fpTxtGp(ticks = gpar(cex = 1.25))
)

print(forest1)

grid.text("RSV mono-infection group",
          x = unit(0.02, "npc"),
          y = unit(0.95, "npc"),
          just = "left",
          gp = gpar(fontsize = 18, fontface = "bold"))

dev.off()

-----------------------------------------------------------
## RSV-viral co-infection group
forest2<-read_excel("Data/forest_viral.xlsx")
forest2$V1 <- ifelse(is.na(forest2$V3),
                     forest2$V1,
                     paste0("   ", forest2$V1))
forest2$V3 <- ifelse(is.na(forest2$V3), "", forest2$V3)
forest2$V4 <- ifelse(is.na(forest2$V4 ), "", forest2$V4 )
n_row2 <- nrow(forest2)  
txt_col2 <- rep("#353338", n_row2)
txt_col2[c(2, 9)] <- "#1B325F"  
pdf_file <- "Output/Figures/forest_viral_glmm.png"
png(pdf_file,width = 3000, height = 2200, res = 300)

forest2 <- forestplot(
  labeltext = as.matrix(forest2[, 1:5]),
  mean = forest2$V6,
  lower = forest2$V7,
  upper = forest2$V8,
  is.summary = c(T, T, F, F, F, F, T, F, T, F, F, F, T, F, T, F, F, F),
  hrzl_lines = list(
    '1' = gpar(lty = 1, lwd = 2),
    '2' = gpar(lty = 1, lwd = 2),
    '15' = gpar(lty = 1, lwd = 2)
  ),
  zero = NA,
  lty.zero = 2,
  boxsize = 0.4,
  lineheight = unit(8, 'mm'),
  colgap = unit(2, 'mm'),
  lwd.zero = 2,
  lwd.ci = 2,
  col = fpColors(box = '#77a8c5', summary = "#285b90", lines = 'black', text = txt_col2, zero = '#4D4D4D'),
  lwd.xaxis = 2,
  lty.ci = "solid",
  graph.pos = 5,
  xticks = c(0.25, 0.5, 0.75, 1),
  txt_gp = fpTxtGp(ticks = gpar(cex = 1.25))
)

print(forest2)
grid.text("RSV-viral co-infection group",
          x = unit(0.02, "npc"),
          y = unit(0.95, "npc"),
          just = "left",
          gp = gpar(fontsize = 18, fontface = "bold"))

dev.off()

-----------------------------------------------------------
## RSV-bacterial co-infection group
file_path <- "Data/forest_bacterial.xlsx"
forest3 <- read_excel(file_path)

forest3$V1 <- ifelse(is.na(forest3$V3),
                     forest3$V1,
                     paste0("   ", forest3$V1))
forest3$V3 <- ifelse(is.na(forest3$V3), "", forest3$V3)
forest3$V4 <- ifelse(is.na(forest3$V4), "", forest3$V4)

n_row3 <- nrow(forest3)
txt_col3 <- rep("#353338", n_row3)
txt_col3[c(2, 10)] <- "#1B325F"
pdf_file <- "Output/Figures/forest_bacterial_glmm.png"
png(pdf_file,width = 3000, height = 1250, res = 300)

forest3 <- forestplot(
  labeltext = as.matrix(forest3[, 1:5]),
  mean = forest3$V6,
  lower = forest3$V7,
  upper = forest3$V8,
  is.summary = c(T, T, F, F, F, F, F, T, F),
  hrzl_lines = list(
    '1' = gpar(lty = 1, lwd = 2),
    '2' = gpar(lty = 1, lwd = 2),
    '8' = gpar(lty = 1, lwd = 2)
  ),
  zero = NA,
  boxsize = 0.4,
  lineheight = unit(8, 'mm'),
  colgap = unit(2.2, 'mm'),
  lwd.zero = 2,
  lwd.ci = 2,
  col = fpColors(box = '#77a8c5', summary = "#285b90", lines = 'black', text = txt_col3, zero = '#4D4D4D'),
  lwd.xaxis = 2,
  lty.ci = "solid",
  graph.pos = 5,
  xticks = c(0.25, 0.5, 0.75, 1),
  txt_gp = fpTxtGp(ticks = gpar(cex = 1.20))
)

print(forest3)

grid.text("RSV-bacterial co-infection group",
          x = unit(0.02, "npc"),
          y = unit(0.95, "npc"),
          just = "left",
          gp = gpar(fontsize = 18, fontface = "bold"))

dev.off()

-----------------------------------------------------------
## Combine forest plots: mono, viral, bacterial
image_paths <- c(
  "Output/Figures/forest_bacterial_glmm.png",
  "Output/Figures/forest_viral_glmm.png",
  "Output/Figures/forest_mono_infection_glmm.png"
)
images <- lapply(image_paths, image_read)
max_width <- max(sapply(images, function(img) image_info(img)$width))
images_padded <- lapply(images, function(img) {
  info <- image_info(img)
  image_extent(img, 
               geometry = paste0(max_width, "x", info$height), 
               gravity = "center",
               color = "white")
})

combined_image <- image_append(image_join(images_padded), stack = TRUE)
output_path <- "Output/Figures/Figure 3.pdf"
image_write(combined_image, path = output_path, format = "pdf")

-----------------------------------------------------------
#02 Main analysis: forest plot(Figure 2)
file_path <- "Data/forest_any_pathogen.xlsx"
forest4 <- read_excel(file_path)

forest4$V1 <- ifelse(is.na(forest4$V3),
                     forest4$V1,
                     paste0("    ", forest4$V1))
forest4$V3 <- ifelse(is.na(forest4$V3), "", forest4$V3)
forest4$V4 <- ifelse(is.na(forest4$V4), "", forest4$V4)

n_row4 <- nrow(forest4)
txt_col4 <- rep("#353338", n_row4)

is.summary <- rep(FALSE, 111)
is.summary[c(2, 58, 61, 86, 89, 97, 100, 104, 107, 110)] <- TRUE

txt_col4[c(2, 61, 89, 100)] <- "#1B325F"

pdf_file <- "Output/Figures/Figure 1.pdf"
pdf(pdf_file, width = 10, height = 36)

forest4 <- forestplot(
  labeltext = as.matrix(forest4[, 1:5]),
  mean = forest4$V6,
  lower = forest4$V7,
  upper = forest4$V8,
  is.summary = is.summary,
  hrzl_lines = list(
    '1' = gpar(lty = 1, lwd = 2),
    '2' = gpar(lty = 1, lwd = 2),
    '107' = gpar(lty = 1, lwd = 2)
  ),
  zero = NA,
  boxsize = 0.25,
  lineheight = unit(8, 'mm'),
  colgap = unit(1, 'mm'),
  lwd.zero = 2,
  lwd.ci = 2,
  col = fpColors(box = '#77a8c5', summary = "#285b90", lines = 'black', text = txt_col4, zero = '#4D4D4D'),
  lwd.xaxis = 2,
  lty.ci = "solid",
  graph.pos = 5,
  xticks = c(0.25, 0.5, 0.75, 1),
  txt_gp = fpTxtGp(
    label = gpar(cex = 1.2),
    ticks = gpar(cex = 1.1),
    xlab = gpar(cex = 1.3),
    summary = gpar(cex = 1.3, fontface = "bold")
  )
)

print(forest4)
dev.off()
