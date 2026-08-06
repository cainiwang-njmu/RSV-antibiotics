#01 Figure 1: Geographical distribution of included studies
## Load world map shapefile
map.R1 <- readOGR(dsn = "Location Data/R1_shp", stringsAsFactors = TRUE)

## Convert to ggplot2-compatible data frame
map.data2 <- fortify(map.R1)

## Extract map metadata (country IDs)
region.metadata2 <- as.data.frame(map.R1)
region.metadata2$id <- row.names(region.metadata2)

## Load WHO region classification and merge
region.raw2 <- read_excel("Location Data/Region2.xlsx")
region.metadata2 <- left_join(region.metadata2, region.raw2, by = "COUNTRY")
map.data2 <- left_join(map.data2, region.metadata2, by = "id")

## Set WHO region display labels
map.data2$whoreg <- factor(
  map.data2$whoreg,
  levels = c("Afr", "Amr", "Emr", "Eur", "Sear", "Wpr"),
  labels = c("Africa region", "Region of the Americas",
             "Eastern Mediterranean region", "European region",
             "Southeast Asia region", "Western Pacific region")
)

## Exclude Antarctica
map.data2 <- map.data2 %>% filter(COUNTRY != "Antarctica")

## Load study data and extract coordinates
file_path <- "Data/data_analysis.xlsx"
sheet_name <- "1. Data Cataloguing"
data <- read_excel(file_path, sheet = sheet_name, skip = 1, col_names = TRUE)

## Filter by age group
filtered_df <- data %>% filter(Age_group != 0 & !is.na(Age_group))

## Extract columns for mapping (Country, LON, LAT)
data_formap <- filtered_df[, c(1, 10, 11, 12)]

## Generate map
p <- ggplot(
  data = rbind(map.data2[map.data2$COUNTRY != "Sudan", ], 
               map.data2[map.data2$COUNTRY == "Sudan", ]), 
  aes(x = long, y = lat, group = group, fill = whoreg)
) +
  geom_polygon(colour = "grey50", size = 0.05, alpha = 0.75) +
  scale_fill_manual(
    name = NULL, 
    values = hsv(c(210, 0, 60, 150, 288, 20) / 360, c(0.3), c(0.9)), 
    na.translate = FALSE
  ) +
  geom_jitter(
    data = data_formap,
    aes(x = LON, y = LAT, group = Country),
    fill = "red2", colour = "white", shape = 21, size = 2,
    width = 1.8, height = 1.8, alpha = 0.7
  ) +
  scale_shape_manual(name = NULL, values = c(21, 22, 24)) +
  ggthemes::theme_map() +
  theme(text = element_text(size = 12)) +
  guides(
    shape = guide_legend(override.aes = list(size = 5)),
    fill = guide_legend(override.aes = list(size = 3.5, colour = NA))
  )

## Save map
ggsave(filename = "Output/Figure 1/map_2025.pdf", plot = p, width = 10, height = 5)
