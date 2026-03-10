
## This script uses the output from SOM_Fire_20240424.R and evaluates the clusters using the input variables

library(tidyverse)
library(terra)
library(sf)
library(ggh4x)
library(ggpubr)
#library(devtools)
#devtools::install_github("morrowcj/remotePARTS")
library(remotePARTS)

rootDir <- "A:/_BioGeo/baumann/CHACO_fires/"
workDir <- paste0(rootDir, "11_FireSOMs/")

## Load the files and merge them
vars <- rast(paste0(workDir, "Stack_10km.tif")) %>% as.data.frame(xy=TRUE)
clust <- rast(paste0(workDir, "Cluster_v4-NoLandUse_10km.tif")) %>% as.data.frame(xy=TRUE)
all <- vars %>% left_join(clust, by=c("x", "y"))
crs <- crs(rast(paste0(workDir, "Stack_10km.tif")))

## Rename variable names
colnames(all) <- c("x", "y", "FRT", "MaxSize", "MeanBA", "MeanSize", "Clearing", "Management", "Other", "Cluster")
all_vars <- all %>% mutate(across(c(Cluster, Clearing, Management, Other, FRT, MeanSize, MeanBA), ~ifelse(FRT == 0, 0, .))) #%>% select(-c(x,y))

## Add the CHIPRS data, calculate trend
chirps <- rast(paste0(workDir, "CHIRPS_Precipitation_1985-2024_10km.tif")) %>% as.data.frame(xy=TRUE) %>% 
  setNames(c("x", "y", paste0("yr", 1985:2024))) %>% select(-"yr2024")

chirps.cols <- grep("yr", names(chirps), value = TRUE)
precipVals <- as.matrix(chirps[, chirps.cols])
coords <- as.matrix(chirps[, c("x", "y")])

ARfit <- fitAR_map(Y = precipVals, coords = coords)

ARfit$coefficients[, "t"] <- ARfit$coefficients[,"t"]/rowMeans(chirps[, chirps.cols])
chirps$AR_coef <- coefficients(ARfit)[, "t"] # save time trend coefficient
chirps_sel <- chirps %>% select(x, y, AR_coef)

## Do the same with the dry days
dryD <- rast(paste0(workDir, "ERA5_DryDays_p50_1985-2023_10km.tif" )) %>% as.data.frame(xy=TRUE) %>% 
  setNames(c("x", "y", paste0("yr", 1985:2023)))
# There are some NAs in the table, we interpolate them 
for (i in 4:(ncol(dryD) - 1)) {
  na_positions <- which(is.na(dryD[[i]]))
  dryD[na_positions, i] <- rowMeans(dryD[na_positions, c(i - 1, i + 1)], na.rm = TRUE)
}

dryD.cols <- grep("yr", names(dryD), value = TRUE)
dryDVals <- as.matrix(dryD[, dryD.cols])
coords <- as.matrix(dryD[, c("x", "y")])

ARfit_dryD <- fitAR_map(Y = dryDVals, coords = coords)

ARfit_dryD$coefficients[, "t"] <- ARfit_dryD$coefficients[,"t"]/rowMeans(dryD[, dryD.cols])
dryD$AR_coef_dry <- coefficients(ARfit_dryD)[, "t"] # save time trend coefficient
dryD_sel <- dryD %>% select(x, y, AR_coef_dry)


# Merge with the all_vars
all_vars_add <- all %>% left_join(chirps_sel, by=c("x", "y")) %>% left_join(dryD_sel, by=c("x", "y")) %>% 
  filter(Cluster > 0) %>% select(x, y, Cluster, Clearing, Management, Other, AR_coef, AR_coef_dry, FRT, MeanSize, MeanBA) %>% drop_na()

## Make some simple density plots to see if we have clusters that need to be merged

all_long <- all_vars_add %>% filter(Cluster > 0) %>% 
  mutate(across(c(AR_coef, AR_coef_dry, Clearing, Management, Other, FRT, MeanSize, MeanBA), function(x){(x - min(x)) / (max(x) - min(x))})) %>% 
  select(Cluster, AR_coef_dry, Clearing, Management, Other) %>%
  pivot_longer(cols = c(AR_coef_dry, Clearing, Management, Other), names_to = "Variable", values_to = "Value")

ggplot(all_long, aes(x = Value, color = Variable, fill = Variable)) +
  geom_density(position = "identity", alpha = 0.5) +
  facet_wrap(~ Cluster, scales = "free") +#
  labs(title = "Density Curves",
       x = "Value",
       y = "Density") +
  theme_minimal()

## Define some rules for merging --> based on some subset-plots
# 1 -> 1, 3, 7
# 2 -> 2
# 3 -> 4, 5, 6
# 4 --> 8
# 5 -> 9, 10

all_long_recl <- all_long %>% 
  mutate(ClusterNew = ifelse(Cluster == 1 | Cluster == 3 | Cluster == 7, 1, 0)) %>% 
  mutate(ClusterNew = ifelse(Cluster == 2, 2, ClusterNew)) %>% 
  mutate(ClusterNew = ifelse(Cluster == 4 | Cluster == 5 | Cluster == 6, 3, ClusterNew)) %>% 
  mutate(ClusterNew = ifelse(Cluster == 8, 4, ClusterNew)) %>% 
  mutate(ClusterNew = ifelse(Cluster == 9 | Cluster == 10, 5, ClusterNew)) %>% 
  mutate(Cluster = ClusterNew) %>% select(-ClusterNew)
 


ggplot(all_long_recl, aes(x = Value, color = Variable, fill = Variable)) +
  geom_density(position = "identity", alpha = 0.5) +
  facet_wrap(~ Cluster, scales = "free") +#
  labs(title = "Density Curves",
       x = "Value",
       y = "Density") +
  theme_minimal()

# #### Make the full plot

# Define colours
# Define custom diverging colors
map_colors <- c("0" = "#000000",  # Light grey for background
                      "1" = "#8dd3c7",  # Custom colors for clusters
                      "2" = "#ffffb3",
                      "3" = "#bebada",
                      "4" = "#fb8072")
cluster_colors <- setNames(map_colors, c("0", "1", "2", "3", "4"))
# Set axis limits to maintain the aspect ratio and avoid skewed pixels
x_range <- range(clust$x)
y_range <- range(clust$y)

## (1) The map

mapTab <- all %>% select(x, y, Cluster) %>% 
  mutate(ClusterNew = ifelse(Cluster == 1 | Cluster == 3 | Cluster == 7, 1, 0)) %>% 
  mutate(ClusterNew = ifelse(Cluster == 2, 2, ClusterNew)) %>% 
  mutate(ClusterNew = ifelse(Cluster == 4 | Cluster == 5 | Cluster == 6 | Cluster == 8, 3, ClusterNew)) %>% 
  mutate(ClusterNew = ifelse(Cluster == 9 | Cluster == 10, 4, ClusterNew)) %>% 
  mutate(Cluster = ClusterNew) %>% select(-ClusterNew)

mapPlot <- mapTab %>% 
  ggplot(aes(x=x, y=y, fill = factor(Cluster))) +
  geom_raster() +
  scale_fill_manual(
    values = cluster_colors,
    name = "Cluster",
    breaks = c("1", "2", "3", "4", "0"),
    labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "No Fires")) +
  coord_fixed(ratio = 1, xlim = x_range, ylim = y_range) +  # Ensure the aspect ratio is fixed
  theme_light() +
  #labs(color = "Cluster") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "right",
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), # Set the panel background to white
        plot.background = element_rect(fill = "white", color = NA)) +
  guides(fill = guide_legend(
    title.position = "top",
    label.position = "right",
    keyheight = unit(0.5, "cm"),  # Control height of legend keys
    ncol = 1,  # Single column legend
    byrow = FALSE,  # Arrange items by row
    override.aes = list(size = 4)),
    )


## (2) The fire types

facet_colors <- c("1" = "#8dd3c7",  # Custom colors for clusters, here zero is removed
                  "2" = "#ffffb3",
                  "3" = "#bebada",
                  "4" = "#fb8072",
                  "5" = "#80b1d3")

get_cluster_color <- function(cluster) {as.character(facet_colors[as.character(cluster)])}
cluster_names <- c("1" = "Cluster 1", "2" = "Cluster 2", "3" = "Cluster 3", 
                   "4" = "Cluster 4", "5" = "Cluster 5")

facet_colors_only <- as.vector(sapply(levels(factor(all_long_recl$Cluster)), get_cluster_color))


ftPlot <- all_long_recl %>% 
  ggplot(aes(x = Value, color = Variable, fill = Variable)) +
  geom_density(position = "identity", alpha = 0.5) +
  facet_wrap2(~ Cluster, strip = strip_themed(background_x = elem_list_rect(fill=facet_colors_only)), labeller = as_labeller(cluster_names), scales = "free", nrow = 3, ncol = 2) +
  labs(title = "Drivers of fire types",
       x = "Normalized size of burned area",
       y = "Density") +
  theme_light() + 
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    strip.background = element_rect(color = "darkgrey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=12, vjust=1.5),
    axis.text.y = element_text(size=12),
    axis.text.x = element_text(size=12),
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank())


## (4) MErge them all together

clusterPlot_all <- ggarrange(mapPlot, ftPlot, ncol = 2, labels = c("A", "B"),
  common.legend = FALSE, legend = "none", widths = c(3, 2))

ggsave(paste0(workDir, "FullPlot_v1_20241218.png"), clusterPlot_all, dpi=600, width = 15, height = 10, units = "in") 


