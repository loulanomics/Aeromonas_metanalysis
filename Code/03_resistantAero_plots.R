#######################################
### Plots of global resistant Aeromonas
### for Jones et al., 2022
### Lou LaMartina, started Apr 1, 2022
#######################################


setwd("~/Desktop/Lab/Projects/Others/Aeromonas/Maps")

library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)
library(reshape2)


# see 02_resistanceData_organize
resist_data <- read.csv("./RData/02_combined_resistance_data.csv")




######################
### global heatmap ###
######################

# load lat/lon data
latlon_data <- map_data("world")


# match countries
resist_data$country[resist_data$country == "United States"] <- "USA"
resist_data$country[resist_data$country == "England"] <- "UK"
resist_data$country[resist_data$country == "Korea"] <- "South Korea"
resist_data$country[resist_data$country == "Trinidad and Tobago"] <- "Trinidad"


# remove groups
globe_data <- resist_data[-grep(",", resist_data$country),]
globe_data <- globe_data[-grep(" and ", globe_data$country),]


# means of each country
stats <- cbind(aggregate(resistance ~ country, mean, data = globe_data),
               aggregate(resistance ~ country, sd, data = globe_data))
colnames(stats)[c(2,4)] <- c("mean", "sd")


# add latlon
globe_data <- dplyr::left_join(latlon_data, stats[-3], by = c("region" = "country")) 


# remove antarctica
globe_data <- globe_data[globe_data$region != "Antarctica",]


# specify colors
heatfill <- c(rev(brewer.pal(11, "RdYlGn")), "black")


# sequence of numbers with mean in middle
heatvals <- c(seq(min(globe_data$mean, na.rm = T), 
                  mean(globe_data$mean, na.rm = T), length.out = 6)[-6], 
              mean(globe_data$mean, na.rm = T), 
              seq(mean(globe_data$mean, na.rm = T), 1, length.out = 7)[-1])

# plot
globe.plot <-
  ggplot(globe_data) +
  
  geom_polygon(aes(x = long, y = lat, group = group, fill = mean, color = mean),
               size = 0.25) +
  
  scale_fill_gradientn(na.value = "grey90", limits = c(0,1),
                       colors = heatfill, values = heatvals) +
  scale_color_gradientn(na.value = "grey90", limits = c(0,1),
                        colors = heatfill, values = heatvals) +
  
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10)) +
  
  guides(color = "none", fill = guide_colorbar(title.vjust = 2)) +
  labs(fill = "Mean\nresistance\nprevalence")
globe.plot

ggsave("./Plots/globe.pdf", plot = globe.plot, device = "pdf", width = 10, height = 5, units = "in")




#########################
### category boxplots ###
#########################

#############
### economics

# relevant data
heat_data <- resist_data


# reorder for plot
heat_data$major_source[heat_data$major_source == "clinical"] <- "1 clinical"
heat_data$major_source[heat_data$major_source == "wastewater"] <- "2 wastewater"
heat_data$major_source[heat_data$major_source == "surfacewater"] <- "3 surfacewater"
heat_data$major_source[heat_data$major_source == "drinking"] <- "4 drinking"
heat_data$major_source[heat_data$major_source == "seafood"] <- "5 seafood"
heat_data$major_source[heat_data$major_source == "agriculture"] <- "6 agriculture"


# facet labels
facets <- c("6 agriculture" = "Agriculture",
            "1 clinical" = "Clinical",
            "4 drinking" = "Drinking water",
            "5 seafood" = "Seafood",
            "3 surfacewater" = "Surface water",
            "2 wastewater" = "Wastewater")


# colors
sourcecols <- c("#08519C", "#4EB3D3", "#A50F15", "#41AB5D", "#6A51A3", "#F16913")
names(sourcecols) <- c("blue", "lightblue", "red", "green", "purple", "orange")
sourcecols <- sourcecols[c(5,6,2,1,3,4)]


# plot
ses.plot <- 
  ggplot(heat_data, aes(x = as.factor(tri_econ), y = resistance)) +
  
  geom_hline(yintercept = c(0.25,0.5,0.75), size = 0.25, color = "grey80", linetype = "dotted") +
  
  geom_boxplot(alpha = 0.25, width = 0.25, outlier.alpha = 0,
               position = position_dodge(preserve = "single", width = 0.75),
               aes(fill = major_source, color = major_source),
               show.legend = F) +
  
  geom_jitter(width = 0.1, size = 0.1, show.legend = F,
              aes(fill = major_source, color = major_source)) +
  
  facet_grid(~ major_source, labeller = labeller(major_source = facets)) +
  
  scale_fill_manual(values = unname(sourcecols)) +
  scale_color_manual(values = unname(sourcecols)) +
  
  scale_x_discrete(labels = c("L-LMI", "MHI", "HI")) +
  
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        strip.text = element_text(size = 10, color = "black", face = "bold"),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(size = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  
  guides(fill = "none") + 
  labs(x = "Economic rank", y = "Pooled resistance prevalence", 
       color = "Resistance\nsource")
ses.plot

ggsave("./Plots/ses_box.pdf", plot = ses.plot, device = "pdf", width = 10, height = 3.5, units = "in")



###############
### antibiotics

heat_data$ab[heat_data$ab == "sulfamethoxazole_trimethoprim_sxt"] <- "sulfamethoxazole\ntrimethoprim (sxt)"


# plot
ab.plot <-
  ggplot(heat_data, aes(x = ab, y = resistance, 
                      color = major_source, fill = major_source)) +
  
  geom_hline(yintercept = c(0.25,0.5,0.75), size = 0.25, color = "grey80", linetype = "dotted") +
  
  geom_boxplot(alpha = 0.25, width = 0.5, outlier.size = 0.1, outlier.alpha = 1,
               position = position_dodge(preserve = "single", width = 0.75)) +
  
  facet_grid(~ ab, scales = "free", space = "free") +
  
  scale_fill_manual(values = unname(sourcecols), labels = facets) +
  scale_color_manual(values = unname(sourcecols), labels = facets) +
  
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "top",
        strip.text = element_blank(),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(size = 0.75, color = "grey80", fill = NA),
        axis.ticks = element_line(size = 0.25),
        axis.line = element_line(size = 0.25)) +
  
  guides(fill = guide_legend(nrow = 1, title.position ="top", title.hjust = 0.5),
         color = guide_legend(nrow = 1, title.position ="top", title.hjust = 0.5)) +
  labs(x = "Antibiotic", y = "Pooled resistance prevalence", 
       color = "Resistance source", fill = "Resistance source")
ab.plot

ggsave("./Plots/ab_box.pdf", plot = ab.plot, device = "pdf", width = 10, height = 4, units = "in")

