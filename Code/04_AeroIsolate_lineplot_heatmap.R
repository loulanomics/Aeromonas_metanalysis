##############################
### Line plot and heatmap
### for Jones et al., 2022
### Lou LaMartina, May 6, 2022
##############################


setwd("~/Desktop/Lab/Projects/Others/Aeromonas")

library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)
library(reshape2)


# see 02_resistanceData_organize
resist_data <- read.csv("./RData/01_resistance_all_data.csv")



#####################
### error barplot ###
#####################

# make variable to test
resist_data$group <- paste0(resist_data$antibioic_abb, "__", tolower(resist_data$major_source))


# get confidence intervals for these groups... there must be a better way of doing this?
ci.ls <- list()
for (i in unique(resist_data$group)) {
  temp <- rma(yi, vi, data = escalc(measure = "PR", xi = num_resistant, ni = num_total,
                                    data = resist_data[resist_data$group == i,]))
  ci.ls[[i]] <- data.frame(upCI = temp$ci.lb, loCI = temp$ci.ub, beta = temp$beta, group = i)
}
ci.df <- do.call(rbind, ci.ls)


# add to data
resist_data <- merge(resist_data, ci.df, by = "group")


# reorder for plot
resist_data$major_source[resist_data$major_source == "clinical"] <- "1 clinical"
resist_data$major_source[resist_data$major_source == "wastewater"] <- "2 wastewater"
resist_data$major_source[resist_data$major_source == "surfacewater"] <- "3 surfacewater"
resist_data$major_source[resist_data$major_source == "drinking"] <- "4 drinking"
resist_data$major_source[resist_data$major_source == "seafood"] <- "5 seafood"
resist_data$major_source[resist_data$major_source == "agriculture"] <- "6 agriculture"


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


# axis labels
unique(resist_data$antibiotic_category)
ab_labs <- c("carbapenem" = "Carbapenem",
             "Chloramphenicol" = "Chloramphenicol",
             "3rd gen" = "3rd gen.\nCephalosporin",
             "Ciprofloxacin" = "Ciprofloxacin",
             "Erythromycin" = "Erythromycin",
             "Gentamicin" = "Gentamicin",
             "SXT" = "Sulfameth-\noxazole-\nTrimethoprim",
             "tetracycline" = "Tetracycline")


# plot
ab.plot <-
  ggplot(resist_data, aes(x = antibiotic_category, y = upCI, 
                          color = major_source, fill = major_source)) +
  geom_hline(yintercept = c(0.25,0.5,0.75), size = 0.25, color = "grey80", linetype = "dotted") +
  geom_errorbar(aes(ymin = loCI, ymax = upCI), position = position_dodge(1), width = 0.5) +
  geom_point(aes(y = beta), position = position_dodge(1), size = 2) +
  
  facet_grid(~ antibiotic_category, scales = "free", space = "free") +
  scale_fill_manual(values = unname(sourcecols), labels = facets) +
  scale_color_manual(values = unname(sourcecols), labels = facets) +
  scale_x_discrete(labels = ab_labs) +
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
  guides(fill = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5),
         color = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5)) +
  labs(x = "Antibiotic", y = "Pooled prevalence", 
       color = "Source", fill = "Source")
ab.plot

ggsave("./Plots/ab_CI.pdf", plot = ab.plot, device = "pdf", width = 10, height = 4, units = "in")




###############
### heatmap ###
###############

# remove nonspecific countries
resist_data <- resist_data[-which(grepl(",", resist_data$country)),]
resist_data <- resist_data[-which(grepl("and", resist_data$country)),]
resist_data <- resist_data[resist_data$country != "Multiple",]


# new labels
abb_labs <- c("carbapenem" = "CAR",
             "Chloramphenicol" = "CHL",
             "3rd gen" = "CEF",
             "Ciprofloxacin" = "CIP",
             "Erythromycin" = "ERY",
             "Gentamicin" = "GEN",
             "SXT" = "SXT",
             "tetracycline" = "TET")


# plot
heat.plot <-
  ggplot(resist_data, aes(x = antibiotic_category, y = country, fill = troy_pooled_prev)) +
  geom_tile() +
  facet_grid(tri_income ~ tri_source, space = "free", scales = "free",
             labeller = labeller(tri_income = c("0-50" = "Low-mid economy", "50-75" = "Mid-hi economy", "75-100" = "Hi economy"))) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), limits = c(0,1)) +
  scale_x_discrete(labels = abb_labs) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10, color = "black", face = "bold"),
        strip.background = element_rect(colour = "grey90", fill = "grey90"),
        panel.border = element_rect(size = 0.75, color = "grey90", fill = NA),
        axis.line = element_line(color = "grey70", size = 0.25),
        axis.ticks = element_line(color = "grey70")) +
  guides(fill = guide_colorbar(barheight = 0.5, barwidth = 12,
                                               title.position = "top", title.hjust = 0.5)) +
  labs(x = "Antibiotic", y = "Country", fill = "Pooled prevalence")
heat.plot

ggsave("./Plots/heatmap.pdf", plot = heat.plot, device = "pdf", width = 10, height = 10, units = "in")

