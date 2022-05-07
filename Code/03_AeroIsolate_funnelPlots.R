###################################
### Funnel plots with trim and fill
### for Jones et al., 2022
### Lou LaMartina, May 6, 2022
###################################


setwd("~/Desktop/Lab/Projects/Others/Aeromonas")
library(metafor)


# load data - see 01_resistanceData_organize.R
all_data <- read.csv("./RData/01_resistance_all_data.csv")



##################
### wastewater ###
##################

# subset relevant groups
#ww_data <- all_data[all_data$sub_WW3 %in% c("Untreated", "Treated"),]
ww_data <- all_data[all_data$major_source == "wastewater",]


# make variable to test
ww_data$group <- paste0(ww_data$antibioic_abb, "__", tolower(ww_data$major_source))


# subcategories
#ww_data$sub_category <- ww_data$sub_WW3
ww_data$sub_category <- ww_data$sub_WW2


# keep relevant data
ww_data <- ww_data[c("group", "sub_category", "study_label", "study", 
                     "author_year", "num_resistant", "num_total", "isolate_label")]


# split into individual data frames
ww_groups.ls <- list()
for (i in unique(ww_data$group)) {
  ww_groups.ls[[i]] <- ww_data[ww_data$group == i,]
  ww_groups.ls[[i]]$order[ww_groups.ls[[i]]$sub_category == "Untreated"] <- 1
  ww_groups.ls[[i]]$order[ww_groups.ls[[i]]$sub_category == "Treated"] <- 2
}



########
### plot

trimfill_L0_left <- list()
trimfill_L0_right <- list()
trimfill_R0_left <- list()
trimfill_R0_right <- list()


# loop to make plots, add text, and save
for (i in unique(ww_data$group)) {

  cat("\n", i, "\n")
  
  # subset data frame
  temp.df <- ww_groups.ls[[i]]
  
  
  # order by subcategory
  temp.df <- temp.df[order(temp.df$sub_category),]
  
  
  # calculate pooled prevalence & variance (see Note 1 below)
  temp.es <- escalc(measure = "PR", xi = num_resistant, ni = num_total, 
                    data = temp.df,
                    slab = study_label)
  
  
  # perform random effects model - all data
  temp.rma <- rma(yi, vi, data = temp.es)

  
  # L0 + left side
  temp.tf.L0.left <- trimfill(temp.rma, estimator = "L0", side = "left")
  pdf(paste0("./Plots/Funnels/Wastewater/L0_left/", i, "_funnel_L0_left.pdf"), width = 10, height = 10)
  funnel(temp.tf.L0.left, legend = FALSE, label = TRUE)
  trimfill_L0_left[[i]] <- data.frame(study = temp.tf.L0.left$slab, poolPrev = temp.tf.L0.left$yi,
                                      variance = temp.tf.L0.left$vi, fill = temp.tf.L0.left$fill,
                                      group = i, stat = "L0, left side")
  dev.off()
  
  # L0 + right side
  temp.tf.L0.right <- trimfill(temp.rma, estimator = "L0", side = "right")
  pdf(paste0("./Plots/Funnels/Wastewater/L0_right/", i, "_funnel_L0_right.pdf"), width = 10, height = 10)
  funnel(temp.tf.L0.right, legend = FALSE, label = TRUE)
  trimfill_L0_right[[i]] <- data.frame(study = temp.tf.L0.left$slab, poolPrev = temp.tf.L0.left$yi,
                                      variance = temp.tf.L0.left$vi, fill = temp.tf.L0.left$fill,
                                      group = i, stat = "L0, right side")
  dev.off()

  # R0 + left side
  temp.tf.R0.left <- trimfill(temp.rma, estimator = "R0", side = "left")
  pdf(paste0("./Plots/Funnels/Wastewater/R0_left/", i, "_funnel_R0_left.pdf"), width = 10, height = 10)
  funnel(temp.tf.R0.left, legend = FALSE, label = TRUE)
  trimfill_R0_left[[i]] <- data.frame(study = temp.tf.L0.left$slab, poolPrev = temp.tf.L0.left$yi,
                                      variance = temp.tf.L0.left$vi, fill = temp.tf.L0.left$fill,
                                      group = i, stat = "R0, left side")
  dev.off()

  # R0 + right side
  temp.tf.R0.right <- trimfill(temp.rma, estimator = "R0", side = "right")
  pdf(paste0("./Plots/Funnels/Wastewater/R0_right/", i, "_funnel_R0_right.pdf"), width = 10, height = 10)
  funnel(temp.tf.R0.right, legend = FALSE, label = TRUE)
  trimfill_R0_right[[i]] <- data.frame(study = temp.tf.L0.left$slab, poolPrev = temp.tf.L0.left$yi,
                                      variance = temp.tf.L0.left$vi, fill = temp.tf.L0.left$fill,
                                      group = i, stat = "R0, right side")
  dev.off()
}


# save results to data frame
wastewater_fill <- rbind(do.call(rbind, trimfill_L0_left), do.call(rbind, trimfill_L0_right),
                         do.call(rbind, trimfill_R0_left), do.call(rbind, trimfill_R0_right))

# save
write.csv(wastewater_fill, na = "", row.names = F, "./Plots/Funnels/Wastewater/wastewater_fill_results.csv")




###############
### seafood ###
###############

# subset relevant groups
fish_data <- all_data[all_data$sub_aquaculture %in% c("Farmed", "Wild"),]


# make variable to test
fish_data$group <- paste0(fish_data$antibioic_abb, "__", tolower(fish_data$major_source))


# subcategories
fish_data$sub_category <- fish_data$sub_aquaculture


# keep relevant data
fish_data <- fish_data[c("group", "sub_category", "study_label", "study", 
                         "author_year", "num_resistant", "num_total", "isolate_label")]


# split into individual data frames
fish_groups.ls <- list()
for (i in unique(fish_data$group)) {
  fish_groups.ls[[i]] <- fish_data[fish_data$group == i,]
  fish_groups.ls[[i]]$order[fish_groups.ls[[i]]$sub_category == "Farmed"] <- 1
  fish_groups.ls[[i]]$order[fish_groups.ls[[i]]$sub_category == "Wild"] <- 2
}



########
### plot

trimfill_L0_left <- list()
trimfill_L0_right <- list()
trimfill_R0_left <- list()
trimfill_R0_right <- list()


# loop to make plots, add text, and save
for (i in unique(fish_data$group)) {
  
  cat("\n", i, "\n")
  
  # subset data frame
  temp.df <- fish_groups.ls[[i]]
  
  
  # order by subcategory
  temp.df <- temp.df[order(temp.df$sub_category),]
  
  
  # calculate pooled prevalence & variance (see Note 1 below)
  temp.es <- escalc(measure = "PR", xi = num_resistant, ni = num_total, 
                    data = temp.df,
                    slab = study_label)
  
  
  # perform random effects model - all data
  temp.rma <- rma(yi, vi, data = temp.es)
  
  
  # L0 + left side
  temp.tf.L0.left <- trimfill(temp.rma, estimator = "L0", side = "left")
  pdf(paste0("./Plots/Funnels/Seafood/L0_left/", i, "_funnel_L0_left.pdf"), width = 10, height = 10)
  funnel(temp.tf.L0.left, legend = FALSE, label = TRUE)
  trimfill_L0_left[[i]] <- data.frame(study = temp.tf.L0.left$slab, poolPrev = temp.tf.L0.left$yi,
                                      variance = temp.tf.L0.left$vi, fill = temp.tf.L0.left$fill,
                                      group = i, stat = "L0, left side")
  dev.off()
  
  # L0 + right side
  temp.tf.L0.right <- trimfill(temp.rma, estimator = "L0", side = "right")
  pdf(paste0("./Plots/Funnels/Seafood/L0_right/", i, "_funnel_L0_right.pdf"), width = 10, height = 10)
  funnel(temp.tf.L0.right, legend = FALSE, label = TRUE)
  trimfill_L0_right[[i]] <- data.frame(study = temp.tf.L0.left$slab, poolPrev = temp.tf.L0.left$yi,
                                       variance = temp.tf.L0.left$vi, fill = temp.tf.L0.left$fill,
                                       group = i, stat = "L0, right side")
  dev.off()
  
  # R0 + left side
  temp.tf.R0.left <- trimfill(temp.rma, estimator = "R0", side = "left")
  pdf(paste0("./Plots/Funnels/Seafood/R0_left/", i, "_funnel_R0_left.pdf"), width = 10, height = 10)
  funnel(temp.tf.R0.left, legend = FALSE, label = TRUE)
  trimfill_R0_left[[i]] <- data.frame(study = temp.tf.L0.left$slab, poolPrev = temp.tf.L0.left$yi,
                                      variance = temp.tf.L0.left$vi, fill = temp.tf.L0.left$fill,
                                      group = i, stat = "R0, left side")
  dev.off()
  
  # R0 + right side
  temp.tf.R0.right <- trimfill(temp.rma, estimator = "R0", side = "right")
  pdf(paste0("./Plots/Funnels/Seafood/R0_right/", i, "_funnel_R0_right.pdf"), width = 10, height = 10)
  funnel(temp.tf.R0.right, legend = FALSE, label = TRUE)
  trimfill_R0_right[[i]] <- data.frame(study = temp.tf.L0.left$slab, poolPrev = temp.tf.L0.left$yi,
                                       variance = temp.tf.L0.left$vi, fill = temp.tf.L0.left$fill,
                                       group = i, stat = "R0, right side")
  dev.off()
}


# save results to data frame
seafood_fill <- rbind(do.call(rbind, trimfill_L0_left), do.call(rbind, trimfill_L0_right),
                         do.call(rbind, trimfill_R0_left), do.call(rbind, trimfill_R0_right))


# save
write.csv(seafood_fill, na = "", row.names = F, "./Plots/Funnels/Seafood/seafood_fill_results.csv")

