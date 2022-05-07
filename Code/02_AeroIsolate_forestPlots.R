#######################################
### Forest plots for Jones et al., 2022
### Lou LaMartina, May 6, 2022
#######################################


setwd("~/Desktop/Lab/Projects/Others/Aeromonas")
library(metafor)


# load data - see 01_resistanceData_organize.R
all_data <- read.csv("./RData/01_resistance_all_data.csv")


# function to extract Q, I^2, tau^2 (see Note 1 at end of script)
mlabfun <- function(text, res) {
  list(bquote(paste(.(text),
                    " (Q = ", .(formatC(res$QE, digits = 2, format = "f")),
                    ", df = ", .(res$k - res$p),
                    ", p ", .(metafor:::.pval(res$QEp, digits = 2, showeq = TRUE, sep = " ")), "; ",
                    I^2, " = ", .(formatC(res$I2, digits = 1, format = "f")), "%, ",
                    tau^2, " = ", .(formatC(res$tau2, digits = 2, format = "f")), ")")))
}


# one for text outside mlab
labfun <- function(res) {
  bquote(italic(paste("Q = ", .(formatC(res$QE, digits = 2, format = "f")),
                      ", df = ", .(res$k - res$p),
                      ", p ", .(metafor:::.pval(res$QEp, digits = 2, showeq = TRUE, sep = " ")), "; ",
                      I^2, " = ", .(formatC(res$I2, digits = 1, format = "f")), "%, ",
                      tau^2, " = ", .(formatC(res$tau2, digits = 2, format = "f")))))
}



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

# loop to make plots, add text, and save
for (i in unique(ww_data$group)) {
  
  cat("\n", i, "\n")
  
  # subset data frame
  temp.df <- ww_groups.ls[[1]]
  
  
  # order by subcategory
  temp.df <- temp.df[order(temp.df$sub_category),]
  
  
  # calculate pooled prevalence & variance (see Note 1 below)
  temp.es <- escalc(measure = "PR", xi = num_resistant, ni = num_total, 
                    data = temp.df,
                    slab = study_label)
  
  
  # perform random effects model - all data
  temp.rma <- rma(yi, vi, data = temp.es)
  
  
  # fit meta regression model to test for subgroup differences
  sub.rma <- rma(yi, vi, mods = ~ sub_category, data = temp.es)
  
  
  # perform random effects model - by subgroup
  temp.rma.1 <- rma(yi, vi, data = temp.es[temp.es$sub_category == "Treated",])
  temp.rma.2 <- rma(yi, vi, data = temp.es[temp.es$sub_category == "Untreated",])
  
  
  # number of resistant and total isolates
  resN = temp.df$num_resistant
  totN = temp.df$num_total
  
  
  # get lengths of subcategories; make sure k1 & k2 are in alphabetical order
  k1 = length(which(temp.df$sub_category == "Treated"))
  k2 = length(which(temp.df$sub_category == "Untreated"))
  
  
  # total number studies = total number of rows in plot
  k = k1 + k2
  
  
  # decide plot size based on number of samples (see Note 2)
  ylo = -3
  yhi = k + 11
  xlo = -1.7
  xhi = 1.8
  
  
  # save to pdfs
  pdf(paste0("./Plots/Forests/Wastewater", i, "_forest.pdf"), width = 11, height = (k / 3) + 5)
  
  
    # plot (see Note 3)
    forest(temp.rma,
         slab = temp.df$study_label,
         xlim = c(xlo, xhi), ylim = c(ylo, yhi),
         ilab = cbind(resN, totN), ilab.xpos = c(-0.75, -0.5),
         rows = c((yhi - 4):(yhi - 3 - k1), (5 + k2 - 1):5),
         xlab = "Pooled prevalence", showweights = T, addpred = T, mlab = NULL)
  
  
    # add summary polygons for the three subgroups
    addpoly(temp.rma.1, row = (yhi - 3), mlab = "", cex = 1)
    addpoly(temp.rma.2, row = (k2 + 5), mlab = "", cex = 1)
  
  
    # subgroup headers
    text(xlo, (yhi - 3), "Treated", font = 4, pos = 4)
    text(xlo, (k2 + 5), "Untreated", font = 4, pos = 4)
  
  
    # subgroup stats
    text(xlo, (yhi - 4 - k1), labfun(temp.rma.1), pos = 4)
    text(xlo, 4, labfun(temp.rma.2), pos = 4)
  
  
    # subgroup diff stats
    text(xlo, 2, "Subgroup differences", font = 4, pos = 4)
    text(xlo, 1, pos = 4,
       bquote(italic(paste(Q[M], " = ", .(formatC(sub.rma$QM, digits = 2, format="f")), 
                           ", df = ", .(sub.rma$p - 1),
                           ", p = ", .(formatC(sub.rma$QMp, digits = 2, format = "f"))))))
  
  
    # stats summary - all
    text(xlo, -1, c("All studies"), font = 4, pos = 4)
    text(xlo, -2, labfun(temp.rma), pos = 4)
  
  
    # study header
    text(xlo, (yhi - 1), c("Author(s), year"), font = 2, pos = 4)
  
  
    # resistant columns headers
    text(-0.35, (yhi - 1), c("Resistant / Total"), font = 2, pos = 2)
  
  
    # stats header
    text(xhi, (yhi - 1), c("Weight / Pooled prevalence [95% CI]"), font = 2, pos = 2)
  
  
    # title
    text(0, yhi, unique(temp.df$isolate_label), cex = 1.25, font = 2)
  
  
    # close
  dev.off()
}




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
  
  
  # fit meta regression model to test for subgroup differences
  sub.rma <- rma(yi, vi, mods = ~ sub_category, data = temp.es)
  
  
  # perform random effects model - by subgroup
  temp.rma.1 <- rma(yi, vi, data = temp.es[temp.es$sub_category == "Farmed",])
  temp.rma.2 <- rma(yi, vi, data = temp.es[temp.es$sub_category == "Wild",])
  
  
  # number of resistant and total isolates
  resN = temp.df$num_resistant
  totN = temp.df$num_total
  
  
  # get lengths of subcategories; make sure k1 & k2 are in alphabetical order
  k1 = length(which(temp.df$sub_category == "Farmed"))
  k2 = length(which(temp.df$sub_category == "Wild"))
  
  
  # total number studies = total number of rows in plot
  k = k1 + k2
  
  
  # decide plot size based on number of samples (see Note 2)
  ylo = -3
  yhi = k + 11
  xlo = -1.7
  xhi = 1.8
  
  
  # save to pdfs
  pdf(paste0("./Plots/Forests/Seafood/", i, "_forest.pdf"), width = 11, height = (k / 3) + 5)
  
  
  # plot (see Note 3)
  forest(temp.rma,
         slab = temp.df$study_label,
         xlim = c(xlo, xhi), ylim = c(ylo, yhi),
         ilab = cbind(resN, totN), ilab.xpos = c(-0.75, -0.5),
         rows = c((yhi - 4):(yhi - 3 - k1), (5 + k2 - 1):5),
         xlab = "Pooled prevalence", showweights = T, addpred = T, mlab = NULL)
  
  
  # add summary polygons for the three subgroups
  addpoly(temp.rma.1, row = (yhi - 3), mlab = "", cex = 1)
  addpoly(temp.rma.2, row = (k2 + 5), mlab = "", cex = 1)
  
  
  # subgroup headers
  text(xlo, (yhi - 3), "Farmed", font = 4, pos = 4)
  text(xlo, (k2 + 5), "Wild", font = 4, pos = 4)
  
  
  # subgroup stats
  text(xlo, (yhi - 4 - k1), labfun(temp.rma.1), pos = 4)
  text(xlo, 4, labfun(temp.rma.2), pos = 4)
  
  
  # subgroup diff stats
  text(xlo, 2, "Subgroup differences", font = 4, pos = 4)
  text(xlo, 1, pos = 4,
       bquote(italic(paste(Q[M], " = ", .(formatC(sub.rma$QM, digits = 2, format="f")), 
                           ", df = ", .(sub.rma$p - 1),
                           ", p = ", .(formatC(sub.rma$QMp, digits = 2, format = "f"))))))
  
  
  # stats summary - all
  text(xlo, -1, c("All studies"), font = 4, pos = 4)
  text(xlo, -2, labfun(temp.rma), pos = 4)
  
  
  # study header
  text(xlo, (yhi - 1), c("Author(s), year"), font = 2, pos = 4)
  
  
  # resistant columns headers
  text(-0.35, (yhi - 1), c("Resistant / Total"), font = 2, pos = 2)
  
  
  # stats header
  text(xhi, (yhi - 1), c("Weight / Pooled prevalence [95% CI]"), font = 2, pos = 2)
  
  
  # title
  text(0, yhi, unique(temp.df$isolate_label), cex = 1.25, font = 2)
  
  
  # close
  dev.off()
}





#############
### NOTES ###
#############

### NOTE 1 ###

# helper function to add Q-test, I^2, tau^2
# https://www.metafor-project.org/doku.php/plots:forest_plot_with_subgroups

# measures for dichotomous variables
# input:
# - xi: number of individuals experiencing the event
# - ni: total number of individuals within each study
# - PR: raw proportion
# output:
# - yi: observed effect sizes or outcomes
# - vi: corresponding sampling variances



### NOTE 2 ###

# parameters help - see first figure here:
# https://wviechtb.github.io/metafor/reference/forest.rma.html


### Y AXIS POINTS EXAMPLE: SAMPLE SIZE = 7

# 18    = yhi = top of plot (top axis line + 2)
# 17    = study header
# 16    = axis line
# 15    = subgroup 1 header (yhi - 3)
# 14:10 = length of k1 = 5
# 9     = subgroup 1 stats
# 8     = space
# 7     = subgroup 2 header
# 6:5   = length of k2 = 2 ALWAYS STARTS ON 5
# 4     = subgroup 1 stats
# 3     = space
# 2     = group diff header
# 1     = group diffs
# 0     = axis line



### NOTE 3 ###

# forest(temp.rma,      model
#        slab =         study labels
#        xlim =         expand x axis of figure to fit resN and totN columns
#        ylim =         expand y axis to fit headers
#        showweights =  weights, to left of CIs
#        ilab =         add resN and totN columns
#        ilab.xpos =    position resN and totN columns
#        rows =         position studies - see Note 2
#        xlab =         x axis label
#        addpred =      dotted line on summary polygon = % of true outcomes expected to fall
#        mlab =         remove summary label, will add it manually

