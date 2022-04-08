######################################
### Combine resistance prevalence data
### for Jones et al., 2022
### Lou LaMartina, started Apr 1, 2022
######################################


setwd("~/Desktop/Lab/Projects/Others/Aeromonas/FINAL")

library(ggplot2)
library(RColorBrewer)
library(readxl)


#######################
### read excel file ###
#######################

# load into list of data frames
data.ls <- list()
for (i in excel_sheets("./RData/compiled_resistance_data.xlsx")){
  
  # read in sheets of excel file
  data.ls[[i]] <- read_excel("./RData/compiled_resistance_data.xlsx", i)
  
}


# combine
data <- do.call(rbind, data.ls)


# remove special characters
colnames(data) <- gsub("-", "_", gsub("#", "num", gsub(" ", "_", colnames(data))))




##########################
### unique study names ###
##########################

# new variable
data$orig <- paste0(data$Study_name, "__", data$Year_study)


# make everything lowercase to avoid redundancies,
# except authors & countries
for(i in colnames(data)[-c(1,10)]) {
  data[[i]] <- tolower(data[[i]])
}


# simplify antibiotic names
data$ab <- gsub("-", "_", gsub(" ", "_", gsub("\\(", "", gsub(")", "", data$antibiotic))))


# simplify and sources
data$major_source <- gsub(" ", "", data$major_source)


# and wastewater
data$sub_Wwtri <- gsub("-", "", data$sub_Wwtri)
data$sub_WW <- gsub("-", "", data$sub_WW)


# new study names
studies <- unique(data[c("Study_name", "Year_study", "orig")])
studies$new <- gsub(" ", "", gsub("\\.", "", gsub("et al", "__", studies$Study_name)))


# which have years?
studies$nameyear <- sapply(strsplit(studies$new, "__"), '[', 2)


# for those that didn't, was it in the years column?
studies$nameyear[is.na(studies$nameyear)] <- studies$Year_study[is.na(studies$nameyear)]


# cool! get author name
studies$author <- sapply(strsplit(studies$new, "__"), '[', 1)


# not sure how R deals with special letters so i'm going to change that too
sort(unique(unlist(strsplit(studies$author, ""))))
# "á", "ã", "é", "í", "ñ", "ú"
studies$author <- gsub("-", "",
                       gsub("ú", "u",
                       gsub("ñ", "n",
                       gsub("í", "i", 
                       gsub("é", "e", 
                       gsub("ã", "a", 
                       gsub("á", "a", studies$author)))))))


# combine author and year
studies$newName <- paste0(studies$author, "__", studies$nameyear)


# combine with data
data <- merge(studies[c("orig", "newName")], data, by = "orig")


# new variable
data$uniq <- paste0(data$newName, "__", data$ab, "__", data$major_source)


# get replicates
gooddata1 <- subset(data, uniq %in% names(which(table(data$uniq) == 1)))
gooddata1$rep <- FALSE


# what major sources are duplicated?
unique(subset(data, uniq %in% names(which(table(data$uniq) > 1)))$major_source)
# [1] "seafood"       "surfacewater" "wastewater"    "clinical" 




###################
### dereplicate ###
###################
# need unique names for study+year+source, using other info

###########
### seafood

# subset
rep.sea <- subset(data, uniq %in% names(which(table(data$uniq) > 1)) & major_source == "seafood")


# unique var with fish type
rep.sea$uniq <- paste0(rep.sea$newName, "__", rep.sea$ab, "__", rep.sea$sub_fishfarm, rep.sea$sub_fish)


# add number - nothing else is different
for (i in names(which(table(rep.sea$uniq) > 1))) {
  for (j in 1:nrow(rep.sea[rep.sea$uniq == i,])) {
    rep.sea$uniq2[rep.sea$uniq == i][j] <- paste0(rep.sea$uniq[rep.sea$uniq == i][j], "__", j)
  }
}


# indicate replicate
rep.sea$rep <- FALSE
rep.sea$rep[is.na(rep.sea$uniq2) == F] <- TRUE


# replace
rep.sea$uniq[is.na(rep.sea$uniq2) == F] <- rep.sea$uniq2[is.na(rep.sea$uniq2) == F]


# add
gooddata2 <- rbind(gooddata1, rep.sea[-25])



#################
### surface water

# subset
rep.surf <- subset(data, uniq %in% names(which(table(data$uniq) > 1)) & major_source == "surfacewater")


# unique var with water type
rep.surf$uniq <- paste0(rep.surf$newName, "__", rep.surf$ab, "__", rep.surf$sub_surface)


# add number - nothing else is different
for (i in names(which(table(rep.surf$uniq) > 1))) {
  for (j in 1:nrow(rep.surf[rep.surf$uniq == i,])) {
    rep.surf$uniq2[rep.surf$uniq == i][j] <- paste0(rep.surf$uniq[rep.surf$uniq == i][j], "__", j)
  }
}


# indicate replicate
rep.surf$rep <- FALSE
rep.surf$rep[is.na(rep.surf$uniq2) == F] <- TRUE


# replace
rep.surf$uniq[is.na(rep.surf$uniq2) == F] <- rep.surf$uniq2[is.na(rep.surf$uniq2) == F]


# add
gooddata3 <- rbind(gooddata2, rep.surf[-25])



##############
### wastewater

# subset
rep.ww <- subset(data, uniq %in% names(which(table(data$uniq) > 1)) & major_source == "wastewater")


# unique var with water type
rep.ww$uniq <- paste0(rep.ww$newName, "__", rep.ww$ab, "__", rep.ww$sub_WW, "_", rep.ww$sub_Wwtri)


# add number - nothing else is different
for (i in names(which(table(rep.ww$uniq) > 1))) {
  for (j in 1:nrow(rep.ww[rep.ww$uniq == i,])) {
    rep.ww$uniq2[rep.ww$uniq == i][j] <- paste0(rep.ww$uniq[rep.ww$uniq == i][j], "__", j)
  }
}


# indicate replicate
rep.ww$rep <- FALSE
rep.ww$rep[is.na(rep.ww$uniq2) == F] <- TRUE


# replace
rep.ww$uniq[is.na(rep.ww$uniq2) == F] <- rep.ww$uniq2[is.na(rep.ww$uniq2) == F]


# add
gooddata4 <- rbind(gooddata3, rep.ww[-25])



############
### clinical

# subset
rep.clin <- subset(data, uniq %in% names(which(table(data$uniq) > 1)) & major_source == "clinical")


# add number - nothing else is different
for (i in names(which(table(rep.clin$uniq) > 1))) {
  for (j in 1:nrow(rep.clin[rep.clin$uniq == i,])) {
    rep.clin$uniq2[rep.clin$uniq == i][j] <- paste0(rep.clin$uniq[rep.clin$uniq == i][j], "__", j)
  }
}


# indicate replicate
rep.clin$rep <- FALSE
rep.clin$rep[is.na(rep.clin$uniq2) == F] <- TRUE


# replace
rep.clin$uniq[is.na(rep.clin$uniq2) == F] <- rep.clin$uniq2[is.na(rep.clin$uniq2) == F]


# add
gooddata5 <- rbind(gooddata4, rep.clin[-25])



############
### organize

# save new data frame
data_new <- gooddata5[c(24, 12, 23, 8:11, 7)]
data_new$resistance <- as.numeric(data_new$resistance)


# fix typos
data_new$ab[data_new$ab == "ceftiraxone"] <- "ceftriaxone"
data_new$ab[data_new$ab == "erythormycin"] <- "erythromycin"
data_new$ab[data_new$ab == "sulfamethoxazone_trimethoprim_sxt"] <- "sulfamethoxazole_trimethoprim_sxt"
data_new$ab[data_new$ab == "imipenem_also_meropenem"] <- "imipenem"




#################
### normalize ###
#################

# 1. get highest resistance value for each antibiotic
multiplier_vals <- list()
for (i in unique(data_new$ab)) {
  multiplier_vals[[i]] <- max(data_new$resistance[data_new$ab == i], na.rm = TRUE)
}
multiplier_vals <- data.frame(ab = names(multiplier_vals), top_resistance = unlist(multiplier_vals))


# 2. divide 100 by that resistance value
multiplier_vals$multiplier <- 100 / multiplier_vals$top_resistance
write.csv(multiplier_vals, "./RData/multiplier_calcs.csv", row.names = F)


# 3. for each antibiotic, multiply resistance values by that ratio
norm_vals <- list()
for (i in unique(data_new$ab)) {
  temp <- data_new[data_new$ab == i,]
  temp$resistance[temp$ab == i] <- temp$resistance[temp$ab == i] * 
    multiplier_vals$multiplier[multiplier_vals$ab == i]
  norm_vals[[i]] <- temp
}
norm_vals <- do.call(rbind, norm_vals)


# add to data frame
colnames(norm_vals)[8] <- "norm_resistance"
data_new <- merge(data_new, norm_vals[c(1,8)], by = "uniq")




############
### save ###
############

write.csv(data_new, "./RData/02_combined_resistance_data.csv", row.names = F, na = "")

