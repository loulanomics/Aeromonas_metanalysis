#######################################
### Organize resistance prevalence data
### for Jones et al., 2022
### Lou LaMartina, May 4, 2022
#######################################


setwd("~/Desktop/Lab/Projects/Others/Aeromonas")


# load
all_data <- read.csv("./RData/Aero_resistance_allData.csv")


#############################
### combine author + year ###
#############################

# create original variable
all_data <- data.frame(orig_author_year = paste0(all_data$author, "__", all_data$year), all_data)


# new data frame
studies <- unique(all_data[c("orig_author_year", "author", "year")])
colnames(studies)[2:3] <- c("orig_author", "orig_year")


# simplify author names
studies$author <- gsub("[[:digit:]]", "", studies$orig_author)
studies$author <- gsub("et al", "", studies$author)
studies$author <- gsub("-", "", gsub(",", "", gsub(" ", "", gsub("\\.", "", 
                                                               studies$author))))


# extract years
studies$year <- gsub("[[:alpha:]]", "", studies$orig_author)
studies$year <- gsub("-", "", gsub(",", "", gsub(" ", "", gsub("\\.", "", 
                                                                studies$year))))


# for empties in year column, add year extracted from name
studies$year[studies$year == ""] <- studies$orig_year[studies$year == ""]


# not sure how R deals with special letters so i'm going to change that too
sort(unique(unlist(strsplit(studies$author, ""))))
studies$author <- gsub("-", "", gsub("ú", "u", gsub("ñ", "n", gsub("í", "i", 
                       gsub("é", "e", gsub("ã", "a", gsub("á", "a", studies$author)))))))


# lowercase, just to prevent redundancies from typos
studies$author <- tolower(studies$author)


# combine author and year
studies <- data.frame(author_year = paste0(studies$author, "__", studies$year), studies)


# add to data
all_data <- merge(studies, all_data[! colnames(all_data) %in% c("author", "year")], by = "orig_author_year")




###############################
### make unique study names ###
###############################

# simplify antibiotic names
all_data$antibioic_abb[all_data$antibioic_abb == "3rd CEF"] <- "CEF_3rd"
all_data$antibiotic <- tolower(gsub("-", "", gsub(" ", "", gsub("\\(", "", gsub(")", "", all_data$antibiotic)))))


# simplify sources
all_data$major_source <- tolower(gsub(" ", "", all_data$major_source))


# add those to new study variable
all_data <- data.frame(study = paste0(all_data$author_year, "__",
                                      all_data$antibiotic, "__", all_data$major_source),
                       all_data)


# get replicates
gooddata1 <- all_data[all_data$study %in% names(which(table(all_data$study) == 1)),]
rep_data <- all_data[! all_data$study %in% gooddata1$study,]


# what major sources are duplicated?
unique(rep_data$major_source)
# [1] "seafood"    "wastewater" "clinical"  



###########
### seafood

# subset
rep.sea <- subset(rep_data, major_source == "seafood")


# unique name var with fish type
rep.sea$study <- paste0(rep.sea$study, "__", tolower(rep.sea$sub_aquaculture), 
                        "_", tolower(rep.sea$sub_seafood))


# still replicates?
names(which(table(rep.sea$study) > 1))


# add number - nothing else is different
rep.sea$study2 <- rep.sea$study
for (i in names(which(table(rep.sea$study) > 1))) {
  for (j in 1:nrow(rep.sea[rep.sea$study == i,])) {
    rep.sea$study2[rep.sea$study == i][j] <- paste0(rep.sea$study[rep.sea$study == i][j], "__", j)
  }
}


# still replicates?
names(which(table(rep.sea$study2) > 1))


# replace
rep.sea$study <- rep.sea$study2


# add
gooddata2 <- rbind(gooddata1, rep.sea[-which(colnames(rep.sea) == "study2")])



##############
### wastewater

# subset
rep.ww <- subset(rep_data, major_source == "wastewater")


# unique name var with water type
rep.ww$study <- gsub("-", "", paste0(rep.ww$study, "__", tolower(rep.ww$sub_WW), "_", tolower(rep.ww$sub_WW3)))


# still replicates?
names(which(table(rep.ww$study) > 1))


# add
gooddata3 <- rbind(gooddata2, rep.ww)



############
### clinical

# subset
rep.clin <- subset(rep_data, major_source == "clinical")


# add number - nothing else is different
rep.clin$study <- paste0(rep.clin$study, "__", 1:2)


# add
new_data <- rbind(gooddata3, rep.clin)
names(which(table(new_data$study) > 1))




############
### save ###
############

# study labels
new_data$study_label <- gsub("[[:digit:]]", "", new_data$orig_author)
new_data$study_label <- gsub(" et al. ", "", new_data$study_label)
new_data$study_label <- gsub(" et al.", "", new_data$study_label)
new_data$study_label <- gsub(" et al", "", new_data$study_label)
new_data$study_label <- gsub(", ", "", new_data$study_label)
new_data$study_label <- paste0(new_data$study_label, " et al., ", new_data$year)


# isolate labels
new_data$isolate_label <- NA
new_data$isolate_label[new_data$major_source == "wastewater"] <- "Wastewater isolates resistant to "
new_data$isolate_label[new_data$major_source == "clinical"] <- "Clinical isolates resistant to "
new_data$isolate_label[new_data$major_source == "drinking"] <- "Drinking water isolates resistant to "
new_data$isolate_label[new_data$major_source == "agriculture"] <- "Agricultural isolates resistant to "
new_data$isolate_label[new_data$major_source == "surfacewater"] <- "Surface water isolates resistant to "
new_data$isolate_label[new_data$major_source == "seafood"] <- "Seafood isolates resistant to "
new_data$isolate_label <- paste0(new_data$isolate_label, new_data$antibiotic_category)

new_data$isolate_label <- gsub("SXT", "Sulfamethoxazole-trimethoprim", new_data$isolate_label)
new_data$isolate_label <- gsub("3rd gen", "Third-generation cephalosporins", new_data$isolate_label)
new_data$isolate_label <- gsub("tetracycline", "Tetracycline", new_data$isolate_label)
new_data$isolate_label <- gsub("carbapenem", "Carbapenem", new_data$isolate_label)


# reorder and save
new_data <- new_data[c(1, 49, 2:48, 50)]
write.csv(new_data, "./RData/01_resistance_all_data.csv", row.names = F, na = "")

