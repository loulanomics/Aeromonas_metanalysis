#################################################
### Subset Aeromonas genus from total communities
### for Jones et al., 2022
### Lou LaMartina, finalized Feb 18, 2022
#################################################


setwd("~/Desktop/Lab/Projects/Others/Aeromonas/FINAL")

library(ggplot2)
library(RColorBrewer)


#################
### load data ###
#################

# # ASV counts
# counts_files.ls <- lapply(list.files("./RData/Counts/Total/", pattern = ".csv", full.names = TRUE), 
#                           function(i) read.csv(i))
# 
# 
# # ASV taxonomy
# taxa_files.ls <- lapply(list.files("./RData/Taxonomy/", pattern = ".csv", full.names = TRUE), 
#                         function(i) read.csv(i))


# shortcut!
all_files.ls <- readRDS("./RData/All_files_list.RData")


# subset
counts_files.ls <- all_files.ls[grep("counts", names(all_files.ls))]
taxa_files.ls <- all_files.ls[grep("tax", names(all_files.ls))]


# simplify dataset names
Datasets <- gsub("_counts", "", names(counts_files.ls))
names(taxa_files.ls) <- Datasets
names(counts_files.ls) <- Datasets


# counts files - change row names to sample names
for (i in Datasets) {
  rownames(counts_files.ls[[i]]) <- counts_files.ls[[i]]$Sample_name
  counts_files.ls[[i]] <- counts_files.ls[[i]][-1]
}


# convert to relative abundance
relabun_files.ls <- list()
for(i in 1:length(counts_files.ls)){
  relabun_files.ls[[i]] <- counts_files.ls[[i]] / rowSums(counts_files.ls[[i]])
}


# taxa files - change row names to FASTA sequences
for (i in Datasets) {
  rownames(taxa_files.ls[[i]]) <- taxa_files.ls[[i]]$FASTA
}


# load all info
info <- read.csv("./RData/AeromonasMetanalysis_16S_metadata.csv")




########################
### subset aeromonas ###
########################

aero_seqs.ls <- list()
aero_relabun.ls <- list()

for (i in 1:length(relabun_files.ls)) {
  aero_seqs.ls[[i]] <- subset(taxa_files.ls[[i]], Genus == "Aeromonas")$FASTA
  
  if (length(aero_seqs.ls[[i]]) > 1) {
    aero_relabun.ls[[i]] <- relabun_files.ls[[i]][, colnames(relabun_files.ls[[i]]) %in% aero_seqs.ls[[i]]]
    
  } else {
    aero_relabun.ls[[i]] <- data.frame(relabun_files.ls[[i]][, colnames(relabun_files.ls[[i]]) %in% aero_seqs.ls[[i]]])
    colnames(aero_relabun.ls[[i]]) <- aero_seqs.ls[[i]]
    rownames(aero_relabun.ls[[i]]) <- rownames(relabun_files.ls[[i]])
    
  }
}


# transpose so FASTA are row names, add source
aero_relabun.t.ls <- list()

for (i in 1:length(aero_relabun.ls)) {
  aero_relabun.t.ls[[i]] <- data.frame(t(aero_relabun.ls[[i]]))
  aero_relabun.t.ls[[i]]$FASTA <- rownames(aero_relabun.t.ls[[i]])
}


# merge by FASTAs
merged_relabun.ls <- list()

for (i in 1:length(aero_relabun.t.ls)) {
  if (i == 1) {
    merged_relabun.ls[[i]] <- merge(aero_relabun.t.ls[[i]], 
                                    aero_relabun.t.ls[[i + 1]], by = "FASTA", all = TRUE)
  } else if (i > 1 & i < length(aero_relabun.t.ls)) {
    merged_relabun.ls[[i]] <- merge(merged_relabun.ls[[i - 1]], 
                                    aero_relabun.t.ls[[i + 1]], by = "FASTA", all = TRUE)
  }
}


# extract the last one - merge() makes you do it one by one
merged_relabun <- data.frame(merged_relabun.ls[[length(aero_relabun.t.ls) - 1]])


# change "NA" to zero
merged_relabun[is.na(merged_relabun)] <- 0


# make FASTA row names, transpose
rownames(merged_relabun) <- merged_relabun$FASTA
merged_relabun <- data.frame(t(merged_relabun[-1]))
dim(merged_relabun) # [1] 333 152




########################
### plot proportions ###
########################

# calculate proportion aeromonas in each sample
propAero <- merged_relabun
propAero$totalAero <- rowSums(propAero)
propAero <- data.frame(File = rownames(propAero), totalAero = propAero$totalAero)


# add sample info
propAero <- merge(propAero, info, by = "File")


# calculate mean & std deviation
stats <- propAero[c(2,7,8)]
stats <- data.frame(mean = aggregate(totalAero ~  SampleSource + SampleType, mean, data = stats),
                    sd = aggregate(totalAero ~ SampleSource + SampleType, sd, data = stats))[-c(4,5)]


# most in each dataset
for(i in Datasets){
  top <- names(sort(rowSums(merged_relabun[rownames(merged_relabun) %in% rownames(counts_files.ls[[i]]),]), decreasing = T)[1])
  cat("\n", i)
  print(info[info$File == top,])
  print(rowSums(merged_relabun[rownames(merged_relabun) %in% top,]))
}


# order facets
facets <- sort(unique(propAero$SampleSource))
names(facets) <- c(3, 5, 1, 4, 2, 6)
facets <- data.frame(facets = paste(names(facets), facets), SampleSource = facets)
propAero <- merge(propAero, facets, by = "SampleSource")


# define colors
colors <- c("#08519C", "#4EB3D3", "#A50F15", "#41AB5D", "#6A51A3", "#F16913")

facets$facets

# facet labels
facslabs <- c("3 Aquaculture" = "Aquaculture",
              "5 Beach" = "Beaches",
              "1 Gut" = "Fish intenstines",         
              "4 River" = "Rivers",       
              "2 Skin" = "Fish exterior",        
              "6 Wastewater" = "Wastewater")


# x axis labels
xaxis <- data.frame(SampleType = unique(propAero$SampleType))
xaxis$lab <- gsub(" ", "\n", xaxis$SampleType)
xaxis$lab <- gsub("-chlorinated", "-\nchlor.", xaxis$lab)
propAero <- merge(propAero, xaxis, by = "SampleType")



# plot
dot.plot <-
  ggplot(propAero, aes(x = lab, y = totalAero)) +
  
  geom_hline(yintercept = c(0.01, 0.1), size = 0.25, color = "grey80", linetype = "dotted") +
  
  geom_boxplot(alpha = 0.25, width = 0.25, outlier.alpha = 0,
               position = position_dodge(preserve = "single", width = 0.75),
               aes(fill = SampleSource, color = SampleSource),
               show.legend = F) +
  
  geom_jitter(width = 0.1, size = 0.1, show.legend = F,
              aes(fill = SampleSource, color = SampleSource)) +
  
  facet_grid(. ~ facets, scales = "free", space = "free",
             labeller = labeller(facets = facslabs)) +
  
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  
  scale_y_continuous(trans = scales::pseudo_log_trans(0.001, 10),
                     breaks = c(0, 0.01, 0.1, 1), 
                     labels = c("0%", "1%", "10%", "100%")) +
  
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
  labs(y = "Percent Aeromonas in samples", x = "Sample source") 
dot.plot

ggsave("./Plots/source_box.pdf", plot = dot.plot, device = "pdf", width = 10, height = 4, units = "in")


