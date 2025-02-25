
# LOAD DATA IN R ----------------------------------------------------------

#set working directory
setwd("N:/SUN-CPR-proteomics_jvo/GIULIA/B_COLLAB_EPIC-XS/UgoC/UgoC_OC_DIA_Prot_SN18")

#load OC proteome DIA data
library(data.table)
mydata <- fread("20231223_113739_20231025_GF_OC_Prot_dDIA_Report.tsv")
mydata <- as.data.frame(mydata)
colnames(mydata)
#6886 ids

newdata <- fread("20241118_135212_20231025_GF_OC_Prot_dDIA_Report.tsv")
newdata <- as.data.frame(newdata)
colnames(newdata)
#6886 ids

#merge data
library(tidyverse)
mydata <- left_join(
  mydata,
  newdata %>% select(
    PG.ProteinGroups,
    PG.Genes,
    PG.ProteinDescriptions,
    PG.ProteinNames,
    PG.FastaFiles,
    `PG.NrOfStrippedSequencesIdentified (Experiment-wide)`,
    `PG.NrOfModifiedSequencesIdentified (Experiment-wide)`,
    `PG.NrOfPrecursorsIdentified (Experiment-wide)`
  ),
  by = c("PG.ProteinGroups", "PG.Genes", "PG.ProteinDescriptions", "PG.FastaFiles")
)

#reorder columns
mydata <- mydata[, c(1:4, 45:48, 5:44)]

# Remove contaminants -----------------------------------------------------
unique(mydata$PG.FastaFiles)
mydata <- mydata[mydata$PG.FastaFiles == "Human_UP_2023-03",]
mydata <- mydata[,-4]
#6782 ids

# Remove rows wo gene name ------------------------------------------------
mydata <- mydata[mydata$PG.Genes != "",]
#6780 ids

# Check for negative values -----------------------------------------------
sum(mydata[,8:47] < 0)

# Rename sample names -----------------------------------------------------
colnames(mydata) <- sub("^.*?\\Prot_", "", colnames(mydata))
colnames(mydata) <- sub("F2_", "", colnames(mydata))
colnames(mydata) <- sub(".raw.PG.Log2Quantity", "", colnames(mydata))

write.table(mydata, "OC_Prot_cleaned.txt", row.names = F, sep = "\t")

mydata <- fread("OC_Prot_cleaned.txt")
mydata <- as.data.frame(mydata)

# id count ----------------------------------------------------------------
id <- as.data.frame(nrow(mydata) - colSums(is.na(mydata[, c(8:47)])))
colnames(id)[1] = "Number_of_IDs"
id$run <- rownames(id)
rownames(id) <- NULL
id$sample <- substr(id$run, 1, nchar(id$run) -3)

library(tidyverse)
Id <- id %>% group_by(sample) %>% summarise(
  mean = mean(Number_of_IDs, na.rm = T),
  st.dev = sd(Number_of_IDs, na.rm = T)
)
Id$mean <- as.numeric(format(round(Id$mean, 0), nsmall = 0))

Id$patient <- sub("\\_.*", "", Id$sample)
Id$patient <-
  factor(Id$patient,
         levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))
Id$group <- sub("^.*?\\_", "", Id$sample)

#average ids per group
Id %>% group_by(group) %>% summarise(mean = mean(mean, na.rm = T))

library(RColorBrewer)
ggplot(Id, aes(x = patient, y = mean, fill = group)) +
  geom_bar(position = position_dodge(),
           stat = "identity",
           colour = "black") +
  scale_fill_brewer(palette = "Set2") +
  geom_errorbar(
    aes(ymin = mean - st.dev, ymax = mean + st.dev),
    width = 0.2,
    position = position_dodge(.9)
  ) +
  geom_text(
    aes(label = mean),
    position = position_dodge(width = .9),
    vjust = -0.5,
    size = 2.5,
    color = "black"
  ) +
  xlab("") +
  ylab("PG") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("PG count - DIA proteome")

# Long format -------------------------------------------------------------
melt <- reshape2::melt(mydata[,c(1,8:47)], id = "PG.ProteinGroups")
melt <- melt[!is.na(melt$value),]
melt$variable <- as.character(melt$variable)

melt$r <- substr(melt$variable, nchar(melt$variable)-1, nchar(melt$variable))
unique(melt$r)

melt$patient <- sub("\\_.*", "", melt$variable)
unique(melt$patient)
melt$patient <-
  factor(melt$patient,
         levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))

melt$patientr <- paste0(melt$patient, "_", melt$r)
unique(melt$patientr)

melt$group <- sub("^.*?\\_", "", melt$variable)
melt$group <- sub("\\_.*", "", melt$group)
unique(melt$group)

#order by patient
melt <- melt[order(melt$patient),]
melt$patientr <- factor(melt$patientr, levels = unique(melt$patientr))

#boxplot
ggplot(melt, aes(x = patientr, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.2, lwd=0.2) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  )) +
  scale_fill_brewer(palette = "Set2") +
  xlab("MS run") +
  ylab("Log2 MS signal")

# normalization by median -------------------------------------------------
sampleMed <- apply(mydata[,8:47], 2, median, na.rm = TRUE) 
meanMed <- mean(sampleMed, na.rm = TRUE)
norm <- t(t(mydata[,8:47]) - sampleMed)
norm <- as.data.frame(norm + meanMed)
norm <- cbind(mydata[,1:7], norm)
colnames(norm)

write.table(norm, "OC_Prot_norm.txt", row.names = F, sep = "\t")

norm <- fread("OC_Prot_norm.txt")
norm <- as.data.frame(norm)

# Long format after normalization ------------------------------------------
melt.norm <- reshape2::melt(norm[,c(1, 8:47)], id = "PG.ProteinGroups")
melt.norm <- melt.norm[!is.na(melt.norm$value),]
melt.norm$variable <- as.character(melt.norm$variable)

melt.norm$r <- substr(melt.norm$variable, nchar(melt.norm$variable)-1, nchar(melt.norm$variable))
unique(melt.norm$r)

melt.norm$patient <- sub("\\_.*", "", melt.norm$variable)
unique(melt.norm$patient)
melt.norm$patient <-
  factor(melt.norm$patient,
         levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))

melt.norm$patientr <- paste0(melt.norm$patient, "_", melt.norm$r)
unique(melt.norm$patientr)

melt.norm$group <- sub("^.*?\\_", "", melt.norm$variable)
melt.norm$group <- sub("\\_.*", "", melt.norm$group)
unique(melt.norm$group)

#order by patient
melt.norm <- melt.norm[order(melt.norm$patient),]
melt.norm$patientr <- factor(melt.norm$patientr, levels = unique(melt.norm$patientr))

#boxplot
ggplot(melt.norm, aes(x = patientr, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.2, lwd=0.2) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  )) +
  scale_fill_brewer(palette = "Set2") +
  xlab("MS run") +
  ylab("Normalized Log2 MS signal")

# Average technical replicates --------------------------------------------
cat <- sub("_01.*", "", colnames(norm[,8:47]))
cat <- sub("_02.*", "", cat)
cat <- factor(cat)
average <-
  as.data.frame(t(apply(norm[, 8:47], 1, function(x) {
    tapply(x, cat, function (y)
      mean(y, na.rm = TRUE))
  })))
average <- cbind(norm[,1:7], average)
colnames(average)

#reorder columns
average <- average[, c(1:9,12:27,10:11)]

write.table(average, "OC_Prot_norm_average.txt", row.names = F, sep = "\t")

average <- fread("OC_Prot_norm_average.txt")
average <- as.data.frame(average)

# Long format after averaging ----------------------------------------------
melt.average <- reshape2::melt(average[,c(1:2,8:27)], id = c("PG.ProteinGroups", "PG.Genes"))
melt.average <- melt.average[!is.na(melt.average$value),]
melt.average$variable <- as.character(melt.average$variable)

melt.average$group <- substr(melt.average$variable, nchar(melt.average$variable)-2, nchar(melt.average$variable))
unique(melt.average$group)

melt.average$patient <- sub("\\_.*", "", melt.average$variable)
unique(melt.average$patient)
melt.average$patient <- factor(melt.average$patient, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))

#boxplot
ggplot(melt.average, aes(x = patient, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.2, lwd=0.2) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  )) +
  scale_fill_brewer(palette = "Set2") +
  xlab("MS run") +
  ylab("Normalized Log2 MS signal")

# PCA ---------------------------------------------------------------------

#define categories
patient <- sub("\\_.*", "", colnames(average[,8:27]))
unique(patient)
patient <-
  factor(patient,
         levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))

group <- sub("^.*?\\_", "", colnames(average[,8:27]))
unique(group)
group <- factor(group, levels = c("Adh", "Sph"))

#remove missing values
nona <- na.omit(average[,8:27])
#5536 proteins

#transpose data
nonat <- t(nona)

#perform PCA
data.pr <- prcomp(nonat, center = T, scale = F)

#plot PCA
library(ggbiplot)
library(ggrepel)

ggbiplot(
  data.pr,
  var.axes = F,
  obs.scale = 2,
  var.scale = 0,
  groups = group,
  ellipse = F
) +
  geom_point(aes(colour = group), size = 4) +
  scale_color_brewer(palette = "Set2") +
  geom_text_repel(aes(label = patient), size = 4, box.padding = 0.5, max.overlaps = Inf) +
  theme_light(base_size = 12)

# Presence/Absence --------------------------------------------------------
col_adh <- seq(8, ncol(average), by = 2)
col_sph <- seq(9, ncol(average), by = 2)

Sph.only <-
  average[rowSums(is.na(average[, col_adh])) == 10 &
            rowSums(is.na(average[, col_sph])) <= 5, ]
Sph.only$presence <- "Sph"

Adh.only <-
  average[rowSums(is.na(average[, col_adh])) <= 5 &
            rowSums(is.na(average[, col_sph])) == 10, ]
Adh.only$presence <- "Adh"

only.in.one.group <- rbind(Sph.only, Adh.only)

#calculate number of values per row (frequency)
only.in.one.group$frequency <- 
  apply(only.in.one.group[, c(8:27)], 1, function(x) {
    sum(!is.na(x))
  })

write.table(only.in.one.group, "only.in.one.group.txt", row.names = F, sep = "\t")

average <-
  merge(average, 
        only.in.one.group[, c(1, 28:29)], 
        by = "PG.ProteinGroups", 
        all.x = T)

write.table(average, "OC_Prot_norm_average_presence-absence.txt", row.names = F, sep = "\t")

average <- fread("OC_Prot_norm_average_presence-absence.txt")
average <- as.data.frame(average)

# Heatmap presence/absence ------------------------------------------------
rownames(only.in.one.group) <- only.in.one.group$PG.Genes
hm <- only.in.one.group[, c(col_adh, col_sph)]
breaks <- seq(min(hm, na.rm=T), max(hm, na.rm=T), length.out = 101)
hm[is.na(hm)] <- 0

#annotate columns
group <- sub(".*\\_", "", colnames(hm))
group <- as.data.frame(group)
rownames(group) <- colnames(hm)

ann_colors <- list(
  group = c(Sph = "#FC8D62", Adh = "#66C2A5")
)

library(pheatmap)
library(RColorBrewer)
pheatmap(hm, 
         color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
         annotation_col = group,
         cluster_cols = F,
         annotation_colors = ann_colors,
         breaks = breaks)

# Proliferation markers ---------------------------------------------------
pi <- fread("proliferation.index.txt")
pi <- as.data.frame(pi)
average$PG <- gsub(";.*", "", average$PG.ProteinGroups)
pi <- left_join(pi, average, by = c("Uniprot" = "PG"))
pindex <- apply(pi[,9:28], 2, median, na.rm = TRUE)

#subtract odd values (Sph) to even values (Adh)
pindex.fc <- pindex[seq(2, length(pindex), by = 2)] - pindex[seq(1, length(pindex), by = 2)]

#one sample t test
t.test(pindex.fc, mu = 0, alternative = "two.sided")
#p = 0.3904

#long format
pi.long <- reshape2::melt(pi[,c(1,9:28)], id = "Uniprot")
pi.long <- na.omit(pi.long)
pi.long$variable <- as.character(pi.long$variable)

pi.long$patient <- sub("\\_.*", "", pi.long$variable)
pi.long$patient <-
  factor(pi.long$patient,
         levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))

pi.long$group <- sub("^.*?\\_", "", pi.long$variable)

#boxplot
ggplot(pi.long, aes(x = patient, y = value, color = group)) +
  geom_boxplot(outlier.size = 0.2, lwd = 0.2, position = position_dodge(width = 0.9)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  )) +
  scale_color_brewer(palette = "Set2") +
  xlab("Patient") +
  ylab("Log2 MS signal") +
  ggtitle("Proliferation markers")

#One-way ANOVA
pi.long.adh <- pi.long[pi.long$group == "Adh",] 
pi.long.sph <- pi.long[pi.long$group == "Sph",]
mod <- aov(value ~ patient, data = pi.long.adh) #0.0452
summary(mod)
mod <- aov(value ~ patient, data = pi.long.sph) #NS
summary(mod)

# Significance analysis ---------------------------------------------------

#filter data
filtered <-
  average[rowSums(!is.na(average[, col_adh])) >= 3 &
            rowSums(!is.na(average[, col_sph])) >= 3, 1:27]
#6558 ids

write.table(filtered, "OC_Prot_norm_average_filtered.txt", row.names = F, sep = "\t")

filtered <- fread("OC_Prot_norm_average_filtered.txt")
filtered <- as.data.frame(filtered)

#perform Limma paired analysis
library(limma)

patient <- sub("\\_.*", "", colnames(filtered[,8:27]))
unique(patient)
patient <-
  factor(patient,
         levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))

group <- sub("^.*?\\_", "", colnames(filtered[,8:27]))
unique(group)
group <- factor(group, levels = c("Adh", "Sph"))

design <- model.matrix(~ patient + group)
colnames(design) <- gsub("patient", "", colnames(design))
colnames(design) <- gsub("group", "", colnames(design))
fit <- lmFit(filtered[,8:27], design)
fit2 <- eBayes(fit)
tab <- topTable(fit2,
                n = Inf,
                coef = "Sph",
                sort.by = "none")
sum(tab$adj.P.Val<=0.05,na.rm = T)
sum(tab$adj.P.Val<=0.01,na.rm = T)
limma <- cbind(filtered[,1:7], tab)

#calculate log2 fc cutoff

2 * sd(limma$logFC[limma$logFC > 0], na.rm = T)
# 0.89
sum(limma$adj.P.Val <= 0.05 & limma$logFC > 0.89, na.rm = T)
sum(limma$adj.P.Val <= 0.01 & limma$logFC > 0.89, na.rm = T)

2 * sd(limma$logFC[limma$logFC < 0], na.rm = T)
# 0.7
sum(limma$adj.P.Val <= 0.05 & limma$logFC < -0.7, na.rm = T)
sum(limma$adj.P.Val <= 0.01 & limma$logFC < -0.7, na.rm = T)

#add categorical column for up and down regulation
limma <- limma %>%
  mutate(
    Sph.vs.Adh = case_when(
      logFC > 0.89 & adj.P.Val <= 0.01 ~ "Up-regulated",
      logFC < -0.7 & adj.P.Val <= 0.01 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  )

limma <- cbind(filtered, limma[,8:14])

write.table(limma, "OC_Prot_norm_averaged_filtered_Limma.txt", row.names = F, sep = "\t")

limma <- fread("OC_Prot_norm_averaged_filtered_Limma.txt")
limma <- as.data.frame(limma)

#density plot of logFC
ggplot(limma, aes(x = logFC)) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 12) +
  xlab("Log2 fold change") +
  ylab("Density")

# Calculate frequency for regulated proteins per group ----------------------
top <- limma[limma$Sph.vs.Adh == "Up-regulated", ]
#count number of values in col_adh
top$Adh <- rowSums(!is.na(top[, col_adh]))
#count number of values in col_sph
top$Sph <- rowSums(!is.na(top[, col_sph]))

bottom <- limma[limma$Sph.vs.Adh == "Down-regulated", ]
#count number of values in col_adh
bottom$Adh <- rowSums(!is.na(bottom[, col_adh]))
#count number of values in col_sph
bottom$Sph <- rowSums(!is.na(bottom[, col_sph]))

freq <- rbind(top, bottom)
#keep first and last 3 columns
freq <- freq[, c(1, 34:36)]
#long format
freq.long <- reshape2::melt(freq, id = c("PG.ProteinGroups", "Sph.vs.Adh"))
freq.long$variable <- factor(freq.long$variable, levels = c("Adh", "Sph"))

#histogram
ggplot(freq.long, aes(x = value, fill = variable)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~Sph.vs.Adh) +
  theme_bw(base_size = 12) +
  xlab("Frequency") +
  ylab("Number of proteins") +
  scale_x_continuous(breaks = seq(3, 10, by = 1))

# Calculate frequency for regulated proteins per patient ---------------------
top <- limma[limma$Sph.vs.Adh == "Up-regulated", ]
bottom <- limma[limma$Sph.vs.Adh == "Down-regulated", ]
reg <- rbind(top, bottom)

#long format
reg.long <- reshape2::melt(reg[,c(1,8:27,34)], id = c("PG.ProteinGroups", "Sph.vs.Adh"))
reg.long$variable <- as.character(reg.long$variable)
reg.long <- reg.long[!is.na(reg.long$value),]
reg.long$patient <- sub("\\_.*", "", reg.long$variable)
reg.long$patient <-
  factor(reg.long$patient,
         levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))
reg.long$Sph.vs.Adh <- factor(reg.long$Sph.vs.Adh, levels = c("Up-regulated", "Down-regulated"))

reg.long %>% group_by(Sph.vs.Adh) %>% 
  summarise(n_distinct(PG.ProteinGroups))

count <- reg.long %>% group_by(patient, Sph.vs.Adh) %>% 
  summarise(count = n_distinct(PG.ProteinGroups))

#plot count
ggplot(count, aes(x = patient, y = count, fill = Sph.vs.Adh)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  )) +
  geom_text(aes(label = count), position = position_dodge(width = 1), vjust = -0.5, size = 3) +
  xlab("Patient") +
  ylab("Number of proteins") +
  ggtitle("Up-regulated proteins = 246 and Down-regulated proteins = 349")

# Fold changes for individual pairs of patients ---------------------------
Adh <- filtered[,col_adh]
Sph <- filtered[,col_sph]
fc <- Sph - Adh
colnames(fc) <- sub("\\_.*", "", colnames(fc))
fc <- cbind(filtered[,1:7], fc)
write.table(fc, "OC_Prot_filtered_fc.txt", row.names = F, sep = "\t")

fc <- fread("OC_Prot_filtered_fc.txt")
fc <- as.data.frame(fc)

# Long format fc ----------------------------------------------------------
fc.long <- reshape2::melt(fc[,c(1,8:17)], id = "PG.ProteinGroups", variable.name = "patient")
fc.long <- fc.long[!is.na(fc.long$value),]

#density plot
ggplot(fc.long, aes(x = value, color = patient)) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 12) +
  xlab("Log2 fold change") +
  ylab("Density")

# Heatmap fold changes ----------------------------------------------------
rownames(fc) <- filtered$PG.ProteinGroups
fc <- na.omit(fc)
breaks <- seq(-2, 2, length.out = 101)

pheatmap(fc, 
         color = colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(100),
         breaks = breaks, 
         show_rownames = F)

# Heatmap fold changes for regulated proteins between groups ---------------
fc <- Sph - Adh
rownames(fc) <- filtered$PG.ProteinGroups
fc_reg <- fc[limma$Sph.vs.Adh != "Unchanged",]
fc_reg <- na.omit(fc_reg)

breaks <- seq(-3, 3, length.out = 101)

pheatmap(fc_reg, 
         color = colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(100),
         breaks = breaks, 
         show_rownames = F)

# Mann-Whitney U Test to compare patient clusters --------------------------
fc_t <- as.data.frame(t(fc[,-c(1:7)]))
colnames(fc_t) <- filtered$PG.ProteinGroups
group <- ifelse(rownames(fc_t) == "P6" | rownames(fc_t) == "P8" | rownames(fc_t) == "P9" | rownames(fc_t) == "P10", "cluster2", "cluster1")
group <- factor(group, levels = c("cluster1", "cluster2"))

# Initialize a data frame to store results
results <- data.frame(
  Protein = colnames(fc_t),
  PValue = numeric(length(colnames(fc_t))),
  MeanDiff = numeric(length(colnames(fc_t)))
)

# Perform Mann-Whitney U test for each protein
for (i in seq_along(colnames(fc_t))) {
  protein <- colnames(fc_t)[i]
  
  # Subset the data for the two groups, removing NAs
  cluster1 <- fc_t[group == "cluster1", protein, drop = FALSE]
  cluster2 <- fc_t[group == "cluster2", protein, drop = FALSE]
  
  cluster1 <- cluster1[!is.na(cluster1)]  # Remove NAs from Group A
  cluster2 <- cluster2[!is.na(cluster2)]  # Remove NAs from Group B
  
  # Check if both groups have sufficient data for the test
  if (length(cluster1) > 1 && length(cluster2) > 1) {
    # Perform the Mann-Whitney U test
    test <- wilcox.test(cluster1, cluster2, exact = FALSE) # Use exact = FALSE for large datasets
    
    # Store the p-value and mean difference
    results$PValue[i] <- test$p.value
    results$MeanDiff[i] <- mean(cluster1) - mean(cluster2)
  } else {
    # If not enough data, assign NA
    results$PValue[i] <- NA
    results$MeanDiff[i] <- NA
  }
}

# Filter significant proteins
significant <- results[results$PValue < 0.02, ]
#remove NA
significant <- na.omit(significant)

# Heatmap fold changes for regulated proteins between patient clusters ----
#keep proteins in fc if they are in significant
fc_reg <- fc[fc$PG.ProteinGroups %in% significant$Protein,]
rownames(fc_reg) <- fc_reg$PG.ProteinGroups
fc_reg <- fc_reg[, -c(1:7)]
fc_reg <- na.omit(fc_reg)

breaks <- seq(-1.5, 1.5, length.out = 101)

set.seed(123)
res <- pheatmap(fc_reg, 
         color = colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(100),
         breaks = breaks, 
         show_rownames = F,
         clustering_distance_rows = "correlation",
         cutree_rows = 2)

fc_reg$PG.ProteinGroups <- rownames(fc_reg)
fc_reg <- left_join(fc_reg, fc[,1:4], by = "PG.ProteinGroups")

#get the proteins in each cluster
cluster <- cutree(res$tree_row, k = 2)
fc_reg$cluster <- cluster
write.table(fc_reg, "clusters.txt", sep = "\t", row.names = F)

fc_reg <- fread("clusters.txt")
fc_reg <- as.data.frame(fc_reg)

# VOLCANO PLOT ------------------------------------------------------------

#keep top 20 up and down
top <- limma[limma$Sph.vs.Adh == "Up-regulated", ]
bottom <- limma[limma$Sph.vs.Adh == "Down-regulated", ]
top <- top[order(top$t, decreasing = T), ]
bottom <- bottom[order(bottom$t, decreasing = F), ]
top30 <- head(top, 30)
bottom30 <- head(bottom, 30)
top_genes <- rbind(top30, bottom30)

limma$Sph.vs.Adh <- factor(limma$Sph.vs.Adh, levels = c("Up-regulated", "Down-regulated", "Unchanged"))

brewer.pal(n = 8, name = "Set2")

ggplot(data = limma, aes(x = logFC, y = -log(adj.P.Val, 10))) +
  geom_point(
    aes(color = Sph.vs.Adh, alpha = Sph.vs.Adh),
    shape = 16,
    size = 2
  ) +
  scale_color_manual(values = c("#FC8D62","#66C2A5","lightgrey")) +
  scale_alpha_manual(values = c(0.9, 0.9, 0.5)) +
  geom_text_repel(
    data = top_genes,
    aes(
      x = logFC,
      y = -log(adj.P.Val, 10),
      label = PG.Genes
    ),
    max.overlaps = Inf,
    size = 2,
    box.padding = 0.4
  ) +
  geom_hline(
    yintercept = -log10(0.01),
    linetype = "dotted",
    col = "black",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = c(-0.62, 0.87),
    linetype = "dotted",
    col = "black",
    linewidth = 1
  ) +
  xlab(label = "Log2FC_Sph.vs.Adh") +
  ylab(label = "-Log10_adj.P.Val") +
  theme_light(base_size = 12)

# hm PDGFR --------------------------------------------------------------
PDGFR <- fc[grepl("PDGFR", fc$PG.Genes),]
#long format
PDGFR <- reshape2::melt(PDGFR[,c(1:2,8:17)], id = c("PG.ProteinGroups", "PG.Genes"))

#order by patient and by group in descending order
PDGFR <- PDGFR[order(PDGFR$variable, decreasing = T),]
PDGFR$variable <- factor(PDGFR$variable, levels = unique(PDGFR$variable))

ggplot(PDGFR, aes(x = PG.Genes, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white", 
    midpoint = 0, 
    name = "Log2 fold-change Sph vs Adh"
  ) +
  geom_text(aes(label = round(value, 2)), size = 2) +
  theme_minimal(base_size = 14) +
  xlab("") +
  ylab("") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# PLOT PDGFRA --------------------------------------------------------------
PDGFRA <- melt.average[grepl("PDGFRA", melt.average$PG.Genes),]

#dot plot
ggplot(PDGFRA, aes(x = patient, y = value, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size=4) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("PDGFRA")

# PLOT PDGFRB --------------------------------------------------------------
PDGFRB <- melt.average[grepl("PDGFRB", melt.average$PG.Genes),]

#dot plot
ggplot(PDGFRB, aes(x = patient, y = value, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size=4) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("PDGFRB")

# PLOT GPNMB --------------------------------------------------------------
GPNMB <- melt.average[grepl("GPNMB", melt.average$PG.Genes),]

#dot plot
ggplot(GPNMB, aes(x = patient, y = value, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size=4) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("GPNMB")

# PLOT CHI3L1 --------------------------------------------------------------
CHI3L1 <- melt.average[grepl("CHI3L1", melt.average$PG.Genes),]

#dot plot
ggplot(CHI3L1, aes(x = patient, y = value, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size=4) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("CHI3L1")

# GO ENRICHMENT ANALYSIS --------------------------------------------------
library(topGO)
library(org.Hs.eg.db)

#remove everything after ";" in PG.Genes
mydata$Gene <- gsub(";.*", "", mydata$PG.Genes)
average$Gene <- gsub(";.*", "", average$PG.Genes)
limma$Gene <- gsub(";.*", "", limma$PG.Genes)

#define the background gene list
geneUniverse <- unique(mydata$Gene)

#define the list of genes of interest
up <- limma %>% filter(Sph.vs.Adh == "Up-regulated")
up <- as.character(up$Gene)
presence <- average %>% filter(presence == "Sph")
presence <- as.character(presence$Gene)
up <- c(up, presence)

down <- limma %>% filter(Sph.vs.Adh == "Down-regulated")
down <- as.character(down$Gene)
absence <- average %>% filter(presence == "Adh")
absence <- as.character(absence$Gene)
down <- c(down, absence)

#tell TopGO where the interesting genes appear in the 'geneUniverse' vector
geneListu <- factor(as.integer(geneUniverse %in% up))
names(geneListu) <- geneUniverse

geneListd <- factor(as.integer(geneUniverse %in% down))
names(geneListd) <- geneUniverse

# Define ontologies and parameters
ontologies <- c("BP", "CC", "MF")
geneLists <- list(up = geneListu, down = geneListd)
nodeSizes <- c(BP = 5, CC = 10, MF = 10)

# Create GO annotations
annotations <- lapply(ontologies, function(onto) {
  annFUN.org(
    whichOnto = onto,
    feasibleGenes = NULL,
    mapping = "org.Hs.eg.db",
    ID = "symbol"
  )
})
names(annotations) <- ontologies

# Create topGOdata objects
topGOdata <- list()
for (onto in ontologies) {
  for (direction in names(geneLists)) {
    topGOdata[[paste0(onto, direction)]] <- new(
      "topGOdata",
      description = paste(onto, "enrichment analysis", direction),
      ontology = onto,
      allGenes = geneLists[[direction]],
      annot = annFUN.GO2genes,
      GO2genes = annotations[[onto]],
      nodeSize = nodeSizes[[onto]]
    )
  }
}

# Perform Fisher exact test
fisherResults <- lapply(topGOdata, function(data) {
  runTest(data, algorithm = "weight01", statistic = "fisher")
})

# Generate tables of results
results <- lapply(names(topGOdata), function(name) {
  allGO <- usedGO(topGOdata[[name]])
  res <- GenTable(
    topGOdata[[name]],
    weightFisher = fisherResults[[name]],
    orderBy = fisherResults[[name]],
    ranksOf = "weightFisher",
    topNodes = length(allGO)
  )
  res$weightFisher <- as.numeric(res$weightFisher)
  res
})
names(results) <- names(topGOdata)

results[["CCup"]]$weightFisher[1:2] <- 1e-20

# Perform BH correction and calculate enrichment %
results <- lapply(results, function(res) {
  res$p.adj <- p.adjust(res$weightFisher, method = "BH")
  res <- res %>% mutate(EnrichmentS = Significant / Annotated)
  res
})

# Export topGO output
for (name in names(results)) {
  results[[name]]$table <- name
}

GO <- do.call(rbind, results)
write.table(GO, "GO_enrichment_new.txt", row.names = F, sep = "\t")

#export string input

mygenes <- genesInTerm(CCu, "GO:0062023")
mygenes <- mygenes$`GO:0062023`
#keep genes present in up
mygenesup <- mygenes[mygenes %in% up]
write.table(mygenesup, file = "GOCC_up_for_string.tsv", row.names = F)

myterms <- c("GO:0006027", "GO:0030574")
mygenes <- genesInTerm(BPu, myterms)
mygenes <- unique(c(mygenes$`GO:0006027`, mygenes$`GO:0030574`))
#keep genes present in up
mygenesup <- mygenes[mygenes %in% up]
write.table(mygenesup, file = "GOBP_up_for_string.tsv", row.names = F)

mygenes <- genesInTerm(CCd, "GO:0005925")
mygenes <- mygenes$`GO:0005925`
#keep genes present in down
mygenesdown <- mygenes[mygenes %in% down]
write.table(mygenesdown, file = "GOCC_down_for_string.tsv", row.names = F)

myterms <- c("GO:0007157", "GO:0034113")
mygenes <- genesInTerm(BPd, myterms)
mygenes <- unique(c(mygenes$`GO:0007157`, mygenes$`GO:0034113`))
#keep genes present in down
mygenesdown <- mygenes[mygenes %in% down]
write.table(mygenesdown, file = "GOBP_down_for_string.tsv", row.names = F)

#Plot significant GO terms

#BP Up
results[["BPup"]]$Term <- tools::toTitleCase(results[["BPup"]]$Term)

ggplot(results[["BPup"]][results[["BPup"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term, -p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 7) +
  scale_fill_gradient(low = "#F2F2F2", high = "#FC8D62") +
  ggtitle("Enrichment.biological.processes_Sph.vs.Adh") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

#BP Down
results[["BPdown"]]$Term <- tools::toTitleCase(results[["BPdown"]]$Term)

ggplot(results[["BPdown"]][results[["BPdown"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term,-p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21, color = "black", size=7) +
  scale_fill_gradient(low="#F2F2F2",high="#66C2A5") +
  ggtitle("Enrichment.biological.processes_Adh.vs-Sph") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

#CC Up
results[["CCup"]]$Term <- tools::toTitleCase(results[["CCup"]]$Term)

ggplot(results[["CCup"]][results[["CCup"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term,-p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21, color = "black", size=7) +
  scale_fill_gradient(low = "#F2F2F2", high = "#FC8D62") +
  ggtitle("Enrichment.cellular.components_Sph.vs.Adh") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

#CC Down
results[["CCdown"]]$Term <- tools::toTitleCase(results[["CCdown"]]$Term)

ggplot(results[["CCdown"]][results[["CCdown"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term, -p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 7) +
  scale_fill_gradient(low="#F2F2F2",high="#66C2A5") +
  ggtitle("Enrichment.cellular.components_Adh.vs.Sph") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

#MF Up
results[["MFup"]]$Term <- tools::toTitleCase(results[["MFup"]]$Term)

ggplot(results[["MFup"]][results[["MFup"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term,-p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21, color = "black", size=7) +
  scale_fill_gradient(low = "#F2F2F2", high = "#FC8D62") +
  ggtitle("Enrichment.molecular.functions_Sph.vs.Adh") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 12)

#MF Down
results[["MFdown"]]$Term <- tools::toTitleCase(results[["MFdown"]]$Term)

ggplot(results[["MFdown"]][results[["MFdown"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term,-p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21, color = "black", size=7) +
  scale_fill_gradient(low="#F2F2F2",high="#66C2A5") +
  ggtitle("Enrichment.molecular.functions_Adh.vs.Sph") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 12)

# KEGG ENRICHMENT ANALYSIS ------------------------------------------------
library(clusterProfiler)

#remove everything after ";" in PG.ProteinGroups
mydata$PG <- gsub(";.*", "", mydata$PG.ProteinGroups)
average$PG <- gsub(";.*", "", average$PG.ProteinGroups)
limma$PG <- gsub(";.*", "", limma$PG.ProteinGroups)

#define background
geneUniverseuniprot <- as.character(mydata$PG)

#define the list of genes of interest
#define the list of genes of interest
up <- limma %>% filter(Sph.vs.Adh == "Up-regulated")
up <- as.character(up$PG)
presence <- average %>% filter(presence == "Sph")
presence <- as.character(presence$PG)
up <- c(up, presence)

down <- limma %>% filter(Sph.vs.Adh == "Down-regulated")
down <- as.character(down$PG)
absence <- average %>% filter(presence == "Adh")
absence <- as.character(absence$PG)
down <- c(down, absence)

# Define input parameters
geneLists <- list(up = up, down = down)
results <- list()

# Perform enrichment in a loop
for (direction in names(geneLists)) {
  enrichment <- enrichKEGG(
    gene = geneLists[[direction]],
    organism = "hsa",
    keyType = "uniprot",
    pvalueCutoff = 0.1,
    pAdjustMethod = "BH",
    universe = geneUniverseuniprot,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
  )
  result <- as.data.frame(enrichment)
  result$direction <- direction
  results[[direction]] <- result
}

# Combine results into a single data frame
keggResults <- do.call(rbind, results)

#calculate count & enrichment %
kegg2 <- separate(keggResults, GeneRatio, c("Count", "Count.Sign"))
kegg2$Count <- as.numeric(kegg2$Count)
kegg2 <- separate(kegg2, BgRatio, c("CountGO", "Count.Bg"))
kegg2$CountGO <- as.numeric(kegg2$CountGO)
kegg2 <- kegg2 %>% mutate(EnrichmentS = Count / CountGO)

#export KEGG enrichment output
write.table(kegg2, "OC_Prot_KEGG_enrichment.txt", row.names = F, sep = "\t")

#plot data
ggplot(kegg2[kegg2$direction == "up",], aes(
  x = -log(p.adjust, 10),
  y = reorder(Description, -p.adjust),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 7) +
  scale_fill_gradient(low = "#F2F2F2", high = "#FC8D62") +
  ggtitle("Enrichment.KEGG_Sph.vs.Adh") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

ggplot(kegg2[kegg2$direction == "down",], aes(
  x = -log(p.adjust, 10),
  y = reorder(Description, -p.adjust),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 7) +
  scale_fill_gradient(low="#F2F2F2",high="#66C2A5") +
  ggtitle("Enrichment.KEGG_Adh.vs.Sph") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

# REACTOME ENRICHMENT ANALYSIS --------------------------------------------

#convert symbol (gene name) to entreid
ids <-
  bitr(geneUniverse,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = "org.Hs.eg.db")

#define background genes 
geneUni <- as.character(ids$ENTREZID)

#define genes of interest after converting symbol to entreid
up.e <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
up.e <- as.character(up.e$ENTREZID)

down.e <- bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
down.e <- as.character(down.e$ENTREZID)

library(ReactomePA)

# Define input parameters
geneLists <- list(up = up.e, down = down.e)
results <- list()

# Perform enrichment in a loop
for (direction in names(geneLists)) {
  enrichment <- enrichPathway(
    gene = geneLists[[direction]],
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    universe = geneUni,
    minGSSize = 10,
    maxGSSize = 500,
    readable = F
  )
  result <- as.data.frame(enrichment)
  result$direction <- direction
  results[[direction]] <- result
}

# Combine results into a single data frame
reactResults <- do.call(rbind, results)

#calculate count & enrichment %
react2 <- separate(reactResults, GeneRatio, c("Count", "Count.Sign"))
react2$Count <- as.numeric(react2$Count)
react2 <- separate(react2, BgRatio, c("CountGO", "Count.Bg"))
react2$CountGO <- as.numeric(react2$CountGO)
react2 <- react2 %>% mutate(EnrichmentS = Count / CountGO)

#export Reactome enrichment output
write.table(react2, "OC_Prot_React_enrichment.txt", row.names = F, sep = "\t")

#plot data
ggplot(react2[react2$direction == "up",], aes(
  x = -log(p.adjust, 10),
  y = reorder(Description, -p.adjust),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 5) +
  scale_fill_gradient(low = "#F2F2F2", high = "#FC8D62") +
  ggtitle("Enrichment.Reactome_Sph.vs.Adh") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 6)

ggplot(react2[react2$direction == "down",], aes(
  x = -log(p.adjust, 10),
  y = reorder(Description, -p.adjust),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 5) +
  scale_fill_gradient(low="#F2F2F2",high="#66C2A5") +
  ggtitle("Enrichment.Reactome_Adh.vs.Sph") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 4)

# Prepare data for Cytoscape ------------------------------------------------
#string
stringup <- fread("string-up.txt")
stringup <- left_join(stringup, mydata, by = c("id" = "Gene"))
write.table(stringup, file = "STRING_up_with_data.txt", row.names = F, sep = "\t")

stringdown <- fread("string-down.txt")
stringdown <- left_join(stringdown, mydata, by = c("id" = "Gene"))
write.table(stringdown, file = "STRING_down_with_data.txt", row.names = F, sep = "\t")

# GSEA --------------------------------------------------------------------

#retrieve annotiations FROM MSigDb
library(msigdbr)

H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol) #hallmark gene sets

C2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::filter(gs_subcat == "CGP") %>%
  dplyr::select(gs_name, human_gene_symbol) #hallmark curated gene sets CGP

#prepare GSEA input
geneList <- mydata$logFC
names(geneList) <- mydata$Gene
geneList <- geneList[!is.na(geneList)]
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

#GSEA
gseaH <- GSEA(
  geneList = geneList,
  TERM2GENE = H
)
gsea.H <- as.data.frame(gseaH)

gseaC2 <- GSEA(
  geneList = geneList,
  TERM2GENE = C2
)
gsea.C2 <- as.data.frame(gseaC2)

#export data
write.table(gsea.H, file = "gsea_hallmark.txt", sep = "\t", row.names = F)
write.table(gsea.C2, file = "gsea_C2_CGP.txt", sep = "\t", row.names = F)

#plot data
library(enrichplot)
gseaplot2(gseaH, geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", title = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", pvalue_table = T)
gseaplot2(gseaH, geneSetID = "HALLMARK_GLYCOLYSIS", title = "HALLMARK_GLYCOLYSIS", pvalue_table = T)
gseaplot2(gseaH, geneSetID = "HALLMARK_HYPOXIA", title = "HALLMARK_HYPOXIA", pvalue_table = T)
gseaplot2(gseaC2, geneSetID = "BOQUEST_STEM_CELL_UP", title = "BOQUEST_STEM_CELL_UP", pvalue_table = T)
gseaplot2(gseaC2, geneSetID = "BOQUEST_STEM_CELL_DN", title = "BOQUEST_STEM_CELL_DN", pvalue_table = T)
gseaplot2(gseaC2, geneSetID = "FOROUTAN_TGFB_EMT_DN", title = "FOROUTAN_TGFB_EMT_DN", pvalue_table = T)

# GO ENRICHMENT ANALYSIS FOR PROTEINS IN CLUSTERS 1 AND 2 ------------------

#remove everything after ";" in PG.Genes
mydata$Gene <- gsub(";.*", "", mydata$PG.Genes)
fc_reg$Gene <- gsub(";.*", "", fc_reg$PG.Genes)

#define the background gene list
geneUniverse <- unique(mydata$Gene)

#define the list of genes of interest
c1 <- fc_reg %>% filter(cluster == 1)
c1 <- as.character(c1$Gene)

c2 <- fc_reg %>% filter(cluster == 2)
c2 <- as.character(c2$Gene)

#tell TopGO where the interesting genes appear in the 'geneUniverse' vector
geneList1 <- factor(as.integer(geneUniverse %in% c1))
names(geneList1) <- geneUniverse

geneList2 <- factor(as.integer(geneUniverse %in% c2))
names(geneList2) <- geneUniverse

# Define ontologies and parameters
ontologies <- c("BP", "CC", "MF")
geneLists <- list(c1 = geneList1, c2 = geneList2)
nodeSizes <- c(BP = 5, CC = 10, MF = 10)

# Create GO annotations
annotations <- lapply(ontologies, function(onto) {
  annFUN.org(
    whichOnto = onto,
    feasibleGenes = NULL,
    mapping = "org.Hs.eg.db",
    ID = "symbol"
  )
})
names(annotations) <- ontologies

# Create topGOdata objects
topGOdata <- list()
for (onto in ontologies) {
  for (direction in names(geneLists)) {
    topGOdata[[paste0(onto, direction)]] <- new(
      "topGOdata",
      description = paste(onto, "enrichment analysis", direction),
      ontology = onto,
      allGenes = geneLists[[direction]],
      annot = annFUN.GO2genes,
      GO2genes = annotations[[onto]],
      nodeSize = nodeSizes[[onto]]
    )
  }
}

# Perform Fisher exact test
fisherResults <- lapply(topGOdata, function(data) {
  runTest(data, algorithm = "weight01", statistic = "fisher")
})

# Generate tables of results
results <- lapply(names(topGOdata), function(name) {
  allGO <- usedGO(topGOdata[[name]])
  res <- GenTable(
    topGOdata[[name]],
    weightFisher = fisherResults[[name]],
    orderBy = fisherResults[[name]],
    ranksOf = "weightFisher",
    topNodes = length(allGO)
  )
  res$weightFisher <- as.numeric(res$weightFisher)
  res
})
names(results) <- names(topGOdata)

# Perform BH correction and calculate enrichment %
results <- lapply(results, function(res) {
  res$p.adj <- p.adjust(res$weightFisher, method = "BH")
  res <- res %>% mutate(EnrichmentS = Significant / Annotated)
  res
})

# Export topGO output
for (name in names(results)) {
  results[[name]]$table <- name
}
GO <- do.call(rbind, results)
write.table(GO, "GO_enrichment_clusters.txt", row.names = F, sep = "\t")

#Plot significant GO terms

#BP c1
results[["BPc1"]]$Term <- tools::toTitleCase(results[["BPc1"]]$Term)

ggplot(results[["BPc1"]][results[["BPc1"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term, -p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 7) +
  scale_fill_gradient(low = "#F2F2F2", high = "#e7298a") +
  ggtitle("Biological processes Cluster1") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

#BP c2
results[["BPc2"]]$Term <- tools::toTitleCase(results[["BPc2"]]$Term)

ggplot(results[["BPc2"]][results[["BPc2"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term,-p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21, color = "black", size=7) +
  scale_fill_gradient(low="#F2F2F2",high="#7570b3") +
  ggtitle("Biological processes Cluster2") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

#CC c1
results[["CCc1"]]$Term <- tools::toTitleCase(results[["CCc1"]]$Term)

#replace NA in results[["CCc1"]]$p.adj with the smallest value
results[["CCc1"]]$p.adj[is.na(results[["CCc1"]]$p.adj)] <- min(results[["CCc1"]]$p.adj, na.rm = TRUE)

ggplot(results[["CCc1"]][results[["CCc1"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term,-p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21, color = "black", size=7) +
  scale_fill_gradient(low = "#F2F2F2", high = "#e7298a") +
  ggtitle("Cellular components Cluster1") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

#CC c2
results[["CCc2"]]$Term <- tools::toTitleCase(results[["CCc2"]]$Term)

ggplot(results[["CCc2"]][results[["CCc2"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term, -p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 7) +
  scale_fill_gradient(low="#F2F2F2",high="#7570b3") +
  ggtitle("Cellular components Cluster2") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

#MF c1
results[["MFc1"]]$Term <- tools::toTitleCase(results[["MFc1"]]$Term)

ggplot(results[["MFc1"]][results[["MFc1"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term,-p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21, color = "black", size=7) +
  scale_fill_gradient(low = "#F2F2F2", high = "#e7298a") +
  ggtitle("Molecular functions Cluster1") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 12)

#MF c2
results[["MFc2"]]$Term <- tools::toTitleCase(results[["MFc2"]]$Term)

ggplot(results[["MFc2"]][results[["MFc2"]]$p.adj<=0.05,], aes(
  x = -log(p.adj, 10),
  y = reorder(Term,-p.adj),
  fill = EnrichmentS
)) +
  geom_point(shape = 21, color = "black", size=7) +
  scale_fill_gradient(low="#F2F2F2",high="#7570b3") +
  ggtitle("Molecular functions Cluster2") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 12)

# KEGG ENRICHMENT ANALYSIS CLUSTER1 and CLUSTER2 ---------------------------
library(clusterProfiler)

#remove everything after ";" in PG.ProteinGroups
mydata$PG <- gsub(";.*", "", mydata$PG.ProteinGroups)
fc_reg$PG <- gsub(";.*", "", fc_reg$PG.ProteinGroups)

#define background
geneUniverseuniprot <- as.character(mydata$PG)

#define the list of genes of interest
c1 <- fc_reg %>% filter(cluster == 1)
c1 <- as.character(c1$PG)

c2 <- fc_reg %>% filter(cluster == 2)
c2 <- as.character(c2$PG)

# Define input parameters
geneLists <- list(c1 = c1, c2 = c2)
results <- list()

# Perform enrichment in a loop
for (cluster in names(geneLists)) {
  enrichment <- enrichKEGG(
    gene = geneLists[[cluster]],
    organism = "hsa",
    keyType = "uniprot",
    pvalueCutoff = 0.1,
    pAdjustMethod = "BH",
    universe = geneUniverseuniprot,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    use_internal_data = FALSE
  )
  result <- as.data.frame(enrichment)
  result$cluster <- cluster
  results[[cluster]] <- result
}

# Combine results into a single data frame
keggResults <- do.call(rbind, results)

#calculate count & enrichment %
kegg2 <- separate(keggResults, GeneRatio, c("Count", "Count.Sign"))
kegg2$Count <- as.numeric(kegg2$Count)
kegg2 <- separate(kegg2, BgRatio, c("CountGO", "Count.Bg"))
kegg2$CountGO <- as.numeric(kegg2$CountGO)
kegg2 <- kegg2 %>% mutate(EnrichmentS = Count / CountGO)

#export KEGG enrichment output
write.table(keggResults, "KEGG_enrichment_clusters.txt", row.names = F, sep = "\t")

#plot data
ggplot(keggResults %>% filter(cluster == "c1" & p.adjust <= 0.001), aes(
  x = -log(p.adjust, 10),
  y = reorder(Description, -p.adjust),
  fill = RichFactor
)) +
  geom_point(shape = 21,
             color = "black",
             size = 7) +
  scale_fill_gradient(low = "#F2F2F2", high = "#e7298a") +
  ggtitle("KEGG Cluster1") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

ggplot(keggResults %>% filter(cluster == "c2"), aes(
  x = -log(p.adjust, 10),
  y = reorder(Description, -p.adjust),
  fill = RichFactor
)) +
  geom_point(shape = 21,
             color = "black",
             size = 7) +
  scale_fill_gradient(low="#F2F2F2",high="#7570b3") +
  ggtitle("KEGG Cluster2") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 8)

# REACTOME ENRICHMENT ANALYSIS FOR PROTEINS IN CLUSTERS 1 AND 2 ------------

#convert symbol (gene name) to entreid
ids <-
  bitr(geneUniverse,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = "org.Hs.eg.db")

#define background genes 
geneUni <- as.character(ids$ENTREZID)

#define genes of interest after converting symbol to entreid
c1.e <- bitr(c1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
c1.e <- as.character(c1.e$ENTREZID)

c2.e <- bitr(c2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
c2.e <- as.character(c2.e$ENTREZID)

library(ReactomePA)

# Define input parameters
geneLists <- list(c1 = c1.e, c2 = c2.e)
results <- list()

# Perform enrichment in a loop
for (cluster in names(geneLists)) {
  enrichment <- enrichPathway(
    gene = geneLists[[cluster]],
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    universe = geneUni,
    minGSSize = 10,
    maxGSSize = 500,
    readable = F
  )
  result <- as.data.frame(enrichment)
  result$cluster <- cluster
  results[[cluster]] <- result
}

# Combine results into a single data frame
reactResults <- do.call(rbind, results)

#calculate count & enrichment %
react2 <- separate(reactResults, GeneRatio, c("Count", "Count.Sign"))
react2$Count <- as.numeric(react2$Count)
react2 <- separate(react2, BgRatio, c("CountGO", "Count.Bg"))
react2$CountGO <- as.numeric(react2$CountGO)
react2 <- react2 %>% mutate(EnrichmentS = Count / CountGO)

#export Reactome enrichment output
write.table(react2, "React_enrichment_clusters.txt", row.names = F, sep = "\t")

react2 <- fread("React_enrichment_clusters.txt")
react2 <- as.data.frame(react2)

#plot data
ggplot(react2 %>% filter(cluster == "c1") %>% slice_head(n = 5), aes(
  x = -log(p.adjust, 10),
  y = reorder(Description, -p.adjust),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 5) +
  scale_fill_gradient(low = "#F2F2F2", high = "#e7298a") +
  ggtitle("Reactome Cluster1") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 10)

ggplot(react2 %>% filter(cluster == "c2") %>% slice_head(n = 5), aes(
  x = -log(p.adjust, 10),
  y = reorder(Description, -p.adjust),
  fill = EnrichmentS
)) +
  geom_point(shape = 21,
             color = "black",
             size = 5) +
  scale_fill_gradient(low="#F2F2F2",high="#7570b3") +
  ggtitle("Reactome Cluster2") +
  xlab("-Log10 p.adj.") +
  ylab("") +
  theme_bw(base_size = 10)





