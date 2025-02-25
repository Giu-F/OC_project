
# Set working directory -----------------------------------------------
setwd("N:/SUN-CPR-proteomics_jvo/GIULIA/B_COLLAB_EPIC-XS/UgoC/UgoC_OC_TMT-Phos_Mq/Mq_TMT.phos_reanalysis_MBR_2.6.6/combined/txt")

# Load OC phosphoproteome TMT data ------------------------------------
library(data.table)
mydata <- fread("Log2_exp_Phospho(STY)Sites.txt") #In Perseus: Log2, expand site table, filter min 1 valid value per row
mydata <- as.data.frame(mydata)

# Add m/z -------------------------------------------------------------
evidence <- fread("evidence.txt")
evidence <- as.data.frame(evidence)

#keep id and m/z columns
library(tidyverse)
evidence_f <- evidence %>% select("id", "m/z")

#left join evidence to mydata
mydata <- left_join(mydata, evidence_f, by = c("Best PEP evidence ID" = "id"))

# Filter data ----------------------------------------------------------

#remove rows without a gene name
mydata <- mydata[!mydata$`Gene names` == "",]

#remove contaminants
mydata <- mydata[!mydata$`Potential contaminant` == "+",]

#remove reverse hits
mydata <- mydata[!mydata$Reverse == "+",]

#remove columns that contain 11 
mydata <- mydata[,!grepl("11", colnames(mydata))]

#remove empty columns
mydata <- mydata[,colSums(is.na(mydata)) < nrow(mydata)]

#remove columns with "" values
mydata <- mydata[,colSums(mydata == "") == 0]

# Rename columns -------------------------------------------------------

#make a condition vector
condition <- rep(c("Adh", "Sph"), 10)

#make a subject vector
subject <- c(rep("P1", 2), rep("P2", 2), rep("P3", 2), rep("P4", 2), rep("P5", 2), rep("P6", 2), rep("P7", 2), rep("P8", 2), rep("P9", 2), rep("P10", 2))

#rename sample columns
colnames(mydata)[1:20] <- paste(subject, condition, sep = "_")

#invert names of columns 19 and 20 (Adh and Sph from P10)
colnames(mydata)[19:20] <- colnames(mydata)[20:19]

#invert position of columns 19 and 20
mydata <- mydata[,c(1:18, 20, 19, 21:49)]

write.table(mydata, file = "OC_Phos_cleaned.txt",row.names=F, sep="\t")

mydata <- fread("OC_Phos_cleaned.txt")
mydata <- as.data.frame(mydata)

# Count number of proteins ---------------------------------------------
length(unique(mydata$`Gene names`))

TMT1 <- mydata[,c(1:10, 27)]
#remove rows with all NA
TMT1 <- TMT1[rowSums(is.na(TMT1[,1:10])) < ncol(TMT1[,1:10]), ]
TMT1 <- TMT1[TMT1$`Localization prob` >= 0.75,]
#count how many rows have NA
sum(rowSums(is.na(TMT1[,1:10])) > 0)

TMT2 <- mydata[,c(11:20, 27)]
#remove rows with all NA
TMT2 <- TMT2[rowSums(is.na(TMT2[,1:10])) < ncol(TMT2[,1:10]), ]
TMT2 <- TMT2[TMT2$`Localization prob` >= 0.75,]
#count how many rows have NA
sum(rowSums(is.na(TMT2[,1:10])) > 0)

# Count ids per TMT ----------------------------------------------------

#sum first 10 columns and then the next 10
TMT1 <- rowSums(mydata[,1:10], na.rm = TRUE)
TMT2 <- rowSums(mydata[,11:20], na.rm = TRUE)

#remove elements = 0
TMT1 <- TMT1[TMT1 != 0]
TMT2 <- TMT2[TMT2 != 0]

TMT <- data.frame("TMT" = c("TMT1", "TMT2"), "Number.of.phosphosites" = c(length(TMT1), length(TMT2)), "Localized" = "<0.75")

#keep localized phosphosites
mydata <- mydata[mydata$`Localization prob` >= 0.75,]

write.table(mydata, file = "OC_Phos_localized.txt",row.names=F, sep="\t")

TMT1 <- rowSums(mydata[,1:10], na.rm = TRUE)
TMT2 <- rowSums(mydata[,11:20], na.rm = TRUE)

#remove elements = 0
TMT1 <- TMT1[TMT1 != 0]
TMT2 <- TMT2[TMT2 != 0]

TMT2 <- data.frame("TMT" = c("TMT1", "TMT2"), "Number.of.phosphosites" = c(length(TMT1), length(TMT2)), "Localized" = ">=0.75")
TMT1 <- TMT
TMT1$Number.of.phosphosites <- TMT1$Number.of.phosphosites - TMT2$Number.of.phosphosites

TMT3 <- rbind(TMT1, TMT2)

write.table(TMT3, file = "OC_Phos_count.txt",row.names=F, sep="\t")

#stacked bar plot
library(RColorBrewer)
ggplot(TMT3, aes(x = TMT, y = Number.of.phosphosites, fill = Localized)) +
  geom_bar(stat = "identity", color = "darkgrey") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw(base_size = 12) +
  scale_y_continuous(breaks = seq(0, 10000, by = 2000))

# Venn diagram between TMT1 and TMT2 -----------------------------------
TMT1 <- mydata[,c(1:10,48)]
TMT2 <- mydata[,c(11:20,48)]

#remove rows with all NA
TMT1 <- TMT1[rowSums(is.na(TMT1[,1:10])) < ncol(TMT1[,1:10]), ]
TMT2 <- TMT2[rowSums(is.na(TMT2[,1:10])) < ncol(TMT2[,1:10]), ]

#keep only the last column
TMT1 <- TMT1$`Unique identifier`
TMT2 <- TMT2$`Unique identifier`

library(eulerr)
overlap <- list(
  TMT1 = TMT1,
  TMT2 = TMT2
)

plot(euler(overlap),fills = c("lightblue", "lightpink"), quantities = TRUE)

# Boxplot --------------------------------------------------------------
melt <- reshape2::melt(mydata[,c(1:20, 48)], id = "Unique identifier")
melt <- melt[!is.na(melt$value),]

ggplot(melt, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 1,
    hjust = 1
  )) +
  xlab("MS run") +
  ylab("Log2 MS signal") +
  ggtitle("Before normalization")

# Density plot ---------------------------------------------------------
ggplot(melt, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 12) +
  xlab("Log2 MS signal") +
  ylab("Density") +
  ggtitle("Before normalization") +
  theme(legend.position = "none")

# Quantile normalization ------------------------------------------------
library(limma)
qb <- limma::normalizeQuantiles(mydata[,1:20])
qb <- cbind(qb, mydata[,21:49])

write.table(qb, file = "OC_Phos_Qb-normalized.txt",row.names=F, sep="\t")

qb <- fread("OC_Phos_Qb-normalized.txt")
qb <- as.data.frame(qb)

# Boxplot after normalization ------------------------------------------
melt.qb <- reshape2::melt(qb[,c(1:20, 48)], id = "Unique identifier")
melt.qb <- melt.qb[!is.na(melt.qb$value),]

ggplot(melt.qb, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 1,
    hjust = 1
  )) +
  xlab("MS run") +
  ylab("Log2 normalized MS signal") +
  ggtitle("After normalization")

# Density plot after normalization -------------------------------------
ggplot(melt.qb, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  theme_bw(base_size = 12) +
  xlab("Log2 normalized MS signal") +
  ylab("Density") +
  ggtitle("After normalization") +
  theme(legend.position = "none")

# Normalize to proteome ------------------------------------------------
prot <- fread("OC_Prot_norm_average.txt")
prot <- as.data.frame(prot)
prot$uniprot <- gsub(";.*", "", prot$PG.ProteinGroups)

qb <- left_join(qb, prot[,8:28], by = c("Protein" = "uniprot"), suffix = c("", ".prot"))

norm <- qb[,1:20] - qb[,50:69]
norm <- cbind(norm, qb[,21:49])
#remove rows with all NA
norm <- norm[rowSums(is.na(norm[,1:20])) < ncol(norm[,1:20]), ]

write.table(norm, file = "OC_Phos_norm-on-prot.txt",row.names=F, sep="\t")

# PCA ------------------------------------------------------------------
#make a condition vector
condition <- rep(c("Adh", "Sph"), 10)

#make a subject vector
subject <- c(rep("P1", 2), rep("P2", 2), rep("P3", 2), rep("P4", 2), rep("P5", 2), rep("P6", 2), rep("P7", 2), rep("P8", 2), rep("P9", 2), rep("P10", 2))

#remove missing values
nona <- na.omit(qb[,1:20])
#nona <- na.omit(norm[,1:20])

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
  obs.scale = 1,
  var.scale = 0,
  groups = condition,
  ellipse = F
) +
  geom_point(aes(colour = condition), size = 4) +
  scale_color_brewer(palette = "Set2") +
  geom_text_repel(aes(label = subject), size = 4, box.padding = 0.5, max.overlaps = Inf) +
  theme_light(base_size = 12)

# Limma paired analysis ------------------------------------------------

#filter data
col_adh <- seq(1, ncol(qb[,1:20]), by = 2)
col_sph <- seq(2, ncol(qb[,1:20]), by = 2)

#filtered <- qb[rowSums(!is.na(qb[, col_adh])) >= 3 & rowSums(!is.na(qb[, col_sph])) >= 3,]

#perform a Limma paired analysis
design <- model.matrix(~ subject + condition)
colnames(design) <- gsub("subject", "", colnames(design))
colnames(design) <- gsub("condition", "", colnames(design))

fit <- lmFit(qb[,1:20], design)
#fit <- lmFit(norm[,1:20], design)

fit2 <- eBayes(fit)
tab <- topTable(fit2, n = Inf, coef = "Sph", sort.by = "none")
limma <- cbind(qb,tab)
#limma <- cbind(norm,tab)

#calculate log2 fc cutoff

2 * sd(limma$logFC [limma$logFC > 0], na.rm = T)
# 0.83
# 0.87
sum(limma$adj.P.Val < 0.05 & limma$logFC > 0.83, na.rm = T)
sum(limma$adj.P.Val < 0.01 & limma$logFC > 0.83, na.rm = T)

2 * sd(limma$logFC [limma$logFC < 0], na.rm = T)
# 0.76
# 0.74
sum(limma$adj.P.Val < 0.05 & limma$logFC < -0.76, na.rm = T)
sum(limma$adj.P.Val < 0.01 & limma$logFC < -0.76, na.rm = T)

#add categorical column for up and down regulation
limma <- limma %>%
  mutate(
    Sph.vs.Adh = case_when(
      logFC >= 0.83 & adj.P.Val <= 0.01 ~ "Up-regulated",
      logFC <= -0.76 & adj.P.Val <= 0.01 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  )

limma <- limma %>%
  mutate(
    Sph.vs.Adh = case_when(
      logFC >= 0.87 & adj.P.Val <= 0.01 ~ "Up-regulated",
      logFC <= -0.74 & adj.P.Val <= 0.01 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  )

#count up and down regulated sites
sum(limma$Sph.vs.Adh == "Up-regulated")
sum(limma$Sph.vs.Adh == "Down-regulated")

aap <- paste(limma$`Amino acid`, limma$Position, sep = "")
limma$Gene <- gsub(";.*", "", limma$`Gene names`)
limma$gene_site <- paste(limma$Gene, aap, sep = "_")

#export Limma table
write.table(limma, file = "OC_Phos_Qb-norm_Limma.txt", row.names=F, sep="\t")
#write.table(limma, file = "OC_Phos_norm-on-prot_Limma.txt", row.names=F, sep="\t")

limma <- fread("OC_Phos_Qb-norm_Limma.txt")
limma <- as.data.frame(limma)

# Fold changes for individual pairs of patients ---------------------------
Adh <- qb[,col_adh]
Sph <- qb[,col_sph]

#Adh <- norm[,col_adh]
#Sph <- norm[,col_sph]

fc <- Sph - Adh
colnames(fc) <- sub("\\_.*", "", colnames(fc))
fc <- cbind(fc, limma[,21:58])

write.table(fc, "OC_Phos_Qb-norm_fc.txt", row.names = F, sep = "\t")
#write.table(fc, "OC_Phos_norm-on-prot_fc.txt", row.names = F, sep = "\t")

fc <- fread("OC_Phos_Qb-norm_fc.txt")
fc <- as.data.frame(fc)

# VOLCANO PLOT ------------------------------------------------------------

#keep top 20 up and down
top <- limma[limma$Sph.vs.Adh == "Up-regulated",]
bottom <- limma[limma$Sph.vs.Adh == "Down-regulated",]
top <- top[order(top$t, decreasing = T),]
bottom <- bottom[order(bottom$t, decreasing = F),]
top30 <- head(top, 30)
bottom30 <- head(bottom, 30)
top_sites <- rbind(top20, bottom20)

#keep PDGFR sites
PDGFR <- limma[grepl("PDGFR", limma$'Gene names'),]
PDGFR <- PDGFR[PDGFR$Sph.vs.Adh == "Up-regulated",]

limma$Sph.vs.Adh <-
  factor(limma$Sph.vs.Adh,
         levels = c("Up-regulated", "Down-regulated", "Unchanged"))
#brewer.pal(n = 8, name = "Set2")

library(ggrepel)
ggplot(data = limma, aes(x = logFC, y = -log(adj.P.Val, 10))) +
  geom_point(
    aes(color = Sph.vs.Adh, alpha = Sph.vs.Adh),
    shape = 16,
    size = 2
  ) +
  scale_color_manual(values = c("#FC8D62","#66C2A5","#B3B3B3")) +
  scale_alpha_manual(values = c(0.9, 0.9, 0.5)) +
  geom_point(
    data = PDGFR,
    aes(x = logFC, y = -log(adj.P.Val, 10)),
    shape = 15, # Square shape
    size = 3,
    color = "black"
  ) +
  geom_text_repel(
    #data = top_sites,
    data = PDGFR,
    aes(
      x = logFC,
      y = -log(adj.P.Val, 10),
      label = gene_site
    ),
    max.overlaps = Inf,
    size = 4,
    box.padding = 0.4
  ) +
  geom_hline(
    yintercept = -log10(0.01),
    linetype = "dotted",
    col = "black",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = c(-0.71, 0.84),
    linetype = "dotted",
    col = "black",
    linewidth = 1
  ) +
  xlab(label = "Log2FC_Sph.vs.Adh") +
  ylab(label = "-Log10_adj.P.Val") +
  theme_light(base_size = 12)

# PDGFR plot ------------------------------------------------------------
PDGFR <- fc[grepl("PDGFR", fc$'Gene names'),]
#PDGFR <- PDGFR[PDGFR$Sph.vs.Adh == "Up-regulated",]
#long format
PDGFR <- reshape2::melt(PDGFR[,c(1:10,48,16)], id = c("gene_site", "Multiplicity"))
PDGFR$Multiplicity <- gsub("_", "", PDGFR$Multiplicity)
PDGFR$id <- paste0(PDGFR$gene_site, "_", PDGFR$Multiplicity)

#order by patient and by group in descending order
PDGFR <- PDGFR[order(PDGFR$variable, decreasing = T),]
PDGFR$variable <- factor(PDGFR$variable, levels = unique(PDGFR$variable))

ggplot(PDGFR, aes(x = id, y = variable, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white", 
    midpoint = 0, 
    name = "Log2 fold-change Sph vs Adh"
  ) +
  geom_text(aes(label = round(value, 2)), size = 3) +
  theme_minimal(base_size = 14) +
  xlab("") +
  ylab("") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# RoKAI analysis -------------------------------------------------------

#generate RoKAI input with logFC

#average the same site
limma_av <- limma %>%
  group_by(gene_site, Protein, `Amino acid`, Position) %>%
  summarise(Quantification = mean(logFC, na.rm = T)) %>%
  ungroup()

aap <- paste(limma_av$`Amino acid`, limma_av$Position, sep = "")
df <- as.data.frame(cbind(limma_av$Protein, aap, limma_av$Quantification))
colnames(df) <- c("Protein", "Position", "Quantification")

write.csv(df, file = "RoKAI_input_logFC.csv", row.names=TRUE)

#plot RoKAI output
RoKAI <- fread("kinase_table.csv")
RoKAI <- RoKAI[RoKAI$PValue<=0.05,]
RoKAI <- RoKAI[RoKAI$NumSubs>=3,]

library(khroma)
ggplot(data = RoKAI, aes(x = reorder (Name,-ZScore), y = ZScore)) +
  geom_bar(stat = "identity", aes(fill = ZScore), color = "black") +
  geom_text(aes(label = NumSubs),
            vjust = 1.3,
            size = 3) +
  scale_fill_sunset() +
  theme_bw(base_size = 14) +
  scale_y_continuous(limits = c(-6, 6)) +
  xlab("Kinase name") +
  theme(axis.text.x = element_text(
    angle = 70,
    vjust = 1,
    hjust = 1,
    color = "black"
  ))

# Plot housekeeping proteins and proteins in the annotated network -------
limma$site_pvalue <- paste0(limma$`Amino acid`, limma$Position, " ", limma$P.Value)
melt.qb <- reshape2::melt(limma[,c(1:20, 42, 48, 59)], id = c("Unique identifier", "Gene names", "site_pvalue"))
melt.qb <- melt.qb[!is.na(melt.qb$value),]
#remove text after _
melt.qb$patient <- gsub("_.*", "", melt.qb$variable)
melt.qb$patient <- factor(melt.qb$patient, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))
# remove text before _
melt.qb$group <- gsub(".*_", "", melt.qb$variable)

min(melt.qb$value)
max(melt.qb$value)

#MYH9
MYH9 <- melt.qb[melt.qb$`Gene names` == "MYH9",]

#dot plot
ggplot(MYH9, aes(x = patient, y = value, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size=4) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("MYH9") +
  ylim(-4, 3) +
  facet_wrap(~site_pvalue, ncol = 3, nrow = 3)

#PKM
PKM <- melt.qb[melt.qb$`Gene names` == "PKM",]

#dot plot
ggplot(PKM, aes(x = patient, y = value, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size=4) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("PKM") +
  ylim(-4, 3) +
  facet_wrap(~site_pvalue)

#GAPDH
GAPDH <- melt.qb[melt.qb$`Gene names` == "GAPDH",]

#dot plot
ggplot(GAPDH, aes(x = patient, y = value, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size=4) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("GAPDH") +
  ylim(-4, 3) +
  facet_wrap(~site_pvalue)

#ENO1
ENO1 <- melt.qb[melt.qb$`Gene names` == "ENO1",]

#dot plot
ggplot(ENO1, aes(x = patient, y = value, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size=4) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("ENO1") +
  ylim(-4, 3) +
  facet_grid (~site_pvalue, ncol = 3, nrow = 3)

#LMNB1
LMNB1 <- melt.qb[melt.qb$`Gene names` == "LMNB1",]

#dot plot
ggplot(LMNB1, aes(x = patient, y = value, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size=4) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("LMNB1") +
  ylim(-4, 3) +
  facet_wrap(~site_pvalue)

# Plot kinase substrates ------------------------------------------------
limma$site <- paste0(limma$`Amino acid`, limma$Position)

setwd("N:/SUN-CPR-proteomics_jvo/GIULIA/PSP_datasets")
ks <- fread("Kinase_Substrate_Dataset.gz")

# List of kinases
kinases <- c("ERK1", "ERK2")

# Initialize a list to store results for each kinase
results <- list()

for (kinase in kinases) {
  # Filter the ks dataset for the current kinase
  ks_filtered <- ks %>% filter(KIN_ORGANISM == "human", SUB_ORGANISM == "human", KINASE == kinase)
  
  # Perform the join with the limma dataset
  matched_data <- limma %>% left_join(ks_filtered, by = c("Gene" = "SUB_GENE", "site" = "SUB_MOD_RSD"))
  
  # Keep rows with a match
  matched_data <- matched_data[!is.na(matched_data$GENE),]
  
  # Melt the data
  melted_data <- reshape2::melt(matched_data[, c(1:20, 48)], id = "Unique identifier")
  
  # Remove rows with NA values
  melted_data <- melted_data[!is.na(melted_data$value),]
  
  # Extract patient information
  melted_data$patient <- gsub("_.*", "", melted_data$variable)
  melted_data$patient <- factor(melted_data$patient, 
                                levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))
  
  # Extract group information
  melted_data$group <- gsub(".*_", "", melted_data$variable)
  
  # Store the processed data in the results list
  results[[kinase]] <- melted_data
}

# Access the results for ERK1 and ERK2
ERK1_melt <- results[["ERK1"]]
ERK2_melt <- results[["ERK2"]]

#boxplot
ggplot(ERK1_melt, aes(x = patient, y = value, color = group)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("ERK1 substrates")

ggplot(ERK2_melt, aes(x = patient, y = value, color = group)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("ERK2 substrates")

#CSK
ks_filtered <- ks %>% filter(KIN_ORGANISM == "human",
                        SUB_ORGANISM == "human",
                        grepl("CSNK", GENE)) %>%
  distinct(SUB_GENE, SUB_MOD_RSD, .keep_all = T)

# Perform the join with the limma dataset
matched_data <- limma %>% left_join(ks_filtered, by = c("Gene" = "SUB_GENE", "site" = "SUB_MOD_RSD"))

# Keep rows with a match
matched_data <- matched_data[!is.na(matched_data$GENE),]

# Melt the data
melted_data <- reshape2::melt(matched_data[, c(1:20, 48)], id = "Unique identifier")

# Remove rows with NA values
melted_data <- melted_data[!is.na(melted_data$value),]

# Extract patient information
melted_data$patient <- gsub("_.*", "", melted_data$variable)
melted_data$patient <- factor(melted_data$patient, 
                              levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))

# Extract group information
melted_data$group <- gsub(".*_", "", melted_data$variable)

ggplot(melted_data, aes(x = patient, y = value, color = group)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("Log2 MS intensity") +
  ggtitle("CK substrates")






  
