
# Set working directory -----------------------------------------------
setwd("N:/SUN-CPR-proteomics_jvo/GIULIA/B_COLLAB_EPIC-XS/UgoC/UgoC_OC_TMT-Phos_Mq/MQ_Txt_Folder")

# Load OC phosphoproteome TMT data ------------------------------------
library(data.table)
mydata <- fread("pSTY_log2_from.Perseus.txt") #In Perseus: Log2, expand site table, filter min 1 valid value per row
mydata <- as.data.frame(mydata)

# Rename columns -------------------------------------------------------

#remove everything before ": " in column names
colnames(mydata) <- gsub(".*: ", "", colnames(mydata))

# Filter data ----------------------------------------------------------

#remove rows without a gene name
mydata <- mydata[!mydata$`Gene names` == "",]

#remove contaminants
mydata <- mydata[!mydata$`Potential contaminant` == "+",]

#remove reverse hits
mydata <- mydata[!mydata$Reverse == "+",]

length(unique(mydata$`Gene names`))

# Count ids per TMT ----------------------------------------------------

#sum first 11 columns and then the next 11
TMT1 <- rowSums(mydata[,1:11], na.rm = TRUE)
TMT2 <- rowSums(mydata[,12:22], na.rm = TRUE)

#remove elements = 0
TMT1 <- TMT1[TMT1 != 0]
TMT2 <- TMT2[TMT2 != 0]

TMT <- data.frame("TMT" = c("TMT1", "TMT2"), "Number.of.phosphosites" = c(length(TMT1), length(TMT2)), "Localized" = "<0.75")

#keep localized phosphosites
mydata <- mydata[mydata$`Localization prob` >= 0.75,]

TMT1 <- rowSums(mydata[,1:11], na.rm = TRUE)
TMT2 <- rowSums(mydata[,12:22], na.rm = TRUE)

#remove elements = 0
TMT1 <- TMT1[TMT1 != 0]
TMT2 <- TMT2[TMT2 != 0]

TMT2 <- data.frame("TMT" = c("TMT1", "TMT2"), "Number.of.phosphosites" = c(length(TMT1), length(TMT2)), "Localized" = ">=0.75")
TMT1 <- TMT
TMT1$Number.of.phosphosites <- TMT1$Number.of.phosphosites - TMT2$Number.of.phosphosites

TMT3 <- rbind(TMT1, TMT2)

#stacked bar plot
library(tidyverse)
library(RColorBrewer)
ggplot(TMT3, aes(x = TMT, y = Number.of.phosphosites, fill = Localized)) +
  geom_bar(stat = "identity", color = "darkgrey") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw(base_size = 12)

# Boxplot --------------------------------------------------------------
melt <- reshape2::melt(mydata[,c(1:22, 46)], id = "Unique identifier")
melt <- na.omit(melt)

library(tidyverse)
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

# Pool normalization ---------------------------------------------------
norm1 <- mydata[,1:10] - mydata$`Reporter intensity 11 TMT1`
norm2 <- mydata[,12:21] - mydata$`Reporter intensity 11 TMT2`
meanPOOL <- rowMeans(cbind(mydata$`Reporter intensity 11 TMT1`,mydata$`Reporter intensity 11 TMT2`), na.rm=TRUE)
norm1 <- norm1 + meanPOOL
norm2 <- norm2 + meanPOOL
norm <- cbind(norm1, norm2,mydata[, 23:47])

#remove rows without numeric values
norm <- norm[rowSums(is.na(norm[,1:20])) < ncol(norm[,1:20]), ] #11660

# Rename columns -------------------------------------------------------

#make a condition vector
condition <- rep(c("Adh", "Sph"), 10)

#make a subject vector
subject <- c(rep("P1", 2), rep("P2", 2), rep("P3", 2), rep("P4", 2), rep("P5", 2), rep("P6", 2), rep("P7", 2), rep("P8", 2), rep("P9", 2), rep("P10", 2))

#rename sample columns
colnames(norm)[1:20] <- paste(subject, condition, sep = "_")

#invert names of columns 19 and 20
colnames(norm)[19:20] <- colnames(norm)[20:19]

#invert position of columns 19 and 20
norm <- norm[,c(1:18, 20, 19, 21:45)]

# Quantile normalization ------------------------------------------------
library(limma)
qb <- limma::normalizeQuantiles(norm[,1:20])
qb <- cbind(qb, norm[, 21:45])

# Boxplot after normalization ------------------------------------------
melt.qb <- reshape2::melt(qb[,c(1:20, 44)], id = "Unique identifier")
melt.qb <- na.omit(melt.qb)

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

# PCA ------------------------------------------------------------------

#remove missing values
nona <- na.omit(qb[,1:20])

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
  groups = condition,
  ellipse = F
) +
  geom_point(aes(colour = condition), size = 4) +
  scale_color_brewer(palette = "Set2") +
  geom_text_repel(aes(label = subject), size = 4, box.padding = 0.5, max.overlaps = Inf) +
  theme_light(base_size = 12)

# Export data ----------------------------------------------------------
write.table(qb, file = "OC_Phos_normalized.txt",row.names=F, sep="\t")

# Limma paired analysis ------------------------------------------------

#perform a Limma paired analysis
design <- model.matrix(~ subject + condition)
colnames(design) <- gsub("subject", "", colnames(design))
colnames(design) <- gsub("condition", "", colnames(design))
fit <- lmFit(qb[,1:20], design)
fit2 <- eBayes(fit)
tab <- topTable(fit2, n = Inf, coef = "Sph", sort.by = "none")
limma <- cbind(qb, tab)

#calculate log2 fc cutoff

2 * sd(limma$logFC [limma$logFC > 0], na.rm = T)
# 0.84
sum(limma$adj.P.Val < 0.05 & limma$logFC > 0.84, na.rm = T)
sum(limma$adj.P.Val < 0.01 & limma$logFC > 0.84, na.rm = T)

2 * sd(limma$logFC [limma$logFC < 0], na.rm = T)
# 0.71
sum(limma$adj.P.Val < 0.05 & limma$logFC < -0.71, na.rm = T)
sum(limma$adj.P.Val < 0.01 & limma$logFC < -0.71, na.rm = T)

#add categorical column for up and down regulation
limma <- limma %>%
  mutate(
    Sph.vs.Adh = case_when(
      logFC >= 0.84 & adj.P.Val <= 0.01 ~ "Up-regulated",
      logFC <= -0.71 & adj.P.Val <= 0.01 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  )

aap <- paste(limma$`Amino acid`, limma$Position, sep = "")
limma$Gene <- gsub(";.*", "", limma$`Gene names`)
limma$gene_site <- paste(limma$Gene, aap, sep = "_")

#export Limma table
write.table(table, file = "OC_Phos_Qb-norm_Limma.txt",row.names=F, sep="\t")

# VOLCANO PLOT ------------------------------------------------------------

#keep top 20 up and down
top <- limma[limma$Sph.vs.Adh == "Up-regulated",]
bottom <- limma[limma$Sph.vs.Adh == "Down-regulated",]
top <- top[order(top$t, decreasing = T),]
bottom <- bottom[order(bottom$t, decreasing = F),]
top20 <- head(top, 20)
bottom20 <- head(bottom, 20)
top_sites <- rbind(top20, bottom20)

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
  geom_text_repel(
    data = top_sites,
    aes(
      x = logFC,
      y = -log(adj.P.Val, 10),
      label = gene_site
    ),
    max.overlaps = Inf,
    size = 3,
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

#PDGFR plot
PDGFR <- limma[grepl("PDGFR", limma$Gene, fixed = T) & limma$Sph.vs.Adh == "Up-regulated", ]
PDGFR <- PDGFR[-4,]
PDGFR_long <- reshape2::melt(PDGFR[,c(1:20, 54)], id = "gene_site")
PDGFR_long <- na.omit(PDGFR_long)
PDGFR_long$condition <- gsub(".*_", "", PDGFR_long$variable)

ggplot(PDGFR_long, aes(x = gene_site, y = value, fill = condition)) +
  geom_point(
    shape = 21,
    color = "black",
    size = 5,
    position = position_jitterdodge(jitter.width = 0.2)
  ) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  stat_summary(
    fun.data = "mean_cl_boot",
    geom = "errorbar",
    aes(ymax = after_stat(y), ymin = after_stat(y)),
    position = position_dodge(0.7),
    width = 0.15,
    color = "black"
  ) +
  xlab("") +
  ylab("Log2 MS signal") +
  theme_minimal(base_size = 16)

library(ggrepel)
ggplot(data = limma, aes(x = logFC, y = -log(adj.P.Val, 10))) +
  geom_point(
    aes(color = Sph.vs.Adh, alpha = Sph.vs.Adh),
    shape = 16,
    size = 2
  ) +
  scale_color_manual(values = c("#FC8D62","#66C2A5","#B3B3B3")) +
  scale_alpha_manual(values = c(0.9, 0.9, 0.5)) +
  geom_text_repel(
    data = PDGFR,
    aes(
      x = logFC,
      y = -log(adj.P.Val, 10),
      label = gene_site
    ),
    max.overlaps = Inf,
    size = 5,
    box.padding = 1, 
    min.segment.length = 0
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
  theme_light(base_size = 16)

# RoKAI analysis -------------------------------------------------------

#generate RoKAI input with logFC
df <- as.data.frame(cbind(limma$Protein, aap, limma$logFC))
colnames(df) <- c("Protein", "Position", "Quantification")
write.csv(df, file = "RoKAI_input_logFC.csv", row.names=TRUE)

#plot RoKAI output
RoKAI <- fread("kinase_table.csv")
RoKAI <- RoKAI[RoKAI$FDR<=0.05,]
RoKAI <- RoKAI[RoKAI$NumSubs>=3,]

ggplot(data = RoKAI, aes(x = reorder (Name,-ZScore), y = ZScore)) +
  geom_bar(stat = "identity", aes(fill = Activity), color = "black") +
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
  
