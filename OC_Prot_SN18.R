
# LOAD DATA IN R ----------------------------------------------------------

#set working directory
setwd("N:/SUN-CPR-proteomics_jvo/GIULIA/B_COLLAB_EPIC-XS/UgoC/UgoC_OC_DIA_Prot_SN18")

#load OC proteome DIA data
library(data.table)
mydata <- fread("20231223_113739_20231025_GF_OC_Prot_dDIA_Report.tsv")
mydata <- as.data.frame(mydata)
#6886 ids

# Remove contaminants -----------------------------------------------------
unique(mydata$PG.FastaFiles)
mydata <- mydata[mydata$PG.FastaFiles == "Human_UP_2023-03",]
mydata <- mydata[,-4]
#6782 ids

# Remove rows wo gene name ------------------------------------------------
mydata <- mydata[mydata$PG.Genes != "",]
#6780 ids

# Rename sample names -----------------------------------------------------
colnames(mydata) <- sub("^.*?\\Prot_", "", colnames(mydata))
colnames(mydata) <- sub("F2_", "", colnames(mydata))
colnames(mydata) <- sub(".raw.PG.Log2Quantity", "", colnames(mydata))

# id count ----------------------------------------------------------------
id <- as.data.frame(nrow(mydata) - colSums(is.na(mydata[, c(4:43)])))
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

# boxplot -----------------------------------------------------------------
melt <- reshape2::melt(mydata, id = c("PG.ProteinGroups", "PG.Genes", "PG.ProteinDescriptions"))
melt <- na.omit(melt)
melt$variable <- as.character(melt$variable)

melt$r <- substr(melt$variable, nchar(melt$variable)-1, nchar(melt$variable))

melt$patient <- sub("\\_.*", "", melt$variable)
melt$patient <-
  factor(melt$patient,
         levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))
melt$patientr <- paste0(melt$patient, "_", melt$r)

melt$group <- sub("^.*?\\_", "", melt$variable)
melt$group <- sub("\\_.*", "", melt$group)

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
sampleMed <- apply(mydata[,4:43], 2, median, na.rm = TRUE) 
meanMed <- mean(sampleMed, na.rm = TRUE)
norm <- t(t(mydata[,4:43]) - sampleMed)
norm <- as.data.frame(norm + meanMed)
norm <- cbind(mydata[,1:3], norm)

# boxplot after normalization ---------------------------------------------
melt.norm <- reshape2::melt(norm, id = c("PG.ProteinGroups", "PG.Genes", "PG.ProteinDescriptions"))
melt.norm <- na.omit(melt.norm)
melt.norm$variable <- as.character(melt.norm$variable)

melt.norm$r <- substr(melt.norm$variable, nchar(melt.norm$variable)-1, nchar(melt.norm$variable))

melt.norm$patient <- sub("\\_.*", "", melt.norm$variable)
melt.norm$patient <-
  factor(melt.norm$patient,
         levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))
melt.norm$patientr <- paste0(melt.norm$patient, "_", melt.norm$r)

melt.norm$group <- sub("^.*?\\_", "", melt.norm$variable)
melt.norm$group <- sub("\\_.*", "", melt.norm$group)

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
cat <- sub("_0.*", "", colnames(norm[,4:43]))
cat <- factor(cat)
average <-
  as.data.frame(t(apply(norm[, 4:43], 1, function(x) {
    tapply(x, cat, function (y)
      mean(y, na.rm = TRUE))
  })))
average <- cbind(mydata[,1:3], average)

#long format
melt.average <- reshape2::melt(average, id = c("PG.ProteinGroups", "PG.Genes", "PG.ProteinDescriptions"))
melt.average <- na.omit(melt.average)
melt.average$variable <- as.character(melt.average$variable)

melt.average$r <- substr(melt.average$variable, nchar(melt.average$variable)-1, nchar(melt.average$variable))

melt.average$patient <- sub("\\_.*", "", melt.average$variable)
melt.average$patient <- factor(melt.average$patient, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"))
melt.average$patientr <- paste0(melt.average$patient, "_", melt.average$r)

melt.average$group <- sub("^.*?\\_", "", melt.average$variable)
melt.average$group <- sub("\\_.*", "", melt.average$group)

# PCA ---------------------------------------------------------------------

#define categories
patient <- Id$patient
group <- factor(Id$group)

#remove missing values
rownames(average) <- average$PG.ProteinGroups
nona <- na.omit(average[,4:23])

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
col_adh <- seq(4, ncol(average), by = 2)
col_sph <- seq(5, ncol(average), by = 2)

Sph.only <-
  average[rowSums(is.na(average[, col_adh])) == 10 &
            rowSums(is.na(average[, col_sph])) <=5, ]
Sph.only$presence <- "Sph"

Adh.only <-
  average[rowSums(is.na(average[, col_adh])) <= 5 &
            rowSums(is.na(average[, col_sph])) == 10, ]
Adh.only$presence <- "Adh"

only.in.one.group <- rbind(Sph.only, Adh.only)

average <-
  merge(average, 
        only.in.one.group[, c(1, 24)], 
        by = "PG.ProteinGroups", 
        all.x = T)

# Proliferation markers ---------------------------------------------------
pi <- fread("proliferation.index.txt")
pi <- as.data.frame(pi)
average$PG <- gsub(";.*", "", average$PG.ProteinGroups)
pi <- left_join(pi, average, by = c("Uniprot" = "PG"))
pindex <- apply(pi[,5:24], 2, median, na.rm = TRUE)
#subtract odd values to even values
pindex.fc <- pindex[seq(2, length(pindex), by = 2)] - pindex[seq(1, length(pindex), by = 2)]

#one sample t test
t.test(pindex.fc, mu = 0, alternative = "two.sided")
#p = 0.3904

#long format
pi.long <- reshape2::melt(pi[,c(1,5:24)], id = "Uniprot")
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
  average[rowSums(is.na(average[, col_adh])) <= 7 &
            rowSums(is.na(average[, col_sph])) <= 7, ]
#6558 ids

#perform Limma paired analysis
library(limma)
design <- model.matrix(~ patient + group)
colnames(design) <- gsub("patient", "", colnames(design))
colnames(design) <- gsub("group", "", colnames(design))
fit <- lmFit(filtered[,4:23], design)
fit2 <- eBayes(fit)
tab <- topTable(fit2,
                n = Inf,
                coef = "Sph",
                sort.by = "none")
sum(tab$adj.P.Val<=0.05,na.rm = T)
sum(tab$adj.P.Val<=0.01,na.rm = T)
limma <- cbind(filtered[,1:3], tab)

#calculate log2 fc cutoff

2 * sd(limma$logFC [limma$logFC > 0], na.rm = T)
# 0.8939356
sum(limma$adj.P.Val < 0.05 & limma$logFC > 0.89, na.rm = T)
sum(limma$adj.P.Val < 0.01 & limma$logFC > 0.89, na.rm = T)

2 * sd(limma$logFC [limma$logFC < 0], na.rm = T)
# 0.698736
sum(limma$adj.P.Val < 0.05 & limma$logFC < -0.7, na.rm = T)
sum(limma$adj.P.Val < 0.01 & limma$logFC < -0.7, na.rm = T)

#add categorical column for up and down regulation
limma <- limma %>%
  mutate(
    Sph.vs.Adh = case_when(
      logFC > 0.89 & adj.P.Val <= 0.01 ~ "Up-regulated",
      logFC < -0.7 & adj.P.Val <= 0.01 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  )

#export limma output
nor <- norm[,4:43]
colnames(nor) <- paste0(colnames(nor), "_", "norm")

mydata <- merge(
  x = cbind(mydata, nor, average[, c(4:24)]),
  y = limma[-c(2:3)],
  by = "PG.ProteinGroups",
  all.x = TRUE
)

write.table(mydata,
            file = "OC_Prot_Limma.txt",
            row.names = F,
            sep = "\t")

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

# CORRELATION RNA-PROTEIN -------------------------------------------------
rna <- fread("RNAdataOC.txt")

rp <- left_join(rna, tab, by=c("Symbol"="Genes"))
rp <- na.omit(rp)
rp$`STEM Fold-change` <- log(rp$`STEM Fold-change`,2)
cor(rp$`STEM Fold-change`, rp$logFC)
cor.test(rp$`STEM Fold-change`, rp$logFC, method=c("pearson"))

ggplot(rp, aes(x = `STEM Fold-change`, y = logFC)) +
  geom_point(
    shape = 21,
    colour = "black",
    fill = "bisque2",
    size = 5
  ) +
  geom_smooth(method = lm,
              linetype = "dashed",
              color = "black") +
  labs(title = "Pearson correlation = 0.36, p-value = 2.277e-13, n = 384") +
  xlab("Affymetrix_Log2FC") +
  ylab("Protein_Log2FC") +
  theme_light(base_size = 16)

# PLOT PDGFR --------------------------------------------------------------

PDGFR <- melt.average[grepl("PDGFR", melt.average$PG.Genes, fixed = T),]

ggplot(PDGFR, aes(x = PG.Genes, y = value, fill = group)) +
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

# GO ENRICHMENT ANALYSIS --------------------------------------------------
library(topGO)
library(org.Hs.eg.db)

#remove everything after ";" in mydata$PG.Genes and mydata$PG.ProteinGroups
mydata$Gene <- gsub(";.*", "", mydata$PG.Genes)
mydata$PG <- gsub(";.*", "", mydata$PG.ProteinGroups)

#get GO annotations
BP <-
  annFUN.org(
    whichOnto = "BP",
    feasibleGenes = NULL,
    mapping = "org.Hs.eg.db",
    ID = "symbol"
  )

CC <-
  annFUN.org(
    whichOnto = "CC",
    feasibleGenes = NULL,
    mapping = "org.Hs.eg.db",
    ID = "symbol"
  )

MF <-
  annFUN.org(
    whichOnto = "MF",
    feasibleGenes = NULL,
    mapping = "org.Hs.eg.db",
    ID = "symbol"
  )

#define the background gene list
geneUniverse <- unique(mydata$Gene)

#define the list of genes of interest
up <- filter(mydata, Sph.vs.Adh == "Up-regulated" | presence == "Sph")
up <- as.character(up$Gene)

down <- filter(mydata, Sph.vs.Adh == "Down-regulated" | presence == "Adh")
down <- as.character(down$Gene)

#tell TopGO where the interesting genes appear in the 'geneUniverse' vector
geneListu <- factor(as.integer(geneUniverse %in% up))
names(geneListu) <- geneUniverse

geneListd <- factor(as.integer(geneUniverse %in% down))
names(geneListd) <- geneUniverse

#create topGOdata object
BPu <- new(
  "topGOdata",
  description = "BP enrichment analysis up",
  ontology = "BP",
  allGenes = geneListu,
  annot = annFUN.GO2genes,
  GO2genes = BP,
  nodeSize = 5
)

BPd <- new(
  "topGOdata",
  description = "BP enrichment analysis down",
  ontology = "BP",
  allGenes = geneListd,
  annot = annFUN.GO2genes,
  GO2genes = BP,
  nodeSize = 5
)

CCu <- new(
  "topGOdata",
  description = "CC enrichment analysis up",
  ontology = "CC",
  allGenes = geneListu,
  annot = annFUN.GO2genes,
  GO2genes = CC,
  nodeSize = 10
)

CCd <- new(
  "topGOdata",
  description = "CC enrichment analysis down",
  ontology = "CC",
  allGenes = geneListd,
  annot = annFUN.GO2genes,
  GO2genes = CC,
  nodeSize = 10
)

MFu <- new(
  "topGOdata",
  description = "MF enrichment analysis up",
  ontology = "MF",
  allGenes = geneListu,
  annot = annFUN.GO2genes,
  GO2genes = MF,
  nodeSize = 10
)

MFd <- new(
  "topGOdata",
  description = "MF enrichment analysis down",
  ontology = "MF",
  allGenes = geneListd,
  annot = annFUN.GO2genes,
  GO2genes = MF,
  nodeSize = 10
)

#perform Fisher exact test
FisherBPu <- runTest(BPu, algorithm="weight01", statistic="fisher")
FisherBPd <- runTest(BPd, algorithm="weight01", statistic="fisher")
FisherCCu <- runTest(CCu, algorithm="weight01", statistic="fisher")
FisherCCd <- runTest(CCd, algorithm="weight01", statistic="fisher")
FisherMFu <- runTest(MFu, algorithm="weight01", statistic="fisher")
FisherMFd <- runTest(MFd, algorithm="weight01", statistic="fisher")

#generate a table of results
allGOBPu <- usedGO(BPu)
allResBPu <-
  GenTable(
    BPu,
    weightFisher = FisherBPu,
    orderBy = "FisherBPu",
    ranksOf = "weightFisher",
    topNodes = length(allGOBPu)
  )
allResBPu$weightFisher <- as.numeric(allResBPu$weightFisher)

allGOBPd <- usedGO(BPd)
allResBPd <-
  GenTable(
    BPd,
    weightFisher = FisherBPd,
    orderBy = "FisherBPd",
    ranksOf = "weightFisher",
    topNodes = length(allGOBPd)
  )
allResBPd$weightFisher <- as.numeric(allResBPd$weightFisher)

allGOCCu <- usedGO(CCu)
allResCCu <-
  GenTable(
    CCu,
    weightFisher = FisherCCu,
    orderBy = "FisherCCu",
    ranksOf = "weightFisher",
    topNodes = length(allGOCCu)
  )
allResCCu$weightFisher <- as.numeric(allResCCu$weightFisher)
allResCCu$weightFisher[1:2] <- 1e-20

allGOCCd <- usedGO(CCd)
allResCCd <-
  GenTable(
    CCd,
    weightFisher = FisherCCd,
    orderBy = "FisherCCd",
    ranksOf = "weightFisher",
    topNodes = length(allGOCCd)
  )
allResCCd$weightFisher <- as.numeric(allResCCd$weightFisher)

allGOMFu <- usedGO(MFu)
allResMFu <-
  GenTable(
    MFu,
    weightFisher = FisherMFu,
    orderBy = "FisherMFu",
    ranksOf = "weightFisher",
    topNodes = length(allGOMFu)
  )
allResMFu$weightFisher <- as.numeric(allResMFu$weightFisher)

allGOMFd <- usedGO(MFd)
allResMFd <-
  GenTable(
    MFd,
    weightFisher = FisherMFd,
    orderBy = "FisherMFd",
    ranksOf = "weightFisher",
    topNodes = length(allGOMFd)
  )
allResMFd$weightFisher <- as.numeric(allResMFd$weightFisher)

#performing BH correction on Fisher's p values
allResBPu$p.adj <- p.adjust(as.numeric(allResBPu$weightFisher), method = "BH")
allResBPd$p.adj <- p.adjust(allResBPd$weightFisher, method = "BH")
allResCCu$p.adj <- p.adjust(allResCCu$weightFisher, method = "BH")
allResCCd$p.adj <- p.adjust(allResCCd$weightFisher, method = "BH")
allResMFu$p.adj <- p.adjust(allResMFu$weightFisher, method = "BH")
allResMFd$p.adj <- p.adjust(allResMFd$weightFisher, method = "BH")

#calculate enrichment %
allResBPu <- allResBPu %>% 
  mutate(EnrichmentS = Significant / Annotated)

allResBPd <- allResBPd %>% 
  mutate(EnrichmentS = Significant / Annotated)

allResCCu <- allResCCu %>% 
  mutate(EnrichmentS = Significant / Annotated)

allResCCd <- allResCCd %>% 
  mutate(EnrichmentS = Significant / Annotated)

allResMFu <- allResMFu %>% 
  mutate(EnrichmentS = Significant / Annotated)

allResMFd <- allResMFd %>% 
  mutate(EnrichmentS = Significant / Annotated)

#export topGO output
allResBPu$table <- "BPu"
allResBPd$table <- "BPd"
allResCCu$table <- "CCu"
allResCCd$table <- "CCd"
allResMFu$table <- "MFu"
allResMFd$table <- "MFd"

GO <- rbind(allResBPu,allResBPd,allResCCu,allResCCd,allResMFu,allResMFd)
write.table(GO, file = "GO_enrichment.txt",row.names=F, sep="\t")

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
resBPu <- allResBPu[allResBPu$p.adj<=0.05,]
resBPu$Term <- tools::toTitleCase(resBPu$Term)

ggplot(resBPu, aes(
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
resBPd <- allResBPd[allResBPd$p.adj<=0.05,]
resBPd$Term <- tools::toTitleCase(resBPd$Term)

ggplot(resBPd, aes(
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
resCCu <- allResCCu [allResCCu$p.adj<=0.05,]
resCCu$Term <- tools::toTitleCase(resCCu$Term)

ggplot(resCCu, aes(
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
resCCd <- allResCCd [allResCCd$p.adj<=0.05,]
resCCd$Term <- tools::toTitleCase(resCCd$Term)

ggplot(resCCd, aes(
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
resMFu <- allResMFu [allResMFu$p.adj<=0.05,]
resMFu$Term <- tools::toTitleCase(resMFu$Term)

ggplot(resMFu, aes(
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
resMFd <- allResMFd [allResMFd$p.adj<=0.05,]
resMFd$Term <- tools::toTitleCase(resMFd$Term)

ggplot(resMFd, aes(
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

#define background
geneUniverseuniprot <- as.character(mydata$PG)

#define the list of genes of interest
upPG <- filter(mydata, Sph.vs.Adh == "Up-regulated" | presence == "Sph")
upPG <- as.character(upPG$PG)

downPG <- filter(mydata, Sph.vs.Adh == "Down-regulated" | presence == "Adh")
downPG <- as.character(downPG$PG)

#perform enrichment
kegg.u <- enrichKEGG(
  gene = upPG,
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
keggu <- as.data.frame(kegg.u)
keggu$updown <- "up"

kegg.d <- enrichKEGG(
  gene = downPG,
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
keggd <- as.data.frame(kegg.d)
keggd$updown <- "down"

keggd_sel <- as.data.frame(keggd[2:4,8])
colnames(keggd_sel) <- "id"
#split string in different columns when there is "/" without knowing how many columns you need
#create a vector with letter from a to z
letters <- LETTERS[1:26]
keggd_sel <- separate(keggd_sel, id, LETTERS[1:26], sep = "/")
keggd_sel <- as.data.frame(t(keggd_sel))
#create a vector by concatenating the 3 different columns
keggd_sel <- na.omit(as.data.frame(c(keggd_sel$V1, keggd_sel$V2, keggd_sel$V3)))
colnames(keggd_sel) <- "id"
keggd_sel <- unique(keggd_sel)
keggd_sel <- left_join(keggd_sel, mydata[,104:113], by = c("id" = "PG"))
write.table(keggd_sel, file = "KEGG_down_for_string.txt", sep = "\t", row.names = F)

kegg <- rbind(keggu, keggd)

#calculate count & enrichment %
kegg2 <- separate(kegg, GeneRatio, c("Count", "Count.Sign"))
kegg2$Count <- as.numeric(kegg2$Count)
kegg2 <- separate(kegg2, BgRatio, c("CountGO", "Count.Bg"))
kegg2$CountGO <- as.numeric(kegg2$CountGO)
kegg2 <- kegg2 %>% mutate(EnrichmentS = Count / CountGO)

#export KEGG enrichment output
write.table(kegg2,
            file = "OC_Prot_KEGG_enrichment.txt",
            row.names = F,
            sep = "\t")

#plot data
ggplot(kegg2[kegg2$updown == "up",], aes(
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

ggplot(kegg2[kegg2$updown == "down",], aes(
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

#define background genes after converting symbol to entreid
ids <-
  bitr(geneUniverse,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = "org.Hs.eg.db")
geneUni <- as.character(ids$ENTREZID)

#define genes of interest after converting symbol to entreid
upe <- bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
upe <- as.character(upe$ENTREZID)

downe <- bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
downe <- as.character(downe$ENTREZID)

#perform enrichment
library(ReactomePA)
react.u <- enrichPathway(
  gene = upe,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  universe = geneUni,
  minGSSize = 10, 
  maxGSSize = 500,
  readable = F)
reactu <- as.data.frame(react.u)
reactu$updown <- "up"

react.d <- enrichPathway(
  gene = downe,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  universe = geneUni,
  minGSSize = 10, 
  maxGSSize = 500,
  readable = F)
reactd <- as.data.frame(react.d)
reactd$updown <- "down"

react <- rbind(reactu, reactd)

#calculate count & enrichment %
react2 <- separate(react, GeneRatio, c("Count", "Count.Sign"))
react2$Count <- as.numeric(react2$Count)
react2 <- separate(react2, BgRatio, c("CountGO", "Count.Bg"))
react2$CountGO <- as.numeric(react2$CountGO)
react2 <- react2 %>% mutate(EnrichmentS = Count / CountGO)

#export Reactome enrichment output
write.table(react2,
            file = "OC_Prot_React_enrichment.txt",
            row.names = F,
            sep = "\t")

#plot data
ggplot(react2[react2$updown == "up",], aes(
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

ggplot(react2[react2$updown == "down",], aes(
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

#plot React data
#dotplot(react.u, showCategory = 20) + ggtitle("Upregulated REACTOME pathways Sph vs Adh")
#dotplot(rd) + ggtitle("Downregulated REACTOME pathways Sph vs Adh")

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









