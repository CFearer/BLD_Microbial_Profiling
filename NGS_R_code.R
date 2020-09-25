#Code obtained and modified from Antonino Malacrinò - The Ohio STate University, Department of EEOB

library("ape")
library("car")
library("data.table")
library("DESeq2")
library("dplyr")
library("emmeans")
library("ggplot2")
library("ggpmisc")
library("ggpubr")
library("grid")
library("lme4")
library("phyloseq")
library("picante")
library("plyr")
library("Rmisc")
library("scales")
library("tibble")
library("vegan")
library("RAM")
library("splitstackshape")


#Import Metadata from QIIME2
data <- read.table("feature-table.txt", header = T, sep = "\t", row.names = 1, check.names=FALSE)
metadata <- read.table("MANIFEST.txt", header = T, sep = "\t", row.names = 1)
taxa <- read.table("taxonomy.tsv", header = T, sep = "\t", row.names = 1)
tree_Q <- read.tree("tree.nwk")

#Rooting the phlogenetic tree
tree_Q <- root(tree_Q, 1, resolve.root = T)

#Tests to see if all the column data from the OTU table are in the metadata; All TRUE=yes
colnames(data) %in% rownames(metadata)
 
#Splitting the Taxon data frame 
taxa$Taxon <- cSplit(data.frame(taxa), 'Taxon', ";")
taxa <- as.data.frame(as.matrix(taxa))
names(taxa) <- c("Kingdom", "Phylum",       "Class",           "Order",          "Family",         "Genus",            "Species")

#Can take any number of phyloseq objects and/or phyloseq components, and attempts to combine them into one larger phyloseq object. 
#This is most-useful for adding separately-imported components to an already-created phyloseq object.
GM <- merge_phyloseq(otu_table(data, taxa_are_rows=T), 
                     sample_data(metadata),
                     tree_Q,
                     tax_table(as.matrix(taxa))
)
GM

#It is intended to allow subsetting complex experimental objects with one function call.
GM <- subset_samples(GM, symptom_type == "asymptomatic" | symptom_type == "banding" | symptom_type == "crinkling" | symptom_type == "naive") 



##----------FILTERS-------------------------
#Filters mitochondria (only do this for ITS library ) and chloroplast DNA
GM <- subset_taxa(GM, Order !="o__Chloroplast")
GM <- subset_taxa(GM, Order !="o__Rickettsiales")
#Filtering samples without at least 1000 reads
GM <- prune_samples(sample_sums(GM)>=1000, GM)
GM <- filter_taxa(GM, function (x) {sum(x > 0) > 1}, prune=TRUE)
GM

##----PERMANOVA--------------------------------------------------------------------------
sampledf <- data.frame(sample_data(GM))
#Calculating the UniFrac Distance for all sample pairs
dist.mat <- phyloseq::distance(GM, method = "unifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~ symptom_type, data = sampledf, permutations = perm, strata = (location + time))
pmv

##----NMDS--------------------------------------------------------------------------
cap_ord <- ordinate(physeq = GM, method = "NMDS", distance = "unifrac", 
                    formula = ~ symptom_type)

cap_plot_source <- plot_ordination(physeq = GM, ordination = cap_ord, axes = c(1,2), color = "symptom_type") +
  stat_ellipse(mapping = aes(linetype = symptom_type),
               level = 0.95,
               type = "norm",
               show.legend=T) +
  theme_bw(base_size = 15) +
  geom_point(size = 5) +
  theme(legend.background = element_rect(color = "black", size = 0.3, linetype = "solid"),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank())
cap_plot_source


##----Alpha diversity ----
otus <- as.data.frame(otu_table(GM))
otus <- t(otus)
comm.pd <- pd(otus, tree_Q, include.root = TRUE)
diversity <- estimate_richness(GM, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
div2 <- cbind(sample_data(GM), diversity)
div2 <- cbind(div2, comm.pd)
div2$Evenness <- div2$Shannon/log(div2$Observed)

#Shannon Diversity 
Shannon_plot1 <- ggplot(div2, aes(x = time, y = Shannon, fill = time)) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  theme_bw(base_size = 12) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position="none",
        legend.background = element_rect(color = "black", size = 0.3, linetype = "solid"),
        legend.text = element_text(size = 12),
        panel.background = element_rect(size = 0.5, linetype = "solid"),
        strip.background = element_rect(fill = "white")) +
  ylim(0, 3)
Shannon_plot1


#Phylogenetic Diversity
Phylogenetic_plot1 <- ggplot(div2, aes(x = time, y = PD, fill = time)) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  theme_bw(base_size = 12) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position="none",
        legend.background = element_rect(color = "black", size = 0.3, linetype = "solid"),
        legend.text = element_text(size = 12),
        panel.background = element_rect(size = 0.5, linetype = "solid"),
        strip.background = element_rect(fill = "white")) +
  ylim(0, 6)
Phylogenetic_plot1

#Simpson Diversity
Simpson_plot1 <- ggplot(div2, aes(x = time, y = Simpson, fill = time)) +
  geom_boxplot(outlier.colour="black", notch=F, outlier.shape=NA) +
  theme_bw(base_size = 12) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.grid = element_blank(),
        legend.position="none",
        legend.background = element_rect(color = "black", size = 0.3, linetype = "solid"),
        legend.text = element_text(size = 12),
        panel.background = element_rect(size = 0.5, linetype = "solid"),
        strip.background = element_rect(fill = "white")) +
  ylim(0, 1)
Simpson_plot1

#Merge all plots
myplots <- list(Shannon_plot1, Phylogenetic_plot1, Simpson_plot1) #let's assume you have 4 plots, you can expand or reduce this list

p <- ggarrange(plotlist = myplots, #this is the list you defined above
               ncol = 3, nrow = 1, #this defines the number of columns and rows. For example, if you have 4 plots and you want to arrange them in a 2x2 grid, you set this number to 2 and 2
               align = "hv", #this is the alignment method. Leave it this way
               widths = 1, heights = 1, #this is the proportion of each plot, you need to modify this if you want your plots to have different proportion. In you case, leave it as it is.
               #labels= c("A)","B)", "C)"), #here you can specify a label for each plot, you can write whatever you like. The order will follow the one you specify in myplots
               common.legend = F) #this can be TRUE or FALSE depending if you want a common lenged for you graphs or not
p

ggsave(p, #this is the graph you just generated
       filename = "GinkgoBox.pdf", #this is the filename
       dpi = 600, device = cairo_pdf, family="Arial Unicode MS", #do not change these
       units = "in", width = 9, height = 6) #here you can set the plot dimension in inches, you can try until you do not get something you like


#Repeated Measures Test
model <- lmer(Shannon ~ time + (1|age) + (1|sex), data = div2)
summary(model)

model <- lmer(Shannon ~ location * (1|symptom.type) * (1|time), data = div2)
Anova(model)
emmeans(model, list(pairwise ~ location), adjust = "tukey")

model <- lmer(PD ~ symptom.type * (1|location) * (1|time), data = div2)
Anova(model)
emmeans(model, list(pairwise ~ location), adjust = "tukey")

model <- lmer(PD ~ location + (1|symptom.type) + (1|time), data = div2)
summary(model)

model <- lmer(Simpson ~ symptom.type * (1|location) * (1|time), data = div2)
Anova(model)
emmeans(model, list(pairwise ~ location), adjust = "tukey")

model <- lmer(Simpson ~ location * (1|symptom.type) * (1|time), data = div2)
Anova(model)
emmeans(model, list(pairwise ~ location), adjust = "tukey")




#DESeq2----Used to perform differential gene expression analyses
#Used to tell use which genes are differentially expressed between groupA and groupB; which OTUs  are differentially abundant between the two groups

#Transforms phyloseq object into deseq object
cds = phyloseq_to_deseq2(GM, ~ symptom.type)
gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(cds), 1, gm_mean)
cds$group <- factor(paste0(cds$symptom.type))
#explains how the dataset is structured
design(cds) <- ~ group
diagdds = estimateSizeFactors(cds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE, contrast = c("group", "crinkling", "naive"))
#log2foldchange tells us how much that OTU changes in which direction
#pajd tells us if is a signficant change
#Double filter to determine which samples are more abundant in symptomatic samples
sigtab = res[which(res$padj < 0.01 & res$log2FoldChange > 0), ]
sigtab1 = as(sigtab, "data.frame")

res = results(diagdds, cooksCutoff = FALSE, contrast = c("group", "banding", "naive"))
sigtab = res[which(res$padj < 0.01 & res$log2FoldChange > 0), ]
sigtab2 = as(sigtab, "data.frame")

l1 <- list(crinkling = row.names(sigtab1),  
           banding = row.names(sigtab2))
v1 <- Venn(l1)
p <- plot(v1, doWeights = FALSE, type = "circles")
p

sigtab1.b <- setDT(sigtab1, keep.rownames = TRUE)[]
sigtab2.b <- setDT(sigtab2, keep.rownames = TRUE)[]
taxa2 <- setDT(taxa, keep.rownames = TRUE)[]


diff.ASVs1 <- merge(sigtab1.b, taxa2, by = "rn")
diff.ASVs2 <- merge(sigtab2.b, taxa2, by = "rn")



#Create a Venn Diagram to see the overrepresented OTUs; Vennerable package currently not updated so is not working
#l1 <- list(ohio = row.names(sigtab1),  
           #columbus = row.names(sigtab2))
#v1 <- Venn(l1)
#p <- plot(v1, doWeights = FALSE, type = "circles")
#p



#Testing Normality 
ggqqplot(div2$PD)
ggdensity(div2$PD,
          main = "Density plot of Ginkgo",
          xlab = "Phylogenetic Diversity")
shapiro.test(div2$PD)



