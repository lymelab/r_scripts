###########
#LIBRARIES
###########
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(compositions)
library(ape)
library(factoextra)
library(RColorBrewer)
library(philr)
library(phyloseq)
library(UpSetR)
library(reshape2)
library(vegan)
library(phylofactor)
library(ggtree)

##########
#ENV SETUP
##########
PATH="/Users/mann/github/triatomes"
setwd(PATH)
cols <- brewer.pal(6, "Set2")

###########
#LOAD DATA
###########
#load sample data
rawmetadata <- read_delim(file = file.path(PATH, "map.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries
seqtab.filtered <- read.table("sequence_table.16s.filtered.txt", header=T, row.names=1)
taxa <- read.table("assigntax/rep_set_fix_tax_assignments.txt", header=F, sep="\t", row.names=1)
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$SampleID) #record samples absent in either metadata or OTU table
notinraw <- setdiff(rawmetadata$SampleID, row.names(seqtab.filtered))
tree <- read.tree("rep_set.filt.tre") #load representative tree

#clean up rownames in filtered sequence table - ONLY RUN THIS ONCE IF NEEDED
# new.rownames <- row.names(seqtab.filtered)
# new.rownames <- gsub("_S.*.fastq.gz", "", new.rownames)
# new.rownames <- paste("Tri", new.rownames, sep="") # append characters to numerical samples names -- these can mess things up in the long run
# rownames(seqtab.filtered) <- new.rownames
# write.table(data.frame("row_names"=rownames(seqtab.filtered),seqtab.filtered),"sequence_table.16s.filtered.txt", row.names=F, quote=F, sep="\t")
# rownames(rawmetadata) <- rawmetadata$SampleID 
# new.samplenames <- paste("Tri", rownames(rawmetadata), sep="") #add character to sample names in metadata file
# rawmetadata$SampleID <- new.samplenames
# rownames(rawmetadata) <- rawmetadata$SampleID #set row names as the sample IDs for the metadata

################
#NORMALIZE DATA
################
#Centered log-ratio transformation
seqtab.clr <- clr(seqtab.filtered)
write.table(as.data.frame(seqtab.clr), "sequence_table.16s.clr.txt", quote=F, sep="\t", col.names=NA)

#############################
#AITCHISON DISTANCE (LINEAR)
#############################
# Heirarchical cluster dendrogram
hc <- hclust(dist(seqtab.clr), method="complete")
df2 <- data.frame(cluster=cutree(hc,6), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$cluster))) + geom_tile() + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("figs/aitchison_dendrogram.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

# Biplot
#have to read back in clr sequence table
comp <- read.table("sequence_table.16s.clr.txt", header=T, row.names=1)
pca <- prcomp(as.matrix(comp))
pdf("figs/aitchison_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/aitchison_pca.pdf")
fviz_pca_biplot(pca, habillage=merge$cluster, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F) + labs(title="Aitchison Distance", ylab="PC1 (9.3%)", xlab="PC2 (11.3%)") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-12.5,11.5) + ylim(-12.5,11.5)
dev.off()

##############################
#PHILR DISTANCE (PHYLOGENETIC)
##############################
#create phyloseq object from "seqtab.filtered", "rawmetadata", "taxa", and "tree"
tree <- read.tree("rep_set.root.tre")
ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa[1])), tree)

philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
is.rooted(phy_tree(philr.dat)) #check that tree is rooted
is.binary.tree(phy_tree(philr.dat)) #check that multichotomies are resolved in tree
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
otu.table <- otu_table(philr.dat)
tree <- phy_tree(philr.dat)
metadata <- sample_data(philr.dat)
tax <- tax_table(philr.dat)
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")

# Heirarchical cluster dendrogram
hc <- hclust(dist(philr.t), method="complete")
df2 <- data.frame(cluster=cutree(hc,5), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
write.table(df2, "philr_cluster.txt", quote=F, sep="\t", col.names=NA)

hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$cluster))) + geom_tile() + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("figs/philr_dendrogram.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

# PCA
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))
pdf("figs/philr_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/philr_pca.pdf")
fviz_pca_ind(pca, habillage=merge$cluster, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F) + labs(title="PhIRL Distance") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-50,50) + ylim(-50,50)
dev.off()

# PCA for supplemental figure
pdf("figs/philr_bloodmeal_pca.pdf")
fviz_pca_ind(pca, habillage=merge$Blood_meal, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F) + labs(title="PhIRL Distance") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-50,50) + ylim(-50,50)
dev.off()

pdf("figs/philr_machres1_pca.pdf")
fviz_pca_ind(pca, habillage=merge$MachRes1, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F) + labs(title="PhIRL Distance") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-50,50) + ylim(-50,50)
dev.off()

pdf("figs/philr_sex_pca.pdf")
fviz_pca_ind(pca, habillage=merge$Sex, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F) + labs(title="PhIRL Distance") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-50,50) + ylim(-50,50)
dev.off()

pdf("figs/philr_straintype_pca.pdf")
fviz_pca_ind(pca, habillage=merge$Strain_type, geom="point", addEllipses=T, ellipse.type="convex", pointshape=19, point.size=1, mean.point=F) + labs(title="PhIRL Distance") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-50,50) + ylim(-50,50)
dev.off()

############
#Upset plot
############
map <- as.matrix(read.table("map.txt", header=T, sep="\t", row.names=1))
merged <- merge(seqtab.filtered, map, by="row.names")
n <- ncol(seqtab.filtered) + 1
agg <- aggregate(merged[,2:n], by=list(merged$Sp_by_key), FUN=sum)
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
#convert to presence absence table -- ignore warnining message, still works
agg[agg>1] <- 1
#transpose again
agg <- setNames(data.frame(t(agg[,-1])), agg[,1])
#upsetR 
pdf("figs/upset.pdf", onefile=F)
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs per species", mb.ratio = c(0.55, 0.45))
dev.off()

####################
#Taxonomy barchart
####################
orderlist <- hc$labels
dat <- read.table("collapsed_simplified_taxonomy_for_plot.txt", header=T, sep="\t")
datmelt <- melt(dat)
datmelt$variable <- factor(datmelt$variable, levels=orderlist)
pdf("figs/taxonomy_barchart.pdf")
ggplot(datmelt, aes(fill=datmelt$taxonomy, x=datmelt$variable, y=datmelt$value)) + geom_bar(stat="identity", position="fill") + theme_classic() + theme(axis.text.x=element_text(angle=90))
dev.off()

#########
#Adonis 
#########
# is there a difference in microbial diversity across samples by some metadata category?
metadata <- as(sample_data(ps.dat), "data.frame")
adonis(philr.dist ~ Sex, data=metadata)
# Call:
# adonis(formula = philr.dist ~ Sex, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Sex        1     141.5 141.452  2.2226 0.02671   0.07 .
# Residuals 81    5154.9  63.641         0.97329
# Total     82    5296.4                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ Sp_by_key, data=metadata)

# Call:
# adonis(formula = philr.dist ~ Sp_by_key, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Sp_by_key  4     510.1 127.535  2.0784 0.09632  0.016 *
# Residuals 78    4786.3  61.362         0.90368
# Total     82    5296.4                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ MachRes1, data=metadata)

# Call:
# adonis(formula = philr.dist ~ MachRes1, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# MachRes1   1     249.3  249.28  4.0006 0.04707  0.006 **
# Residuals 81    5047.1   62.31         0.95293
# Total     82    5296.4                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ County, data=metadata)

# Call:
# adonis(formula = philr.dist ~ County, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# County    27    1743.3  64.568 0.99948 0.32915  0.496
# Residuals 55    3553.1  64.601         0.67085
# Total     82    5296.4                 1.00000

##############
#Phylofactor
##############
OTUTable <- as.matrix(t(seqtab.filtered))
filt.list <- colnames(OTUTable)
filtmap <- rawmetadata[rawmetadata$SampleID %in% filt.list,]
filtmap <- filtmap[match(filt.list, filtmap$SampleID),]
x <- as.factor(filtmap$MachRes1) # CHANGE ME to your variable of interest
tree <- phy_tree(philr.dat)
tax <- read.table("assigntax/rep_set_tax_assignments_phylofactor.txt", sep="\t", header=T)
common.otus <- which(rowSums(OTUTable>0)>10)
OTUTable <- OTUTable[common.otus,]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(OTUTable)))
PF <- PhyloFactor(OTUTable, tree, x, nfactors=3)
PF$Data <- PF$Data[PF$tree$tip.label,]
gtree <- pf.tree(PF,layout="rectangular")
pdf("figs/phylofactor_tree.pdf")
gtree$ggplot + geom_tiplab()
dev.off()

y <- t(PF$basis[,1]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor1_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1')
dev.off()

wilcox.test(dat[dat$V1 == "Negative",]$V2, dat[dat$V1 == "Positive",]$V2)

# 	Wilcoxon rank sum test with continuity correction

# data:  dat[dat$V1 == "Negative", ]$V2 and dat[dat$V1 == "Positive", ]$V2
# W = 553, p-value = 0.01745
# alternative hypothesis: true location shift is not equal to 0

y <- t(PF$basis[,2]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor2_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor2')
dev.off()

wilcox.test(dat[dat$V1 == "Negative",]$V2, dat[dat$V1 == "Positive",]$V2)

# 	Wilcoxon rank sum test with continuity correction

# data:  dat[dat$V1 == "Negative", ]$V2 and dat[dat$V1 == "Positive", ]$V2
# W = 569, p-value = 0.02599
# alternative hypothesis: true location shift is not equal to 0

y <- t(PF$basis[,3]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor3_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[3]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor3')
dev.off()

wilcox.test(dat[dat$V1 == "Negative",]$V2, dat[dat$V1 == "Positive",]$V2)

# 	Wilcoxon rank sum test with continuity correction

# data:  dat[dat$V1 == "Negative", ]$V2 and dat[dat$V1 == "Positive", ]$V2
# W = 545, p-value = 0.01419
# alternative hypothesis: true location shift is not equal to 0

PF$factors
#                               Group1                       Group2      ExpVar
# Factor 1 2 member Monophyletic clade 30 member Monophyletic clade 0.006977436
# Factor 2 3 member Monophyletic clade 27 member Paraphyletic clade 0.006067390
# Factor 3                         tip 26 member Paraphyletic clade 0.004715486
#                 F      Pr(>F)
# Factor 1 7.072276 0.009432920
# Factor 2 5.632511 0.019998950
# Factor 3 9.177399 0.003287377
