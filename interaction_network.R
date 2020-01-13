#run: env R_MAX_VSPACE=10000Gb Rscript interaction_network.r
###########
#LIBRARIES
###########
library(SpiecEasi)
library(Matrix)
library(igraph)

##########
#ENV SETUP
##########
PATH="/home/lymelab/aem/tortoise"
setwd(PATH)

###########
#LOAD DATA
###########
#load sample data
seqtab.filtered <- read.table("sequence_table.16s.filtered.txt", header=T, row.names=1)

###############################
#MICROBIAL CO-OCCURANCE NETWORK
###############################
# solar samples only
solar.seqtab <- seqtab.filtered[grepl("S", rownames(seqtab.filtered)),]
solar.seqtab <- as.matrix(solar.seqtab[,colSums(solar.seqtab != 0) > 0])

se.mb <- spiec.easi(solar.seqtab, method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
se.gl <- spiec.easi(solar.seqtab, method='glasso', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
se.sparcc <- sparcc(solar.seqtab)
## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(se.sparcc$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
## Create igraph objects
ig.mb     <- adj2igraph(getRefit(se.mb))
ig.gl     <- adj2igraph(getRefit(se.gl))
ig.sparcc <- adj2igraph(sparcc.graph)

vsize    <- rowMeans(clr(solar.seqtab, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

pdf("test.pdf")
par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
dev.off()



