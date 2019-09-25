################
#Random forest 
################
library(plyr)
library(randomForest)
library(rfUtilities)

otu_table <- read.table("sequence_table.16s.filtered.txt", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
otu_table <- t(otu_table)
metadata <- read.table("map.txt", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

otu_nonzero_counts <- apply(otu_table, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_rare_removed <- remove_rare(table=otu_table, cutoff_pro=0.1)
otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)  

# can model predict tortoise group? random subset of 31 each (all solar, random 31 necropsy)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))  
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "tortoise_group"]
otu_table_nec <- subset(otu_table_scaled_var, var=="necropsy")
otu_table_sol <- subset(otu_table_scaled_var, var=="solar")
randomRows = function(df,n){
   return(df[sample(nrow(df),n),])
}
otu_table_nec31 <- randomRows(otu_table_nec, 31)
#bind together
subset_otu_table <- rbind(otu_table_nec31, otu_table_sol)
#reset factors
subset_otu_table$var <- factor(subset_otu_table$var)
# run random forest
set.seed(151)  
rf_tgroup <- randomForest(x=subset_otu_table[,1:(ncol(subset_otu_table)-1)] , y=subset_otu_table[ , ncol(subset_otu_table)] , ntree=10000, importance=TRUE, proximities=TRUE )  
rf_tgroup

# Call:
#  randomForest(x = subset_otu_table[, 1:(ncol(subset_otu_table) -      1)], y = subset_otu_table[, ncol(subset_otu_table)], ntree = 10000,      importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 18

#         OOB estimate of  error rate: 8.06%
# Confusion matrix:
#          necropsy solar class.error
# necropsy       30     1  0.03225806
# solar           4    27  0.12903226

rf.significance(x=rf_tgroup, xdata=subset_otu_table[,1:(ncol(subset_otu_table)-1)], nperm=1000, ntree=1000)

# Number of permutations:  1000
# p-value:  0
# Model signifiant at p = 0
#    Model OOB error:  0.08064516
#    Random OOB error:  0.516129
#    min random global error: 0.2580645
#    max random global error:  0.7580645
#    min random within class error: 0.5
#    max random within class error:  0.5

# what about tortoise sex?
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "sex"]
otu_table_f <- subset(otu_table_scaled_var, var=="female")
otu_table_m <- subset(otu_table_scaled_var, var=="male")
dim(otu_table_f)
# [1]  27 356
dim(otu_table_m)
# [1]  50 356
# randomly subset males to 27 individuals
otu_table_m27 <- randomRows(otu_table_m, 27)
subset_otu_table <- rbind(otu_table_f, otu_table_m27)
subset_otu_table$var <- factor(subset_otu_table$var)

set.seed(151)  
rf_sex <- randomForest(x=subset_otu_table[,1:(ncol(subset_otu_table)-1)] , y=subset_otu_table[ , ncol(subset_otu_table)] , ntree=10000, importance=TRUE, proximities=TRUE )  
rf_sex

# Call:
#  randomForest(x = subset_otu_table[, 1:(ncol(subset_otu_table) -      1)], y = subset_otu_table[, ncol(subset_otu_table)], ntree = 10000,      importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 18

#         OOB estimate of  error rate: 35.19%
# Confusion matrix:
#        female male class.error
# female     16   11   0.4074074
# male        8   19   0.2962963

rf.significance(x=rf_sex, xdata=subset_otu_table[,1:(ncol(subset_otu_table)-1)], nperm=1000, ntree=1000)
# Number of permutations:  1000
# p-value:  0.015
# Model signifiant at p = 0.015
#    Model OOB error:  0.3518519
#    Random OOB error:  0.5185185
#    min random global error: 0.2407407
#    max random global error:  0.8148148
#    min random within class error: 0.4444444
#    max random within class error:  0.4444444

