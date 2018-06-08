library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
# args <- c("human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt")
input <- fread(args[1], sep = "\t", header = T)
scores <- data.matrix(input[,-1])
rownames(scores) <- as.matrix(input[,1])
nullcount <- apply(scores, 1, function(x) {length(x[x==0])})
scores.cons <- scores[nullcount < 500, ]
write.table(data.frame(rownames(scores.cons), scores.cons), 
	file = paste(args[1], "notnull.tsv", sep = "."), 
	sep = "\t", quote = F, row.names = F)