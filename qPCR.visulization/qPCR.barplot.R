# setwd("~/Documents/GitHub/BioinformaticsDaily/qPCR.visulization")
# install packages before the first use
# install.packages("ggpubr")
# install.packages("data.table")

library(data.table)
library(ggpubr)

############  your inputs ############
inputfile <- "example.data.txt"
replicates <- 3
normalizedTo <- "GAPDH"
control <- "MCF7 Control gRNA"
# the width and height for pdf file
width <- 8
height <- 8
######################################

######################################
#######  DO NOT NEED TO CHANGE  ######
######################################
datax <- fread(inputfile, header = T, sep = "\t")
colnames(datax)[1:2] <- c("targets", "treatments")
mytargets <- unique(datax$targets)
mytargets <- mytargets[mytargets!=normalizedTo]
mytreats <- unique(datax$treatments)

xa <- c()
xb <- c()
xc <- c()
xd <- c()
for(i in 1:length(mytargets)){
    for(j in 1:length(mytreats)){
        for(k in 1:replicates){
            vals <- 2^(datax[targets==normalizedTo & treatments==mytreats[j], .SD, .SDcols = 2 + k] - 
                           datax[targets==mytargets[i] & treatments==mytreats[j], .SD, .SDcols = 2 + k])
            xa <- c(xa, mytargets[i])
            xb <- c(xb, mytreats[j])
            xc <- c(xc, paste("rep", k))
            xd <- c(xd, as.numeric(vals))
        }
    }
}
myresults <- data.table(targets = xa, treatments = xb, replicates = xc, results = xd)
mycompare <- compare_means(results ~ treatments, data = myresults, ref.group = control,
                           group.by = "targets", paired = TRUE, method = "t.test")
mymax <- max(myresults$results, na.rm = T)

p <- ggbarplot(myresults, x = "treatments", y = "results", add = c("dotplot"),
          color = "treatments", palette = "jco", facet.by = "targets",
          fill = "treatments", position = position_dodge(0.9),
          add.params = list(fill = "red", size = 0.3)) +
    stat_pvalue_manual(data = mycompare, label = "p.signif",
                      xmin = "group1", xmax = "group2",
                      y.position = seq(mymax, mymax * (1 + length(mytreats) / 20),
                                       mymax * length(mytreats) / 20))
# stat_compare_means(label = "p.signif", method = "t.test",
#                    ref.group = control)
myplot <- ggpar(p, xlab = "", ylab = paste("Normalized to", normalizedTo, sep = " "), x.text.angle = 90)
pdf(file = paste("ggpubr.barplot", inputfile, "pdf", sep = "."), width = width, height = height)
print(myplot)
dev.off()
