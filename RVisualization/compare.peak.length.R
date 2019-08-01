# args <- commandArgs(TRUE)
library(data.table)
library(ggpubr)
mydir <- "../bed_file/"
myfiles <- dir(path = mydir, pattern = ".bed")
factors <- unique(matrix(unlist(strsplit(myfiles, split = "-")), byrow = T, ncol = 4)[,2])
for(i in 1:length(factors)){
    mydata <- data.table()
    factorbeds <- myfiles[grep(factors[i], myfiles)]
    for(j in 1:length(factorbeds)){
        input <- fread(paste(mydir, factorbeds[j], sep = "/"), header = F)[, c(1:3)]
        input[, Files := factorbeds[j]]
        mydata <- rbindlist(list(mydata, input))
    }
    mydata[, PeakLength := V3 - V2]
    myplot <- gghistogram(data = mydata, x = "PeakLength", y = "..density..", 
                          color = "royalblue", fill = "royalblue", alpha = 0.5,
                          bins = 500, title = factors[i], 
                          xlim = c(0, quantile(mydata$PeakLength, probs = 0.95)),
                          facet.by = "Files")
    ggsave(filename = paste("gghistogram", factors[i], "facet.pdf", sep = "."), 
           plot = myplot, width = 20, height = 12)
    
    myplot <- gghistogram(data = mydata, x = "PeakLength", y = "..density..",
                          color = "Files", fill = "Files", alpha = 0.2,
                          bins = 500, title = factors[i], 
                          xlim = c(0, quantile(mydata$PeakLength, probs = 0.95)))
    ggsave(filename = paste("gghistogram", factors[i], "combined.pdf", sep = "."), 
           plot = myplot, width = 10, height = 6)
    
    myplot <- ggdensity(data = mydata, x = "PeakLength", y = "..density..",
                        color = "Files", title = factors[i], 
                        xlim = c(0, quantile(mydata$PeakLength, probs = 0.95)))
    ggsave(filename = paste("ggdensity", factors[i], "combined.pdf", sep = "."), 
           plot = myplot, width = 10, height = 6)
}