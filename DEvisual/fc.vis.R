# rclone sync /work/04935/shaojf/stampede2/JooPlots/ mygoogle:JooPlots/
library(data.table)
library(ggplot2)
# args <- commandArgs(TRUE)
args <- c("Joo.TFF1e.MS.noNA.txt")
input <- fread(args[1], sep = "\t", header = T, na.strings = "")
datax <- data.matrix(input[,lapply(.SD,function(x){ifelse(is.na(x),0.01,x)})][,3:6])
rep1fc <- log2((datax[,3] + 0.01) / (datax[,1] + 0.01))
rep2fc <- log2((datax[,4] + 0.01) / (datax[,2] + 0.01))
rep1fc[rep1fc > 4] <- 4
rep2fc[rep2fc > 4] <- 4
rep1fc[rep1fc < -4] <- -4
rep2fc[rep2fc < -4] <- -4
low.count <- apply(datax, 1, function(x){length(x[x < 0.5])})

res <- data.frame(rep1fc, rep2fc, datax)
rownames(res) <- as.matrix(input[,2])
flag <- rep("ND", nrow(res))
# flag[(res[,1] > 0.5 & res[,5] > 1) & (res[,2] > 0.5 & res[,6] > 1)] <- "up"
# flag[(res[,1] < -0.5 & res[,3] > 1) & (res[,2] < -0.5 & res[,4] > 1)] <- "down"
flag[res[,1] > 0.5 & res[,2] > 0.5] <- "up"
flag[res[,1] < -0.5 & res[,2] < -0.5] <- "down"
res <- data.frame(res, flag)
res$flag <- factor(res$flag, levels = c("up", "down", "ND"))
res <- res[low.count < 2,]

write.table(data.frame(rownames(res), res), file = paste(args[1],"labeled.xls", sep = "."), sep = "\t", row.names = F)

myplot <- ggplot(data = res, aes(x = rep1fc, y = rep2fc, colour = flag)) + 
	geom_point() +
	labs(title = args[1], caption = date(), x = "log2(Methyl/Non-methyl) rep1", y = "log2(Methyl/Non-methyl) rep2") + 
	theme(legend.position = "right", axis.text.x = element_text(angle = 60, hjust = 1))
png(filename = paste(args[1], "scatterplot.png", sep = "."), width = 600, height = 600)
print(myplot)
dev.off()

# final <- res[(res[,1] > 0.5 & res[,5] > 1) | (res[,2] > 0.5 & res[,6] > 1) | (res[,1] < -0.5 & res[,3] > 1) | (res[,2] < -0.5 & res[,4] > 1), ]
# myplot <- ggplot(data = final, aes(x = rep1fc, y = rep2fc, colour = flag)) + 
# 	geom_point() +
# 	labs(title = args[1], caption = date()) + 
# 	theme(legend.position = "right", axis.text.x = element_text(angle = 60, hjust = 1))
# png(filename = paste(args[1], "sim.scatterplot.png", sep = "."), width = 800, height = 800)
# print(myplot)
# dev.off()
