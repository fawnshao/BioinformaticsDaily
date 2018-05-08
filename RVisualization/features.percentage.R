library(data.table)
library(ggplot2)
args <- commandArgs(TRUE)
input <- fread(args[1], sep = "\t", header = T)
count <- table(input[,2])
class <- length(count)

for (i in 3:ncol(input)){
	png(filename = paste(args[1], colnames(input)[i], "barplot.png", sep = "."), width = 400, height = 400)
	datax <- input[,c(2,i),with = FALSE]
	colnames(datax) <- c("type", "value")
	a <- table(datax[datax$value > 0,1])
	if(length(a) < class){
		a <- count - table(datax[datax$value == 0,1])
	}
	b <- a/count
	datay <- as.data.frame(b)
	myplot <- ggplot(data = datay, aes(x = Var1, y = Freq, fill = Var1)) + 
		geom_bar(stat = "identity") +
		ggtitle(colnames(input)[i]) + theme(legend.position = "none")
	print(myplot)
	dev.off()
}

