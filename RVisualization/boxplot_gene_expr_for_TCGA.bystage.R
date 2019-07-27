# args <- commandArgs(TRUE)
library(ggpubr)
library(data.table)
library(gridExtra)
# args <- c("xena.mygene.txt", "xena.stage.info.txt")
args <- c("ICOSLG.txt", "xena.stage.info.txt")
exprmat <- fread(args[1], header = T)
samples <- fread(args[2], header = T)

genes <- unique(exprmat$sample)
samples[, ajcc_pathologic_tumor_stage := gsub("Stage ", replacement = "", ajcc_pathologic_tumor_stage)]
stages <- names(sort(table(samples$ajcc_pathologic_tumor_stage), decreasing=T))[c(2:16)]
cancers <- c("BRCA", "KIRC", "LUAD", "THCA", "HNSC", "LUSC", "UCEC", "COAD", "BLCA", "STAD",
	"SKCM", "LIHC", "OV", "KIRP", "CESC", "ESCA", "PAAD", "READ", "KICH", "MESO", 
	"TGCT", "UVM", "ACC", "THYM", "CHOL", "UCS", "DLBC")

exprmatsim <- data.matrix(exprmat[,-1])
rownames(exprmatsim) <- exprmat$sample
exprmatsim <- t(exprmatsim)
cancersample <- matrix(unlist(strsplit(rownames(exprmatsim), split = "-")), ncol = 4, byrow = T)[,4]
cancersample <- as.numeric(cancersample)
sampledesc <- data.frame(SampleTypeCodes = cancersample, 
  CancerType = samples[match(rownames(exprmatsim), samples$sample)]$`cancer type abbreviation`, 
  TumorStage = samples[match(rownames(exprmatsim), samples$sample)]$`ajcc_pathologic_tumor_stage`)
sampledesc <- data.frame(sampledesc, 
  SampleCodeTumorStage = paste(sampledesc$SampleTypeCodes, sampledesc$TumorStage))

for (i in 1:length(genes)) {
	datasim <- data.table(Expression = exprmatsim[,i], sampledesc)
	colnames(datasim)[1] <- "Expression"
	colnames(datasim)[3] <- "Cancer"
	# datasim <- datasim[ajcc_pathologic_tumor_stage!="" & ajcc_pathologic_tumor_stage!="[Discrepancy]" & ajcc_pathologic_tumor_stage!="/"]
	datasim <- datasim[TumorStage %in% stages]
	datasim <- datasim[Cancer %in% cancers]
	datasim$SampleCodeTumorStage <- factor(datasim$SampleCodeTumorStage, levels = sort(unique(datasim$SampleCodeTumorStage)))
	datasim$TumorStage <- factor(datasim$TumorStage, levels = sort(unique(datasim$TumorStage)))
	datasim$SampleTypeCodes <- factor(datasim$SampleTypeCodes, levels = sort(unique(datasim$SampleTypeCodes)))
	for(j in 1:length(cancers)){
		p1 <- ggboxplot(data = datasim[Cancer == cancers[j]], x = "SampleTypeCodes", y = "Expression", 
			fill = "TumorStage", color = "TumorStage", alpha = 0.3,
			title = paste(genes[i], cancers[j]), 
			outlier.shape = NA, notch = F, add = "median")
		ggsave(paste("ggboxplot", genes[i], cancers[j], "png", sep = "."), p1,
			width = 18, height = 10)
	}
}


data <- data.frame(info[,2], log2(as.numeric(x[-1])+1))
genename <- x[1]
# data <- melt(expr)
colnames(data) <- c("tissue", "expression")
png(filename = paste(genename, "violin.png", sep = "."), width = 1500, height = 600)
ggplot(data, aes(x = tissue, y = expression, fill = tissue)) + 
	geom_violin() + 
	stat_summary(fun.data = "mean_sdl", geom = "pointrange", color = "black") +
	theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1)) 
dev.off()


p1 <- ggboxplot(data = mysigs[!is.na(headtails)], x = "experiment", y = "CGILength", 
                fill = "headtails", color = "headtails", palette = "jco", 
                outlier.shape = NA, notch = T)
p11 <- ggpar(p1, x.text.angle = 60, ylim = c(0, quantile(unique(mysigs$CGILength), probs = 0.95)))
ggsave("ggboxplot.H3K4me3.peak.toCGI.headtails10.sig.png", 
       grid.arrange(p11, p22, p33, ncol = 3),
       width = 25, height = 10)