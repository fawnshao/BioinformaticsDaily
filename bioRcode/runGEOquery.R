library(Biobase)
library(GEOquery)
library(limma)
# library(GEOmetadb)
# load series and platform data from GEO

gset <- getGEO("GSE50588", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("XXXXX0XXXXXXXXXXXX0X0XXXXXXXXXXXXXXXX0XXXXX0XXXXXX",
               "XXXXXXXXX0XX0XXXXXXXXXXXXXXXXXXXXXX00XXXXXX0XXXXXX",
               "XXXXXXX0XXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX111",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXX000000XXXXXXXXXXXXXXXXX",
               "X")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


################################## get the expression
mydata <- getGEO("GSE50588", GSEMatrix = FALSE)
mygse <- getGEO("GSE50588", GSEMatrix = TRUE)
myExpr <- exprs(mygse[[1]])
dim(myExpr)
probesets <- Table(GPLList(mydata)[[1]])
sampledesc <- pData(phenoData(mygse[[1]]))[,c(1,8)]
myExprAnnos <- data.frame(ID = rownames(myExpr),
                          Symbol = probesets[match(rownames(myExpr),probesets$ID),]$Symbol
                          )
mySampleAnnos <- data.frame(ACC = colnames(myExpr),
                            Desc = sampledesc[match(colnames(myExpr), rownames(sampledesc)),]$title,
                            Sample = sampledesc[match(colnames(myExpr), rownames(sampledesc)),]$source_name_ch1
                            )
ctrl <- myExpr[,match(mySampleAnnos[grep("NS",mySampleAnnos$Sample),1], colnames(myExpr))]
# "siSP1","siE2F4","siYY1","siEP300","siRAD21","siGTF2B"
a <- myExpr[,match(mySampleAnnos[grep("SP1",mySampleAnnos$Sample),1], colnames(myExpr))]
b <- myExpr[,match(mySampleAnnos[grep("E2F4",mySampleAnnos$Sample),1], colnames(myExpr))]
c <- myExpr[,match(mySampleAnnos[grep("YY1",mySampleAnnos$Sample),1], colnames(myExpr))]
d <- myExpr[,match(mySampleAnnos[grep("EP300",mySampleAnnos$Sample),1], colnames(myExpr))]
e <- myExpr[,match(mySampleAnnos[grep("RAD21",mySampleAnnos$Sample),1], colnames(myExpr))]
f <- myExpr[,match(mySampleAnnos[grep("GTF2B",mySampleAnnos$Sample),1], colnames(myExpr))]
# platforms <- annotation(data1[[1]])
# sqlfile <- getSQLiteFile()
# con <- dbConnect("SQLite",sqlfile)
# dbGetQuery(con,"select gpl,title,bioc_package from gpl where gpl=platforms")
