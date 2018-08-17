library(circlize)
library(reshape2)
args <- commandArgs(TRUE)
# args <- c("regions","segments","links")
regions <- read.table(args[1], header = F, sep = "\t")
segs <- read.table(args[2], header = F, sep = "\t")
links <- read.table(args[3], header = F, sep = "\t")

df <- melt(regions, id.vars = "V1")
df <- df[,-2]
colnames(df) <- c("factors", "x")

circos.clear()
png(filename = paste(args[1], "circlize.png", sep = "."), width = 1000, height = 1000)

circos.par("track.height" = 0.1)
circos.initialize(factors = df$factors, x = df$x)

circos.track(factors = df$factors, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index, cex = 0.8)
               circos.axis(labels.cex = 0.6, 
                           major.at = df$x[1] + 5e5 * seq(1:floor((df$x[2] - df$x[1]) / 5e6)))
             })

for(i in 1:nrow(segs)){
    shift <- rnorm(1, mean = 0, sd = 0.1)
    circos.segments(segs[i,2], 0.5 + shift, segs[i,3], 0.5 + shift, col = as.character(segs[i,5]), lwd = 2)
    # circos.text(segs[i,3], 1 + shift, segs[i,4], col = as.character(segs[i,5]), cex = 0.6, niceFacing = F, facing = "bending.inside")
}

# links[links[,3] < 5,] <- 0
for(i in 1:nrow(links)){
    shift <- rnorm(1, mean = 0, sd = 0.1)
    circos.link(levels(df$factors), links[i,1], levels(df$factors), links[i,2], h = 0.4 + shift, lwd = min(links[i,3], 10))
}

dev.off()
# circos.clear()
