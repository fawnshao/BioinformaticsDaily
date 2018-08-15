library(circlize)
args <- commandArgs(TRUE)
# args <- c("chr1","11166592", "11802263", "subnet",
#           "11321564", "11323564", "ENSG00000198793.8|MTOR", "HKG", 
#           "11798028","11802263","enh1", "enh")
# chr1	11166592	11322564	ENSG00000198793.8|MTOR	.	-
# chr1:11798028-11802263

df = data.frame(factors = args[1],
                x = c(as.numeric(args[2:3])))

circos.clear()
circos.par("track.height" = 0.1)
circos.initialize(factors = df$factors, x = df$x)

circos.track(factors = df$factors, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index, cex = 0.8)
               circos.axis(labels.cex = 0.6, 
                           major.at = df$x[1] + 1e5 * seq(1:floor((df$x[2] - df$x[1]) / 1e5)))
             })

for(i in 1:(length(args)/4)){
    coords <- c(as.numeric(args[(i*4+1):(i*4+2)]))
    circos.segments(coords[1] ,0.5, coords[2], 0.5, col = factor(args[(i*4+4)]), lwd = 2)
    circos.text(coords[1], -1, args[(i*4+3)], col = "blue", niceFacing = T, facing = "bending.inside")
}
circos.link("seg", 10501, "seg", 16001, h = 0.4)

circos.segments(10201,0.5,10850,0.5, col = factor(), lwd = 2)
circos.text(10550, -1, "gene1", col = "blue", niceFacing = T, facing = "bending.inside")

circos.segments(16001,0.5,17900,0.5, col = "red", lwd = 2)
circos.text(16550, -1, "gene2", col = "blue", niceFacing = T, facing = "bending.inside")

circos.link("seg", 10501, "seg", 16001, h = 0.4)
# circos.link("seg", c(10501,10520), "seg", c(15001,15201),
#             col = "green", border = "black", h = 0.3, lwd = 2, lty = 1)

circos.clear()
