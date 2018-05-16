args <- commandArgs(TRUE)
# args <- c("shNUP53-2_Dox_E2", "shNUP53-2_E2", "3562", "30895", "51030")
library(VennDiagram)
pdf(file = paste(args[1], args[2], "venn.pdf", sep = "."))
venn.plot <- draw.pairwise.venn(
    area1 = as.numeric(args[3]) + as.numeric(args[5]),
    area2 = as.numeric(args[4]) + as.numeric(args[5]),
    cross.area = as.numeric(args[5]),
    category = args[1:2],
    fill = c("blue", "red"),
    lty = "blank",
    scaled = TRUE,
    # cex = 2,
    alpha = rep(0.3, 2),
    # cat.cex = 1,
    # cat.pos = c(285, 105),
    # cat.dist = 0.05,
    # cat.just = list(c(-1, -1), c(1, 1)),
    # ext.pos = 30,
    # ext.dist = -0.05,
    # ext.length = 0.5,
    # ext.line.lwd = 2,
    ext.line.lty = "dashed"
)
dev.off()
