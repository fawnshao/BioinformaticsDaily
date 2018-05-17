args <- commandArgs(TRUE)
# args <- c("shNUP53-2_Dox_E2", "shNUP53-2_E2", "3562", "30895", "51030")
library(VennDiagram)
if(length(args)==5){
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
}
if(length(args)==10){
    pdf(file = paste(args[1], args[2], args[3], "venn.pdf", sep = "."))
    venn.plot <- draw.triple.venn(
        area1 = as.numeric(args[4]) + as.numeric(args[7] + as.numeric(args[8] + as.numeric(args[10]),
        area2 = as.numeric(args[5]) + as.numeric(args[8] + as.numeric(args[9] + as.numeric(args[10]),
        area3 = as.numeric(args[6]) + as.numeric(args[7] + as.numeric(args[9] + as.numeric(args[10]),
        n12 = as.numeric(args[7]), 
        n23 = as.numeric(args[8]),
        n13 = as.numeric(args[9]),
        n123 = as.numeric(args[10]),
        category = args[1:3],
        fill = c("blue", "red", "green"),
        lty = "blank",
        scaled = TRUE,
        alpha = rep(0.3, 2),
        ext.line.lty = "dashed"
    )
    dev.off()
}
