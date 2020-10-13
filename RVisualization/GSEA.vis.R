setwd("/Users/fawnshao/gsea_home/output/jan22/")
library(data.table)
library(ggpubr)
library(gridExtra)
hallmarks.pos <- fread("E14K.NT.GseaPreranked.1579754044842/gsea_report_for_na_pos_1579754044842.xls", header = T)
hallmarks.neg <- fread("E14K.NT.GseaPreranked.1579754044842/gsea_report_for_na_neg_1579754044842.xls", header = T)
gobp.pos <- fread("E14K.NT.GOBP.GseaPreranked.1579754821363/gsea_report_for_na_pos_1579754821363.xls", header = T)
gobp.neg <- fread("E14K.NT.GOBP.GseaPreranked.1579754821363/gsea_report_for_na_neg_1579754821363.xls", header = T)
hallmarks <- data.table(rbindlist(list(hallmarks.pos, hallmarks.neg)), 
                        Type = c(rep("UP in Tumor", nrow(hallmarks.pos)), rep("DN in Tumor", nrow(hallmarks.neg))))
hallmarks$Type <- factor(hallmarks$Type, levels = rev(sort(unique(hallmarks$Type))))
gobps <- data.table(rbindlist(list(gobp.pos, gobp.neg)), 
                    Type = c(rep("UP in Tumor", nrow(gobp.pos)), rep("DN in Tumor", nrow(gobp.neg))))
gobps$Type <- factor(gobps$Type, levels = rev(sort(unique(gobps$Type))))

p1 <- ggbarplot(hallmarks[`FDR q-val` < 0.05], x = "NAME", y = "NES", 
                fill = "Type", color = "Type", palette = "npg", title = "Hallmarks FDR < 0.05")
p11 <- ggpar(p1, rotate = T)
# p11
ggsave("/Users/fawnshao/Data/myworkData/pancancer.data/NUP93.E14K/ggbarplot.hallmarks.pdf", width = 10, height = 8, p11)
p2 <- ggbarplot(gobps[`FDR q-val` < 0.001], x = "NAME", y = "NES", 
                fill = "Type", color = "Type", palette = "npg", title = "GO BP FDR < 0.001")
p22 <- ggpar(p2, rotate = T)
ggsave("/Users/fawnshao/Data/myworkData/pancancer.data/NUP93.E14K/ggbarplot.gobp.pdf", width = 18, height = 15, p22)

