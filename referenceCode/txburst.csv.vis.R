library(ggpubr)
library(dplyr)
library(data.table)
library(gridExtra)

# txburstML.py outputs a pickled pandas dataframe which for each gene contains an array with three entries [k_on, k_off, k_syn] 
# where the burst frequency = k_on and burst size = k_syn/k_off, 
# and a boolean indicating whether gene passed a rudimentary filtering step based on the quality of the inference procedure.

# txburstTEST.py outputs a pickled pandas dataframe with the two point estimates from txburstML.py and p-values for the significance tests for differential bursting kinetics. For example:
#     
#     gene name	point estimates 1	point estimates 2	bf_pvalue	bs_pvalue
# gene1	[k_on,k_off,k_syn]_gene1_sample1	[k_on,k_off,k_syn]_gene1_sample2	burst frequency pvalue gene1	burst size pvalue gene1
# gene2	[k_on,k_off,k_syn]_gene2_sample1	[k_on,k_off,k_syn]_gene2_sample2	burst frequency pvalue gene2	burst size pvalue gene2
setwd("/Volumes/GoogleDrive/My Drive/NUP_project/NUP_related_others/NS27.NUP93.scRNA/gene.eRNA/seurat/txburst.res/")
nup93 <- fread(input = "fmt.filtered.X1.mat2csv_vs_filtered.X2.mat2csv_TEST.pkl.tsv", header = F, sep = "\t")
rad21 <- fread(input = "fmt.filtered.R1.mat2csv_vs_filtered.R2.mat2csv_TEST.pkl.tsv", header = F, sep = "\t")
nup93[, burstsize1 := V4/V3]
nup93[, burstsize2 := V7/V6]
rad21[, burstsize1 := V4/V3]
rad21[, burstsize2 := V7/V6]

# nup93.sim <- nup93[!is.na(V8), .SD, .SDcols = c(1,2,5,10:11)]
# rad21.sim <- rad21[!is.na(V8), .SD, .SDcols = c(1,2,5,10:11)]
nup93.sim <- nup93[, .SD, .SDcols = c(1,2,5,10:11)]
rad21.sim <- rad21[, .SD, .SDcols = c(1,2,5,10:11)]
colnames(nup93.sim) <- c("Gene", "burst.freq.UT", "burst.freq.IAA", "burst.size.UT", "burst.size.IAA")
colnames(rad21.sim) <- c("Gene", "burst.freq.UT", "burst.freq.IAA", "burst.size.UT", "burst.size.IAA")

nup93[V1=="ENSG00000049130"]
rad21[V1=="ENSG00000049130"]

nup93[V1=="ENSG00000160888"]
rad21[V1=="ENSG00000160888"]

data.merge <- data.table(rbindlist(list(nup93.sim, rad21.sim)), 
                         Target = c(rep("NUP93", nrow(nup93.sim)), rep("RAD21", nrow(rad21.sim))))
data.merge.m <- melt(data.merge, id.vars = c("Gene", "Target"))
data.merge.m[grep("burst.freq", variable), Type := "burst.freq"]
data.merge.m[grep("burst.size", variable), Type := "burst.size"]
data.merge.m[grep("UT", variable), Treatment := "UT"]
data.merge.m[grep("IAA", variable), Treatment := "IAA"]
data.merge.m$Treatment <- factor(data.merge.m$Treatment, levels = c("UT","IAA"))

p1 <- ggboxplot(data = data.merge.m, x = "Target", y = "value", 
                fill = "Treatment", color = "Treatment", 
                alpha = 0.8, palette = "npg", facet.by = "Type",
                outlier.shape = NA, ylim = c(0,20))
# p1 <- ggdotplot(data = data.merge.m, x = "Target", y = "value", 
#                 fill = "Treatment", color = "Treatment", 
#                 alpha = 0.8, palette = "npg", facet.by = "Type",
#                 binwidth = 0.1, size = 0.3, ylim = c(0,20))
# + stat_compare_means()
# wilcoxtest1 <- wilcox.test(x = data.merge.m[Type=="burst.freq" & Target=="NUP93" & Treatment=="UT"]$value, 
#                           y = data.merge.m[Type=="burst.freq" & Target=="NUP93" & Treatment=="IAA"]$value, 
#                           alternative = "two.sided")
# wilcoxtest2 <- wilcox.test(x = data.merge.m[Type=="burst.freq" & Target=="RAD21" & Treatment=="UT"]$value, 
#                            y = data.merge.m[Type=="burst.freq" & Target=="RAD21" & Treatment=="IAA"]$value, 
#                            alternative = "two.sided")
stat.test1 <- compare_means(value ~ Treatment, data = data.merge.m[Type=="burst.freq"], 
                           group.by = "Target",method = "wilcox.test")
stat.test1 <- stat.test1 %>% mutate(y.position = quantile(data.merge.m[Type=="burst.freq"]$value, probs = 1))
# stat.test
p1 <- ggboxplot(data = data.merge.m[Type=="burst.freq"], x = "Target", y = "value", 
                fill = "Treatment", color = "Treatment", 
                alpha = 0.8, palette = "npg", width = 0.6,
                ylab = "burst.freq", xlab = "", outlier.size = 0.1,
                ylim = quantile(data.merge.m[Type=="burst.freq"]$value, probs = c(0,1))) +
    stat_pvalue_manual(stat.test1, label = "p.adj = {p.adj}", x = "Target", position = position_dodge(0.8))
stat.test2 <- compare_means(value ~ Treatment, data = data.merge.m[Type=="burst.size"], 
                           group.by = "Target",method = "wilcox.test")
stat.test2 <- stat.test2 %>% mutate(y.position = quantile(data.merge.m[Type=="burst.size"]$value, probs = 0.98))
p2 <- ggboxplot(data = data.merge.m[Type=="burst.size"], x = "Target", y = "value", 
                fill = "Treatment", color = "Treatment", 
                alpha = 0.8, palette = "npg", width = 0.6,
                ylab = "burst.size", xlab = "", outlier.size = 0.1,
                ylim = quantile(data.merge.m[Type=="burst.size"]$value, probs = c(0,0.98))) +
    stat_pvalue_manual(stat.test2, label = "p.adj = {p.adj}", x = "Target", position = position_dodge(0.8))

ggsave(filename = "ggboxplot.burst.freq.pdf", width = 8, height = 5,
       plot = grid.arrange(p1, p2, ncol = 2))
####
p1 <- ggdotplot(data = data.merge.m[Type=="burst.freq"], x = "Target", y = "value", 
                fill = "Treatment", color = "Treatment", 
                alpha = 0.8, palette = "npg", add = "mean_sd",
                binwidth = 0.05, size = 0.3, ylab = "burst.freq", xlab = "",
                ylim = quantile(data.merge.m[Type=="burst.freq"]$value, probs = c(0,1))) +
    stat_pvalue_manual(stat.test1, label = "p.adj = {p.adj}", x = "Target", position = position_dodge(0.8))
p2 <- ggdotplot(data = data.merge.m[Type=="burst.size"], x = "Target", y = "value", 
                fill = "Treatment", color = "Treatment", 
                alpha = 0.8, palette = "npg", add = "mean",
                binwidth = 0.2, size = 0.3, ylab = "burst.size", xlab = "",
                ylim = quantile(data.merge.m[Type=="burst.size"]$value, probs = c(0,0.99))) +
    stat_pvalue_manual(stat.test2, label = "p.adj = {p.adj}", x = "Target", position = position_dodge(0.8))
ggsave(filename = "ggdotplot.burst.freq.pdf", width = 10, height = 5,
       plot = grid.arrange(p1, p2, ncol = 2))


####
# p1 <- ggscatter(data = data.merge[Target=="NUP93"], x = "burst.freq.UT", y = "burst.size.UT")
geneannotation <- fread("../genes.annotation.txt", header = F)
data.merge[, burst.freq.lab := ""]
data.merge[, SignificantDots := "NotSig"]
cutoff <- quantile(data.merge.m[Type=="burst.freq"]$value, probs = 0.99)
data.merge[((burst.freq.IAA+0.1)/(burst.freq.UT+0.1) > 3 | (burst.freq.UT+0.1)/(burst.freq.IAA+0.1) > 3) 
           & (burst.freq.UT > cutoff | burst.freq.IAA > cutoff), 
           burst.freq.lab := geneannotation[match(Gene, geneannotation$V1)]$V2]
data.merge[((burst.freq.IAA+0.1)/(burst.freq.UT+0.1) > 3 | (burst.freq.UT+0.1)/(burst.freq.IAA+0.1) > 3) 
           & (burst.freq.UT > cutoff | burst.freq.IAA > cutoff), 
           SignificantDots := "Sig"]
# rug = T, add = c("reg.line"), 
p1 <- ggscatterhist(data = data.merge[Target=="NUP93"], x = "burst.freq.UT", y = "burst.freq.IAA", 
                    cor.coef = T, conf.int = T,
                    color = "SignificantDots",
                    xlim = quantile(data.merge.m[Type=="burst.freq"]$value, probs = c(0,1)), 
                    ylim = quantile(data.merge.m[Type=="burst.freq"]$value, probs = c(0,1)), 
                    bins = 100, margin.plot = c("histogram"), 
                    margin.params = list(fill = "skyblue", color = "royalblue"), 
                    repel = T, label = "burst.freq.lab")
p2 <- ggscatterhist(data = data.merge[Target=="RAD21"], x = "burst.freq.UT", y = "burst.freq.IAA", 
                    cor.coef = T, conf.int = T,
                    color = "SignificantDots",
                    xlim = quantile(data.merge.m[Type=="burst.freq"]$value, probs = c(0,1)), 
                    ylim = quantile(data.merge.m[Type=="burst.freq"]$value, probs = c(0,1)), 
                    bins = 100, margin.plot = c("histogram"), 
                    margin.params = list(fill = "skyblue", color = "royalblue"), 
                    repel = T, label = "burst.freq.lab")
ggsave(filename = "ggscatterhist.burst.freq.pdf", width = 16, height = 8,
       plot = grid.arrange(p1, p2, ncol = 2))


data.merge[, burst.size.lab := ""]
data.merge[, SignificantDots := "NotSig"]
cutoff <- quantile(data.merge.m[Type=="burst.size"]$value, probs = 0.98)
data.merge[((burst.size.IAA+0.1)/(burst.size.UT+0.1) > 3 | (burst.size.UT+0.1)/(burst.size.IAA+0.1) > 3) 
           & (burst.size.UT >cutoff | burst.size.IAA > cutoff), 
           burst.size.lab := geneannotation[match(Gene, geneannotation$V1)]$V2]
data.merge[((burst.size.IAA+0.1)/(burst.size.UT+0.1) > 3 | (burst.size.UT+0.1)/(burst.size.IAA+0.1) > 3) 
           & (burst.size.UT > cutoff | burst.size.IAA > cutoff), 
           SignificantDots := "Sig"]
# data.merge[SignificantDots=="Sig"]
# rug = T, add = c("reg.line"), 
p1 <- ggscatterhist(data = data.merge[Target=="NUP93"], x = "burst.size.UT", y = "burst.size.IAA", 
                    cor.coef = T, conf.int = T,
                    color = "SignificantDots",
                    xlim = quantile(data.merge.m[Type=="burst.size"]$value, probs = c(0,0.995)), 
                    ylim = quantile(data.merge.m[Type=="burst.size"]$value, probs = c(0,0.995)), 
                    bins = 100, margin.plot = c("histogram"), 
                    margin.params = list(fill = "royalblue", color = "royalblue"), 
                    repel = T, label = "burst.size.lab")
p2 <- ggscatterhist(data = data.merge[Target=="RAD21"], x = "burst.size.UT", y = "burst.size.IAA", 
                    cor.coef = T, conf.int = T,
                    color = "SignificantDots",
                    xlim = quantile(data.merge.m[Type=="burst.size"]$value, probs = c(0,0.995)), 
                    ylim = quantile(data.merge.m[Type=="burst.size"]$value, probs = c(0,0.995)), 
                    bins = 100, margin.plot = c("histogram"), 
                    margin.params = list(fill = "royalblue", color = "royalblue"), 
                    repel = T, label = "burst.size.lab")
ggsave(filename = "ggscatterhist.burst.size.pdf", width = 16, height = 8,
       plot = grid.arrange(p1, p2, ncol = 2))

saveRDS(data.merge, file = "burst.data.merge.RDS")
########################
data.comb <- data.merge[, .SD, .SDcols = c(1:6)]
data.comb[, GeneName := geneannotation[match(Gene, geneannotation$V1)]$V2]

# nup93[!is.na(V8) & V8 < 5e-2 & V9 < 5e-2]
nup93.sim1 <- nup93[!is.na(V8) & V8 < 5e-2 & V9 < 5e-2, .SD, .SDcols = c(1,2,5,10:11)]
rad21.sim1 <- rad21[!is.na(V8) & V8 < 5e-2 & V9 < 5e-2, .SD, .SDcols = c(1,2,5,10:11)]
colnames(nup93.sim1) <- c("Gene", "burst.freq.UT", "burst.freq.IAA", "burst.size.UT", "burst.size.IAA")
colnames(rad21.sim1) <- c("Gene", "burst.freq.UT", "burst.freq.IAA", "burst.size.UT", "burst.size.IAA")
data.merge1 <- data.table(rbindlist(list(nup93.sim1, rad21.sim1)), 
                          Target = c(rep("NUP93", nrow(nup93.sim1)), rep("RAD21", nrow(rad21.sim1))))
data.merge1[, GeneName := geneannotation[match(Gene, geneannotation$V1)]$V2]
data.merge1[Target=="NUP93" & 
                ((burst.freq.IAA+0.1)/(burst.freq.UT+0.1) > 3 | (burst.freq.UT+0.1)/(burst.freq.IAA+0.1) > 3) &
                ((burst.size.IAA+0.1)/(burst.size.UT+0.1) > 3 | (burst.size.UT+0.1)/(burst.size.IAA+0.1) > 3) ]

data.merge1[Target=="NUP93" & 
                ((burst.freq.IAA+0.1)/(burst.freq.UT+0.1) > 5 | (burst.freq.UT+0.1)/(burst.freq.IAA+0.1) > 5)]$GeneName
data.merge1[Target=="NUP93" & 
                ((burst.size.IAA+0.1)/(burst.size.UT+0.1) > 5 | (burst.size.UT+0.1)/(burst.size.IAA+0.1) > 5)]$GeneName

paste(data.comb[Target=="NUP93" & 
              ((burst.freq.IAA+0.1)/(burst.freq.UT+0.1) > 3 | (burst.freq.UT+0.1)/(burst.freq.IAA+0.1) > 3) &
              ((burst.size.IAA+0.1)/(burst.size.UT+0.1) > 3 | (burst.size.UT+0.1)/(burst.size.IAA+0.1) > 3) ]$GeneName, 
      collapse = ",")

# > data.comb[GeneName=="MALAT1"]
# Gene burst.freq.UT burst.freq.IAA burst.size.UT burst.size.IAA Target GeneName
# 1: ENSG00000251562      10.00018       1.565488      1.000182        641.076  NUP93   MALAT1
# > data.merge1[GeneName=="MALAT1"]