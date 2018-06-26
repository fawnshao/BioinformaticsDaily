library("pheatmap")
library("colorspace")
library("RColorBrewer")
library("clusterProfiler")
##repeat scale in row and column for n times.Make result scaled both in row and column 
##Number of rows and columns must >2
repeat_scale=function(mat, n){
  mat.norm=mat
  for (i in 1:n) {
    mat.norm=scale(mat) ##scale in column
    mat.norm=t(scale(t(mat))) ##scale in row
  }
  return(mat.norm)
}


##Description: 
##  Select diff gene from DESeq results. return a vector of diff gene name.
##Args:
##  defiles, DESeq2 result files. c(file1, file2,...)
##  set, "union" or "intersect".
##  fc, threshold of log2 Fold change.
##  padj, threshold of adjusted p values.
select_diff_gene=function(defiles, set="OR", fc=1, padj=0.01){
  all.diff.gene=c()
  for (i in 1:length(defiles)) {
    de.df=read.table(defiles[i], header = T, sep="\t", quote = "\"", stringsAsFactors = F )
    ind=(de.df$padj < padj) & (abs(de.df$log2FoldChange) > fc)
    diff.gene=unique(rownames(de.df[which(ind), ]))
    cat(length(diff.gene), "genes in", defiles[i], "\n")
    if(i==1){
      all.diff.gene=diff.gene
    }else{
      if(grepl("OR", set, ignore.case = T)){
        all.diff.gene=union(all.diff.gene, diff.gene)
      }else if(grepl("AND", set, ignore.case = T)){
        all.diff.gene=intersect(all.diff.gene, diff.gene)
      }else{
        stop("unknow operrations!\n")
      }
    }
  }
  all.diff.gene=unique(all.diff.gene)
  cat("Totally ", length(all.diff.gene), "diff genes\n")
  
  return(all.diff.gene)
}


##Description:
##  Extract a subset matrix from whole FPKM matrix by selecting diff genes and samples. Return a matrix.
##Args:
##  
extract_mat=function(FPKM.mat, diff.gene, samples=colnames(FPKM.mat)[4:length(colnames(FPKM.mat))]){
  diff.gene=unique(intersect(diff.gene, FPKM.mat$tracking_id))
  cat(length(diff.gene), "diff genes have same ID with FPKM matrix...\n")
  diff.mat=subset(FPKM.mat, tracking_id %in% diff.gene)
  low.gene=diff.mat$tracking_id[rowSums(diff.mat[, c(-1,-2,-3)])<10]
  cat(length(unique(low.gene)), "genes expressed lower than 10 in all samples\n")
  diff.gene=setdiff(diff.gene, low.gene)
  diff.mat=subset(diff.mat, tracking_id %in% diff.gene)
  cat(length(which(duplicated(diff.mat$tracking_id))), "duplicated genes are removed\n")
  diff.mat=diff.mat[!duplicated(diff.mat$tracking_id), ]    ##remove duplicate genes
  
  plot.mat=subset(diff.mat, select=samples)
  rownames(plot.mat)=diff.mat$tracking_id
  plot.mat=log2(plot.mat+1)
  name.table=diff.mat[, 1:2]
  cat("Totally ", nrow(plot.mat), "genes are used for heatmap\n")
  
  return(list(plot.mat, name.table))
}

##Description:
##  Draw heatmap with different clustering methods.
##Args:
##  plot.mat, log2(FPKM) matrix
##  graph.dir, output dir
##  my.clust, a list of clustering object
draw_heatmap=function(mat.file, graph.dir, my.clust, ntree){
  plot.mat=read.table(mat.file, header = T, sep="\t", quote = "\"", stringsAsFactors = F)
  ncolor=20
  #mycol=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(ncolor)
  mycol=colorRampPalette(c("green", "black", "red"))(ncolor)
  my.range=round(quantile(unlist(plot.mat), c(0.1, 0.9)), 2)
  my.break=c(seq(min(plot.mat), my.range[1], length.out = ncolor*0.1),
               seq(my.range[1], my.range[2], length.out=ncolor*0.8+1)[-1], 
               seq(my.range[2],max(plot.mat), length.out = ncolor*0.1)[-1])
  cmp.name=sub("FPKM_matrix_(.+).txt", "\\1", basename(mat.file))
  
  for (i in 1:length(my.clust)) {
    graph=paste(graph.dir, "/", names(my.clust)[i], "_heatmap_", cmp.name , ".pdf", sep="")
    groups=cutree(my.clust[[i]], ntree)
    gene.group=paste("group", groups, sep="")
    gene.anno=data.frame(group=gene.group, stringsAsFactors = F)
    rownames(gene.anno)=names(groups)
    pheatmap(plot.mat, 
             cluster_rows=T, 
             show_rownames=F,
               cluster_cols=F, 
               cutree_rows = ntree,
               clustering_distance_rows = "euclidean",
               clustering_method=names(my.clust)[i],
               color=mycol,
               #breaks = my.break,
               scale = "row",
               annotation_row=gene.anno, 
               #annotation_colors = list(group=rainbow(4)),
               width = 6, height = 8, 
               main = paste(nrow(plot.mat), "differential genes for", cmp.name, sep=" "),
               filename = graph)
  }
}


##Description:
##  Gene function enrichment
function_enrichment=function(groups, species, id.type="ENSEMBL", data.dir, graph.dir){
  cmp.name=basename(data.dir)
  p.cut=0.05
  show.num=10
  if(species=="mmu"){##mouse
    orgdb="org.Mm.eg.db"
  }else if(species=="hsa"){##human
    orgdb="org.Hs.eg.db"
  }else if(species=="sce"){##yeast
    orgdb="org.Sc.sgd.db"
  }
  else{
    stop(paste("unknown species", species, "!\n"))
  }
  gene.id=sub("\\.\\d+$", "", names(groups)) ##remove version number of gene id
  gene.group=paste("group", groups, sep="")
  gene.df=data.frame(gene.id, gene.group, stringsAsFactors = F)
  colnames(gene.df)=c(id.type, "group")
  id.df <- bitr(gene.id, fromType=id.type, toType=c("ENTREZID"), OrgDb=orgdb)
  id.df=merge(gene.df, id.df, by=id.type)
  gene.list=lapply(unique(id.df$group), function(x){id.df$ENTREZID[id.df$group==x]})
  names(gene.list)=unique(id.df$group)
  ##KEGG enrichment
  kegg <- compareCluster(geneCluster = gene.list, fun = "enrichKEGG", 
                         organism=species, 
                         pAdjustMethod = "BH",
                         pvalueCutoff = p.cut,
                         #qvalueCutoff  = 0.05,
                         minGSSize=10
                         )
  write.table(kegg@compareClusterResult, file=paste(data.dir, "/KEGG_enrichment_", cmp.name, ".txt", sep=""), quote = F, sep="\t", row.names = F, col.names = T)
  kegg.plot=dotplot(kegg, showCategory=show.num, by="count")+
    theme(axis.text.x=element_text(angle=90, hjust=0))
  ggsave(paste(graph.dir, "/KEGG_enrichment_", cmp.name, ".pdf", sep=""), kegg.plot, width=16, height=16, units = "cm")
  
  ##MKEGG enrichment
  mkegg=compareCluster(geneClusters = gene.list, fun="enrichMKEGG",
                       organism=species, 
                       pAdjustMethod = "BH",
                       pvalueCutoff = p.cut,
                       #qvalueCutoff  = 0.05,
                       minGSSize=10
                    )
  write.table(mkegg@compareClusterResult, file=paste(data.dir, "/MKEGG_enrichment_", cmp.name, ".txt", sep=""), quote = F, sep="\t", row.names = F, col.names = T)
  mkegg.plot=dotplot(mkegg, showCategory=show.num, by="count")+
    theme(axis.text.x=element_text(angle=90, hjust=0))
  ggsave(paste(graph.dir, "/MKEGG_enrichment_", cmp.name, ".pdf", sep=""), mkegg.plot, width=16, height=16, units = "cm")
  
  ##GO-MF enrichment
  go.mf=compareCluster(geneClusters = gene.list, fun="enrichGO",
                    OrgDb=orgdb, 
                    ont="MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = p.cut,
                    #qvalueCutoff  = 0.05,
                    minGSSize=10,
                    readable=TRUE
  )
  go.mf=gofilter(go.mf, level = 4)
  write.table(go.mf@compareClusterResult, file=paste(data.dir, "/GO-MF_enrichment_", cmp.name, ".txt", sep=""), quote = F, sep="\t", row.names = F, col.names = T)
  go.mf.plot=dotplot(go.mf, showCategory=show.num, by="count")+
    theme(axis.text.x=element_text(angle=90, hjust=0), axis.text.y = element_text(size=8))
  ggsave(paste(graph.dir, "/GO-MF_enrichment_", cmp.name, ".pdf", sep=""), go.mf.plot, width=20, height=16, units = "cm")
  
  ##GO-BP enrichment
  go.bp=compareCluster(geneClusters = gene.list, fun="enrichGO",
                    OrgDb=orgdb, 
                    ont="BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = p.cut,
                    #qvalueCutoff  = 0.05,
                    minGSSize=10,
                    readable=TRUE
  )
  go.bp=gofilter(go.bp, level = 4)
  write.table(go.bp@compareClusterResult, file=paste(data.dir, "/GO-BP_enrichment_", cmp.name, ".txt", sep=""), quote = F, sep="\t", row.names = F, col.names = T)
  go.bp.plot=dotplot(go.bp, showCategory=show.num, by="count")+
    theme(axis.text.x=element_text(angle=90, hjust=0), axis.text.y = element_text(size=8))
  ggsave(paste(graph.dir, "/GO-BP_enrichment_", cmp.name, ".pdf", sep=""), go.bp.plot, width=20, height=16, units = "cm")
  
  if(species=="hsa"){
  ##MSigDB-H enrichment
  hallmark=read.gmt("/export/data2/big_data/MSigDB/msigdb_v6.1_GMTs/h.all.v6.1.entrez.gmt")
  h=enricher(gene = gene.list[[1]],
                          pAdjustMethod = "BH",
                          pvalueCutoff = p.cut,
                          minGSSize=10,
                          TERM2GENE=hallmark
                          )
  }
}



##input
FPKM.file="../data/RNA-SEQ/Expression/FPKM_gene_matrix.txt"
deseq.files=c("../data/RNA-SEQ/DESeq2_result/EB_WT_VS_ES_WT.txt"
             #"../data/RNA-SEQ/DESeq2_result/EB_KO_VS_ES_KO.txt"
             )
mysample=c("ES_WT", "ES_WT2", "EB_WT", "EB_WT2", "ES_KO", "ES_KO2", "EB_KO", "EB_KO2")
SPEC="mmu" ##mouse, "mmu"; human, "hsa"; yeast, "sce". https://www.genome.jp/kegg/catalog/org_list.html

##output
diff.dir="../data/RNA_comparison"
if(!dir.exists(diff.dir)){dir.create(diff.dir)}
diff.graph.dir="../graph/RNA_comparison"
if(!dir.exists(diff.graph.dir)){dir.create(diff.graph.dir)}



cat("Step1: Select diff gene...\n")
set="OR"
diff.gene=select_diff_gene(deseq.files, set, fc=1, padj=0.01)


cat("Step2: Extract heatmap matrix...\n")
cmp.name=sub(".txt$", "", basename(deseq.files))
cmp.name=paste(cmp.name, collapse = paste("-", set, "-", sep=""))
cmp.dir=paste(diff.dir, "/", cmp.name, sep="")
if(!dir.exists(cmp.dir)){dir.create(cmp.dir)}
mat.file=paste(cmp.dir, "/", "FPKM_matrix_", cmp.name, ".txt", sep="")
gene.file=paste(cmp.dir, "/", "Gene_table_", cmp.name, ".txt", sep="")
FPKM.mat=read.table(FPKM.file, header = T, sep="\t", quote = "", stringsAsFactors = F)
plot.list=extract_mat(FPKM.mat, diff.gene, samples=mysample) ##select your wanted samples in FPKM matrix
plot.mat=plot.list[[1]]
gene.df=plot.list[[2]]
write.table(plot.mat, file=mat.file, quote = F, sep="\t", row.names = T, col.names = T)
write.table(gene.df, file=gene.file, quote = F, sep="\t", row.names = F, col.names = T)

cat("Step3: Plot heatmap clustering...\n")
md=c("ward.D", "ward.D2", "average", "median")
ntree=4
plot.mat=read.table(mat.file, header = T, sep="\t", quote = "\"", stringsAsFactors = F)
mat.norm=t(scale(t(plot.mat))) ##scale in row
my.clust=lapply(md, function(x){result=hclust(dist(as.matrix(mat.norm),method = "euclidean"),method = x)})
names(my.clust)=md
draw_heatmap(mat.file, diff.graph.dir, my.clust, ntree)

cat("Step4: Compare groups on one cluster method...\n")
one.clust=my.clust[["ward.D2"]]
groups=cutree(one.clust, 4)
function_enrichment(groups, SPEC, id.type="ENSEMBL", data.dir=cmp.dir, graph.dir=diff.graph.dir)



