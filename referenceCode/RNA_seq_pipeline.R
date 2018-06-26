trim_fastq=function(input, output){
  tool="/opt/Trimmomatic-0.36/trimmomatic-0.36.jar"
  N_thread=8
  
  ##single end
  cmd=paste("nohup java -jar", tool, "SE",
            "-threads ", N_thread, 
            input, output,
            "ILLUMINACLIP:./adapters/Pool-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
            ">>run_trim.log 2>&1",
            "&", sep=" ")
  system(cmd)
}

trim_fastq_PE=function(input1, input2, adp.file, outdir){
  tool="/opt/Trimmomatic-0.36/trimmomatic-0.36.jar"
  N_thread=8
  if(!dir.exists(outdir)){dir.create(outdir)}
  
  ##pair end
  fq.name=string_overlap(basename(input1), basename(input2)) ##overlap string of two files
  fq.name=sub("[._]*(\\.fq\\.gz)?$", "", fq.name) ##the common prefix of forward and reverse files
  out.1p=paste(outdir, "/", fq.name, "_1P.fq.gz", sep="")
  out.1u=paste(outdir, "/", fq.name, "_1U.fq.gz", sep="")
  out.2p=paste(outdir, "/", fq.name, "_2P.fq.gz", sep="")
  out.2u=paste(outdir, "/", fq.name, "_2U.fq.gz", sep="")
  cmd=paste("java -jar", tool, "PE",
            "-threads ", N_thread, 
            input1, input2, 
            out.1p, out.1u, out.2p, out.2u,
            paste("ILLUMINACLIP:",adp.file, ":2:30:10", sep=""),
            "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
            ">>run_trim_PE.log 2>&1",
            sep=" ")
  system(cmd)
  
}

##Need FastQC installed
FastQC_control=function(indir, outdir){
  tool="fastqc"
  N_thread=8
  if(!dir.exists(outdir)){dir.create(outdir)}
  
  cmd=paste("nohup", tool, 
            paste(indir, "/*.fq.gz", sep=""),
            "-o", outdir, "-t", N_thread, "-q",
            "&", sep=" " )
  system(cmd)
  
}


##Download Gencode GTF annotation and generate genome index for STAR
generate_STAR_index=function(GTFfile, FASTAfile, outdir){
  N_thread=8
  
  cmd=paste("nohup STAR",
            "--runMode genomeGenerate",
            "--genomeDir", outdir,
            "--genomeFastaFiles", FASTAfile,
            "--sjdbOverhang 100",
            "--sjdbGTFfile", GTFfile,
            "--runThreadN", N_thread,
            ">>run_STAR_index.log 2>&1",
            "&", sep=" ")
  system(cmd)
  
}


##Need STAR installed. Genome index should be generated previously
STAR_align=function(indir, outdir, ref){
  tool="STAR"
  N_thread=8
  if(!dir.exists(outdir)){dir.create(outdir)}
  genome.index=ref
  infiles=list.files(trim.dir, pattern = ".fq.gz$", full.names = T)
  
  if(length(infiles)==1){
  ##single end
  cat("Single end alignment...\n")
  cmd=paste(tool,
    "--genomeDir", genome.index,
    "--readFilesIn", infiles,
    "--outFileNamePrefix",  paste(outdir, "/", sep=""),
    "--runThreadN", N_thread,
    "--outFilterMultimapScoreRange 1",
    "--outFilterMultimapNmax 20", 
    "--outFilterMismatchNmax 10", 
    "--alignIntronMax 500000",
    "--alignMatesGapMax 1000000",
    "--sjdbScore 2", 
    "--alignSJDBoverhangMin 1",
    "--genomeLoad NoSharedMemory", 
    "--limitBAMsortRAM 0", 
    "--readFilesCommand zcat", 
    "--outFilterMatchNminOverLread 0.33", 
    "--outFilterScoreMinOverLread 0.33", 
    "--sjdbOverhang 100", 
    "--outSAMstrandField intronMotif", 
    "--outSAMattributes NH HI NM MD AS XS", 
    "--outSAMunmapped Within", 
    "--outSAMtype BAM SortedByCoordinate", 
    "--outSAMheaderHD @HD VN:1.4",
    ">>run_STAR.log 2>&1",
    sep=" ")
  system(cmd)
  }else if(length(infiles)>1){
    cat("Pair end alignment...\n")
    file1=list.files(trim.dir, pattern = "_1P.fq.gz$", full.names = T)
    files1=paste(file1, collapse = ",")
    file2=list.files(trim.dir, pattern = "_2P.fq.gz$", full.names = T)
    files2=paste(file2, collapse = ",")
    cmd=paste(tool,
              "--genomeDir", genome.index,
              "--readFilesIn", paste(files1, files2, sep=" "),
              "--outFileNamePrefix",  paste(outdir, "/", sep=""),
              "--runThreadN", N_thread,
              "--outFilterMultimapScoreRange 1",
              "--outFilterMultimapNmax 20", 
              "--outFilterMismatchNmax 10", 
              "--alignIntronMax 500000",
              "--alignMatesGapMax 1000000",
              "--sjdbScore 2", 
              "--alignSJDBoverhangMin 1",
              "--genomeLoad NoSharedMemory", 
              "--limitBAMsortRAM 0", 
              "--readFilesCommand zcat", 
              "--outFilterMatchNminOverLread 0.33", 
              "--outFilterScoreMinOverLread 0.33", 
              "--sjdbOverhang 100", 
              "--outSAMstrandField intronMotif", 
              "--outSAMattributes NH HI NM MD AS XS", 
              "--outSAMunmapped Within", 
              "--outSAMtype BAM SortedByCoordinate", 
              "--outSAMheaderHD @HD VN:1.4",
              ">>run_STAR.log 2>&1",
              sep=" ")
    system(cmd)
  }else{
    cat("No files!\n")
  }
  
}


HTSeq=function(infile, outfile, gtf.file){
  cmd=paste("samtools view -F 4", infile, ##Filter the low qulity read. Used in TCGA pipeline
            "|", 
            "htseq-count", 
            "-", gtf.file,
            "-m intersection-nonempty", ##Used in TCGA pipeline
            "-i gene_id",
            "-r pos", ##order by position
            "-s no", ##non-stranded sequence data
            "-q",
            ">", outfile,
            sep=" ")
  system(cmd)
  
}


Cufflinks=function(infile, outdir, gtf.file){
  N_thread=8
  if(!dir.exists(outdir)){dir.create(outdir)}
  cmd=paste("cufflinks", infile,
            "-G", gtf.file,
            "-p", N_thread,
            "-o", outdir,
            "-q",
            ">>run_cufflinks.log 2>&1",
            sep=" ")
  system(cmd)
}


##return the overlap between two strings
string_overlap=function(s1, s2){
  diff.index=which(unlist(strsplit(s1, ""))!=unlist(strsplit(s2, ""))) ##split the string and find out the position with difference
  same=substr(s1, 1,diff.index[1]-1)
  return(paste(same, collapse = ""))
}


##combine FPKM of each sample into FPKM matrix
FPKM_matrix=function(fpkm.file, fpkm.mat.file){
  fpkm.mat=c()
  for (i in 1:length(fpkm.file)) {
    fpkm=read.table(fpkm.file[i], header = T, sep="\t", quote = "", stringsAsFactors = F  )
    if(i==1){
      fpkm.mat=subset(fpkm, select = c("tracking_id", "gene_short_name", "locus", "FPKM"))
      colnames(fpkm.mat)[4]=names(fpkm.file[i]) ##sample name
      fpkm.mat=fpkm.mat[order(fpkm.mat[, 1]), ] ##sort by tracking id
      next
    }else{
      f=subset(fpkm, select=c("tracking_id", "gene_short_name", "locus", "FPKM"))
      colnames(f)[4]=names(fpkm.file[i])
      fpkm.mat=merge(fpkm.mat, f, by=c("tracking_id", "gene_short_name", "locus"))
    }
  }
  write.table(fpkm.mat, file=fpkm.mat.file, quote = F, sep="\t", row.names = F, col.names = T)
}

##combine read count of each sample into Count matrix
Count_matrix=function(htseq.file, count.mat.file){
  count.mat=c()
  for (i in 1:length(htseq.file)) {
    htseq=read.table(htseq.file[i], header = F, sep="\t", quote = "", stringsAsFactors = F)
    colnames(htseq)=c("gene_id", names(htseq.file[i]))
    htseq=htseq[!grepl("^_", htseq$gene_id), ]
    if(i==1){
      count.mat=htseq
      next
    }else{
      count.mat=merge(count.mat, htseq, by="gene_id")
    }
  }
  write.table(count.mat, file=count.mat.file, quote = F, sep="\t", row.names = F, col.names = T)
}

##merge sample information with fastq files
sample_info=function(sample.anno.file, filelist1, filelist2){
  sample.anno=read.table(sample.anno.file, header = T, sep="\t", quote = "", stringsAsFactors = F)
  filelist=c(filelist1, filelist2)
  file.df=data.frame(file_path=filelist, file_name=basename(filelist))
  if(all(file.df$file_name %in% sample.anno$file_name)){
    cat("Good! All fastq files have sample annotations.\n")
    sample.info=merge(file.df, sample.anno)
  }else{
    cat("Not all fastq files have sample annotations!\n")
    stop()
  }
  return(sample.info)
  
}

##DESeq2, differential analysis
make_diff_analysis=function(count.mat.file, sample.df, DEseq.dir){
  library("DESeq2")
  library("BiocParallel")
  library("pheatmap")
  library("ggplot2")
  N_thread=2
  
  count.mat=read.table(count.mat.file, header = T, sep="\t", quote = "",row.names = 1, stringsAsFactors = F)
  coldata=sample.df[, -1, drop=FALSE]
  rownames(coldata)=sample.df[, 1]
  if(all(rownames(coldata) %in% colnames(count.mat))){
    cat("Good! sample name matches with sample annotatio file\n")
  }else{
    cat("sample name doesn't match with sample annotatio file!\n")
    stop()
  }
  count.mat <- count.mat[, rownames(coldata)]
  if(all(rownames(coldata) == colnames(count.mat))){
    cat("Good! Sample order is same in count matrix and sample annotation\n")
  }else{
    cat("Sample order is wrong!\n")
    stop()
  }
  
  dds <- DESeqDataSetFromMatrix(countData = count.mat,
                                colData = coldata,
                                design = ~ sample_classification)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds
               , BPPARAM=MulticoreParam(N_thread), parallel = F
  )
  
  
  # this gives log2(n + 1)
  ntd <- normTransform(dds)
  
  pcaData <- plotPCA(ntd, intgroup=c("sample_classification"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pca.plot=ggplot(pcaData, aes(PC1, PC2, color=sample_classification, shape=sample_classification)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  ggsave(paste(DEseq.dir, "/PCA_plot.pdf", sep=""), pca.plot, width=8, height=8)
  
  comp=unique(coldata$sample_classification)
  all.diff=c()
  for (i in 1:(length(comp)-1)) {
    for (j in (i+1):length(comp)) {
      res=results(dds, contrast = c("sample_classification", comp[i], comp[j]))
      res.df=as.data.frame(res)
      out.file=paste(DEseq.dir, "/", comp[i],"_VS_", comp[j],".txt", sep="")
      write.table(res.df, out.file, quote = F, sep="\t", row.names = T, col.names = T)
      
      diff.sel=res.df$padj<0.05 & (abs(res.df$log2FoldChange) > 1)
      diff.sel[is.na(diff.sel)]=FALSE
      if(i==1){
        all.diff=diff.sel
      }else{
        all.diff= all.diff | diff.sel
      }
      
      i.anno=subset(coldata, sample_classification==comp[i])
      j.anno=subset(coldata, sample_classification==comp[j])
      one.dds=assay(ntd)[, c(rownames(i.anno), rownames(j.anno))]
      one.dds=one.dds[diff.sel, ]
      pheatmap(one.dds, cluster_rows=T, show_rownames=F,
               cluster_cols=T, annotation_col=rbind(i.anno, j.anno), width = 8, height = 8, 
               main = paste(length(which(diff.sel)), "differential genes for", comp[i], "VS", comp[j], sep=" "),
               filename = paste(DEseq.dir, "/diff_gene_heatmap_of", comp[i], "_VS_", comp[j], ".pdf", sep=""))
      
    }
  }
  
  all.diff.num=length(which(all.diff))
  ntd.df=assay(ntd)
  ntd.df=ntd.df[all.diff, ]
  pheatmap(ntd.df, cluster_rows=T, show_rownames=F,
           cluster_cols=T, annotation_col=coldata, width = 8, height = 8, 
           main = paste(all.diff.num, "differential genes among samples", sep=" "),
           filename = paste(DEseq.dir, "/All_diff_gene_heatmap.pdf", sep=""))
  
}


##input
pre.dir="../raw/RNA-SEQ/C202SC18031099"
fq.file.f=c("ES_WT/ES_WT_USR17003586L_HKY5KBBXX_L4_1.fq.gz",
            "ES_KO/ES_KO_USR17003587L_HKY5KBBXX_L4_1.fq.gz", 
            "EBKO2/EBKO2_USR17003593L_HKY5KBBXX_L4_1.fq.gz",
            "ESWT2/ESWT2_USR17003590L_HKY5KBBXX_L4_1.fq.gz",
            "EBWT2/EBWT2_USR17003592L_HKY5KBBXX_L4_1.fq.gz",
            "EBKO/EBKO_USR17003589L_HKY5KBBXX_L4_1.fq.gz",
            "ESKO2/ESKO2_USR17003591L_HKY5KBBXX_L4_1.fq.gz",
            "EBWT/EBWT_USR17003588L_HKY5KBBXX_L4_1.fq.gz"
            )
fq.file.f=paste(pre.dir, "/", fq.file.f, sep="")
fq.file.r=c("ES_WT/ES_WT_USR17003586L_HKY5KBBXX_L4_2.fq.gz",
            "ES_KO/ES_KO_USR17003587L_HKY5KBBXX_L4_2.fq.gz", 
            "EBKO2/EBKO2_USR17003593L_HKY5KBBXX_L4_2.fq.gz",
            "ESWT2/ESWT2_USR17003590L_HKY5KBBXX_L4_2.fq.gz",
            "EBWT2/EBWT2_USR17003592L_HKY5KBBXX_L4_2.fq.gz",
            "EBKO/EBKO_USR17003589L_HKY5KBBXX_L4_2.fq.gz",
            "ESKO2/ESKO2_USR17003591L_HKY5KBBXX_L4_2.fq.gz",
            "EBWT/EBWT_USR17003588L_HKY5KBBXX_L4_2.fq.gz"
            )
fq.file.r=paste(pre.dir, "/", fq.file.r, sep="")
adp.file="./adapters/liang-PE.fa"
sample.anno.file="../raw/RNA-SEQ/sample_annotation.txt"
ref="/export/data2/big_data/STAR_ref_genome/GRCm38"
gtf.file="/export/data2/big_data/Genome_annotation/GRCm38/gencode.vM17.annotation.gtf"
mode="PE"

##output
Trim.dir="../data/trimmed_fastq"
if(!dir.exists(Trim.dir)){dir.create(Trim.dir)}
Fastqc.dir="../data/FastQC"
if(!dir.exists(Fastqc.dir)){dir.create(Fastqc.dir)}
Align.dir="../data/alignment"
if(!dir.exists(Align.dir)){dir.create(Align.dir)}
Htseq.dir="../data/HTSeq"
if(!dir.exists(Htseq.dir)){dir.create(Htseq.dir)}
Cufflinks.dir="../data/Cufflinks"
if(!dir.exists(Cufflinks.dir)){dir.create(Cufflinks.dir)}
Exp.dir="../data/Expression"
if(!dir.exists(Exp.dir)){dir.create(Exp.dir)}
DEseq.dir="../data/DESeq2_result"
if(!dir.exists(DEseq.dir)){dir.create(DEseq.dir)}



if(FALSE){
generate_STAR_index("/export/data2/big_data/Genome_annotation/GRCm38/gencode.vM17.annotation.gtf",
                    "/export/data2/big_data/Genome_annotation/GRCm38/GRCm38.primary_assembly.genome.fa",
                    "/export/data2/big_data/STAR_ref_genome/GRCm38")
generate_STAR_index("/export/data2/big_data/Genome_annotation/NCBIM37/gencode.vM1.annotation.gtf",
                    "/export/data2/big_data/Genome_annotation/NCBIM37/NCBIM37.genome.fa",
                    "/export/data2/big_data/STAR_ref_genome/NCBIM37")
}


if(FALSE){
##guess sample name from fastq file
fq.name=c()
for (i in 1:length(fq.file.f)) {
  fq.name[i]=string_overlap(basename(fq.file.f[i]), basename(fq.file.r[i])) ##overlap string of two files
  fq.name[i]=sub("[._]*(\\.fq\\.gz)?$", "", fq.name[i]) ##the common prefix of forward and reverse files
}
}

sample.anno=sample_info(sample.anno.file, fq.file.f, fq.file.r)
fq.name=sample.anno$sample[sapply(fq.file.f, function(x){which(sample.anno$file_path==x)})]


for(i in 1:length(fq.name)){
cat(i, fq.name[i], "...\n")

cat("Step 1: Remove adaptors from fastq...\n")
trim.dir=paste(Trim.dir, "/", fq.name[i], sep="")
if(mode=="PE"){
  ##input the two fastq file, the adapter fasta file, and outdir
  #trim_fastq_PE(fq.file.f[i], fq.file.r[i], adp.file, trim.dir)
  
}else if(mode == "SE"){
  #trim_fastq(fq.file[1], trim.file[1])
}


cat("Step 2: FastQC quality control...\n")
qc.dir=paste(Fastqc.dir, "/", fq.name[i], sep="")
#FastQC_control(trim.dir, qc.dir)


cat("Step 3: Genome alignment using STAR...\n")
star.dir=paste(Align.dir, "/", fq.name[i], sep="")
#STAR_align(trim.dir, star.dir, ref)


cat("Step 4: Gene quantification using HTSeq-count...\n")
bam.file=paste(star.dir, "/Aligned.sortedByCoord.out.bam", sep="")
count.file=paste(Htseq.dir, "/", fq.name[i], ".txt", sep="")
#HTSeq(bam.file, count.file, gtf.file)

cat("Step 5: Cal FPKM using cufflinks...\n")
cuff.dir=paste(Cufflinks.dir, "/", fq.name[i], sep = "")
#Cufflinks(bam.file, cuff.dir, gtf.file)

}


##Generate expression matrix
fpkm.gene.file=paste(Cufflinks.dir, "/", fq.name, "/genes.fpkm_tracking", sep="") ##all FPKM files
names(fpkm.gene.file)=fq.name
fpkm.transcript.file=paste(Cufflinks.dir, "/", fq.name, "/isoforms.fpkm_tracking", sep="") ##all FPKM files of transcripts
names(fpkm.transcript.file)=fq.name
htseq.file=paste(Htseq.dir, "/", fq.name, ".txt", sep="") ##all read count files
names(htseq.file)=fq.name
if(all(file.exists(fpkm.gene.file)) & all(file.exists(fpkm.transcript.file)) & all(file.exists(htseq.file)) ){
  cat("All sample expression is ready. Make expression matrix...\n")
  
  ##FPKM gene matrix
  fpkm.gene.mat.file=paste(Exp.dir, "/FPKM_gene_matrix.txt", sep="")
  FPKM_matrix(fpkm.gene.file, fpkm.gene.mat.file)
  
  ##FPKM transcript matrix
  fpkm.transcript.mat.file=paste(Exp.dir, "/FPKM_transcript_matrix.txt", sep="")
  FPKM_matrix(fpkm.transcript.file, fpkm.transcript.mat.file)
  
  ##HTSeq-count matrix
  count.mat.file=paste(Exp.dir, "/HTSeq_count_matrix.txt", sep="")
  Count_matrix(htseq.file, count.mat.file)
}else{
  cat("Not all samples are done.")
}


cat("Step 6: Differential gene analysis...\n")
sample.df=unique(subset(sample.anno, select = c("sample", "sample_classification")))
make_diff_analysis(count.mat.file, sample.df, DEseq.dir)





