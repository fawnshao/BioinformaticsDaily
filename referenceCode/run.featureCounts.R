args <- commandArgs(TRUE)
library(Rsubread)
mybams <- c("mapping_results1.bam","mapping_results2.bam")
featureCounts(files = mybams, isPairedEnd = TRUE, nthreads = 272)
# # Summarize single-end reads using built-in RefSeq annotation for mouse genome ‘mm10’ (‘mm10’is the default inbuilt genome annotation):
# featureCounts(files="mapping_results_SE.sam")
# # Summarize single-end reads using a user-provided GTF annotation file:
# featureCounts(files="mapping_results_SE.sam",annot.ext="annotation.gtf",isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id")
# # Summarize single-end reads using 5 threads:
# featureCounts(files="mapping_results_SE.sam",nthreads=5)
# # Summarize BAM format single-end read data:
# featureCounts(files="mapping_results_SE.bam")
# # Summarize multiple libraries at the same time:
# featureCounts(files=c("mapping_results1.bam","mapping_results2.bam"))
# # Summarize paired-end reads and counting fragments (instead of reads):
# featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE)
# # Count fragments satisfying the fragment length criteria, eg. [50bp, 600bp]:
# featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE,checkFragLength=TRUE,minFragLength=50,maxFragLength=600)
# # Count fragments which have both ends successfully aligned without considering the fragment length constraint:
# featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE,requireBothEndsMapped=TRUE)
# # Exclude chimeric fragments from fragment counting:
# featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE,countChimericFragments=FALSE)