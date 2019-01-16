sam-dump SRR1636861 | samtools view -bS - > SRR1636861.bam

# Use examples:
sam-dump SRR390728
# Output SAM format data to standard out. Alignment information is not required to output in this format.
sam-dump --aligned-region 1:6484848-6521430 --output-file SRR390728.sam SRR390728
# Store output in the file SRR390728.sam for only the region 6484848-6521430 on chromosome 1. The sequence name as submitted for the alignment (@SQ SN in SAM/BAM files) or the reference sequence accession must be used.
sam-dump SRR390728 | samtools view -bS - > SRR390728.bam
# With "samtools" installed the above command pipes (|) the sam-dump output directly to samtools for conversion directly into .bam format (view -bS; the "-" following -bS allow samtools to read the streaming data from sam-dump).
sam-dump -r --gzip --output-file SRR390728.sam.gz SRR390728
# Produces gzip’d (--gzip) output file (--output-file) SRR390728.sam.gz that has a reconstructed header(-r). Will include reference accessions but may not include some header info from the submitted data.