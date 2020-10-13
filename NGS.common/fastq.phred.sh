# if the quality scores contain character 0 it is either Sanger phred+33 or Illumina 1.8+ phred+33. 
# When they also contain the character J, it is Illumina 1.8+ phred 33, otherwise it is Sanger phred + 33.
# When the quality scores do not contain 0, it is either Solexa +64, Illumina 1.3+ Phred+64, Illumina 1.5+ Phred+64.
# Then it is Solexa +64 when it contains character =
# It is Illumina 1.3 phred + 64 when it contains A
# It is Illumina 1.5 phred +64 when it contains no A or =
# Take a look at the wiki and try to understand the table

 #  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
 #  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
 #  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
 #  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ.....................
 #  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
 #  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 #  |                         |    |        |                              |                     |
 # 33                        59   64       73                            104                   126
 #  0........................26...31.......40                                
 #                           -5....0........9.............................40 
 #                                 0........9.............................40 
 #                                    3.....9..............................41 
 #  0.2......................26...31........41                              

 # S - Sanger        Phred+33,  raw reads typically (0, 40)
 # X - Solexa        Solexa+64, raw reads typically (-5, 40)
 # I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 # J - Illumina 1.5+ Phred+64,  raw reads typically (3, 41)
 #     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
 #     (Note: See discussion above).
 # L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

###############################
head -n 40 file.fastq | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding\!";}'
###############################

###############################
# As noted by medhat above, GNU od or hexdump can be used to convert the quality scores to their decimal value, so
cat file.fq | awk 'NR%4==0' | tr -d '\n' | hexdump -v -e'/1 "%u\n"' | sort -nu
# will display which (decimal) quality scores exist in your file.
# According to brentp's "guess-encoding.py" script the possible ranges are 33-93 (Sanger/Illumina1.8), 
# 64-104 (Illumina1.3 or Illumina1.5) and 59-104 (Solexa). 
# Similarly FastQC assumes that anything with some scores in the 33-63 range is Sanger 
# and that the rest is Illumina1.3-1.5 (it doesn't know about Solexa scores).
###############################
