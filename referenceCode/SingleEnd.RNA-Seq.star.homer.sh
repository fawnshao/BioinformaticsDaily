#/bin/bash
hg19indexdir=/data/shaojf/myReference/star_index/hg19.star
WorkDir=/data/shaojf/pubseq/ThoreenCC2017ATF4
FastqDir=$WorkDir/FastqDir
BamDir=$WorkDir/BamDir
BwDir=$WorkDir/BwDir
HomerTagDir=$WorkDir/HomerTagDir
hubname="ThoreenCC2017.293T.RNASeq"
# cd RawData
# tail -n +6 runinfo.txt | cut -f 5 | xargs -n 1 prefetch
# mv ~/ncbi/public/*.sra .
# tail -n +6 runinfo.txt | cut -f 5,7 | \
# 	paste - gse.description | sed 's?,?\t?g;s? ??g' | \
# 	awk -F"\t" '{print "mv "$1".sra "$1"."$5"."$6".sra"}'
# for f in *.sra; do fastq-dump -O ../FastqDir/ $f & done
cd $WorkDir/
# --readFilesCommand zcat \
for f in $FastqDir/*.fastq
do
	date
	pre=`echo $f | awk -F"/" '{print $NF}' | sed 's/.fastq//'`
	echo $pre

	STAR --runMode alignReads --runThreadN 48 \
	--genomeDir $hg19indexdir --genomeLoad LoadAndRemove \
	--readFilesIn $FastqDir/${pre}.fastq \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 1 \
	--outFilterScoreMinOverLread  0.5 \
	--outFilterMatchNminOverLread 0.5 \
	--outFileNamePrefix $BamDir/${pre}.toGenome. \
	--outSAMattributes All \
	--outSAMtype BAM Unsorted \
	--outFilterType BySJout --outReadsUnmapped Fastx \
	--outFilterScoreMin 10 --outSAMattrRGline ID:$pre \
	--alignEndsType Local

	samtools sort -@ 12 -o $BamDir/${pre}.toGenome.srt.bam $BamDir/${pre}.toGenome.Aligned.out.bam &
done
wait
for f in $BamDir/*.toGenome.srt.bam
do
	samtools index -@ 16 $f &
done
wait
for f in $BamDir/*.toGenome.srt.bam
do
	pre=`echo $f | awk -F"/" '{print $NF}' | sed 's/.toGenome.srt.bam//'`
	makeTagDirectory $HomerTagDir/$pre.mTD -fragLength given $f &
done
wait
makeMultiWigHub.pl $hubname hg19 -url $BwDir -webdir $BwDir -d $HomerTagDir/*.mTD &
wait

# move to /data/webshare/
mkdir -p /data/webshare/$hubname/
mv $BwDir/$hubname/hg19/*.mTD.ucsc.bigWig /data/webshare/$hubname/

# generate index.html for listing files
echo "<!DOCTYPE html>" > /data/webshare/$hubname/index.html
echo "<html lang=\"en\">" >> /data/webshare/$hubname/index.html
echo "<head>" >> /data/webshare/$hubname/index.html
echo "    <meta charset=\"UTF-8\">" >> /data/webshare/$hubname/index.html
echo "    <title>DataSharing</title>" >> /data/webshare/$hubname/index.html
echo "</head>" >> /data/webshare/$hubname/index.html
echo "<body>" >> /data/webshare/$hubname/index.html
echo "" >> /data/webshare/$hubname/index.html
echo "<h1>$hubname</h1>" >> /data/webshare/$hubname/index.html
echo "<ul>" >> /data/webshare/$hubname/index.html
for f in `ls /data/webshare/$hubname/`
do
	echo "	<li><a href=$f>$f</a></li>" >> /data/webshare/$hubname/index.html
done
f=igv_session.xml
echo "	<li><a href=$f>$f</a></li>" >> /data/webshare/$hubname/index.html
echo "</ul>" >> /data/webshare/$hubname/index.html
echo "</body>" >> /data/webshare/$hubname/index.html

# generate igv.session.xml for IGV
echo "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" > /data/webshare/$hubname/igv_session.xml
echo "<Session genome=\"hg19\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"chr22:39914815-39920942\" path=\"http://139.52.17.130/webshare/$hubname/igv_session.xml\" version=\"8\">" >> /data/webshare/$hubname/igv_session.xml
echo "    <Resources>" >> /data/webshare/$hubname/igv_session.xml
for f in `ls /data/webshare/$hubname/ | grep bigWig`
do
    echo "    <Resource path=\"http://139.52.17.130/webshare/$hubname/$f\"/>" >> /data/webshare/$hubname/igv_session.xml
done
echo "    </Resources>" >> /data/webshare/$hubname/igv_session.xml
echo "    <Panel height=\"435\" name=\"DataPanel\" width=\"1606\">" >> /data/webshare/$hubname/igv_session.xml
for f in `ls /data/webshare/$hubname/ | grep bigWig`
do
    echo "    <Track altColor=\"0,0,178\" autoScale=\"true\" clazz=\"org.broad.igv.track.DataSourceTrack\" color=\"0,0,178\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" id=\"http://139.52.17.130/webshare/$hubname/$f\" name=\"$f\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">" >> /data/webshare/$hubname/igv_session.xml
    echo "    <DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"880.0367\" minimum=\"0.0\" type=\"LINEAR\"/>" >> /data/webshare/$hubname/igv_session.xml
    echo "	</Track>" >> /data/webshare/$hubname/igv_session.xml
done
echo "	</Panel>" >> /data/webshare/$hubname/igv_session.xml
echo "	<Panel height=\"214\" name=\"FeaturePanel\" width=\"1606\">" >> /data/webshare/$hubname/igv_session.xml
echo "	<Track altColor=\"0,0,178\" autoScale=\"false\" color=\"0,0,178\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" id=\"Reference sequence\" name=\"Reference sequence\" sortable=\"false\" visible=\"true\"/>" >> /data/webshare/$hubname/igv_session.xml
echo "	<Track altColor=\"0,0,178\" autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"0,0,178\" colorScale=\"ContinuousColorScale;0.0;423.0;255,255,255;0,0,178\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" height=\"35\" id=\"hg19_genes\" name=\"RefSeq Genes\" renderer=\"BASIC_FEATURE\" sortable=\"false\" visible=\"true\" windowFunction=\"count\">" >> /data/webshare/$hubname/igv_session.xml
echo "	<DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"423.0\" minimum=\"0.0\" type=\"LINEAR\"/>" >> /data/webshare/$hubname/igv_session.xml
echo "	</Track>" >> /data/webshare/$hubname/igv_session.xml
echo "	</Panel>" >> /data/webshare/$hubname/igv_session.xml
echo "	<PanelLayout dividerFractions=\"0.6661585365853658\"/>" >> /data/webshare/$hubname/igv_session.xml
echo "	<HiddenAttributes>" >> /data/webshare/$hubname/igv_session.xml
echo "	<Attribute name=\"DATA FILE\"/>" >> /data/webshare/$hubname/igv_session.xml
echo "	<Attribute name=\"DATA TYPE\"/>" >> /data/webshare/$hubname/igv_session.xml
echo "	<Attribute name=\"NAME\"/>" >> /data/webshare/$hubname/igv_session.xml
echo "	</HiddenAttributes>" >> /data/webshare/$hubname/igv_session.xml
echo "</Session>" >> /data/webshare/$hubname/igv_session.xml
### make folder accessible
chcon -R -t httpd_sys_content_t /data/webshare/$hubname/
chmod -R 777 /data/webshare/$hubname/


