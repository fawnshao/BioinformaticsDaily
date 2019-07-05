# args <- commandArgs(TRUE)
# args <- c("htseq_6sample.csv")
# for bamfile in NS23*.bam
# do
# 	samtools view $bamfile | awk '$7=="=" && $4 > $8 {a=substr($6,1,3);b=split(a,arr,"M");print $4-$8+arr[1]}' > is.add.$bamfile.txt &
# done

library(ggpubr)
library(data.table)
files <- dir(pattern = "is.NS")
allsize <- data.table()
for(i in 1:length(files)){
	insertsize <- fread(files[i], header = F)
	insertsize[, expr:=files[i]]
	allsize <- rbindlist(list(allsize, insertsize))
}
p1 <- gghistogram(data = allsize, x = "V1", y = "..density..", 
	color = "coral3", fill = "coral3", bins = 1000, facet.by = "expr")
ggsave("gghistogram.insertsize.pdf", p1, width = 20, height = 10)

p1 <- gghistogram(data = allsize[V1 < 1e5], x = "V1", y = "..density..", 
	color = "coral3", fill = "coral3", bins = 1000, facet.by = "expr")
ggsave("gghistogram.insertsize.lt10k.pdf", p1, width = 20, height = 10)

p1 <- gghistogram(data = allsize[V1 < 5e3], x = "V1", y = "..density..", 
	color = "coral3", fill = "coral3", bins = 1000, facet.by = "expr")
ggsave("gghistogram.insertsize.lt5k.pdf", p1, width = 20, height = 10)

p1 <- gghistogram(data = allsize[V1 < 1e3], x = "V1", y = "..density..", 
	color = "coral3", fill = "coral3", bins = 1000, facet.by = "expr")
ggsave("gghistogram.insertsize.lt1k.pdf", p1, width = 20, height = 10)

allsize.man <- allsize
allsize.man[V1 > 1000, V1:=1000]
p1 <- gghistogram(data = allsize.man, x = "V1", y = "..density..", 
	color = "coral3", fill = "coral3", bins = 1000, facet.by = "expr")
ggsave("gghistogram.insertsize.man1k.pdf", p1, width = 20, height = 10)

allsize <- data.table()
for(i in 1:length(files)){
	insertsize <- fread(files[i], header = F)
	insertsize[, expr:=files[i]]
	allsize <- rbindlist(list(allsize, insertsize))
}
allsize[V1 > 5000, V1:=5000]
p1 <- gghistogram(data = allsize, x = "V1", y = "..density..", 
	color = "coral3", fill = "coral3", bins = 1000, facet.by = "expr")
ggsave("gghistogram.insertsize.man5k.pdf", p1, width = 20, height = 10)

allsize <- data.table()
for(i in 1:length(files)){
	insertsize <- fread(files[i], header = F)
	insertsize[, expr:=files[i]]
	allsize <- rbindlist(list(allsize, insertsize))
}
allsize[V1 > 1e5, V1:=1e5]
p1 <- gghistogram(data = allsize, x = "V1", y = "..density..", 
	color = "coral3", fill = "coral3", bins = 1000, facet.by = "expr")
ggsave("gghistogram.insertsize.man100k.pdf", p1, width = 20, height = 10)


allsize <- data.table()
for(i in 1:length(files)){
	insertsize <- fread(files[i], header = F)
	insertsize[, expr:=files[i]]
	allsize <- rbindlist(list(allsize, insertsize))
}
allsize[, log10IS:=log10(V1)]
p1 <- gghistogram(data = allsize, x = "log10IS", y = "..density..", 
	color = "coral3", fill = "coral3", bins = 1000, facet.by = "expr")
ggsave("gghistogram.insertsize.log10IS.pdf", p1, width = 20, height = 10)

for(i in 1:length(files)){
	print(files[i])
	print(median(allsize[expr==files[i]]$V1))
}

for(i in 1:length(files)){
	print(files[i])
	print(median(allsize[expr==files[i] & V1 < 1000]$V1))
}

