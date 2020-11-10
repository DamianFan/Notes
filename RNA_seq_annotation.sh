#This script is used for paralleled Data analysis for ATAC-seq!

1.软链接
指令：ln -s /home/disk/diatre_data/*/*.fq.gz /home/disk/*/data/

2.fastqc质量控制
指令：/home/software/bin/fastqc *.fq.gz -t 20（线程） -o fastqc_out（文件夹）

3.物种参考基因组下载，NCBI上下载fna和gff文件
ftp://ftp.ncbi.nlm.nih.gov/genomes
#法一
#==================================================================================================
4.建立索引
4.1.双端索引
/home/software/bin/STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir /home/disk/caijingtao/index/*（目标文件夹）--genomeFastaFiles *.fna \
--sjdbGTFfile *.gff --sjdbOverhang 149  （双端：149；单端：49）


5.STAR比对输出BAM文件
/home/software/bin/STAR \
--runThreadN 20 \
--genomeDir $STAR_index_path \ #STAR索引
--readFilesCommand gunzip -c \
--readFilesIn ${DATA_PATH}/$i\_1.fq.gz ${DATA_PATH}/$i\_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${RESULTS_PATH}/star_out/$i-


6.samtools进行sort
#/home/software/bin/samtools sort -@ 20 ${RESULTS_PATH}/star_out/$i-Aligned.sortedByCoord.out.bam\
#-o ${RESULTS_PATH}/sort_out/$i-sort.bam

7.cuffdiff输出差异基因
/home/software/bin/cuffdiff -p 10  -o /home/disk/caijingtao/sunyu-RNA/results/diff_out/* \
-b /home/disk/*/*.fa -u /home/disk/*/*.gtf -L BLEO(对照组),BLEO_M(实验组) \
BLEO-sort1.bam,BLEO-sort2.bam,BLEO-sort3.bam BLEO_M-sort1.bam,BLEO_M-sort2.bam,BLEO_M-sort3.bam \
--library-type fr-firststrand &
#法二
#==================================================================================================
8.HISAT2输出差异基因(基因组)
程序路径：/root/anaconda3/bin
8.1建立基因索引：
hisat2-build –p 20  genome.fa genome
包括转录组：
extract_exons.py hg19.gtf > hg19.exon     
extract_splice_sites.py hg19.gtf > hg19.ss
hisat2-build -p 30 hg19.fa --ss hg19.ss --exon hg19.exon hg19
8.2比对：
hisat2 -p 40 -x /home/zhimin/ref/hisat2/rnaseqref/mm10/genome -1 *_1.fq -2  *_2.fq -S *.sam
8.3转换
/home/software/bin/samtools view -bS  *.sam  > *.bam
8.4双端排序，# -n 按reads的ID排序
/home/software/bin/samtools sort -n -@ 20 *.bam > *-sort.bam
8.5建索引
ls *.bam |xargs -i samtools index {} 
8.6计数,gff转gtf
gffread test.gff -T -o test.gtf
#注:这边使用STAR产生的*-sort.bam
/home/disk/wangxuelong/miniconda2/bin/htseq-count -f bam *_sort.bam -s no -i gene_name /home/zhimin/ref/rnaseqref/*.gtf > *.txt
                       
#==================================================================================================
9.R包分析，对照样csc11和csc12，实验样asc11和asc12
9.1数据采集 
#CONTROL
rm(list=ls())
BM1=read.table('multicov-BM1_distal-enhancer.out',sep = '\t', row.names=1)
BM2=read.table('BM2.txt',sep = '\t', row.names=1
BM=cbind(BM1,BM2)
colnames(BM)=c("BM1", "BM2")

#TREAT
L1=read.table('multicov-BM1.out',sep = '\t', row.names=1)
L2=read.table('multicov-BM2.out',sep = '\t', row.names=1)
L=cbind(L1, L2)
colnames(L)=c("L1", "L2")
L_vs_BM=cbind(L, BM)

9.2DESeq2
write.table(L_vs_BM, file="L_vs_BM.txt",row.names =T,sep = "\t")
countdata=read.table("L_vs_BM.txt",header=T,row.names=1)

#countdata <- read.csv("BM_LI-couints.csv", sep=",")
#row.names(countdata) <- countdata$Gene_name
#data <- data[,-1]
#样本数修改
nctrl =2
ntreat=2
condition = c(rep('treat',ntreat), rep('control',nctrl))
#pattern = c("one","two","one","two")  注意个数 前后
#paired_end = pattern
library(DESeq2)
coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
#dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition+pattern)
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition", "treat", "control"))
res <- res[order(res$pvalue),]
head(res)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
diff<-resdata[! is.na(resdata[,6]),]
#需整理处理杂数据,注意样本数
diff1<-subset(diff,(padj > 0 | padj < 0) & BM1 > 0 & BM2 > 0 & L1 > 0 & L2 > 0)
BM1 = log10(diff1$BM1)
BM2 = log10(diff1$BM2)
L1 = log10(diff1$L1)
L2 = log10(diff1$L2)

gene = diff1$Row.name
baseMean = diff1$baseMean
log2FoldChange = diff1$log2FoldChange
lfcSE = diff1$lfcSE
stat = diff1$stat
pvalue = diff1$pvalue
padj = diff1$padj
diff1a <- data.frame(gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,BM1,BM2,L1,L2)
#对应几个样本进行计算取平均值,放入结果
BM_mean <- apply(diff1a[,8:9],1,mean)
L_mean <- apply(diff1a[,10:11],1,mean)
diff1b <- data.frame(gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,BM1,BM2,L1,L2,BM_mean,L_mean)

write.table(diff1b, file = "BM_vs_L_allgene.csv", sep = ",", row.names = FALSE)
#截取up和down数据

diff2 <-subset(diff1b,pvalue < 0.05 )
diff4 <- subset(diff2,abs(log2FoldChange)>0.5 )  #注：一般为1
write.table(diff4, file = "BM_vs_L-diffGene.csv", sep = ",", row.names = FALSE)
up <- subset(diff4,log2FoldChange>0 )
write.table(up,'BM_vs_L-up-regulate.xls',sep='\t',quote = F,row.names = F)
down <- subset(diff4,log2FoldChange<0 )
write.table(down,'BM_vs_L-down-regulate.xls',sep='\t',quote = F,row.names = F)
#截取non数据
diff5<-subset(diff1b,pvalue > 0.05 )
write.table(diff5,'BM_vs_L-none-regulate.xls',sep='\t',quote = F,row.names = F)

9.1热图pheatmap
#处理*-diffGene.xls数据，计算lgfpkm平均值,取三列数据：gene，样品，对照
library(pheatmap)
data1 <- read.csv("BM_vs_L-diffGene.csv", sep=",")
Gene_Name = data1$gene
#BM_mean = data1$BM_mean
#L_mean = data1$L_mean
BM1 = data1$BM1
BM2 = data1$BM2
L1 = data1$L1
L2 = data1$L2
#data2 <- data.frame(Gene_Name,BM_mean,L_mean)
data2 <- data.frame(Gene_Name,BM1,BM2,L1,L2)
write.table(data2, file = "BM_vs_L_heatmap.csv", sep = ",", row.names = FALSE)
data <- read.csv("BM_vs_L_heatmap.csv", sep=",")
row.names(data) <- data$gene
data <- data[,-1]
data_matrix <- data.matrix(data)
colnames(data_matrix)=paste(c("BM1","BM1","L1","L2"),sep="")
pheatmap(data_matrix, Rowv=NA, Colv=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
revC=FALSE, show_rownames=F,scale='column',filename="BM_vs_L-heatmap.pdf",width=4,
height=8,main = "BM_vs_L-Heatmap",cluster_cols=F,cluster_rows = T)






































