# RNA-Seq analysis
# http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2

library(DESeq2)
library(DEFormats)
library(openxlsx)
library(apeglm)
library(vsn)
library(pheatmap)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
#--- Load DESeq2 files
##########
##### Child Trio 5 analysis
directory = file.path("HTSeq_output")
sampleFiles <- list.files(directory, pattern="*.counts*")
print(sampleFiles)

sampleNames <- c("T1_P","T4_M","T1_F","T5_P","T1_M","T5_F","T3_P","T5_M","T3_M","T7_P","T3_F","T7_F","T4_P","T7_M","T4_F","T5_M_2")
sampleCondition <- c("control","control","control","proband","control","control","control","control","control","control","control","control","control","control","control","control")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments =c("control","proband")

# Set the prefix for each output file name
outputPrefix <- "T5_Gmax_DESeq2"

snv <- read.csv("snv/trio5_mendelian.avinput.exonic_variant_function",sep="\t", header=FALSE)
sv <- read.csv("pacbio+novaseq_DNA17-06166.DEL.sorted.annotated.tsv",sep="\t")

##########
###### Child Trio 1 analysis
directory = file.path("HTSeq_output")
sampleFiles <- list.files(directory, pattern="*.counts*")
print(sampleFiles)

sampleNames <- c("T1_P","T4_M","T1_F","T5_P","T1_M","T5_F","T3_P","T5_M","T3_M","T7_P","T3_F","T7_F","T4_P","T7_M","T4_F","T5_M_2")
sampleCondition <- c("proband","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments =c("control","proband")

# Set the prefix for each output file name
outputPrefix <- "T1_Gmax_DESeq2"

snv <- read.csv("snv/trio1_mendelian.avinput.exonic_variant_function",sep="\t", header=FALSE)
sv <- read.csv("sv/hg38.DNA17-07724.pbsv2.20180823.sorted.annotated.tsv",sep="\t")


##########
###### Child Trio 3 analysis
directory = file.path("HTSeq_output")
sampleFiles <- list.files(directory, pattern="*.counts*")
print(sampleFiles)

sampleNames <- c("T1_P","T4_M","T1_F","T5_P","T1_M","T5_F","T3_P","T5_M","T3_M","T7_P","T3_F","T7_F","T4_P","T7_M","T4_F","T5_M_2")
sampleCondition <- c("control","control","control","control","control","control","proband","control","control","control","control","control","control","control","control","control")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments =c("control","proband")

# Set the prefix for each output file name
outputPrefix <- "T3_Gmax_DESeq2"

snv <- read.csv("snv/trio3_mendelian.avinput.exonic_variant_function",sep="\t", header=FALSE)
sv <- read.csv("sv/hg38.DNA17-06463.pbsv2.20180823.sorted.annotated.tsv",sep="\t")

##########
###### Child Trio 4 analysis
directory = file.path("HTSeq_output")
sampleFiles <- list.files(directory, pattern="*.counts*")
print(sampleFiles)

sampleNames <- c("T1_P","T4_M","T1_F","T5_P","T1_M","T5_F","T3_P","T5_M","T3_M","T7_P","T3_F","T7_F","T4_P","T7_M","T4_F","T5_M_2")
sampleCondition <- c("control","control","control","control","control","control","control","control","control","control","control","control","proband","control","control","control")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments =c("control","proband")

# Set the prefix for each output file name
outputPrefix <- "T4_Gmax_DESeq2"


##########
###### Child Trio 7 analysis
directory = file.path("HTSeq_output")
sampleFiles <- list.files(directory, pattern="*.counts*")
print(sampleFiles)

sampleNames <- c("T1_P","T4_M","T1_F","T5_P","T1_M","T5_F","T3_P","T5_M","T3_M","T7_P","T3_F","T7_F","T4_P","T7_M","T4_F","T5_M_2")
sampleCondition <- c("control","control","control","control","control","control","control","control","control","proband","control","control","control","control","control","control")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments =c("control","proband")

# Set the prefix for each output file name
outputPrefix <- "T7_Gmax_DESeq2"

snv <- read.csv("snv/trio7_mendelian.avinput.exonic_variant_function",sep="\t", header=FALSE)
sv <- read.csv("sv/hg38.DNA17-06575.pbsv2.20180823.sorted.annotated.tsv",sep="\t")

#the denovo longshot SNV already found on CG.
#Uploaded_variation	Location	Allele	Consequence	IMPACT	SYMBOL	Gene	Feature_type	Feature	BIOTYPE	EXON	INTRON	HGVSc	HGVSp	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	DISTANCE	STRAND	FLAGS	SYMBOL_SOURCE	HGNC_ID	TSL	APPRIS	REFSEQ_MATCH	GIVEN_REF	USED_REF	BAM_EDIT	SIFT	PolyPhen	AF	gnomAD_AF	gnomAD_AFR_AF	gnomAD_AMR_AF	gnomAD_ASJ_AF	gnomAD_EAS_AF	gnomAD_FIN_AF	gnomAD_NFE_AF	gnomAD_OTH_AF	gnomAD_SAS_AF	CLIN_SIG	SOMATIC	PHENO	PUBMED	MOTIF_NAME	MOTIF_POS	HIGH_INF_POS	MOTIF_SCORE_CHANGE	ada_score	rf_score
#70790653	7:70790653-70790653	A	missense_variant	MODERATE	AUTS2	26053	Transcript	NM_001127231.2	protein_coding	18/18	-	-	-	4100	3365	1122	G/E	gGa/gAa	rs150926322	-	1	-	EntrezGene	-	-	-	-	G	G	OK	tolerated(0.71)	probably_damaging(0.918)	-	3.742e-05	0	2.987e-05	0.0001029	5.847e-05	4.518e-05	2.797e-05	0	6.518e-05	likely_benign	-	-	25741868	-	-	-	-	-	-
#hg19 7       70255639        .       G       A       .       .       SYMBOL=AUTS2;STRAND=+;TYPE=E;DIST=-905;DS_AG=0.0000;DS_AL=0.0000;DS_DG=0.0000;DS_DL=0.0000;DP_AG=3;DP_AL=45;DP_DG=-3;DP_DL=-33
#other splicing variant in CG : 7       70252185        .       A       G       .       .       SYMBOL=AUTS2;STRAND=+;TYPE=I;DIST=10;DS_AG=0.0005;DS_AL=0.0000;DS_DG=0.0000;DS_DL=0.0000;DP_AG=-3;DP_AL=36;DP_DG=-5;DP_DL=2


##################################################
##################################################

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)
dds <- DESeq(ddsHTSeq)
res <- results(dds)

# order results by padj value (most significant to least)
res <- subset(res, padj<0.05)
res <- res[order(res$padj),]
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'

# add gene name
ensembl <- read.csv("mart_export.txt",sep="\t")
colnames(ensembl)[1] <- "gene"
resdata <- separate(resdata,gene,into=c("gene","version"), sep = "[:punct:]")
ens_resdata <- merge(resdata,ensembl[,c(1,4)],by="gene",all.x=TRUE)
ens_resdata <- unique(ens_resdata)
# add omim
omim <- read.csv("genemap2_20190103.txt",sep="\t")
colnames(omim)[11] <- "gene"

# T5
mim_ens_resdata <- merge(ens_resdata,omim[,c(11,13)],by="gene",all.x=TRUE)
mim_ens_resdata$min_except_P <- apply(mim_ens_resdata[,c(9:11,13:24)],1,min)
mim_ens_resdata$max_except_P <- apply(mim_ens_resdata[,c(9:11,13:24)],1,max)

mim_ens_resdata <- mim_ens_resdata[order(mim_ens_resdata$padj),]
ordered_mim_ens_resdata <- mim_ens_resdata[,c(1,25,26,4:8,12,14,16,24,3,27:28,9:11,13,15,17:23)]
ordered_mim_ens_resdata_1.6 <- subset(ordered_mim_ens_resdata, log2FoldChange <= -1.6)

# T7
mim_ens_resdata <- merge(ens_resdata,omim[,c(11,13)],by="gene",all.x=TRUE)
mim_ens_resdata$min_except_P <- apply(mim_ens_resdata[,c(9:17,19:24)],1,min)
mim_ens_resdata$max_except_P <- apply(mim_ens_resdata[,c(9:17,19:24)],1,max)

mim_ens_resdata <- mim_ens_resdata[order(mim_ens_resdata$padj),]
ordered_mim_ens_resdata <- mim_ens_resdata[,c(1,25,26,4:8,18,20,22,3,27:28,9:17,19,21,23)]
ordered_mim_ens_resdata_1.6 <- subset(ordered_mim_ens_resdata, log2FoldChange <= -1.6)


# graph for Alex Hoischen
mim_ens_resdata_alex <- mim_ens_resdata[mim_ens_resdata$HGNC.symbol == "AUTS2", c(9:24)]
mim_ens_resdata_alex <- t(mim_ens_resdata_alex)
mim_ens_resdata_alex <- as.data.frame(mim_ens_resdata_alex)
library(data.table)
mim_ens_resdata_alex <- setDT(mim_ens_resdata_alex, keep.rownames = TRUE)
colnames(mim_ens_resdata_alex)[1] <- "Sample"
colnames(mim_ens_resdata_alex)[2] <- "AUTS2_gene_expression"
mim_ens_resdata_alex$colors <- 1
mim_ens_resdata_alex[which(mim_ens_resdata_alex$AUTS2_gene_expression < 20) ,"colors"] <- 0

library("ggplot2")
library("ggrepel")
ggplot(mim_ens_resdata_alex, aes(x=Sample,y=AUTS2_gene_expression))  + 
  geom_point(size = 2,aes(color=as.factor(colors))) +
  geom_hline(yintercept = mean(mim_ens_resdata_alex$AUTS2_gene_expression), linetype ="dashed", color = "#FF4500") +
  geom_hline(yintercept = mean(mim_ens_resdata_alex$AUTS2_gene_expression) + (mean(mim_ens_resdata_alex$AUTS2_gene_expression)/(2**0.7010039)), linetype ="dashed", color = "darkgray") +
  geom_hline(yintercept = mean(mim_ens_resdata_alex$AUTS2_gene_expression) - (mean(mim_ens_resdata_alex$AUTS2_gene_expression)/(2**0.7010039)), linetype ="dashed", color = "darkgray") +
  geom_text(aes(0,mean(mim_ens_resdata_alex$AUTS2_gene_expression),label = "mean", vjust = -1, hjust=-0.5),colour="darkgray",size=3) +
  geom_text(aes(0,mean(mim_ens_resdata_alex$AUTS2_gene_expression) - (mean(mim_ens_resdata_alex$AUTS2_gene_expression)/(2**0.7010039)),label = "lfc standard error", vjust = +1, hjust=-0.5),colour="darkgray",size=3) +
  geom_text(aes(0,mean(mim_ens_resdata_alex$AUTS2_gene_expression) + (mean(mim_ens_resdata_alex$AUTS2_gene_expression)/(2**0.7010039)),label = "logarithm 2 fold change (lfc) standard error", vjust = -1, hjust=-0.2),colour="darkgray",size=3) +
  annotate("text", x = 15, y = 30, label = "padj= 0.00528" ,colour="#FF4500") +
  annotate("text", x = 15, y = 45, label = "lfc= -3.389",colour="#FF4500") + 
  theme(legend.position="none") +
  scale_color_manual(values=c("#FF4500","#C0C0C0"))+
  labs(x="Samples", y="RNA-seq AUTS2 normalized gene counts")

(mean(mim_ens_resdata_alex$AUTS2_gene_expression)/(2**0.7010039))

# add denovo SNV
snv <- separate(snv,V3,into=c("HGNC.symbol","NM","exon","c.","p."), sep = ":", perl = TRUE) 
snv <- unite(snv,snv.nomenclature,NM,exon,c.,p.,sep=":")
snv <- unite(snv,pos,V5,V6,sep="-")
snv <- unite(snv,subst,V7,V8,sep=">")
snv <- unite(snv,snv.genomicloc,V4,pos,subst,sep=":")

colnames(snv)[2] <- "snv.type"
ordered_mim_ens_resdata_1.6_snv <- merge(ordered_mim_ens_resdata_1.6,snv[,c(2:5)],by="HGNC.symbol",all.x=TRUE)

# add SV coding T7
sv <- subset(sv,AnnotSV.type == "split")
sv <- subset(sv,SV.type == "DEL")
colnames(sv)[16] <- "HGNC.symbol"
colnames(sv)[colnames(sv) == "X.hom.DNA_06575."] <- "Hom_SNV"
colnames(sv)[colnames(sv) == 'X.htz.DNA_06575.']<- "Het_SNV"
sv_sum <- sv %>%
  dplyr::group_by(HGNC.symbol) %>%
  dplyr::summarize(SV.PacBio=paste(unique(DNA17.06575),collapse=";"), SV.CDS=paste(unique(CDS.length),collapse=";"), SV.snv_hom_number=paste(unique(Hom_SNV),collapse=";"),SV.snv_het_number=paste(unique(Het_SNV),collapse=";") )
#ordered_mim_ens_resdata_1.6_snv_sv <- merge(ordered_mim_ens_resdata_1.6_snv,sv[,c(1:3,8:9,11:46)],by="HGNC.symbol",all.x=TRUE)
ordered_mim_ens_resdata_1.6_snv_sv <- merge(ordered_mim_ens_resdata_1.6_snv,sv_sum,by="HGNC.symbol",all.x=TRUE)
colnames(sv_sum)[1] <- "HGNC.symbol"
sv_rna <- merge(sv_sum,ordered_mim_ens_resdata_1.6_snv,by="HGNC.symbol",all.x=TRUE)


# add SV coding T5
sv <- subset(sv,AnnotSV.type == "split")
colnames(sv)[11] <- "HGNC.symbol"
colnames(sv)[45] <- "Hom_SNV"
colnames(sv)[46] <- "Het_SNV"
sv_sum <- sv %>%
  dplyr::group_by(HGNC.symbol) %>%
  dplyr::summarize(SV.PacBio=paste(unique(DNA17.06166),collapse=";"), SV.Illumina=paste(unique(BvB41_child_manta),collapse=";"), SV.CDS=paste(unique(CDS.length),collapse=";"), SV.snv_hom_number=paste(unique(Hom_SNV),collapse=";"),SV.snv_het_number=paste(unique(Het_SNV),collapse=";") )
#ordered_mim_ens_resdata_1.6_snv_sv <- merge(ordered_mim_ens_resdata_1.6_snv,sv[,c(1:3,8:9,11:46)],by="HGNC.symbol",all.x=TRUE)
ordered_mim_ens_resdata_1.6_snv_sv <- merge(ordered_mim_ens_resdata_1.6_snv,sv_sum,by="HGNC.symbol",all.x=TRUE)
colnames(sv_sum)[1] <- "HGNC.symbol"
sv_rna <- merge(sv_sum,ordered_mim_ens_resdata_1.6_snv,by="HGNC.symbol",all.x=TRUE)

# save table 
write.table(sv_rna, file = paste0(outputPrefix, "_gene_SV_rna.txt"), sep = '\t', row.names = FALSE, na = "")
write.table(ordered_mim_ens_resdata, file = paste0(outputPrefix, "-results-with-normalized.txt"), sep = '\t', row.names = FALSE)
write.table(ordered_mim_ens_resdata_1.6_snv_sv, file = paste0(outputPrefix, "-results-with-normalized-filterlog-1.6-snv-sv.txt"), sep = '\t', row.names = FALSE,  na = "")

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.table(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.txt"), sep = '\t')

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
dev.copy(png, paste0(outputPrefix, "-MAplot_initial_analysis.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld), file = paste0(outputPrefix, "-rlog-transformed-counts.txt"), sep = '\t'))
write.table(as.data.frame(assay(vsd), file = paste0(outputPrefix, "-vst-transformed-counts.txt"), sep = '\t'))


# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,paste0(outputPrefix, "-variance_stabilizing.png"))
dev.off()

# clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library("RColorBrewer")
library("gplots")
library(SummarizedExperiment)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
dev.copy(png, paste0(outputPrefix, "-clustering.png"))
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
dev.off()

#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

# set condition
condition <- treatments
scores <- data.frame(pc$x, condition)

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.9,.2),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

ggsave(pcaplot,file=paste0(outputPrefix, "-ggplot2.pdf"))

# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
head(assay(rld))
# plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,6:5],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

