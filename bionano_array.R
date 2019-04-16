#!/usr/bin/env Rscript

##### Library 
#tools is useful to write output with filename input.
#GenomicsRanges is used to find overlaps between SNP array and bionano calls.
#Tidyr and dplyr are used for data manipulation.
#karyoploteR and diffloop are used for data visualization.
#Xlsx is used to produce a excel output file.

library(tools)
library(tidyr)
library(dplyr)
library(GenomicRanges)
library(karyoploteR)
library(diffloop)
library(xlsx)


##### OS command line 
args = commandArgs(trailingOnly=TRUE)

# test if there is at least three argument: if not, return an error
if (length(args)!=3) {
  stop("You need 3 arguments, first is SNP array file, second in CNV bionano CSV file and third is SV bionano CSV file.", call.=FALSE)
}

arrayfile <- args[1]
cnv_nanofile <- args[2]
sv_nanofile <- args[3]

##### Function

### Pre processing raw data

# Get sample name

pureBasename <- function(file){
      require(tools)
  return(file_path_sans_ext(basename(file)))
}

# Process CNV output from ChAS
arraycnv <- function(arrayfile,hg)
{
  require("tidyr")
  arraycnv <- read.delim(arrayfile,sep="\t", header=TRUE, na.strings = "", fill=TRUE)
  arraycnv_clean <- arraycnv %>%
    separate(col= Microarray.Nomenclature..ISCN.2016.,sep="\\(",into=c("Cytoband","coordinate")) %>%
    separate(col= coordinate,sep="\\)",into=c("startend","Copy.Number")) %>%
    separate(col= startend,sep="\\_",into=c("Start","End")) %>%
    dplyr::select(Chromosome,Start,End, Type, CN.State, Size..kbp.,Marker.Count,Gene.Count,Genes,File)
  arraycnv_clean <- tibble::rowid_to_column(arraycnv_clean, "CNV.array.ID")
  return(arraycnv_clean)
}


# Process SV output from bionano access 

nanosv <- function(sv_nanofile,hg)
{
  require("tidyr")
  nanosv <- read.delim(sv_nanofile,sep=",", header=TRUE, na.strings = "", fill=TRUE)
  nanosv_clean <- nanosv %>%
    dplyr::select(RefcontigID1, RefStartPos, RefEndPos, Confidence, Type, Zygosity, Size, Present_in_._of_BNG_control_samples_with_the_same_enzyme, Fail_assembly_chimeric_score )
  nanosv_clean <- dplyr::rename(nanosv_clean, Chromosome = RefcontigID1, Start = RefStartPos, End = RefEndPos, bionano_sv_Type = Type)
  nanosv_clean$Start <- as.integer(nanosv_clean$Start)
  nanosv_clean$End <- as.integer(nanosv_clean$End)
  nanosv_clean <- nanosv_clean[nanosv_clean$End - nanosv_clean$Start >= 0,]
  nanosv_clean$Chromosome[nanosv_clean$Chromosome == 23] <- "X"
  nanosv_clean$Chromosome[nanosv_clean$Chromosome == 24] <- "Y"
  return(nanosv_clean)
}

# Process CNV output from bionano access

nanocnv <- function(cnv_nanofile,hg)
{
  nanocnv <- read.csv(cnv_nanofile,sep=",", header=TRUE)
  colnames(nanocnv)[6] <- "bionano_cnv_Type"
  nanocnv$Chromosome[nanocnv$Chromosome == 23] <- "X"
  nanocnv$Chromosome[nanocnv$Chromosome == 24] <- "Y"
  #
  nanocnv$Start <- trunc(nanocnv$Start)
  nanocnv$End <- trunc(nanocnv$End)
  nanocnv$Width <- trunc(nanocnv$Width)
  #
  nanocnv$Start <- as.integer(nanocnv$Start)
  nanocnv$End <- as.integer(nanocnv$End)
  nanocnv$Width <- as.integer(nanocnv$Width)
  return(nanocnv)
} 
 
### Comparative function for SNP array calls <- SV & CNV bionano calls

# Merge the CNV bionano output with SNP array calls

nanoarraycnv <- function(cnv_nanofile,arrayfile,hg)
{
  require("tidyr")
  require("GenomicRanges")
  # Import file and transform into GenomicsRange 
  arraycnv <- arraycnv(arrayfile,hg)
  arraycnv.GenomicRanges <- makeGRangesFromDataFrame(arraycnv, keep.extra.columns=TRUE,ignore.strand=TRUE)
  
  nanocnv <- nanocnv(cnv_nanofile,hg)
  nanocnv.GenomicRanges <- makeGRangesFromDataFrame(nanocnv, keep.extra.columns=TRUE,ignore.strand=TRUE)
  
  # Find overlaps : SNP array calls <- CNV bionano calls 
  hits <- findOverlaps(nanocnv.GenomicRanges, arraycnv.GenomicRanges, type="any")
  nanoarraycnv.Ann <- cbind(arraycnv[subjectHits(hits),], nanocnv.GenomicRanges[queryHits(hits),])
  
  # Get complete CNV array file with and without bionano CNV overlap
  nanoarraycnv <- dplyr::full_join(x=arraycnv,y=nanoarraycnv.Ann)
  nanoarraycnv <- unique(nanoarraycnv)
  
  # Get overlap and size ratio between common CNV calls
  nanoarraycnv <- dplyr::mutate(nanoarraycnv, overlap_ratio = pmin(((as.numeric(pmin(End,end)) - as.numeric(pmax(Start,start)))) / (Size..kbp. * 1000), 1), size_ratio =  as.numeric(width) / (as.numeric(Size..kbp.) * 1000) )
  
  colnames(nanoarraycnv)[12:ncol(nanoarraycnv)] <- paste0("cnv.bng_", colnames(nanoarraycnv)[12:ncol(nanoarraycnv)])
  return(nanoarraycnv)
}

# Merge the SV bionano output with SNP array calls

nanoarraysv <- function(sv_nanofile,arrayfile,hg)
{
  require("tidyr")
  require("GenomicRanges")
  # Import file and transform into GenomicsRange 
  arraycnv <- arraycnv(arrayfile,hg)
  arraycnv.GenomicRanges <- makeGRangesFromDataFrame(arraycnv, keep.extra.columns=TRUE,ignore.strand=TRUE)
  
  nanosv <- nanosv(sv_nanofile,hg)
  nanosv.GenomicRanges <- makeGRangesFromDataFrame(nanosv, keep.extra.columns=TRUE,ignore.strand=TRUE)

  # Find overlaps : SNP array calls <- SV bionano calls 
  hits <- findOverlaps(nanosv.GenomicRanges, arraycnv.GenomicRanges, type="any")
  nanoarraysv.Ann <- cbind(arraycnv[subjectHits(hits),], nanosv.GenomicRanges[queryHits(hits),])
  
  # Get complete CNV array file with and without bionano SV overlap
  nanoarraysv <- dplyr::full_join(arraycnv,nanoarraysv.Ann)
  nanoarraysv <- unique(nanoarraysv)
  
  # Get overlap and size ratio between common CNV calls
  nanoarraysv <- dplyr::mutate(nanoarraysv, overlap_ratio = pmin(((as.numeric(pmin(End,end)) - as.numeric(pmax(Start,start)))) / (Size..kbp. * 1000),1), size_ratio =  as.numeric(Size) / (as.numeric(Size..kbp.) * 1000) )
  colnames(nanoarraysv)[12:ncol(nanoarraysv)] <- paste0("sv.bng_", colnames(nanoarraysv)[12:ncol(nanoarraysv)])
  return(nanoarraysv)
}

# Merge the results from the 2 output file with SNP Array

nanoarray_cnv_sv <- function(cnv_nanofile, sv_nanofile, arrayfile, hg)
{
  require("tidyr")
  # get input file names for output file name
  prefix <- pureBasename(sv_nanofile)
  # merge technologies 
  nanoarraycnv <- nanoarraycnv(cnv_nanofile,arrayfile,hg)
  nanoarraysv <- nanoarraysv(sv_nanofile,arrayfile,hg)
  nanoarray_cnv_sv  <- dplyr::full_join(nanoarraycnv,nanoarraysv)
  nanoarray_cnv_sv  <- unique(nanoarray_cnv_sv)
  return(nanoarray_cnv_sv)
}

#### Comparative function for SV & CNV bionano calls <- SNP array calls

# CNV bionano output vs SNP array calls 

nano_noarraycnv <- function(cnv_nanofile,arrayfile,hg)
{
  require("tidyr")
  require("GenomicRanges")
  # Import file and transform into GenomicsRange
  arraycnv <- arraycnv(arrayfile,hg)
  arraycnv.GenomicRanges <- makeGRangesFromDataFrame(arraycnv, keep.extra.columns=TRUE,ignore.strand=TRUE)
  
  nanocnv <- nanocnv(cnv_nanofile,hg)
  colnames(nanocnv)[2] <- "chromosome"
  colnames(nanocnv)[3] <- "start"
  colnames(nanocnv)[4] <- "end"
  colnames(nanocnv)[5] <- "width"
  nanocnv.GenomicRanges <- makeGRangesFromDataFrame(nanocnv, keep.extra.columns=TRUE,ignore.strand=TRUE)
  
  # Find overlaps between bionano CNV calls and SNP array
  hits <- findOverlaps(arraycnv.GenomicRanges, nanocnv.GenomicRanges, type="any")
  nano_noarraycnv.Ann <- cbind(nanocnv.GenomicRanges[subjectHits(hits),], arraycnv[queryHits(hits),])
  colnames(nano_noarraycnv.Ann)[1] <- "chromosome"
  # remove these 2 columns for merging since width computation is 1bp different from each pipeline.
  nano_noarraycnv.Ann$strand <- NULL
  nano_noarraycnv.Ann$width <- NULL
  
  # Get complete bionano CNV file with and without CNV from SNP array overlap
  nano_noarraycnv <- dplyr::full_join(x=nanocnv,y= nano_noarraycnv.Ann)
  nano_noarraycnv <- unique(nano_noarraycnv)

  # Get overlap and size ratio between common CNV calls
  nano_noarraycnv <- dplyr::mutate(nano_noarraycnv, overlap_ratio = pmin(((as.numeric(pmin(End,end)) - as.numeric(pmax(Start,start)))) / (Size..kbp. * 1000), 1), size_ratio =  as.numeric(width) / (as.numeric(Size..kbp.) * 1000) )

  return(nano_noarraycnv)
}

# SV bionano output vs SNP array calls 
 
nano_noarraysv <- function(sv_nanofile,arrayfile,hg)
{
  require("tidyr")
  require("GenomicRanges")
  # Import file and transform into GenomicsRange 
  arraycnv <- arraycnv(arrayfile,hg)
  arraycnv.GenomicRanges <- makeGRangesFromDataFrame(arraycnv, keep.extra.columns=TRUE,ignore.strand=TRUE)
  
  nanosv <- nanosv(sv_nanofile,hg)
  nanosv$Chromosome[nanosv$Chromosome == 23] <- "X"
  nanosv$Chromosome[nanosv$Chromosome == 24] <- "Y"
  colnames(nanosv)[1] <- "chromosome"
  colnames(nanosv)[2] <- "start"
  colnames(nanosv)[3] <- "end"
  nanosv.GenomicRanges <- makeGRangesFromDataFrame(nanosv, keep.extra.columns=TRUE,ignore.strand=TRUE)
  
  # Find overlaps between bionano CNV calls and SNP array
  hits <- findOverlaps(arraycnv.GenomicRanges, nanosv.GenomicRanges, type="any")
  nano_noarraysv.Ann <- cbind(nanosv.GenomicRanges[subjectHits(hits),], arraycnv[queryHits(hits),])
  colnames(nano_noarraysv.Ann)[1] <- "chromosome"  
  nano_noarraysv.Ann$strand <- NULL
  
  # Get complete bionano CNV file with and without SNP overlap
  nano_noarraysv <- dplyr::full_join(x=nanosv,y= nano_noarraysv.Ann)
  nano_noarraysv <- unique(nano_noarraysv)

  # Get overlap and size ratio between common CNV calls
  nano_noarraysv <- dplyr::mutate(nano_noarraysv, overlap_ratio = pmin(((as.numeric(pmin(End,end)) - as.numeric(pmax(Start,start)))) / (Size..kbp. * 1000), 1), size_ratio =  as.numeric(Size) / (as.numeric(Size..kbp.) * 1000) )

  return(nano_noarraysv)
}

#### Visualization


print_karyoplot <- function(arrayfile, cnv_nanofile, sv_nanofile,hg) {
  require(karyoploteR)
  require(diffloop)
  
  # process array data
  arraycnv <- arraycnv(arrayfile,hg)
  arraycnv.gain <- arraycnv[arraycnv$Type== "Gain", ]
  arraycnv.loss <- arraycnv[arraycnv$Type== "Loss", ]
  
  arraycnv.gain.GenomicRanges <- makeGRangesFromDataFrame(arraycnv.gain, keep.extra.columns=TRUE,ignore.strand=TRUE,starts.in.df.are.0based=TRUE)
  arraycnv.gain.GenomicRanges <- addchr(arraycnv.gain.GenomicRanges)
  arraycnv.loss.GenomicRanges <- makeGRangesFromDataFrame(arraycnv.loss, keep.extra.columns=TRUE,ignore.strand=TRUE,starts.in.df.are.0based=TRUE)
  arraycnv.loss.GenomicRanges <- addchr(arraycnv.loss.GenomicRanges)
  
  # process sv bionano data
  nanosv <- nanosv(sv_nanofile,hg)
  nanosv$Chromosome[nanosv$Chromosome == 23] <- "X"
  nanosv$Chromosome[nanosv$Chromosome == 24] <- "Y"
  nanosv <- unique(nanosv)
  nanosv.gain <- nanosv[nanosv$bionano_sv_Type != "deletion",]
  nanosv.loss <- nanosv[nanosv$bionano_sv_Type == "deletion",]

  nanosv.gain.GenomicRanges <- makeGRangesFromDataFrame(nanosv.gain, keep.extra.columns=TRUE,ignore.strand=TRUE,starts.in.df.are.0based=TRUE)
  nanosv.gain.GenomicRanges <- addchr(nanosv.gain.GenomicRanges)
  nanosv.loss.GenomicRanges <- makeGRangesFromDataFrame(nanosv.loss, keep.extra.columns=TRUE,ignore.strand=TRUE,starts.in.df.are.0based=TRUE)
  nanosv.loss.GenomicRanges <- addchr(nanosv.loss.GenomicRanges)
  
  # process cnv bionano data
  nanocnv <- read.csv(cnv_nanofile,sep=",", header=TRUE)
  colnames(nanocnv)[6] <- "bionano_cnv_Type"
  nanocnv$Chromosome[nanocnv$Chromosome == 23] <- "X"
  nanocnv$Chromosome[nanocnv$Chromosome == 24] <- "Y"
  nanocnv.GenomicRanges <- makeGRangesFromDataFrame(nanocnv, keep.extra.columns=TRUE,ignore.strand=TRUE,starts.in.df.are.0based=TRUE)
  nanocnv.gain <- nanocnv[nanocnv$bionano_cnv_Type != "deletion",]
  nanocnv.loss <- nanocnv[nanocnv$bionano_cnv_Type == "deletion",]

  nanocnv.gain.GenomicRanges <- makeGRangesFromDataFrame(nanocnv.gain, keep.extra.columns=TRUE,ignore.strand=TRUE,starts.in.df.are.0based=TRUE)
  nanocnv.gain.GenomicRanges <- addchr(nanocnv.gain.GenomicRanges)
  nanocnv.loss.GenomicRanges <- makeGRangesFromDataFrame(nanocnv.loss, keep.extra.columns=TRUE,ignore.strand=TRUE,starts.in.df.are.0based=TRUE)
  nanocnv.loss.GenomicRanges <- addchr(nanocnv.loss.GenomicRanges)
  
  # print plot as file
  prefix <- pureBasename(sv_nanofile)
  png(paste0(prefix,"_karyoplotR.png"), width = 3440, height = 4000,res=450)
  
  kp <- plotKaryotype(genome=hg)
  kpPlotRegions(kp, data=arraycnv.gain.GenomicRanges ,col="#005073", layer.margin = 0.01, border=NA, r0=0, r1=0.25)
  kpPlotRegions(kp, data=arraycnv.loss.GenomicRanges ,col="#d50000", layer.margin = 0.01, border=NA, r0=0, r1=0.25)
  kpAddLabels(kp, labels="array",  r0=0, r1=0.25, cex=0.2, col="#C0C0C0")
  
  kpDataBackground(kp, r0=0.35, r1=0.6, col= "#fafafa")
  kpPlotRegions(kp, data=nanocnv.gain.GenomicRanges ,col="#189ad3", layer.margin = 0.01, border=NA, r0=0.35, r1=0.6)
  kpPlotRegions(kp, data=nanocnv.loss.GenomicRanges ,col="#ff5252", layer.margin = 0.01, border=NA, r0=0.35, r1=0.6)
  kpAddLabels(kp, labels="bng_cnv", r0=0.35, r1=0.6, cex=0.2, col="#C0C0C0")

  kpDataBackground(kp, r0=0.75, r1=1, col= "#e0e0e0")
  kpPlotRegions(kp, data=nanosv.gain.GenomicRanges ,col="#107dac", layer.margin = 0.01, border=NA, r0=0.75, r1=1, avoid.overlapping = FALSE)
  kpPlotRegions(kp, data=nanosv.loss.GenomicRanges ,col="#ff1744", layer.margin = 0.01, border=NA, r0=0.75, r1=1, avoid.overlapping = FALSE)
  kpAddLabels(kp, labels="bng_sv", r0=0.75, r1=1, cex=0.2, col="#C0C0C0")
  
  kpAddCytobandLabels(kp, col="gray", cex=0.2)
  dev.off()
}

#### Statistics

# Function to create a unique Excel file output

xlsx.writeMultipleData <- function (file, ...)
  {
    require(xlsx, quietly = TRUE)
    objects <- list(...)
    fargs <- as.list(match.call(expand.dots = TRUE))
    objnames <- as.character(fargs)[-c(1, 2)]
    nobjects <- length(objects)
    for (i in 1:nobjects) {
        if (i == 1)
            write.xlsx2(objects[[i]], file, sheetName = objnames[i],showNA=FALSE, row.names=FALSE)
        else write.xlsx2(objects[[i]], file, sheetName = objnames[i], row.names=FALSE, showNA=FALSE,
            append = TRUE)
    }
  }

nanoarray_cnv_sv_statistics <- function(cnv_nanofile, sv_nanofile, arrayfile, hg) {  
  
  # get input file names to get Sample name
  prefix <- pureBasename(sv_nanofile)
  
  # import data for descriptive statistics
  array_calls <- nanoarray_cnv_sv(cnv_nanofile, sv_nanofile, arrayfile, hg)
  sv_bionano_calls <- nano_noarraysv(sv_nanofile,arrayfile,hg)
  cnv_bionano_calls <- nano_noarraycnv(cnv_nanofile,arrayfile,hg)
  
  # Count CNV event on differents filters
  only_array <- nrow(unique(array_calls[is.na(array_calls$cnv.bng_seqnames) == TRUE & is.na(array_calls$sv.bng_seqnames) ==TRUE,c(1:12)]))
  common_nanoarray <- nrow(unique(array_calls[is.na(array_calls$cnv.bng_seqnames) == FALSE | is.na(array_calls$sv.bng_seqnames) == FALSE ,c(1:12)]))
  total_array <- nrow(unique(arraycnv(arrayfile)))
  
  only_nanosv <- nrow(unique(sv_bionano_calls[is.na(sv_bionano_calls$CNV.array.ID) == TRUE,]))
  total_nanosv <- nrow(unique(nanosv(sv_nanofile)))
  common_nanosv <- nrow(unique(sv_bionano_calls[is.na(sv_bionano_calls$CNV.array.ID) == FALSE ,c(1:9)]))
  
  only_nanocnv <- nrow(unique(cnv_bionano_calls[is.na(cnv_bionano_calls$CNV.array.ID) == TRUE,]))
  total_nanocnv <- nrow(unique(nanocnv(cnv_nanofile,hg)))
  common_nanocnv <- nrow(unique(cnv_bionano_calls[is.na(cnv_bionano_calls$CNV.array.ID) == FALSE ,c(1:9)]))
  
  # Gather all counts into a data frame
  results <- data.frame(Sample = prefix, Array_total_calls = total_array, Array_common_calls = common_nanoarray, Array_only_calls = only_array, SV_bionano_total_calls = total_nanosv, SV_bionano_common_calls = common_nanosv , SV_bionano_only_calls = only_nanosv, CNV_bionano_total_calls = total_nanocnv , CNV_bionano_common_calls = common_nanocnv, CNV_bionano_only_calls = only_nanocnv)
  return(results)
}

#### Complete function with comparison and plot

nanoarray_cnv_sv_plot_statistics <- function(cnv_nanofile, sv_nanofile, arrayfile, hg) {

   # get input file names for output file name
  prefix <- pureBasename(sv_nanofile)
  nanocnv <- pureBasename(cnv_nanofile)
  arrayname <- pureBasename(arrayfile)
  
  # load processed data with previous functions
  array_calls <- nanoarray_cnv_sv(cnv_nanofile, sv_nanofile, arrayfile, hg)
  sv_bionano_calls <- nano_noarraysv(sv_nanofile,arrayfile,hg)
  cnv_bionano_calls <- nano_noarraycnv(cnv_nanofile,arrayfile,hg)
  array_bionano_statistics <- nanoarray_cnv_sv_statistics(cnv_nanofile, sv_nanofile, arrayfile, hg)
  
  # Write data into one Excel file
  xlsx.writeMultipleData(paste0(prefix,"_SNParray_statistics.xlsx"), array_bionano_statistics, array_calls, sv_bionano_calls, cnv_bionano_calls)
  print_karyoplot(arrayfile, cnv_nanofile, sv_nanofile,hg)
  print(paste0("CNV comparaison with ",arrayname, " SNP array is done in ",hg," with ",prefix," and ",nanocnv,"."))
}

### Run

nanoarray_cnv_sv_plot_statistics(cnv_nanofile, sv_nanofile ,arrayfile, hg="hg19")




