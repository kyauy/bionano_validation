Deletion\_accuracy
================
Kevin Yauy

  - [Library](#library)
  - [Experiments](#experiments)
      - [pbsv2 with longshot](#pbsv2-with-longshot)

Radboudumc

Date of publication : 2019-04-17

## Library

Tidyr and dplyr are used for data manipulation.

``` r
library(tidyr)
library(dplyr)
```

## Experiments

### pbsv2 with longshot

#### Function

``` r
valid_SV <- function(file,sample,tech)
{
  require("dplyr")
  valid_sv <- read.delim(file,sep="\t", header=TRUE)
  valid_sv_del  <-  valid_sv[valid_sv$SV.type == "DEL",]
  valid_sv_del_full <-  valid_sv_del[valid_sv_del$AnnotSV.type == "full",]
  
  valid_sv_del_full <- dplyr::select(valid_sv_del_full, grep(sample,colnames(valid_sv_del_full)) )
  
  colnames(valid_sv_del_full)[1] <- paste("longshot.hom.DNA-",sample,sep="")
  colnames(valid_sv_del_full)[2] <- paste("longshot.htz.DNA-",sample,sep="")
  
  hom <- paste("longshot.hom.DNA-",sample,sep="")
  het <- paste("longshot.htz.DNA-",sample,sep="")

  valid_sv_del_full$longshot.diff <- valid_sv_del_full[,hom] - valid_sv_del_full[,het]

  #longshot_valid_sv_snv_ratio <- sum(valid_sv_del_full[,hom]) / (sum(valid_sv_del_full[,hom]) + sum(valid_sv_del_full[,het]) )
  #longshot_valid_sv_snv_hethomratio <- sum(valid_sv_del_full[,het]) / sum(valid_sv_del_full[,hom])
  #print(valid_sv_del_full)
  soft <- nrow(valid_sv_del_full[valid_sv_del_full$longshot.diff>0, ])
  noinfo <- nrow(valid_sv_del_full[valid_sv_del_full$longshot.diff == 0, ])
  filtered <- nrow(valid_sv_del_full[valid_sv_del_full$longshot.diff<0, ])
  total <- nrow(valid_sv_del_full)
  hard <- nrow(valid_sv_del_full[valid_sv_del_full$longshot.diff>0 & valid_sv_del_full[,2] <= 1 , ])
  #hard <- valid_sv_del_full %>% filter(longshot.diff > 0 & grep(het,colnames(valid_sv_del_full)) <= 1) %>% count()

  results <- data.frame(Sample=sample, DEL.total = total, ls.DEL.PASS.soft = soft, ls.DEL.PASS.hard = hard, ls.DEL.no.info = noinfo, ls.DEL.filtered= filtered )
  colnames(results)[2:ncol(results)] <- paste(tech,colnames(results)[2:ncol(results)],sep="_")
  return(results)
}

valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06166.pbsv2.20180823.sorted.annotated.tsv",sample ="06166.",tech="pbsv2")
```

    ##   Sample pbsv2_DEL.total pbsv2_ls.DEL.PASS.soft pbsv2_ls.DEL.PASS.hard
    ## 1 06166.           10579                   1942                   1699
    ##   pbsv2_ls.DEL.no.info pbsv2_ls.DEL.filtered
    ## 1                 6686                  1951

#### Run

``` r
s06166 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06166.pbsv2.20180823.sorted.annotated.tsv",sample ="06166.",tech="pbsv2")
s06167 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06167.pbsv2.20180823.sorted.annotated.tsv",sample ="06167.",tech="pbsv2")
s06168 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06168.pbsv2.20180823.sorted.annotated.tsv",sample ="06168.",tech="pbsv2")
s06463 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06463.pbsv2.20180823.sorted.annotated.tsv",sample ="06463.",tech="pbsv2")
s06464 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06464.pbsv2.20180823.sorted.annotated.tsv",sample ="06464.",tech="pbsv2")
s06467 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06467.pbsv2.20180823.sorted.annotated.tsv",sample ="06467.",tech="pbsv2")
s06468 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06468.pbsv2.20180823.sorted.annotated.tsv",sample ="06468.",tech="pbsv2")
s06469 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06469.pbsv2.20180823.sorted.annotated.tsv",sample ="06469.",tech="pbsv2")
s06470 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06470.pbsv2.20180823.sorted.annotated.tsv",sample ="06470.",tech="pbsv2")
s06575 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06575.pbsv2.20180823.sorted.annotated.tsv",sample ="06575.",tech="pbsv2")
s06576 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06576.pbsv2.20180823.sorted.annotated.tsv",sample ="06576.",tech="pbsv2")
s06577 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-06577.pbsv2.20180823.sorted.annotated.tsv",sample ="06577.",tech="pbsv2")
s07724 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-07724.pbsv2.20180823.sorted.annotated.tsv",sample ="07724.",tech="pbsv2")
s07725 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-07725.pbsv2.20180823.sorted.annotated.tsv",sample ="07725.",tech="pbsv2")
s07726 <- valid_SV(file= "~/SMRT_SV_Project/bionano/pbsv2/hg38.DNA17-07726.pbsv2.20180823.sorted.annotated.tsv",sample ="07726.",tech="pbsv2")


all_trio <- rbind(s06166,s06167,s06168,s06463,s06464,s06467,s06468,s06469,s06470,s06575,s06576,s06577,s07724,s07725,s07726)

all_trio$soft_perc <- all_trio[,3] *100 / all_trio[,2]
all_trio$hard_perc <- all_trio[,4] *100 / all_trio[,2]
all_trio$noinfo_perc <- all_trio[,5] *100 / all_trio[,2]
all_trio$filtered_perc <- all_trio[,6] *100 / all_trio[,2]

all_trio_order <- all_trio[,c(1,2,3,7,4,8,5,9,6,10)]

write.table(all_trio_order,"15trio_pbsv2_DEL.tsv", na = ".", sep="\t",row.names = FALSE)
```
