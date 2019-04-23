Compound htz and bi allelic coding SV
================
Kevin Yauy

  - [Library](#library)
  - [Function](#function)
  - [Run](#run)

Radboudumc

Date of publication : 2019-04-17

## Library

``` r
library(tidyr)
library(dplyr)
```

## Function

Biallelic SV with coding regions

``` r
biSV <- function(file,sample)
{
  require("dplyr")
  sv <- read.delim(file,sep="\t", header=TRUE)
  sv_split <- sv[sv$AnnotSV.type =="split",]
  bisv <- sv_split %>%  filter(CDS.length > 0) %>% group_by(Gene.name) %>% summarise(n= n_distinct(Gene.name)) %>%  filter(n > 1)
  return(bisv)
}
```

Coding SV with predicted SNV damaging

``` r
SNVSV <- function(file,sample)
{
  require("tidyr")
  require("dplyr")
  sv <- read.delim(file,sep="\t", header=TRUE, na.strings = "", fill=TRUE)
  sv_split <- sv[sv$AnnotSV.type =="split",]
  svsnv <- sv_split %>%  filter(CDS.length > 0)  %>% select(AnnotSV.ID,Gene.name,CDS.length,contains("compound"))  %>% select(AnnotSV.ID,Gene.name, CDS.length, contains(sample)) %>% drop_na()
  write.table(svsnv,paste0(sample,"_codingSV_SNVdamaging.tsv"), sep="\t", row.names = FALSE)
  return(svsnv)
}
```

## Run

Biallelic SV with coding regions

    ## # A tibble: 0 x 2
    ## # ... with 2 variables: Gene.name <fct>, n <int>

    ## # A tibble: 0 x 2
    ## # ... with 2 variables: Gene.name <fct>, n <int>

    ## # A tibble: 0 x 2
    ## # ... with 2 variables: Gene.name <fct>, n <int>

    ## # A tibble: 0 x 2
    ## # ... with 2 variables: Gene.name <fct>, n <int>

    ## # A tibble: 0 x 2
    ## # ... with 2 variables: Gene.name <fct>, n <int>

Coding SV with predicted SNV damaging

    ##                    AnnotSV.ID Gene.name CDS.length compound.htz.DNA.06166.
    ## 77    6_32471455_32528184_DEL  HLA-DRB5        706              6_32530128
    ## 78    6_32485221_32571277_DEL  HLA-DRB5        807              6_32530128
    ## 80    6_32487059_32573147_DEL  HLA-DRB5        807              6_32530128
    ## 82    6_32524834_32584981_DEL  HLA-DRB1        706              6_32584364
    ## 83    6_32524834_32584981_DEL  HLA-DRB5        101              6_32530128
    ## 85    6_32526199_32586366_DEL  HLA-DRB1        706              6_32584364
    ## 86    6_32526199_32586366_DEL  HLA-DRB5        101              6_32530128
    ## 93  7_100958296_100958325_DEL     MUC3A         30             7_100960609
    ## 94  7_100958300_100958301_INS     MUC3A          2             7_100960609
    ## 150  12_10882513_10882574_DEL      PRH1         62             12_11061546
    ## 175  15_20535329_20535370_DEL  GOLGA6L6         42             15_20535014
    ## 176  15_23130316_23130357_DEL  GOLGA6L1         42             15_23129472
    ## 177  15_23130316_23130357_DEL GOLGA6L22         42             15_23129472
    ## 178  15_23136600_23136646_DEL  GOLGA6L1         47             15_23129472
    ## 179  15_23327869_23327870_INS GOLGA6L22          2             15_23129472
    ## 180  15_23334090_23334091_INS GOLGA6L22          2             15_23129472
    ## 254    19_2939220_2939290_DEL     ZNF77         11              19_2936537

    ##                    AnnotSV.ID Gene.name CDS.length compound.htz.DNA_06463.
    ## 78  7_100958296_100958325_DEL     MUC3A         30             7_100959263
    ## 79  7_100958300_100958301_INS     MUC3A          2             7_100959263
    ## 143  12_10882513_10882574_DEL      PRH1         62             12_11061546

    ##                   AnnotSV.ID Gene.name CDS.length compound.htz.DNA_06468.
    ## 36 3_195779684_195779685_INS      MUC4          2             3_195778790
    ## 37 3_195780226_195780227_INS      MUC4          2             3_195778790
    ## 38 3_195781258_195781259_INS      MUC4          2             3_195778790
    ## 39 3_195781275_195781322_DEL      MUC4         48             3_195778790
    ## 40 3_195781570_195781665_DEL      MUC4         96             3_195778790
    ## 41 3_195781618_195781665_DEL      MUC4         48             3_195778790
    ## 42 3_195782080_195782415_DEL      MUC4        336             3_195778790
    ## 43 3_195783304_195783783_DEL      MUC4        480             3_195778790
    ## 44 3_195784486_195784581_DEL      MUC4         96             3_195778790
    ## 45 3_195785035_195785132_DEL      MUC4         98             3_195778790
    ## 46 3_195785085_195785132_DEL      MUC4         48             3_195778790
    ## 47 3_195785729_195785730_INS      MUC4          2             3_195778790
    ## 48 3_195786295_195788165_DEL      MUC4       1873             3_195778790
    ## 49 3_195786527_195786528_INS      MUC4          3             3_195778790
    ## 50 3_195786928_195786975_DEL      MUC4         48             3_195778790
    ## 51 3_195788495_195788542_DEL      MUC4         48             3_195778790
    ## 78 7_100958296_100958325_DEL     MUC3A         30             7_100960609
    ## 79 7_100958300_100958301_INS     MUC3A          2             7_100960609

    ##                   AnnotSV.ID Gene.name CDS.length compound.htz.DNA_06575.
    ## 135 12_10882513_10882574_DEL      PRH1         62             12_11061546
    ## 206   18_9887391_9887435_DEL    TXNDC2         45              18_9887032
    ## 225 19_55358511_55358534_DEL   FAM71E2         24             19_55359721

    ##                    AnnotSV.ID Gene.name CDS.length compound.htz.DNA_07724.
    ## 71    6_32471455_32528184_DEL  HLA-DRB5        706              6_32521954
    ## 73    6_32524834_32584981_DEL  HLA-DRB5        101              6_32521954
    ## 76    6_32526199_32586366_DEL  HLA-DRB5        101              6_32521954
    ## 101 9_137050283_137050284_INS    ENTPD2          2             9_137053898
    ## 136  12_10882513_10882574_DEL      PRH1         62             12_11061546
    ## 137  12_11308619_11308620_INS      PRB4          2             12_11308868
    ## 159  15_20534338_20534505_DEL  GOLGA6L6        168             15_20535014
    ## 160  15_23327869_23327870_INS GOLGA6L22          2             15_23129688
    ## 161  15_23334090_23334091_INS GOLGA6L22          2             15_23129688
    ## 211    19_2939220_2939290_DEL     ZNF77         11              19_2936537
