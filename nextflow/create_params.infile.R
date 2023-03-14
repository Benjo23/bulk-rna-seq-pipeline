#set nextflow working directory
setwd("")

#sample id
#directory that contains fastq files
SAMPLE_ID<-list.files(path="", pattern=NULL, all.files=FALSE, full.names=FALSE)
SAMPLE_ID<-sub("_R1_001.*", "",SAMPLE_ID)  
SAMPLE_ID<-sub("_R2_001.*", "",SAMPLE_ID)
SAMPLE_ID<-unique(SAMPLE_ID)

#full path of R1 & R2 
R1<-paste0("",SAMPLE_ID,"_paired_R1.fastq.gz")
R2<-paste0("",SAMPLE_ID,"_paired_R2.fastq.gz")
  
#create df
df<- data.frame (INDIVIDUAL_ID = c(""),
                 SAMPLE_ID = SAMPLE_ID ,
                 LIBRARY_ID=c(""),
                 RG_ID=c(""),
                 PLATFORM_UNIT=(""),
                 PLATFORM=(""),
                 PLATFORM_MODEL=(""),
                 RUN_DATE=(""),
                 CENTER=(""),
                 R1 = R1,
                 R2= R2
)

df<-df[-1,]

#save tsv file
library(readr)
write_tsv(df, file='')
