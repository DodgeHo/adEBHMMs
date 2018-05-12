library("sleuth")
library("EBSeq")
library("EBSeqHMM")
source('FmtSleuthToEBHMM.R')
source('adEBHMMNBfun.R')
source('adEBHMMNBfunForMulti.R')
source('adEBSeqHMMTest.R')
source('adEBHMMNBMultiEM_2chain.R')
source('adEBTest_ext.R')

#####################
s2c <- read.table(file.path(".", "metadata", "info4.txt"), header = TRUE, stringsAsFactors=FALSE)
order_list <-s2c[,2]
levels<-order_list[!duplicated(order_list)]
sample_id <- s2c[,1]
#######################
kal_dirs <- file.path(".", "kall", sample_id)
s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
so <- sleuth_prep(s2c)
#####################
levels<-order_list[!duplicated(order_list)]
Conditions <- factor(order_list, levels)
#####################
Sizes <- MedianNorm(data, alternative = TRUE)
data <- GetNormalizedMat(data, Sizes)
EBSeqHMMGeneOut <- EBSeqHMMTest(Data=data, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=10)
adEBSeqHMMGeneOut <- adEBSeqHMMTest(Data=data, sizeFactors=Sizes, Conditions=Conditions, UpdateRd=10,s2c=s2c,filter_fold=0)

