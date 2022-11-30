# AdEBHMMs

AdEBHMMs is an ordered RNA-seq differential analysis method with error model.

You should follow the Rscript.R to run the code.

Written by Lang HE in May 2018.

And 3 libraries are needed for the code: "sleuth", "EBSeq"and "EBSeqHMM".



# A difference Analysis Method for Ordered RNA-Seq with Error Model: adEBHMM

With improvements in next-generation sequencing technologies and reductions in price, ordered RNA-seq experiments are becoming common. Of primary interest in these experiments is identifying genes that are changing over time or space, for example, and then characterizing the speciﬁc expression changes. We propose an empirical Bayes mixture modeling approach called advanced-EBSeq-HMM-Sleuth(abbreviation: adEBHMMs), which can utilize bootstrapping in conjunction with response error linear modeling to decouple biological variance from inferential variance, then conduct ordered DE analysis. We conduct analysis on simulation dataset and case, to illustrate that: comparing to old-line method EBSeq-HMM, adEBHMMs sustains its performance in identifying ordered DE genes, and has higher sensitivity and accuracy in specifying gene-speciﬁc expression paths.
