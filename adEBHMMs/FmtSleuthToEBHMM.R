FmtSleuthToEBHMM <- function(so,kal_dirs,filter_fold=0.1){
	newTable <-so$obs_raw$est_counts
	ncol<-length(kal_dirs)
	nrow<-length(newTable)/length(kal_dirs)
	newTable<-matrix(newTable,byrow=T,nrow=nrow,ncol=ncol)
	exp=data.frame(newTable,stringsAsFactors=FALSE)
	rownames<-so$obs_raw$target_id 
	rownames<-rownames[!duplicated(rownames)]
	rownames(exp)<-rownames[]
	exp[is.na(exp)]<-1
	exp<-exp[which(rowSums(exp)>=ncol+20),]
	colnames(exp) <- NULL
	data<-data.matrix(exp)
	data[is.na(data)]<-0
	rm(so,newTable,exp)
	gc()
	data
}
FitToGetSleuthAvgSigma<- function(so){
	so <- invisible(sleuth_fit(so, ~condition, 'full'))
	so <- invisible(sleuth_fit(so, ~1, 'reduced'))
	so <- invisible(sleuth_lrt(so, 'reduced', 'full'))
	#models(so)
	sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE)
	final_sigma_list<-data.matrix(cbind(sleuth_table[,1],as.numeric(sleuth_table[, 9]),as.numeric(sleuth_table[, 12])))
	sig=data.frame(final_sigma_list,stringsAsFactors=FALSE)
	colnames(sig)<-c("target_id","var_obs","final_sigma_sq")
	rm(so)
               gc()
	sig
}