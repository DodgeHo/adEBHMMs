#' @title Extented EBTest function
#' @usage EBTest_ext(Data,NgVector=NULL,Conditions, 
#'	sizeFactors, maxround, Pool=FALSE, NumBin=1000,
#'	ApproxVal=10^-10, Alpha=NULL, Beta=NULL,
#'	PInput=NULL,RInput=NULL,PoolLower=.25, 
#'	PoolUpper=.75,OnlyCalcR=FALSE,Print=TRUE)
#' @param Data Input data, rows are genes/isoforms and columns are samples. Data should come from a two condition experiment
#' @param NgVector Ng vector; NULL for gene level data
#' @param Conditions A factor indicates the condition (time/spatial point) which each sample belongs to. Only two levels are allowed.
#' @param sizeFactors a vector indicates library size factors
#' @param maxround number of iteration
#' @param Pool While working without replicates, user could define the Pool
#'          = TRUE in the EBTest function to enable pooling.
#' @param NumBin By defining NumBin = 1000, EBSeq will group the genes with
#'          similar means together into 1,000 bins.	
#' @param PoolLower,PoolUpper With the assumption that only subset of the genes
#'          are DE in the data set, we take genes whose FC are in the
#'          PoolLower - PoolUpper quantile of the FCs as the candidate
#'					          genes (default is 25%-75%).
#' For each bin, the bin-wise variance estimation is defined as
#'          the median of the cross condition variance estimations of the
#'          candidate genes within that bin.
#'					  We use the cross condition variance estimations for the
#'					          candidate genes and the bin-wise variance estimations of the
#'										          host bin for the non-candidate genes.
#' @param ApproxVal The variances of the transcripts with mean < var will be
#'										          approximated as mean/(1-ApproxVal).
#' @param Alpha,Beta,PInput,RInput If the parameters are known and the user
#'          doesn't want to estimate them from the data, user may
#'					          specify them here.
#' @param Print Whether print the elapsed-time while running the test.					
#' @param OnlyCalcR if OnlyCalcR=TRUE, the function will only return estimation of r's.
#' @author Ning Leng
#' @examples data(GeneExampleData)
#' Data=GeneExampleData[,1:6]
#' CondVector <- rep(paste("t",1:2,sep=""),each=3)
#' Conditions <- factor(CondVector, levels=c("t1","t2"))
#' Sizes <- MedianNorm(Data[1:10,])
#' Out <- EBTest_ext(Data=Data[1:10,], sizeFactors=Sizes, Conditions=Conditions,
#'          maxround=1)
#' @details    
#' EBSeq_ext() function is an extension of EBTest() function, which is used to calculate the conditional probability P(X_g,t | X_g,t-1).
#' In EBSeqHMM, we assume the conditional distribution is Beta-Negative Binomial.
#' @return See \code{\link{EBTest}}

adEBTest_ext <-
function(Data,NgVector=NULL,Conditions, sizeFactors, maxround,SigmaFoldList=NULL,SigInput=FALSE,
Pool=FALSE, NumBin=1000,ApproxVal=10^-10, Alpha=NULL, 
Beta=NULL,PInput=NULL,RInput=NULL,PoolLower=.25, PoolUpper=.75,
OnlyCalcR=FALSE,Print=TRUE){
	
	
	 if (!is.factor(Conditions)) 
		         Conditions = as.factor(Conditions)
		if (!is.matrix(Data)) stop("The input Data is not a matrix")
    if (length(Conditions) != ncol(Data)) 
		         stop("The number of conditions is not the same as the number of samples! ")
   if (nlevels(Conditions) > 2) 
		         stop("More than 2 conditions! Please use EBMultiTest() function")
   if (nlevels(Conditions) < 2) 
		        stop("Less than 2 conditions - Please check your input")
    if (length(sizeFactors) != length(Data) & length(sizeFactors) != 
						         ncol(Data)) 
			     stop("The number of library size factors is not the same as the number of samples!")
		     
	
	Vect5End <- NULL
	Vect3End <- NULL
	CI <- NULL
	CIthre <- NULL
	tau <-NULL
	Dataraw <- Data
	AllZeroNames <- which(rowMeans(Data)==0)
	NotAllZeroNames <- which(rowMeans(Data)>0)
	if(length(AllZeroNames)>0 & Print==TRUE) cat("Remove transcripts with all zero \n")
	Data <- Data[NotAllZeroNames,]
	SigmaFoldList<-SigmaFoldList[NotAllZeroNames]
	if(!is.null(NgVector))NgVector <- NgVector[NotAllZeroNames]
	if(!length(sizeFactors)==ncol(Data))sizeFactors <- sizeFactors[NotAllZeroNames,]

	if(is.null(NgVector))NgVector <- rep(1,nrow(Data))

	if(SigInput==TRUE){
		SigmaFoldList<-as.numeric(SigmaFoldList)
	    if (length(SigmaFoldList) != nrow(Data)) {
	    	cat(paste0("Length(SigmaFoldList): ",length(SigmaFoldList)," nrow(Data):",length(idx2) ," \n"))
		    stop("BeforeCal:The number of SigmaFoldList is not the same as the number of samples! ")
	    }
	    else{
	    	cat(paste0("Using sleuth provided Sigma_sq:",
						" length(SigmaFoldList,nrow(Data)): ",length(SigmaFoldList)," ",nrow(Data) ," \n"))
	    }
	}	
	#Rename Them
	IsoNamesIn <- rownames(Data)
	Names <- paste("I",c(1:dim(Data)[1]),sep="")
	names(IsoNamesIn) <- Names
	rownames(Data) <- paste("I",c(1:dim(Data)[1]),sep="")
	names(NgVector) <- paste("I",c(1:dim(Data)[1]),sep="")
	

	if(!length(sizeFactors)==ncol(Data)){
		rownames(sizeFactors) <- rownames(Data)
		colnames(sizeFactors) <- Conditions
	}
	
	NumOfNg <- nlevels(as.factor(NgVector))
	NameList <- sapply(1:NumOfNg,function(i)Names[NgVector==i],simplify=FALSE)
	names(NameList) <- paste("Ng",c(1:NumOfNg),sep="")
	NotNone <- NULL
	for (i in 1:NumOfNg) {
		if (length(NameList[[i]])!=0) 
			NotNone <- c(NotNone,names(NameList)[i])
		}
	NameList <- NameList[NotNone]
		
	NoneZeroLength <- length(NameList)
	DataList <- vector("list",NoneZeroLength)
	DataList <- sapply(1:NoneZeroLength , function(i) Data[NameList[[i]],],simplify=FALSE)
	names(DataList) <- names(NameList)
    
	NumEachGroup <- sapply(1:NoneZeroLength , function(i)dim(DataList[[i]])[1])
	# Unlist 
	DataList.unlist <- do.call(rbind, DataList)

	# Divide by SampleSize factor
	
	if(length(sizeFactors)==ncol(Data))
	DataList.unlist.dvd <- t(t( DataList.unlist)/sizeFactors)
	
	if(length(sizeFactors)!=ncol(Data))
	DataList.unlist.dvd <- DataList.unlist/sizeFactors
	
	MeanList <- rowMeans(DataList.unlist.dvd)


	# Get FC and VarPool for pooling - Only works on 2 conditions
	if(ncol(Data)==2)
		stop("Data less than two col! Or run with no replicates in each Conditions!")

	#DataListSP Here also unlist.. Only two lists
	DataListSP <- vector("list",nlevels(Conditions))
	DataListSP.dvd <- vector("list",nlevels(Conditions))
	SizeFSP <- DataListSP
	MeanSP <- DataListSP
	VarSP <- DataListSP
	#NewVarSP <- DataListSP
	GetPSP <- DataListSP
	RSP <- DataListSP
	CISP <- DataListSP
	tauSP <- DataListSP
	NumSampleEachCon <- rep(NULL,nlevels(Conditions))

	for (lv in 1:nlevels(Conditions)){
		DataListSP[[lv]] <-  matrix(DataList.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist)[1])
		rownames(DataListSP[[lv]]) <- rownames(DataList.unlist)
		DataListSP.dvd[[lv]] <-  matrix(DataList.unlist.dvd[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist.dvd)[1])
		NumSampleEachCon[lv] <- ncol(DataListSP[[lv]])

	if(ncol(DataListSP[[lv]])==1 & !is.null(CI))
		stop("Data less than two col! Or run with no replicates in each Conditions!")
	# no matter sizeFactors is a vector or a matrix. Matrix should be columns are the normalization factors
	# may input one for each 
	if(length(sizeFactors)==ncol(Data))SizeFSP[[lv]] <- sizeFactors[Conditions==levels(Conditions)[lv]]
	if(length(sizeFactors)!= ncol(Data))SizeFSP[[lv]] <- sizeFactors[,Conditions==levels(Conditions)[lv]]
	
	
	MeanSP[[lv]] <- rowMeans(DataListSP.dvd[[lv]])
	
	if(length(sizeFactors)==ncol(Data))PrePareVar <- sapply(1:ncol( DataListSP[[lv]]),function(i)( DataListSP[[lv]][,i]- SizeFSP[[lv]][i]*MeanSP[[lv]])^2 /SizeFSP[[lv]][i])
	if(length(sizeFactors)!=ncol(Data))PrePareVar <- sapply(1:ncol( DataListSP[[lv]]),function(i)( DataListSP[[lv]][,i]- SizeFSP[[lv]][,i]*MeanSP[[lv]])^2 /SizeFSP[[lv]][,i])

	if(ncol(DataListSP[[lv]])==1 & !is.null(CI))
		stop("Data less than two col! Or run with no replicates in each Conditions!")
	if(ncol(DataListSP[[lv]])!=1){
		if(SigInput==TRUE & nrow(PrePareVar)==length(SigmaFoldList)){
			VarSP[[lv]] <- rowSums(PrePareVar)/ncol( DataListSP[[lv]])*SigmaFoldList
			cat(paste0("From adEBTest_ext: Calling sleuth provided Sigma_sq Succeeded in adjusting VarSP!","\n"))
		}
		else{
			VarSP[[lv]] <- rowSums(PrePareVar)/ncol( DataListSP[[lv]])
			cat(paste0("Warning: Cannot use sleuth provided Sigma_sq:",
						" length(SigmaFoldList),nrow(PrePareVar): ",length(SigmaFoldList)," ",nrow(PrePareVar) ," \n"))
		}
		names(MeanSP[[lv]]) <- rownames(DataList.unlist)
		names(VarSP[[lv]]) <- rownames(DataList.unlist)
		GetPSP[[lv]] <- MeanSP[[lv]]/VarSP[[lv]]
		RSP[[lv]] <- MeanSP[[lv]]*GetPSP[[lv]]/(1-GetPSP[[lv]])
	}
}	
	VarList <- apply(DataList.unlist.dvd, 1, var)
	if(SigInput==TRUE & length(VarList)==length(SigmaFoldList)){
		VarList <- VarList*SigmaFoldList
		cat(paste0("From adEBTest_ext: Calling sleuth provided Sigma_sq Succeeded in adjusting VarList!","\n"))
	}
			
	if(ncol(Data)==2)
		stop("Data less than two col! Or run with no replicates in each Conditions!")
	if(!ncol(Data)==2){
		CondWithRep <- which(NumSampleEachCon>1)
		VarCondWithRep <- do.call(cbind,VarSP[CondWithRep])
		PoolVar <- rowMeans(VarCondWithRep)
	
	}
	GetP <- MeanList/PoolVar
	
    EmpiricalRList <- MeanList*GetP/(1-GetP) 
	EmpiricalRList[EmpiricalRList==Inf]	 <- max(EmpiricalRList[EmpiricalRList!=Inf])
#####################
	if(ncol(Data)!=2){
	Varcbind <- do.call(cbind,VarSP)
	VarrowMin <- apply(Varcbind,1,min)
	}

	if(ncol(Data)==2)
		stop("Data less than two col! Or run with no replicates in each Conditions!")
	
	# 
	# 
	GoodData <- names(MeanList)[EmpiricalRList>0 &  VarrowMin!=0 & EmpiricalRList!=Inf & !is.na(VarrowMin) & !is.na(EmpiricalRList)]
	NotIn <- names(MeanList)[EmpiricalRList<=0 | VarrowMin==0 | EmpiricalRList==Inf |  is.na(VarrowMin) | is.na(EmpiricalRList)]
	#print(paste("ZeroVar",sum(VarrowMin==0), "InfR", length(which(EmpiricalRList==Inf)), "Poi", length(which(EmpiricalRList<0)), ""))
	EmpiricalRList.NotIn <- EmpiricalRList[NotIn]
	EmpiricalRList.Good <- EmpiricalRList[GoodData]
	EmpiricalRList.Good[EmpiricalRList.Good<1] <- 1+EmpiricalRList.Good[EmpiricalRList.Good<1]
	if(length(sizeFactors)==ncol(Data)){
		EmpiricalRList.Good.mat <-  outer(EmpiricalRList.Good, sizeFactors)	
		EmpiricalRList.mat <-  outer(EmpiricalRList, sizeFactors)
	}
	if(!length(sizeFactors)==ncol(Data)){
	EmpiricalRList.Good.mat <- EmpiricalRList.Good* sizeFactors[GoodData,]
	EmpiricalRList.mat <- EmpiricalRList* sizeFactors
	}

	# Only Use Data has Good q's
	DataList.In <- sapply(1:NoneZeroLength, function(i)DataList[[i]][GoodData[GoodData%in%rownames(DataList[[i]])],],simplify=FALSE)
	DataList.NotIn <- sapply(1:NoneZeroLength, function(i)DataList[[i]][NotIn[NotIn%in%rownames(DataList[[i]])],],simplify=FALSE)
	DataListIn.unlist <- do.call(rbind, DataList.In)
	DataListNotIn.unlist <- do.call(rbind, DataList.NotIn)
	
	DataListSPIn <- vector("list",nlevels(Conditions))
	DataListSPNotIn <- vector("list",nlevels(Conditions))
	EmpiricalRList.Good.mat.SP = EmpiricalRList.mat.SP=vector("list",nlevels(Conditions))
	for (lv in 1:nlevels(Conditions)){
		DataListSPIn[[lv]] <-  matrix(DataListIn.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataListIn.unlist)[1])
		if(length(NotIn)>0){	
			DataListSPNotIn[[lv]] <-  matrix(DataListNotIn.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataListNotIn.unlist)[1])
			rownames(DataListSPNotIn[[lv]]) <- rownames(DataListNotIn.unlist)
	}
		rownames(DataListSPIn[[lv]]) <- rownames(DataListIn.unlist)
		EmpiricalRList.Good.mat.SP[[lv]] <- matrix(EmpiricalRList.Good.mat[,Conditions==levels(Conditions)[lv]],nrow=dim(EmpiricalRList.Good.mat)[1])
		EmpiricalRList.mat.SP[[lv]] <- matrix(EmpiricalRList.mat[,Conditions==levels(Conditions)[lv]],nrow=dim(EmpiricalRList.mat)[1])
	}	

	NumOfEachGroupIn <- sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.In[[i]])[1]))
	NumOfEachGroupNotIn <- sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.NotIn[[i]])[1]))
	

#################
# For output
#################
RealName.EmpiricalRList <- sapply(1:NoneZeroLength,function(i)EmpiricalRList[names(EmpiricalRList)%in%NameList[[i]]], simplify=FALSE)
RealName.MeanList <- sapply(1:NoneZeroLength,function(i)MeanList[names(MeanList)%in%NameList[[i]]], simplify=FALSE)
RealName.C1MeanList <- sapply(1:NoneZeroLength,function(i)MeanSP[[1]][names(MeanSP[[1]])%in%NameList[[i]]], simplify=FALSE)
RealName.C2MeanList <- sapply(1:NoneZeroLength,function(i)MeanSP[[2]][names(MeanSP[[2]])%in%NameList[[i]]], simplify=FALSE)
RealName.C1VarList <- sapply(1:NoneZeroLength,function(i)VarSP[[1]][names(VarSP[[1]])%in%NameList[[i]]], simplify=FALSE)
RealName.C2VarList <- sapply(1:NoneZeroLength,function(i)VarSP[[2]][names(VarSP[[2]])%in%NameList[[i]]], simplify=FALSE)
RealName.DataList <- sapply(1:NoneZeroLength,function(i)DataList[[i]][rownames(DataList[[i]])%in%NameList[[i]],], simplify=FALSE)



RealName.VarList <- sapply(1:NoneZeroLength,function(i)VarList[names(VarList)%in%NameList[[i]]], simplify=FALSE)
RealName.PoolVarList <- sapply(1:NoneZeroLength,function(i)PoolVar[names(PoolVar)%in%NameList[[i]]], simplify=FALSE)


RealName.QList1 <- sapply(1:NoneZeroLength,function(i)GetPSP[[1]][names(GetPSP[[1]])%in%NameList[[i]]], simplify=FALSE)
RealName.QList2 <- sapply(1:NoneZeroLength,function(i)GetPSP[[2]][names(GetPSP[[2]])%in%NameList[[i]]], simplify=FALSE)


for (i in 1:NoneZeroLength){
tmp <- NameList[[i]]
names <- IsoNamesIn[tmp]

RealName.MeanList[[i]] <- RealName.MeanList[[i]][NameList[[i]]]
RealName.VarList[[i]] <- RealName.VarList[[i]][NameList[[i]]]
RealName.QList1[[i]] <- RealName.QList1[[i]][NameList[[i]]]
RealName.QList2[[i]] <- RealName.QList2[[i]][NameList[[i]]]
RealName.EmpiricalRList[[i]] <- RealName.EmpiricalRList[[i]][NameList[[i]]]
RealName.C1MeanList[[i]] <- RealName.C1MeanList[[i]][NameList[[i]]]
RealName.C2MeanList[[i]] <- RealName.C2MeanList[[i]][NameList[[i]]]
RealName.PoolVarList[[i]] <- RealName.PoolVarList[[i]][NameList[[i]]]
RealName.C1VarList[[i]] <- RealName.C1VarList[[i]][NameList[[i]]]
RealName.C2VarList[[i]] <- RealName.C2VarList[[i]][NameList[[i]]]
RealName.DataList[[i]] <- RealName.DataList[[i]][NameList[[i]],]

names(RealName.MeanList[[i]]) <- names
names(RealName.VarList[[i]]) <- names
if(ncol(DataListSP[[1]])!=1){
	names(RealName.QList1[[i]]) <- names
	names(RealName.C1VarList[[i]]) <- names
}
if(ncol(DataListSP[[2]])!=1){
	names(RealName.QList2[[i]]) <- names
	names(RealName.C2VarList[[i]]) <- names
}

names(RealName.EmpiricalRList[[i]]) <- names
names(RealName.C1MeanList[[i]]) <- names
names(RealName.C2MeanList[[i]]) <- names
names(RealName.PoolVarList[[i]]) <- names
rownames(RealName.DataList[[i]]) <- names


}
  



	output <- list(

			 RList=RealName.EmpiricalRList, MeanList=RealName.MeanList, 
			 VarList=RealName.VarList, QList1=RealName.QList1, QList2=RealName.QList2, 
			 C1Mean=RealName.C1MeanList, C2Mean=RealName.C2MeanList,
			 C1EstVar=RealName.C1VarList, C2EstVar=RealName.C2VarList)
	#filename<-paste0(format(Sys.time(),format="%H_%M%S"),"_output_adEB")	
	#saveRDS(output,file=paste0(filename,".Rds"))
	return(output)

}