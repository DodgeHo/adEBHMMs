#' @title Baum-Welch algorithm for a single hidden markov chain
#' @usage EBHMMNBfun(Data,NgVector=NULL,Conditions, sizeFactors,
#'	PriorFC=1.5,homo=TRUE, maxround=5,
#'	Pi0=NULL, Tran=NULL,NoTrend=FALSE, NumTranStage=3,
#'	FCParam=NULL, AlphaIn=NULL,BetaIn=NULL,
#'	StateNames=c("Up","NC","Down"),
#'	EM=TRUE, UpdateParam=TRUE, Print=TRUE,
#'	OnlyQ=FALSE,WithinCondR=TRUE, 
#'	PenalizeLowMed=TRUE, PenalizeLowMedQt=.2,PenalizeLowMedVal=10)
#' @param Data input data, rows are genes/isoforms and columns are samples
#' @param NgVector Ng vector; NULL for gene level data
#' @param Conditions A factor indicates the condition (time/spatial point) which each sample belongs to.
#' @param sizeFactors a vector indicates library size factors
#' @param Tran initial values for transition matrices
#' @param Pi0 initial values for starting probabilities
#' @param NumTranStage number of states
#' @param PriorFC target FC for gridient change
#' @param StateNames name of the hidden states
#' @param homo whether the chain is assumed to be homogenious
#' @param maxround max number of iteration
#' @param AlphaIn,BetaIn If the parameters are known and the user
#'          doesn't want to estimate them from the data, user may
#'                    specify them here.
#' @param NoTrend if NoTrend=TRUE, initial transition probabilities will be set to be equal
#' @param FCParam not in use 
#' @param EM Whether estimate the prior parameters of the beta distribution by EM
#' @param UpdateParam Whether update starting probabilities and transition probabilities
#' @param OnlyQ If OnlyQ=TRUE, the function will only return estimated q's.
#' @param WithinCondR By defining WithinCondR=TRUE, estimation of r's are obtained within each condition. (WithinCondR=FALSE is not suggested here)
#' @param Print Whether print the elapsed-time while running the test.
#' @param PenalizeLowMed,PenalizeLowMedQt,PenalizeLowMedVal 
#' Transcripts with median quantile < = PenalizeLowMedQt  will be penalized
#' @author Ning Leng
#' @examples 
#' data(GeneExampleData)
#' CondVector <- rep(paste("t",1:5,sep=""),each=3)
#' Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
#' Sizes <- MedianNorm(GeneExampleData)
#' tmp <- EBHMMNBfun(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
#'           maxround=2, OnlyQ=TRUE)
#' @details EBHMMNBfun() function implements the Balm-Welch algorithm that estimates the starting probabilities and transition probabilities
#' of a single hidden Markov model. Here the emission distribution of each gene is assumed to be a 
#' Beta-Negative Binomial distribution with parameters (r_g, alpha, beta)
#' , in which alpha and beta are shared by all the genes and r_g is gene specific.
#' If not specified, r_g, alpha and beta will be estimated using method of moments.
#' For isoform data, we assume isoforms from the same Ig group share the same beta^Ig. alpha is shared by all the isoforms and r_gi is isoform specific.
#' The user also needs to specify an expected FC.
#' @return MAPTerm: the most likely path of each gene/isoform. MAPTermNum: numeric version of MAPTerm.
#'
#' AllTerm: all possible expression paths considered in the model. PP: posterior probability of being each expression path.
#'
#' WhichMax: index of the most likely path. Allf: prior probability of being each path.
#'
#' Pi0Track: estimated starting probabilities of each iteration.
#'
#' TranTrack: estimated transition probabilities of each iteration.
#'
#' AlphaTrack, BetaTrack: estimated alpha and beta(s).
#'
#' LLAll=PostSumForLL.Sum: log likelihood of the model.



adEBHMMNBfun <- function(
			Data,NgVector=NULL,Conditions, sizeFactors,
			PriorFC=1.5,homo=TRUE, maxround=5,
			Pi0=NULL, Tran=NULL,NoTrend=FALSE, NumTranStage=3,
			FCParam=NULL, AlphaIn=NULL,BetaIn=NULL,
			StateNames=c("Up","NC","Down"),
			EM=TRUE, UpdateParam=TRUE, Print=TRUE, OnlyQ=FALSE,WithinCondR=TRUE,
			PenalizeLowMed=TRUE, PenalizeLowMedQt=.2,PenalizeLowMedVal=10,s2c=NULL,filter_fold=0.1
		){

	Dataraw <- Data
	AllZeroNames <- which(rowMeans(Data)==0)
	NotAllZeroNames <- which(rowMeans(Data)>0)
	if(length(AllZeroNames)>0 & Print==TRUE) cat("Remove transcripts with all zero \n")
	Data <- Data[NotAllZeroNames,]
	Dataraw <- Data
	if(!is.null(NgVector))NgVector <- NgVector[NotAllZeroNames]
	if(is.null(NgVector))NgVector <- rep(1,nrow(Data))
	if(!is.factor(Conditions))Conditions <- as.factor(Conditions)


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
		
	NumCond <- nlevels(Conditions)
	CondLevels <- levels(Conditions)

	NoneZeroLength <- nlevels(as.factor(NgVector))
	NameList <- sapply(1:NoneZeroLength,function(i)names(NgVector)[NgVector==i],simplify=FALSE)
	DataList <- sapply(1:NoneZeroLength , function(i) Data[NameList[[i]],],simplify=FALSE)
	names(DataList) <- names(NameList)
    
	NumEachGroup <- sapply(1:NoneZeroLength , function(i)dim(DataList)[i])
	# Unlist 
	DataList.unlist <- do.call(rbind, DataList)

	# Divide by SampleSize factor
	
	if(length(sizeFactors)==ncol(Data))
	DataList.unlist.dvd <- t(t( DataList.unlist)/sizeFactors)
	
	if(length(sizeFactors)!=ncol(Data))
	DataList.unlist.dvd <- DataList.unlist/sizeFactors

	# Calculate R for all conditions
	# median for conditions	
	PairRList <- vector("list",NumCond-1)
	#NamesInR<-rep("name",NumCond-1)
	MedList <- matrix(0,ncol=NumCond,nrow=nrow(DataList.unlist))
	rownames(MedList) <- rownames(DataList.unlist)
	for (lv in 1:(NumCond-1)){
		idx1 <- which(Conditions==levels(Conditions)[lv])
		idx2 <- which(Conditions==levels(Conditions)[lv+1])	
	 

		t1 <-  matrix(DataList.unlist[,idx1],nrow=dim(DataList.unlist)[1])
		t2 <-  matrix(DataList.unlist[,idx2],nrow=dim(DataList.unlist)[1])
		tmpdata <- cbind(t1,t2)
		rownames(tmpdata) <- rownames(DataList.unlist)
		Cond.use <- as.factor(as.vector(Conditions)[c(idx1,idx2)])
		
		#cat(paste0("NOTMulti idx1: ",length(idx1)," idx2:",length(idx2) ," \n"))
		tmp_s2c<-NULL
		SigmaFoldList<-NULL
		SigInput<-FALSE
		if(length(c(idx1,idx2))>=2){
			if(!is.null(s2c) && (length(idx1)!=0)&&(length(idx2)!=0)){
				tmp_s2c<-s2c[c(idx1,idx2),]
				cat("EBFun Calling sleuth_prep: ----------------------------------------------")
				tmpso <- invisible(sleuth_prep(tmp_s2c, extra_bootstrap_summary = TRUE))
				cat("Calling sleuth_fit: |")
				tmpsig<-invisible(FitToGetSleuthAvgSigma(tmpso))
				cat("Calling sleuth done!----------------------------------------------\n")
				GeneNameList<-rownames(Dataraw)
				SigmaFoldList<-rep(1,length(GeneNameList))
				for (w in 1:length(GeneNameList)){
					if(length(which(tmpsig$target_id==GeneNameList[w]))!=0){
						idx<-which(tmpsig$target_id==GeneNameList[w])
						tmpsigsq<-as.numeric(tmpsig[idx,3])/as.numeric(tmpsig[idx,2])
						if(!is.na(tmpsigsq)){
							if (tmpsigsq==0)
								SigmaFoldList[w]<-0.000001
							if(tmpsigsq>0)
								SigmaFoldList[w]<-tmpsigsq
						}
					}
				}
				rm(tmpso,tmpsig)
					gc()
				names(SigmaFoldList)<-rownames(DataList.unlist)
				SigInput<-TRUE
				#filename=paste0(format(Sys.time(),format="%H_%M%S"),"_SigmaFoldList.csv")
				#write.table(SigmaFoldList ,filename,row.names = TRUE) 
 				#write.table(GeneNameList,"GeneNameList.csv",row.names = TRUE) 
			}
		}
		PairRList[[lv]]<- adEBTest_ext(Data=tmpdata,NgVector=NgVector,
					 Conditions=Cond.use,sizeFactors=sizeFactors[c(idx1,idx2)], 
					 maxround=5, OnlyCalcR=TRUE, Print=FALSE,SigmaFoldList=SigmaFoldList,SigInput=SigInput)

		#tempPairRList<- adEBTest_ext(Data=tmpdata,NgVector=NgVector,
		#			 Conditions=Cond.use,sizeFactors=sizeFactors[c(idx1,idx2)], 
		#			 maxround=5, OnlyCalcR=TRUE, Print=FALSE,SigmaFoldList=SigmaFoldList,SigInput=SigInput)
		#filename=paste0("PairRList_",as.character(lv),".setting.Rds")
		#saveRDS(tempPairRList,file=paste0(filename))
		cat(paste0("Calling adEBTest_ext done!----","-----------------------\n"))
		#rm(tempPairRList)
		#gc()

		tt1 <-  matrix(DataList.unlist.dvd[,idx1],nrow=dim(DataList.unlist)[1])
		tt2 <-  matrix(DataList.unlist.dvd[,idx2],nrow=dim(DataList.unlist)[1])
		MedList[,i] <- apply(tt1,1,median)
		MedList[,i+1] <- apply(tt2,1,median)
	}
	#PairRList <- vector("list",NumCond-1)
	#for (lv in 1:(NumCond-1)){
	#	filename=paste0("PairRList_",as.character(lv),".setting.Rds")
	#	PairRList[[lv]]<-readRDS(filename)
	#}
	#if(is.null(PairRList[[1]])){
	#	PairRList<-PairRList[length(PairRList)]

	#}
	# if only need Q
	if(OnlyQ==TRUE){
		cat("Returning Q!-----------------------\n")
		return(PairRList)
	}
	# Only use genes with R in all conditions
			#NamesInR[lv]<-names(unlist(PairRList[[lv]][[1]]))
	#filename=paste0(format(Sys.time(),format="%H_%M%S"),"_PairRList_")
	#saveRDS(PairRList,file=paste0(filename,".Rds"))

	NamesInR <- sapply(PairRList,function(i)names(unlist(i[[1]])),simplify=FALSE)

	#	NamesInR <- sapply(PairRList,function(i)names(unlist(i)),simplify=FALSE)
	save(NamesInR,file="NamesInR.R")
	#write.table(NamesInR ,"NamesInR.csv",row.names = TRUE) 
	GoodData0 <- Reduce(intersect, NamesInR)
	GoodData <- rownames(DataList.unlist)# consider all
	NotIn <- setdiff(rownames(DataList.unlist), GoodData)
	cat(paste0("GoodData0: ",length(GoodData0)," \n"))
	cat(paste0("NoneZeroLength: ",NoneZeroLength," \n"))
	cat(paste0("NumCond-1: ",NumCond-1," \n"))
	# Only Use Data has Good q's
	DataList.In <- sapply(1:NoneZeroLength, function(i)DataList[[i]][GoodData[GoodData%in%rownames(DataList[[i]])],],simplify=FALSE)
	DataList.NotIn <- sapply(1:NoneZeroLength, function(i)DataList[[i]][NotIn[NotIn%in%rownames(DataList[[i]])],],simplify=FALSE)
	DataListIn.unlist <- do.call(rbind, DataList.In)
	DataListNotIn.unlist <- do.call(rbind, DataList.NotIn)
	MedList.In <- MedList[GoodData,]
	#cat(paste0("DataList.In : ",length(DataList.In)," \n"))
	# R for good ones	
	# pairwise R
	RPairsListIn <- vector("list",NumCond-1)
	for(lv in 1:(NumCond-1)){
		#cat(paste0("Iterations round begin ,lv: ",lv," \n"))
		RPairsListIn[[lv]] <- sapply(1:NoneZeroLength, 
									function(i){
										tmpname <- GoodData0[GoodData0%in%rownames(DataList[[i]])]
										tmp1 <- PairRList[[lv]][[1]][[i]][tmpname]
										tmpMean <- PairRList[[lv]][[2]][[i]][tmpname]
										tmpVar <- PairRList[[lv]][[3]][[i]][tmpname]
										NAOnes <- which(is.na(tmp1)|tmp1<0)
										#print(length(NAOnes))
										PNotIn <- 1-10^-6
										R.NotIn <- tmpMean[NAOnes]*PNotIn/(1-PNotIn)
										R.out <- tmp1
										R.out[NAOnes] <- R.NotIn
										R.out					
										
										R.norm <- rep(100, nrow(DataList[[i]]))
										names(R.norm) <- rownames(DataList[[i]])
										R.norm[names(R.out)] <- R.out
										R.norm

									})
		#cat("Calling NAones done!-----------------------\n")
	# not pairwise
       if(WithinCondR==TRUE){
        # If calculate R from t-1 but not jointly
                        RPairsListIn[[lv]] <- sapply(1:NoneZeroLength,
                                                                        function(i){
                                                                                tmpname <- GoodData0[GoodData0%in%rownames(DataList[[i]])]
                                                                                #cat("Calling tmpname done!-----------------------\n")
                                                                                tmpMean <- PairRList[[lv]][[6]][[i]][tmpname]
                                                                                #cat("Calling tmpMean done!-----------------------\n")
                                                                                tmpVar <- PairRList[[lv]][[8]][[i]][tmpname]
                                                                                #cat("Calling tmpVar done!-----------------------\n")
                                                                                tmpQ <- PairRList[[lv]][[4]][[i]][tmpname]
                                                                                #cat("Calling tmpQ done!-----------------------\n")
                                                                                tmp1 <- tmpMean*tmpQ/(1-tmpQ)
                                                                                NAOnes <- which(is.na(tmp1)|tmp1<0)
                                                                                #print(length(NAOnes))
         
	 									PNotIn <- 1-10^-6
						       	R.NotIn <- tmpMean[NAOnes]*PNotIn/(1-PNotIn)
                   	R.out <- tmp1
                    R.out[NAOnes] <- R.NotIn
                    R.out 
										R.norm <- rep(100, nrow(DataList[[i]]))
										names(R.norm) <- rownames(DataList[[i]])
										R.norm[names(R.out)] <- R.out
										R.norm


                                                                        })

        }
        #cat(paste0("Iterations round Done ,lv: ",lv," \n"))
	}
	#cat("Calling Rout done!-----------------------\n")
	RPairsListIn.unlist.raw <- sapply(1:(NumCond-1),function(lv)unlist(RPairsListIn[[lv]]))
	rownames(RPairsListIn.unlist.raw) <- rownames(DataListIn.unlist)
	#cat("Calling Rpairrowname done!-----------------------\n")
	# Penalty on R
	if(PenalizeLowMed==FALSE)RPairsListIn.unlist <- RPairsListIn.unlist.raw

	if(PenalizeLowMed==TRUE){
		RPairsListIn.unlist <- RPairsListIn.unlist.raw
		MedQt <- apply(MedList.In,2,function(jj)quantile(jj,PenalizeLowMedQt))
		for(kk in 1:(NumCond-1)){
		WhichLowMed <- which(MedList.In[,kk]<=MedQt[kk])
		RPairsListIn.unlist[WhichLowMed,kk] <- RPairsListIn.unlist.raw[WhichLowMed,kk]* PenalizeLowMedVal
		}
	}
	#cat("Calling RpairsListin done!-----------------------\n")
	NumOfEachGroupIn <- sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.In[[i]])[1]))
	NumOfEachGroupNotIn <- sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.NotIn[[i]])[1]))


# Remove any transcripts with Input NA  
NumPPDim <- NumTranStage
#NumGene <- nrow(DataList.In[[1]])
NumGene <- nrow(DataListIn.unlist)
NumPoints <- NumCond-1

NamesGene <- rownames(DataListIn.unlist)
NamesTranStage <- paste("St",1:NumTranStage,sep="")
NamesPoints <- paste("Pt",1:NumPoints,sep="")
#GeneNameMap <- NamesGene
#names(GeneNameMap) <- IsoNamesIn
GeneNameMap <- IsoNamesIn


if(is.null(Pi0)){Pi0 <- rep(1/NumTranStage, NumTranStage)}
# transition matrix.. homo
# define as from rowname to colname
if(is.null(Tran)&&homo==TRUE){
	if(NumTranStage==2){
		if(NoTrend==FALSE)Tran <- matrix(c(.7,.3,.3,.7), ncol=NumTranStage,nrow=NumTranStage)
		if(NoTrend==TRUE)Tran <- matrix(c(.5,.5,.5,.5), ncol=NumTranStage,nrow=NumTranStage)
	}
	if(NumTranStage==3){
		if(NoTrend==FALSE)Tran <- matrix(c(.6,.2,.2,.2,.6,.2,.2,.2,.6), 
																 ncol=NumTranStage,nrow=NumTranStage)	
		if(NoTrend==TRUE)Tran <- matrix(c(.5,.5,.5,.5,.5,.5,.5,.5,.5),
									             ncol=NumTranStage,nrow=NumTranStage)
	}
	}
# non-homo, 3d matrix, third dim is T-1 transitions
if(is.null(Tran)&&homo==FALSE){
	if(NumTranStage==2){
		Tran <- array(1/NumTranStage, dim=c(NumTranStage,NumTranStage,NumPoints-1))
			if(NoTrend==FALSE)	for(i in 1:(NumPoints-1))Tran[,,i] <- c(.7,.3,.3,.7)
			if(NoTrend==TRUE)  for(i in 1:(NumPoints-1))Tran[,,i] <- c(.5,.5,.5,.5)
	}
	if(NumTranStage==3){
		Tran <- array(1/NumTranStage, dim=c(NumTranStage,NumTranStage,NumPoints-1))
			if(NoTrend==FALSE)for(i in 1:(NumPoints-1))Tran[,,i] <- c(.6,.2,.2,.2,.6,.2,.2,.2,.6)
			if(NoTrend==FALSE)for(i in 1:(NumPoints-1))Tran[,,i] <- c(.5,.5,.5,.5,.5,.5,.5,.5,.5)
	}

#cat("Calling homo done!-----------------------\n")
}
# Beta parameter	
# Time points x Dim
if(is.null(FCParam)){
	if(NumTranStage==2)FCParam <- cbind(c(1/sqrt(PriorFC), sqrt(PriorFC)),
						 c(sqrt(PriorFC), 1/sqrt(PriorFC)))
	if(NumTranStage==3)FCParam <- cbind(c(1/sqrt(PriorFC),1, sqrt(PriorFC)),
				                  c(sqrt(PriorFC),1, 1/sqrt(PriorFC)))

	if(WithinCondR==TRUE){
	# If within condition R based on t-1
	if(NumTranStage==2)FCParam <- cbind(c(1, 1),
					 c(PriorFC, 1/PriorFC))
	if(NumTranStage==3)FCParam <- cbind(c(1,1, 1),
		                   c(PriorFC,1, 1/PriorFC))

	}



}
	AlphaIn <- 1
	BetaIn <- rep(1,NoneZeroLength)

# Name them
rownames(FCParam) <- NamesTranStage

if(homo==TRUE)rownames(Tran)=colnames(Tran)=NamesTranStage
if(homo==FALSE){
	dimnames(Tran)[[1]]=dimnames(Tran)[[2]]=NamesTranStage
	dimnames(Tran)[[3]] <- paste("Cp",1:(NumPoints-1),sep="")
}

TranTrack <- list(Tran)
FCParamTrack <- list(FCParam)
Pi0Track <- Pi0
AlphaTrack <- AlphaIn
BetaTrack <- BetaIn
#cat("Calling AlphaBetaTrack done!-----------------------\n")
################################
# Iter
################################

for (iter in 1:maxround){
Time.start <- proc.time()[3]
# Forward

# rows are States; cols are genes
B <- sapply(1:NumPoints,function(pp)
 sapply(1:NumTranStage,function(ss){
	idx1 <- which(Conditions==levels(Conditions)[pp])
	idx2 <- which(Conditions==levels(Conditions)[pp+1])	
	#
	t1 <-  matrix(DataListIn.unlist[,idx1],nrow=dim(DataListIn.unlist)[1])
	t2 <-  matrix(DataListIn.unlist[,idx2],nrow=dim(DataListIn.unlist)[1])
	tmpdata <- cbind(t1,t2)
	rownames(tmpdata) <- rownames(DataListIn.unlist)
	Rtmp1 <- outer(RPairsListIn.unlist[,pp],sizeFactors[idx1])			
	Rtmp2 <- outer(RPairsListIn.unlist[,pp],sizeFactors[idx2])			
	Log <- f0(tmpdata,
					        AlphaIn=AlphaIn, BetaIn=BetaIn, 
									cbind(Rtmp1*FCParam[ss,1],
										Rtmp2*FCParam[ss,2]), NumOfGroups=NumOfEachGroupIn, log=FALSE)	
}),simplify=FALSE)	
#cat("Calling f0 done!-----------------------\n")

for(i in 1:NumPoints){
B[[i]] <- t(B[[i]])
rownames(B[[i]]) <- NamesTranStage
colnames(B[[i]]) <- NamesGene
}
# alpha_g,0 = pi_i * b_i (O_g)
# row gene column stages
Alpha0 <- B[[1]]*Pi0
Alpha.tmp <- Alpha0
AlphaMat <- NULL
AlphaMat <- list(Alpha0)
# step
for (ctr in 2:(NumPoints)){
if(homo==TRUE)
Sumpart <- t(Tran)%*%Alpha.tmp
if(homo==FALSE)
Sumpart <- t(Tran[,,ctr-1])%*%Alpha.tmp

Alpha.upd <- Sumpart*B[[ctr]]
Alpha.tmp <- Alpha.upd
AlphaMat <- c(AlphaMat,list(Alpha.upd))
}

# P (O| lambda)= sum_i alpha_T-1 (i)
#P_O_lmbda=colSums(Alpha.upd)

#Back
Beta_T1 <- matrix(1,ncol=NumGene,nrow=NumTranStage)
Beta.tmp <- Beta_T1
BetaMat <- vector("list",NumPoints)
BetaMat[[NumPoints]] <- Beta_T1
if(is.null(Tran)&&homo==TRUE){Tran <- matrix(1/NumTranStage, ncol=NumTranStage,nrow=NumTranStage)}
# non-homo, 3d matrix, third dim is T-1 transitions

for(ctr in (NumPoints-1):1){
if(homo==TRUE)tranuse <- Tran
if(homo==FALSE)tranuse <- Tran[,,ctr]
tmpupd <- B[[ctr+1]]*Beta.tmp
Beta.upd <- tranuse%*%tmpupd
Beta.tmp <- Beta.upd
BetaMat[[ctr]] <- Beta.tmp
}

# T-1 elements
# col genes row states
Gamma_i_num <- sapply(1:(NumPoints-1),function(i)AlphaMat[[i]]*BetaMat[[i]],simplify=FALSE)
Gamma_i_denom <- sapply(1:(NumPoints-1),function(i)colSums(Gamma_i_num[[i]]))
Gamma_i <- sapply(1:(NumPoints-1),function(i)t(t(Gamma_i_num[[i]])/Gamma_i_denom[,i]),simplify=FALSE)
# Update Pi and Tran
Pi.Up.allgene <- Gamma_i[[1]]
Pi.Up <- rowMeans(Pi.Up.allgene,na.rm=TRUE)

# Sum g first
#Gamma_i_sumg <- sapply(1:(NumPoints),function(i)rowSums(Gamma_i_num[[i]],na.rm=TRUE))
#Pi.Up <- Gamma_i_sumg[,1]/sum(Gamma_i_sumg[,1],na.rm=TRUE)


# dims: i, j, T-1, g
Gamma_ij_array_raw = Eta_ij_array_raw=
array(data=NA,dim=c(NumTranStage, NumTranStage, NumPoints-1, NumGene))

#cat("Calling gamma done!-----------------------\n")
for(i in 1:NumTranStage){
for(j in 1:NumTranStage){
	for(k in 1:(NumPoints-1)){
		for(l in 1:NumGene){
								if(homo==TRUE)tranuse <- Tran
								if(homo==FALSE)tranuse <- Tran[,,k]
			Gamma_ij_array_raw[i,j,k,l] <- tranuse[i,j]*AlphaMat[[k]][i,l]*BetaMat[[k+1]][j,l]*B[[k+1]][j,l]
			
	}}
}
}
# eta first... then Sum g first!
#cat("Calling gamma array done!-----------------------\n")
for(i in 1:NumTranStage){
for(j in 1:NumTranStage){
	for(k in 1:(NumPoints-1)){
		for(l in 1:NumGene){
			Eta_ij_array_raw[i,j,k,l] <- Gamma_ij_array_raw[i,j,k,l]/sum(Gamma_ij_array_raw[,,k,l],na.rm=TRUE)
			
	}}
}
}
#cat("Calling eta raw array done!-----------------------\n")
Eta_ij_array_raw_sumg <- apply(Eta_ij_array_raw,c(1,2,3),function(oo)sum(oo,na.rm=TRUE))
Eta_ij_array_raw_sumgj <- sapply(1:(NumPoints-1),function(i)rowSums(Eta_ij_array_raw_sumg[,,i],
																																	 na.rm=TRUE),simplify=FALSE)
if(homo==FALSE){
	Tran.Up.raw <- sapply(1:(NumPoints-1),function(i)
										 Eta_ij_array_raw_sumg[,,i]/Eta_ij_array_raw_sumgj[[i]],
										 simplify=FALSE)
	Tran.Up <- array(data=NA,dim=c(NumTranStage, NumTranStage, NumPoints-1))
	for(i in 1:(NumPoints-1))Tran.Up[,,i] <- Tran.Up.raw[[i]]
}

if(homo==TRUE){
	Eta_ij_array_raw_sumgt <- apply(Eta_ij_array_raw_sumg,c(1,2),function(oo)sum(oo,na.rm=TRUE))
	Eta_ij_array_raw_sumgjt_raw <- do.call(cbind,Eta_ij_array_raw_sumgj)
	Eta_ij_array_raw_sumgjt <- rowSums(Eta_ij_array_raw_sumgjt_raw, na.rm=TRUE)
	Tran.Up <- Eta_ij_array_raw_sumgt/Eta_ij_array_raw_sumgjt

}

#cat("Calling all eta done!-----------------------\n")
##############################
#Optim
if (EM==TRUE){

StartValue <- c(AlphaIn, BetaIn)

AlphaToEMraw <- do.call(cbind, AlphaMat)* do.call(cbind, BetaMat)
AlphaToEM <- as.vector(t(AlphaToEMraw))
#print(summary(AlphaToEM))
AlphaToEM[which(is.na(AlphaToEM))] <- 10^-6
AlphaToEM[which(AlphaToEM==Inf)] <- 10^6


InputPool <- list(DataListIn.unlist, AlphaToEM, NumPoints, 
							 NumGene, NumTranStage, RPairsListIn.unlist, sizeFactors, 
							 Conditions, NumOfEachGroupIn, FCParam)
###################	
Result<-optim(par=StartValue,fn=LikefunNBHMM,InputPool=InputPool)

out1 <- Result[[1]]
AlphaIn <- out1[1]
BetaIn <- out1[2:(1+length(NumOfEachGroupIn))]
AlphaTrack <- c(AlphaTrack,AlphaIn)
BetaTrack <- rbind(BetaTrack,BetaIn)

}
#cat("Calling em done!-----------------------\n")
if(UpdateParam==TRUE){
	Pi0Track <- rbind(Pi0Track,Pi.Up)
	TranTrack <- c(TranTrack,list(Tran.Up))

	Pi0 <- Pi.Up
	Tran <- Tran.Up
}
Time.end <- proc.time()[3]
Time.use <- Time.end-Time.start
cat(paste("-------------------\n iteration","time",round(Time.use,1), "\n-------------------\n"))
}
#cat("Calling updateparam done!-----------------------\n")
############# 
# R for all
# NA take mean of others
#############
RPairsListAll <- matrix(NA,ncol=(NumCond-1),nrow=nrow(DataList.unlist))
rownames(RPairsListAll) <- rownames(DataList.unlist)
for(lv in 1:(NumCond-1)){
	tmpName <- names(unlist(PairRList[[lv]][[1]]))
	tmp1 <- unlist(PairRList[[lv]][[1]])
	tmpMean <- unlist(PairRList[[lv]][[2]])
	NAOnes <- which(is.na(tmp1)|tmp1<0)
	PNotIn <- 1-10^-6
	R.NotIn <- tmpMean[NAOnes]*PNotIn/(1-PNotIn)
	R.out <- tmp1
	R.out[NAOnes] <- R.NotIn
	RPairsListAll[tmpName,lv] <- R.out
}
#cat("Calling iteralmean in numcond done!-----------------------\n")
RMean <- rowMeans(RPairsListAll,na.rm=TRUE)
WhichNA <- which(is.na(RPairsListAll),arr.ind=TRUE)
if(nrow(WhichNA)>0)for(i in 1:nrow(WhichNA))RPairsListAll[WhichNA[i]] <- RMean[WhichNA[i,1]]


	# Penalty on R
	if(PenalizeLowMed==FALSE)RPairsListAllPen <- RPairsListAll

	if(PenalizeLowMed==TRUE){
		RPairsListAllPen <- RPairsListAll
		MedQt <- apply(MedList,2,function(jj)quantile(jj,PenalizeLowMedQt))
		for(kk in 1:(NumCond-1)){
		WhichLowMed <- which(MedList[,kk]<=MedQt[kk])
		RPairsListAllPen[WhichLowMed,kk] <- RPairsListAll[WhichLowMed,kk]* PenalizeLowMedVal
		}
	}
#cat("Calling iteralmean NA in numcond done!-----------------------\n")
AllStList <- vector("list",NumPoints-1)
for(i in 1: (NumPoints))AllStList[[i]] <- StateNames
AllPosi <- expand.grid(AllStList)


	LogB <- sapply(1:NumPoints,function(pp)
			 sapply(1:NumTranStage,function(ss){
	idx1 <- which(Conditions==levels(Conditions)[pp])
	idx2 <- which(Conditions==levels(Conditions)[pp+1])	

	t1 <-  matrix(DataList.unlist[,idx1],nrow=dim(DataList.unlist)[1])
	t2 <-  matrix(DataList.unlist[,idx2],nrow=dim(DataList.unlist)[1])
	tmpdata <- cbind(t1,t2)
	rownames(tmpdata) <- rownames(DataList.unlist)
	Rtmp1 <- outer(RPairsListAllPen[,pp],sizeFactors[idx1])			
	Rtmp2 <- outer(RPairsListAllPen[,pp],sizeFactors[idx2])			
	Log <- f0(tmpdata,
					        AlphaIn=AlphaIn, BetaIn=BetaIn, 
									cbind(Rtmp1*FCParam[ss,1],
										Rtmp2*FCParam[ss,2]), NumOfGroups=NumOfEachGroupIn+NumOfEachGroupNotIn, log=TRUE)
	
}),simplify=FALSE)	
#cat("Calling ALLpisi done!-----------------------\n")
for(i in 1:NumPoints){
LogB[[i]] <- t(LogB[[i]])
rownames(LogB[[i]]) <- NamesTranStage
colnames(LogB[[i]]) <- rownames(DataList.unlist)
}

AlllogPost <- sapply(1:nrow(AllPosi),function(i){
				 logPi0.use <- log(Pi0[AllPosi[i,1]])
				logB.use <- sapply(1:NumPoints,function(j)LogB[[j]][AllPosi[i,j],])
				# gene by time
				logB.sumtime <- rowSums(logB.use)
				if(homo==TRUE)tran.use <- sapply(1:(NumPoints-1),function(j)Tran[as.numeric(AllPosi[i,j]),as.numeric(AllPosi[i,j+1])])
				if(homo==FALSE)tran.use <- sapply(1:(NumPoints-1),function(j)Tran[as.numeric(AllPosi[i,j]),as.numeric(AllPosi[i,j+1]),j])
				logtran.sum=sum(log(tran.use))
				logPi0.use+logB.sumtime+logtran.sum 
				 })
#cat("Calling ALLOGPOST done!-----------------------\n")
AlllogPost.Max <- apply(AlllogPost,1,max)
PostSum <- apply(exp(AlllogPost),1,sum)
Which0 <- which(PostSum==0)
AlllogPost.Adj0 <- AlllogPost
if(length(Which0)>0)AlllogPost.Adj0[Which0,] <- AlllogPost[Which0,]-AlllogPost.Max[Which0]

PostSum.Adj0 <- rowSums(exp(AlllogPost.Adj0))
PPraw <- exp(AlllogPost.Adj0)/PostSum.Adj0
PP <- PPraw

LogPostSum.Adj0 <- log(PostSum.Adj0)
LogPostSum <- LogPostSum.Adj0
if(length(Which0)>0)LogPostSum[Which0] <- LogPostSum.Adj0[Which0]+AlllogPost.Max[Which0]
#PP[which(is.na(PP),arr.ind=TRUE)]=-Inf
WhichMax <- apply(PP,1,function(oo){
							 tt <- which.max(oo)
							 if(length(tt)>0)out <- tt
								if(length(tt)==0)out <- NA
				 out})
WhichMax.NonNA.ind <- which(!is.na(WhichMax))
WhichMax.NonNA <- WhichMax[WhichMax.NonNA.ind]
MatTerm <- AllPosi[WhichMax.NonNA,]
rownames(MatTerm) <- rownames(Data)[WhichMax.NonNA.ind]	
#cat("Calling wich max done!-----------------------\n")
AllStListNum <- vector("list",NumPoints-1)
for(i in 1: (NumPoints))AllStListNum[[i]] <- 1:2
AllPosiNum <- expand.grid(AllStListNum)
MatTermNum <- AllPosiNum[WhichMax.NonNA,]
rownames(MatTermNum) <- rownames(Data)[WhichMax.NonNA.ind]

PostSumForLL <- log(PostSum)
ISNA <- which(PostSumForLL==-Inf | is.na(PostSumForLL))
#print(ISNA)
if(length(ISNA)>0)ISNOTNA <- c(1:length(PostSumForLL))[-ISNA]
if(length(ISNA)==0)ISNOTNA <- c(1:length(PostSumForLL))
Min.ISNOT <- min(PostSumForLL[ISNOTNA])
PostSumForLL[ISNA] <- Min.ISNOT
PostSumForLL.Sum <- sum(PostSumForLL)
#cat("Calling postsum done!-----------------------\n")
NameOut <- as.character(IsoNamesIn[rownames(PP)])
#message(str(NameOut))
rownames(MatTerm) = rownames(MatTermNum)=
NameOut[WhichMax.NonNA.ind]
rownames(PP)=
names(WhichMax)=names(LogPostSum)=
names(PostSum)=rownames(AlllogPost)=
NameOut
#cat("Calling names done!-----------------------\n")
tmp <- AllPosi
tmpL <- levels(tmp[[1]])
qq <- sapply(1:nrow(tmp),function(i)paste0(tmpL[as.numeric(tmp[i,])],collapse="-"))
colnames(PP) <- qq
cat("adEBHMMNBfun ready to out-----------------------\n")
rm(PairRList)
gc()
out <- list(MAPTerm=MatTerm, MAPTermNum=MatTermNum,
				 AllTerm=AllPosi, PP=PP, WhichMax=WhichMax, Allf=AlllogPost,
					 Pi0Track=Pi0Track,
					 TranTrack=TranTrack,
					 AlphaTrack=AlphaTrack,
					 BetaTrack=BetaTrack,
					 GeneNameMap=GeneNameMap,
					 LogLike=LogPostSum, 
					 LL=PostSum, LLAll=PostSumForLL.Sum
						 )
}
