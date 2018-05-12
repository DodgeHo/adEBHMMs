#' @title Run EBSeqHMM model with a fixed expected FC
#' @usage EBHMMNBMultiEM_2chain(Data, 
#'	NgVector=NULL, Conditions, AllTran=NULL, 
#'	AllPi0=NULL, Terms=NULL,  
#'	sizeFactors, NumTranStage=c(3,2),PriorFC=2, 
#'	StateNames=c("Up","Down"),homo=FALSE, 
#'	UpdateRd=5, PIBound=.9, UpdatePI=FALSE,Print=FALSE,
#'	WithinCondR=TRUE,
#'	PenalizeLowMed=TRUE, PenalizeLowMedQt=.1,PenalizeLowMedVal=10)
#' @param Data input data, rows are genes and columns are samples
#' @param NgVector Ng vector; NULL for gene level data
#' @param Conditions A factor indicates the condition (time/spatial point) which each sample belongs to.
#' @param AllTran initial values for transition matrices
#' @param AllPi0 initial values for starting probabilities
#' @param Terms Terms
#' @param sizeFactors a vector indicates library size factors
#' @param StateNames names of the hidden states
#' @param NumTranStage number of states in two chains
#' @param PriorFC target FC for gridient change
#' @param homo whether the chain is assumed to be homogenious
#' @param UpdateRd max number of iteration
#' @param UpdatePI whether update the mixture proportion of two chains
#' @param PIBound upper bound of the mixture proportion of the two chains
#' @param Print Whether print the elapsed-time while running the test.
#' @param WithinCondR By defining WithinCondR=TRUE, estimation of r's are obtained within each condition. (WithinCondR=FALSE is not suggested here)
#' @param PenalizeLowMed,PenalizeLowMedQt,PenalizeLowMedVal 
#' Transcripts with median quantile < = PenalizeLowMedQt  will be penalized
#' @author Ning Leng
#' @examples 
#' data(GeneExampleData)
#' CondVector <- rep(paste("t",1:5,sep=""),each=3)
#' Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
#' Sizes <- MedianNorm(GeneExampleData)
#' tmp <- EBHMMNBMultiEM_2chain(Data=GeneExampleData,sizeFactors=Sizes, Conditions=Conditions,
#'           UpdateRd=2)
#' @details EBHMMNBMultiEM_2chain() function implements the EBSeqHMM model to perform statistical analysis in an RNA-seq experiment with ordered conditions.
#' EBHMMNBMultiEM_2chain() calls EBHMMNBfunForMulti() function to perform Balm-Welch algorithm that estimates the starting probabilities and transition probabilities.
#' Here the emission distribution of each gene is assumed to be a 
#' Beta-Negative Binomial distribution with parameters (r_g, alpha, beta)
#' , in which alpha and beta are shared by all the genes and r_g is gene specific.
#' If not specified, r_g, alpha and beta will be estimated using method of moments.
#' For isoform data, we assume isoforms from the same Ig group share the same beta^Ig. alpha is shared by all the isoforms and r_gi is isoform specific.
#' The user also needs to specify an expected FC.
#' Function EBSeqHMMTest() runs several models with varying FCs and returns the model with maximum likelihood.
#' @return 
#' Pi0Out: estimated starting probabilities of each iteration.
#'
#' TranOut: estimated transition probabilities of each iteration.
#'
#' Pi: estimated probability of being each chain.
#'
#' Alpha, Beta: estimated alpha and beta(s).
#'
#' LLSum: log likelihood of the model.
#'
#' QList: estimated q's.
#'
#' MgAllPP: marginal PP for all paths.
#'
#' MgAllMAPChar: Most likely path based on MgAllPP.
#'
#' MgAllMaxVal: Highest PP based on MgAllPP.
#'
#' PPMatW: Posterior probabilities of being each of the chains.
#'

adEBHMMNBMultiEM_2chain <- function(Data, 
												NgVector=NULL, Conditions,
												AllTran=NULL, 
												AllPi0=NULL, Terms=NULL,  
												sizeFactors, NumTranStage=c(3,2),PriorFC=2, 
												StateNames=c("Up","Down"),homo=FALSE, 
												UpdateRd=5, PIBound=.9, UpdatePI=FALSE,Print=FALSE
,WithinCondR=TRUE,
PenalizeLowMed=TRUE, PenalizeLowMedQt=.1,PenalizeLowMedVal=10,s2c=NULL,filter_fold=0.01){
		
##########################
# Define transition parameters if they are not defined in inputs
# First chain is noise chain - 3 states Up EE Down
# Second chain is signal chain - 2 states Up Down
##########################
	NumCond <- nlevels(as.factor(Conditions))
	if(is.null(AllTran)){
	AllTran <- vector("list",2)
	AllTran[[1]] <- array(1/3,dim=c(3,3,NumCond-2))
	AllTran[[2]] <- array(1/2,dim=c(2,2,NumCond-2))
	}
	if(is.null(AllPi0)){
	AllPi0 <- vector("list",2)
	AllPi0[[1]] <- rep(1/3,3)
	AllPi0[[2]] <- rep(1/2,2)
	}

##########################
# Initialization and removing non expressed genes
##########################
		PI <- rep(1/(length(AllTran)), length(AllTran))
		PI.up <- PI
		PI.track <- PI
		Tran.Up <- AllTran
		Pi0.Up <- AllPi0
		DataRename <- Data
		Dataraw<-Data
		NamesIn <- rownames(Data)
		#NewNames <- paste("g",1:nrow(Data),sep="")
		NewNames <- rownames(Data)
		#rownames(DataRename)<-rownames(Data)
		#rownames(DataRename) <- NewNames
		Data.NZ <- DataRename[which(rowSums(Data)>0),]
		Data.NZ.OriName <- Data[which(rowSums(Data)>0),]
		if(!is.null(NgVector))names(NgVector) <- NewNames
		NotNA <- c(1:nrow(Data.NZ))
		
		StateNameList <- vector("list",2)
		StateNameList[[1]] <- c("Up","EE","Down")
		StateNameList[[2]] <- c("Up","Down")
		UpParam <- c(FALSE,TRUE)
		if(UpdateRd==0){
		UpParam <- c(FALSE, FALSE)
		UpdateRd <- 1
		}


######################
# MOM alpha beta
######################
progress_rate<-0
cat(paste0("Rate of Calculation progress: ",progress_rate*100,"% \n"))

QList <- suppressWarnings(EBHMMNBfun(DataRename,NgVector=NgVector,Conditions=Conditions, sizeFactors=sizeFactors,
					  PriorFC=PriorFC,homo=homo, maxround=1,
					 Pi0=AllPi0[[1]], Tran=AllTran[[1]], NumTranStage=NumTranStage[1],
					 FCParam=NULL, AlphaIn=NULL,BetaIn=NULL,
               				 StateNames=StateNameList[[i]],
	               			EM=FALSE, UpdateParam=FALSE, Print=Print, OnlyQ=TRUE,WithinCondR=WithinCondR,
					PenalizeLowMed=PenalizeLowMed,
					                        PenalizeLowMedQt=PenalizeLowMedQt,
								                        PenalizeLowMedVal=PenalizeLowMedVal))
saveRDS(QList,"oldQList.rds")  
QList <- suppressWarnings(adEBHMMNBfun(DataRename,NgVector=NgVector,Conditions=Conditions, sizeFactors=sizeFactors,
					  PriorFC=PriorFC,homo=homo, maxround=1,
					 Pi0=AllPi0[[1]], Tran=AllTran[[1]], NumTranStage=NumTranStage[1],
					 FCParam=NULL, AlphaIn=NULL,BetaIn=NULL,
               				 StateNames=StateNameList[[i]],
	               			EM=FALSE, UpdateParam=FALSE, Print=Print, OnlyQ=TRUE,WithinCondR=WithinCondR,
					PenalizeLowMed=PenalizeLowMed,
					                        PenalizeLowMedQt=PenalizeLowMedQt,
								                        PenalizeLowMedVal=PenalizeLowMedVal,s2c=s2c,filter_fold=filter_fold))

saveRDS(QList,"newQList.rds")  
progress_rate<-0.1
cat(paste0("Rate of Calculation progress: ",progress_rate*100,"% \n"))

QUse <- vector("list",length(QList)+1)
for(i in 1:length(QList))
	QUse[[i]] <- QList[[i]]$QList1
QUse[[length(QUse)]] <- QList[[length(QList)]]$QList2
#write.csv(QUse,"Quse.csv",row.names = TRUE) 
# list length N, N=#conditions

QEach <- vector("list",length(QUse[[1]]))

for(i in 1:length(QUse[[1]]))
	QEach[[i]] <- as.vector(unlist(sapply(QUse,function(j)j[[i]])))
#write.csv(QEach,"QEach.csv",row.names = TRUE) 
cat(paste0("Done cal QEach. "," \n"))


QInRange <- sapply(QEach,function(i)i[which(i>0 & i<1)],simplify=FALSE)

#Est <- sapply(QInRange,beta.mom,simplify=FALSE)
Est <- sapply(QInRange,
					function(qs.in){
						xbar <- mean(qs.in)
						s2 <- var(qs.in)
						# if(!is.null(s2c)){
						# 	so <- sleuth_prep(s2c)
						# 	data<-FmtSleuthToEBHMM(so,kal_dirs,filter_fold=filter_fold)
						# 	sig<-FitToGetSleuthAvgSigma(so)
						# 	sigma_fold_list<-as.numeric(sig[,3])/as.numeric(sig[,2])
						# 	sigma_fold_list[which(sigma_fold_list>1)]<-1
						# 	sigma_fold<-mean(sigma_fold_list,na.rm=T)
						# 	if(!is.null(sigma_fold))
						# 		cat(paste0("Total sigma fold= ",sigma_fold ," \n"))
						# 		if(sigma_fold<=1&sigma_fold>0){
						# 			s2 <- s2*sigma_fold
						# 			cat(paste0("Calculate Alpha Beta by sleuth. "," \n"))
						# 		}
						# 	rm(so,sig)
   			# 				gc()
						# }
						term <- (xbar*(1-xbar))/s2
						alpha.hat <- xbar*(term-1)
						beta.hat <- (1-xbar)*(term-1)
						return(c(alpha.hat,beta.hat))
					},
				simplify=FALSE)


AlphaMom <- mean(unlist(sapply(Est,function(i)i[1])))
BetaMom <- unlist(sapply(Est,function(i)i[2]))

# Q for check
NumGp <- length(QUse[[1]])
QOut <- vector("list",NumGp)
for(i in 1:NumGp){
	QOut[[i]] <- vector("list",length(QUse))
	for(j in 1:length(QUse)) QOut[[i]][[j]] <- QUse[[j]][[i]]
}

##########################
# Iterations
##########################



		for(rd in 1:UpdateRd){
		
			
			EachEB = EachLike=vector("list",length(AllTran))
			for(i in 1:length(AllTran)){
				# First run -  P(Noise group) = P(Signal group)
				if(rd==1){
					EachEB[[i]] <- suppressWarnings(adEBHMMNBfun(DataRename,NgVector=NgVector,Conditions=Conditions, sizeFactors=sizeFactors,
									   PriorFC=PriorFC,homo=homo, maxround=1,
										 Pi0=AllPi0[[i]], Tran=AllTran[[i]], NumTranStage=NumTranStage[i],
										 FCParam=NULL, AlphaIn=AlphaMom,BetaIn=BetaMom,
		               							 StateNames=StateNameList[[i]],
		               							EM=FALSE, UpdateParam=FALSE, Print=Print,
										WithinCondR=WithinCondR,
										PenalizeLowMed=PenalizeLowMed,
										                        PenalizeLowMedQt=PenalizeLowMedQt,
													                        PenalizeLowMedVal=PenalizeLowMedVal,s2c=s2c,filter_fold=filter_fold))
				}
				# other runs - P(N group) and P(S group) estimated from last step
				if(rd>1){
					# Isoform level  
					if(!is.null(NgVector))
					EachEB[[i]] <- suppressWarnings(adEBHMMNBfunForMulti(DataRename[NotNA,],PPtmp0[,i],
										NgVector=NgVector[NotNA],Conditions=Conditions, sizeFactors=sizeFactors,
									   PriorFC=PriorFC,homo=homo, maxround=1,
										 Pi0=Pi0.Up[[i]], Tran=Tran.Up[[i]], NumTranStage=NumTranStage[i],
										 FCParam=NULL, AlphaIn=AlphaMom,BetaIn=BetaMom,
		                 						StateNames=StateNameList[[i]],
			              						 EM=FALSE, UpdateParam=UpParam[i],Print=Print,
										 WithinCondR=WithinCondR,PenalizeLowMed=PenalizeLowMed,
										                         PenalizeLowMedQt=PenalizeLowMedQt,
													                         PenalizeLowMedVal=PenalizeLowMedVal,s2c=s2c,filter_fold=filter_fold))
					# gene level
					if(is.null(NgVector))
					EachEB[[i]] <- suppressWarnings(adEBHMMNBfunForMulti(DataRename[NotNA,],PPtmp0[,i],
										NgVector=NULL,Conditions=Conditions, sizeFactors=sizeFactors,
									   PriorFC=PriorFC,homo=homo, maxround=1,
										 Pi0=Pi0.Up[[i]], Tran=Tran.Up[[i]], NumTranStage=NumTranStage[i],
										 FCParam=NULL, AlphaIn=AlphaMom,BetaIn=BetaMom,
		        						         StateNames=StateNameList[[i]],
			              						 EM=FALSE, UpdateParam=UpParam[i],Print=Print,
										 WithinCondR=WithinCondR,
										 PenalizeLowMed=PenalizeLowMed,
										                         PenalizeLowMedQt=PenalizeLowMedQt,
													                         PenalizeLowMedVal=PenalizeLowMedVal,s2c=s2c,filter_fold=filter_fold))
				}
			}
			# Update parameters
			if(rd>1){
				Tran.Up <- sapply(1:length(EachEB),function(i)EachEB[[i]]$TranTrack[[length(EachEB[[i]]$TranTrack)]], simplify=FALSE)
				if(Print==TRUE)print(Tran.Up)
				Pi0.Up <- sapply(1:length(EachEB),function(i)EachEB[[i]]$Pi0Track[length(EachEB[[i]]$TranTrack),], simplify=FALSE)
			}

			# Calculate likelihood 
			AllTerms <- EachEB[[1]][[3]]
			# Allf P(path)*P(data|path)
			EachLike64 <- sapply(1:length(AllTran),function(i)EachEB[[i]]$Allf,
											                  simplify=FALSE)
		
			EachSum64Mat <- matrix(0, ncol=length(AllTran), nrow=nrow(EachEB[[i]]$Allf))
			rownames(EachSum64Mat) <- 
			rownames(EachEB[[i]]$Allf)
			NotNA <- rownames(EachEB[[i]]$Allf)
			for(i in 1:(length(AllTran)))EachSum64Mat[rownames(EachLike64[[i]]),i] <- rowSums(exp(EachLike64[[i]]))
				LogSum64Mat <- log(EachSum64Mat)
			#	PI=rep(1/(length(AllTran)), length(AllTran))
			#	PI.up=PI
			#	PI.track=PI
			#	for(j in 1:maxround){
				LogPI <- log(PI.up)
				LogPILike <- t(t(LogSum64Mat)+LogPI)
				PILike <- exp(LogPILike)
				SumPILike <- rowSums(PILike)
				PPtmp <- PILike/SumPILike
				PI.up <- colMeans(PPtmp,na.rm=TRUE)
				if(!is.null(PIBound)){
					if(max(PI.up)>PIBound){
						PI.up[which.max(PI.up)] <- PIBound
						PI.up[which.min(PI.up)] <- 1-PIBound
					}
				}
				if(UpdatePI==FALSE)PI.up <- PI
					PI.track <- rbind(PI.track, PI.up)
			#	}

			# Calc PP based on PI
			LogPI <- log(PI.up)
			LogPILike <- t(t(LogSum64Mat)+LogPI)
			# likelihood
			PILike <- exp(LogPILike)
			SumPILike <- rowSums(PILike)
			#################
			# P (W=1 | Y)
	    	##################
			# PP
			PPtmp0 <- PILike/SumPILike
			# Max term
			PPMaxT0 <- sapply(1:nrow(PPtmp0),function(i){
										tt=which.max(PPtmp0[i,])
										if(length(tt)==0)tt=NA
										tt
										})
			PPMaxVal0 <- sapply(1:nrow(PPtmp0),function(i)max(PPtmp0[i,]))	
			progress_rate<-(0.1+(rd*0.9/UpdateRd))
			cat(paste0("Rate of Calculation progress: ",progress_rate*100,"% \n"))
		}
	
		#################
		# P (Z=1 | Y)
   		 ##################
		PostLikeLog_BG <- EachLike64[[1]] 
		PostLike_BG <- exp(PostLikeLog_BG) # P(path)P(Data|Path) in background chain
		PostLikeMerge <- PostLike_BG*PI.up[1] # P(backgroud)P(Data|Background, path)
		if(ncol(PostLikeMerge)==length(StateNameList[[1]]))colnames(PostLikeMerge) <- StateNameList[[1]]
		# assume P(Data|path,signal) for EE related pathes are 0
		for(i in colnames(EachLike64[[2]]))PostLikeMerge[,i] <- PostLikeMerge[,i]+PI.up[2]*exp(EachLike64[[2]][,i])
		PPMerge <- PostLikeMerge/rowSums(PostLikeMerge)



		WhichMax <- apply(PPMerge,1,function(oo){
							 tt=which.max(oo)
							 if(length(tt)>0)out=tt
								if(length(tt)==0)out=NA
				 out})
		WhichMax.NonNA.ind <- which(!is.na(WhichMax))
		WhichMax.NonNA <- WhichMax[WhichMax.NonNA.ind]
		MapTerm <- colnames(PPMerge)[WhichMax.NonNA]
		MapVal <- sapply(1:length(WhichMax.NonNA),function(i)PPMerge[i,WhichMax.NonNA[i]])
		names(MapTerm) = names(MapVal)=rownames(Data.NZ.OriName)[WhichMax.NonNA.ind]	


		#########################
		# Rename
		#########################
		rownames(PPtmp0) = names(PPMaxT0)=names(PPMaxVal0)=
		rownames(LogPILike)=rownames(PILike)=names(SumPILike)=
		rownames(PPMerge)<-
		rownames(Data.NZ.OriName)
		
		# log likwlihood
		# -Inf are imputed with smallest readable valyes
		SumPILikeForLL <- log(SumPILike)
		ISNA <- which(SumPILikeForLL==-Inf | is.na(SumPILikeForLL))
		ISNOTNA <- c(1:length(SumPILikeForLL))[-ISNA]
		Min.ISNOT <- min(SumPILikeForLL[ISNOTNA])
		SumPILikeForLL[ISNA] <- Min.ISNOT
		SumPILikeForLL.Sum <- sum(SumPILikeForLL)

		####################
		#Tran
		####################
		TranOut <- sapply(1:length(EachEB),function(i)EachEB[[i]]$TranTrack[[length(EachEB[[i]]$TranTrack)]],simplify=FALSE)
		Pi0Out <- sapply(1:length(EachEB),function(i)EachEB[[i]]$Pi0Track[length(EachEB[[i]]$TranTrack),],simplify=FALSE)

		####################
		#
		####################
		# conditional probabilities
		CondC1PPAll <- EachEB[[1]]$PP
		CondC1PatAll <- EachEB[[1]]$AllTermChar
	
		CondC1MAPChar <- apply(CondC1PPAll,1,function(i){
                                        tt <- CondC1PatAll[which.max(i)]
					 if(length(tt)==0)tt <- NA
						tt})
		CondC1MaxVal <- sapply(1:nrow(CondC1PPAll),function(i){
				    if(is.na(CondC1MAPChar[i]))tt <- NA
				    else tt <-  CondC1PPAll[i,CondC1MAPChar[i]]
				 tt})
		names(CondC1MaxVal)=rownames(CondC1PPAll)=
		names(CondC1MAPChar)=rownames(Data.NZ.OriName)
		
		CondC2PPAll <- EachEB[[2]]$PP
		CondC2PatAll <- EachEB[[2]]$AllTermChar
		CondC2MAPChar <- apply(CondC2PPAll,1,function(i){
                                        tt <- CondC2PatAll[which.max(i)]
					 if(length(tt)==0)tt <- NA
						tt})
	
		CondC2MaxVal <- sapply(1:nrow(CondC2PPAll),function(i){
				    if(is.na(CondC2MAPChar[i]))tt <- NA
				    else tt <-  CondC2PPAll[i,CondC2MAPChar[i]]
				 tt})
		names(CondC2MaxVal)=rownames(CondC2PPAll)=
		names(CondC2MAPChar)=rownames(Data.NZ.OriName)
	
		# Marginal probability of gradient patterns
		PI <- PI.track
		SigPPTerm <- colnames(CondC2PPAll)
		SigPPCB <- (CondC2PPAll*PI[2]+CondC1PPAll[,SigPPTerm]*PI[1])
		SigMAPCBChar <- apply(SigPPCB,1,function(i){
                                        tt <- SigPPTerm[which.max(i)]
					 if(length(tt)==0)tt <- NA
						 tt})
		SigMaxValCB <- sapply(1:nrow(SigPPCB),function(i){
		                         if(is.na(SigMAPCBChar[i]))tt <- NA
                                           else tt <- SigPPCB[i,SigMAPCBChar[i]]
					 tt})
		names(SigMAPCBChar)=names(SigMaxValCB)=rownames(Data.NZ.OriName)

		
		# Marginal probability for all patterns
		SigPPLarge <- cbind(SigPPCB,CondC1PPAll[,colnames(CondC1PPAll)[which(!colnames(CondC1PPAll)%in%colnames(SigPPCB))]]*PI[1])
		SigMAPLargeChar <- apply(SigPPLarge,1,function(i){
				          tt <- colnames(SigPPLarge)[which.max(i)]
					  if(length(tt)==0)tt <- NA
							 tt})
		SigMaxValLarge <- sapply(1:nrow(SigPPLarge),function(i){
		                                        if(is.na(SigMAPLargeChar[i]))tt <- NA
							  else tt <- SigPPLarge[i,SigMAPLargeChar[i]]
					                tt})
		names(SigMAPLargeChar)=names(SigMaxValLarge)=rownames(Data.NZ.OriName)



		out <- list(
						PPMatW=PPtmp0, PPMatZ=PPMerge,MAPZ=MapTerm, MAPVal=MapVal,
						SignalPattern=colnames(EachLike64[[2]]),
						LL=SumPILike,
						PPMaxIdxW=PPMaxT0, PPMaxValW=PPMaxVal0, fW=PILike,
						AllLike=EachLike64, PI=PI.track, TranOut=TranOut, Pi0Out=Pi0Out,
						ResInChain=EachEB,LLEach=SumPILikeForLL,LLSum=SumPILikeForLL.Sum,
					       	CondC1PP=CondC1PPAll, CondC1Pat=CondC1PatAll, CondC1MAPChar=CondC1MAPChar,
						CondC1MaxVal=CondC1MaxVal,
					     	CondC2PP=CondC2PPAll, CondC2Pat=CondC2PatAll, CondC2MAPChar=CondC2MAPChar,
						CondC2MaxVal=CondC2MaxVal,
						MgGrdPP=SigPPCB,MgGrdTerm=SigPPTerm, MgGrdMAPChar=SigMAPCBChar,
						MgGrdMaxVal=SigMaxValCB, MgAllPP=SigPPLarge, MgAllMAPChar=SigMAPLargeChar,
						MgAllMaxVal=SigMaxValLarge, QList=QOut, Alpha=matrix(AlphaMom,nrow=1),
						Beta=matrix(BetaMom,nrow=1), PPpattern="PPpattern"
						)
}
