
#' @title Identify DE genes and classify them into their most likely path in an RNA-seq experiment with ordered conditions
#' @usage EBSeqHMMTest(Data, 
#'	NgVector=NULL, Conditions, AllTran=NULL, 
#'	AllPi0=NULL, Terms=NULL,  
#'	sizeFactors, NumTranStage=c(3,2),FCV=2, 
#'	homo=FALSE, UpdateRd=10, PIBound=.9, UpdatePI=FALSE,
#'	Print=FALSE,WithinCondR=TRUE,Qtrm=.75,QtrmCut=10,
#'	PenalizeLowMed=TRUE, PenalizeLowMedQt=.1,PenalizeLowMedVal=10)
#' @param Data input data, rows are genes and columns are samples
#' @param NgVector Ng vector; NULL for gene level data
#' @param Conditions A factor indicates the condition (time/spatial point) which each sample belongs to.
#' @param AllTran initial values for transition matrices
#' @param AllPi0 initial values for starting probabilities
#' @param Terms Terms
#' @param FCV candidate values for expected FC. Default is 2.  If user
#' wants to search through a list of candidate FCs, he/she may define FCV as a vector. e.g. By defining FCV=seq(1.4,2,.2), the function 
#' will search over (1.4 1.6 1.8 2.0). Note that searching over a number of candidate FCs will increase the running time.
#' @param sizeFactors a vector indicates library size factors
#' @param NumTranStage number of states in two chains
#' @param homo whether the chain is assumed to be homogenious
#' @param UpdateRd max number of iteration
#' @param UpdatePI whether update the mixture proportion of two chains
#' @param PIBound upper bound of the mixture proportion of the two chains
#' @param Qtrm,QtrmCut Transcripts with Qtrm th quantile < = QtrmCut  will be removed before testing. The default value is Qtrm = 0.75 and QtrmCut=10.
#' By default setting, transcripts that have >75\% of the samples with expression less than 10
#' won't be tested.
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
#' EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
#'           UpdateRd=2)
#' @details EBSeqHMMTest() function applies EBSeqHMM model with differentexpected FC's and select the optimal FC that maximizes the log likelohood. 
#' EBSeqHMMTest() calls EBHMMNBMultiEM_2chain() function which implements the EBSeqHMM model to perform statistical analysis in an RNA-seq experiment with ordered conditions based on a fixed expected FC.
#' EBSeqHMMTest() runs EBHMMNBMultiEM_2chain() with varying FCs (default is seq(1.4,2,.2)).
#' And it will return the results of the model with optimal FC.
#' Here the emission distribution of each gene is assumed to be a 
#' Beta-Negative Binomial distribution with parameters (r_g, alpha, beta)
#' , in which alpha and beta are shared by all the genes and r_g is gene specific.
#' If not specified, r_g, alpha and beta will be estimated using method of moments.
#' For isoform data, we assume isoforms from the same Ig group share the same beta^Ig. alpha is shared by all the isoforms and r_gi is isoform specific.
#' The user also needs to specify an expected FC.
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
#' FCLikelihood: log likelihood of each FC.

adEBSeqHMMTest <- function(Data, 
												NgVector=NULL, Conditions,
												AllTran=NULL, 
												AllPi0=NULL, Terms=NULL,  
												sizeFactors, NumTranStage=c(3,2)
												,FCV=2, 
												homo=FALSE, UpdateRd=10,
												 PIBound=.9, UpdatePI=FALSE,Print=FALSE
												,WithinCondR=TRUE,Qtrm=.75,QtrmCut=10,
											PenalizeLowMed=TRUE, PenalizeLowMedQt=.1,PenalizeLowMedVal=10,s2c=NULL,filter_fold=0.1	){

	if(!is.factor(Conditions))stop("Object Conditions is not in factor format")
	if(is.null(rownames(Data)))stop("Please add gene/isoform names to the data matrix")

	if(!is.matrix(Data))stop("The input Data is not a matrix")
	if(length(Conditions)!=ncol(Data))stop("The number of conditions is not the same as the number of samples! ")
	if(nlevels(Conditions)<2)stop("Less than 2 conditions - Please check your input")
	if(ncol(Data)<=2) stop("Less than 3 points - Please check your input")
	if(length(sizeFactors)!=length(Data) &  length(sizeFactors)!=ncol(Data))
	stop("The number of library size factors is not the same as the number of samples!")		

	cat("\n Conditions are ordered as: \n ")
	cat(levels(Conditions))

	Dataraw <- Data

	DataNorm <- GetNormalizedMat(Data, sizeFactors)
	Levels <- levels(as.factor(Conditions))


	QuantileFor0 <- apply(DataNorm,1,function(i)quantile(i,Qtrm))
	AllZeroNames <- which(QuantileFor0<=QtrmCut)
	NotAllZeroNames <- which(QuantileFor0>QtrmCut)
	if(length(AllZeroNames)>0 & Print==TRUE)
	   cat(paste0("Removing transcripts with ",Qtrm*100,
						" th quantile < = ",QtrmCut," \n",
	length(NotAllZeroNames)," transcripts will be tested \n"))
	if(length(NotAllZeroNames)==0)stop("0 transcript passed")
	Data <- Data[NotAllZeroNames,]
	if(!is.null(NgVector))NgVector <- NgVector[NotAllZeroNames]
	if(length(sizeFactors)==length(Data))sizeFactors <- sizeFactors[NotAllZeroNames,]


	List <- vector("list",length(FCV))
	LL <- rep(NA,length(FCV))
	for(i in 1:length(FCV)){
	cat(paste0("\n Estimating parameters when expected FC = ",FCV[i] ," \n"))
	List[[i]] <- suppressWarnings(adEBHMMNBMultiEM_2chain(Data,
            NgVector=NgVector, Conditions=Conditions,
                         AllTran=AllTran,filter_fold=filter_fold,
                        AllPi0=AllPi0,Terms=Terms,			                                         
			sizeFactors=sizeFactors, PriorFC=FCV[i],
			StateNames=c("Up","EE","Down"),homo=homo, UpdateRd=UpdateRd, UpdatePI=UpdatePI,
			Print=Print, WithinCondR=WithinCondR,
			PenalizeLowMed=PenalizeLowMed,
			PenalizeLowMedQt=PenalizeLowMedQt,
			PenalizeLowMedVal=PenalizeLowMedVal,s2c=s2c))
	LL[i] <- List[[i]]$LLSum
}
Max <- which.max(LL)
names(LL) <- FCV
cat(paste0("\n Estimated expected FC ",FCV[Max] ," \n"))
Out <- c(List[[Max]], FCLikelihood=LL)

}
