#' BEAST_R0_Network_Pull_Tracer
#'
#' Summary Statistics are clock rate, growthRate, ePop and TMRCA.
#'
#'
#' @param sampling.reps Sampling replicates total. Default = 20
#' @param rep.name The character name used for each sampling rep, i.e. Unifrom, Random, Oldest_Bar_One
#' @param outgroup Boolean describing if an outgroup was present
#' @param tree.reps number of subdirectories to be scanned
#' @param vars Tracer variables t collect in total. Default = 5
#' @param var.R0 Size parameter for offspring distribtuion. Default = 0
#' @param sigma value Sigma value used
#' @param N.v column vector of N. Default: N.v=c(rep(100000,tree.reps)),
#' @param mu.v column vector of mutation rates in directory. Default: mu.v=c(rep(1.57e-05,tree.reps)),
#' @param Tg.v column vector of generation times in directory. Default: Tg.v=c(rep(2.6,tree.reps)),
#' @param R0.v column vector of R0 in directory. Default: R0.v=c(rep(2,tree.reps)),
#' @param network.root.dir  network root directory
#' @param local.root.dir local root directory
#'
#' @importFrom ape read.nexus.data
#' @return Returns a data frame of collated statstics
#'
#'
#'
#'

BEAST_R0_Network_Pull_Tracer <- function(sampling.reps=20,
                                         rep.name="bs",
                                         outgroup=FALSE,
                                         tree.reps=10,
                                         vars=5,
                                         var.R0=0,
                                         sigma=1,
                                         N.v=c(rep(1e5,tree.reps)),
                                         mu.v=c(rep(1.57e-05,tree.reps)),
                                         Tg.v=c(rep(2.6,tree.reps)),
                                         R0.v=c(rep(2,tree.reps)),
                                         network.root.dir,
                                         local.root.dir){

  ## WARNINGS ##

  ## AUXILLARY FUNCTIONS ##
  nex_read <- function(nex.file, case.tree, tree.age){

    # pull nexus file and isolate the names for meta data
    nex <- ape::read.nexus.data(nex.file)
    nex.names <- names(nex)
    split.names <- strsplit(nex.names,split = "_")

    # set up results to store from each split
    Seq.Res <- vector(mode="character",length=length(split.names))
    Time.Res <- vector(mode="character",length=length(split.names))

    for (i in 1:length(split.names)){

      Seq.Res[i] <- split.names[[i]][2]
      Time.Res[i] <- split.names[[i]][3]

    }

    # calculate needed statistics before returning
    Sample.Size <- length(Seq.Res)
    Sequence.Count <- length(unique(Seq.Res))
    Time.Span <- max(as.numeric(Time.Res))-min(as.numeric(Time.Res))
    Youngest.Sample <- tree.age-max(as.numeric(Time.Res))
    Oldest.Sample <- min(as.numeric(Time.Res))
    Percent.Cover <- Time.Span/max(case.tree$Time)

    res <- c(Sample.Size,Sequence.Count,Time.Span, Youngest.Sample,Oldest.Sample, Percent.Cover)
    return(res)
  }

  RMSE_calculation <- function(Full.Table){

    l <- rep(0,dim(Full.Table)[1])
    Errors <- list("RMedSE_mu"=l,"RMedSE_Gr"=l,"RMedSE_ePop"=l,"RMedSE_TMRCA"=l,"RMeanSE_mu"=l,"RMeanSE_Gr"=l,"RMeanSE_ePop"=l,"RMeanSE_TMRCA"=l)

    Errors$RMedSE_mu <- sqrt((Full.Table$O_Mu_Med - Full.Table$E_Mu)^2)
    Errors$RMedSE_Gr <- sqrt((Full.Table$O_Gr_Med - Full.Table$E_Gr)^2)
    Errors$RMedSE_ePop <- sqrt((Full.Table$O_ePop_Med - Full.Table$E_ePop)^2)
    Errors$RMedSE_TMRCA <- sqrt((Full.Table$O_TMRCA_Med - Full.Table$E_TMRCA)^2)
    Errors$RMeanSE_mu <- sqrt((Full.Table$O_Mu_Mean - Full.Table$E_Mu)^2)
    Errors$RMeanSE_Gr <- sqrt((Full.Table$O_Gr_Mean - Full.Table$E_Gr)^2)
    Errors$RMeanSE_ePop <- sqrt((Full.Table$O_ePop_Mean - Full.Table$E_ePop)^2)
    Errors$RMeanSE_TMRCA <- sqrt((Full.Table$O_TMRCA_Mean - Full.Table$E_TMRCA)^2)
    res <- cbind(Full.Table,data.frame(Errors))

    return(res)
  }

  ## Handle network path

  network.dir <- network.root.dir
  root.dir <- local.root.dir


  for (i in 1:tree.reps){

    ## create file and directory locations from the given root.dir
    case.file <- paste(root.dir,"rep",i,"/pseudo_case_tree.rds",sep="")
    summary.file <- paste(root.dir,"rep",i,"/",rep.name,"_summary.txt",sep="")

    nex.dir <- paste(root.dir,"rep",i,"/",sep="")
    case.tree <- readRDS(case.file)

    tree.age <- max(unlist(case.tree$sim.out$times))
    population.size <- case.tree$sim.out$total.infections

    ## load the summary statistic file exported from TRACER

    if(file.exists(summary.file)==FALSE){
      Tracer_Addon(file.dir=nex.dir,log.dir=paste(network.dir,"rep",i,"\\",sep=""),reps=sampling.reps,replicate.name = rep.name)
    }

    all <- read.delim(summary.file, row.names=1,stringsAsFactors = FALSE)

    ## create results

    Result <- list("Sampling"=rep(rep.name,sampling.reps),"Outgroup"=rep(outgroup,sampling.reps),"Sigma"=rep(sigma,sampling.reps),"Var"=rep(var.R0,sampling.reps),"Rep"=1:sampling.reps,"Tg"=rep(Tg.v[i],sampling.reps),
                   "E_Mu"=rep(mu.v[i],sampling.reps), "O_Mu_Med"=as.numeric(all["median",seq(3,vars*(sampling.reps-1)+3,vars)]),
                   "O_Mu_Mean"=as.numeric(all["mean",seq(3,vars*(sampling.reps-1)+3,vars)]),"O_Mu_StdErr"=as.numeric(all["stdErr",seq(3,vars*(sampling.reps-1)+3,vars)]),
                   "ESS_mu"=as.numeric(all["effective sample size (ESS)",seq(3,vars*(sampling.reps-1)+3,vars)]),
                   "O_Mu_95L"=as.numeric(all["95% HPD Interval lower",seq(3,vars*(sampling.reps-1)+3,vars)]),
                   "O_Mu_95H"=as.numeric(all["95% HPD Interval higher",seq(3,vars*(sampling.reps-1)+3,vars)]),
                   "E_Gr"=rep((log(R0.v[i])/Tg.v[i]),sampling.reps), "O_Gr_Med"=as.numeric(all["median",seq(5,vars*(sampling.reps-1)+5,vars)]),
                   "O_Gr_Mean"=as.numeric(all["mean",seq(5,vars*(sampling.reps-1)+5,vars)]),"O_Gr_StdErr"=as.numeric(all["stdErr",seq(5,vars*(sampling.reps-1)+5,vars)]),
                   "ESS_Gr"=as.numeric(all["effective sample size (ESS)",seq(5,vars*(sampling.reps-1)+5,vars)]),
                   "O_Gr_95L"=as.numeric(all["95% HPD Interval lower",seq(5,vars*(sampling.reps-1)+5,vars)]),
                   "O_Gr_95H"=as.numeric(all["95% HPD Interval higher",seq(5,vars*(sampling.reps-1)+5,vars)]),
                   "E_ePop"=rep(population.size,sampling.reps), "O_ePop_Med"=as.numeric(all["median",seq(4,vars*(sampling.reps-1)+4,vars)]),
                   "O_ePop_Mean"=as.numeric(all["mean",seq(4,vars*(sampling.reps-1)+4,vars)]),"O_ePop_StdErr"=as.numeric(all["stdErr",seq(4,vars*(sampling.reps-1)+4,vars)]),
                   "ESS_ePop"=as.numeric(all["effective sample size (ESS)",seq(4,vars*(sampling.reps-1)+4,vars)]),
                   "O_ePop_95L"=as.numeric(all["95% HPD Interval lower",seq(4,vars*(sampling.reps-1)+4,vars)]),
                   "O_ePop_95H"=as.numeric(all["95% HPD Interval higher",seq(4,vars*(sampling.reps-1)+4,vars)]),
                   "E_TMRCA"=rep(tree.age,sampling.reps), "O_TMRCA_Med"=as.numeric(all["median",seq(2,vars*(sampling.reps-1)+2,vars)]),
                   "O_TMRCA_Mean"=as.numeric(all["mean",seq(2,vars*(sampling.reps-1)+2,vars)]),"O_TMRCA_StdErr"=as.numeric(all["stdErr",seq(2,vars*(sampling.reps-1)+2,vars)]),
                   "ESS_TMRCA"=as.numeric(all["effective sample size (ESS)",seq(2,vars*(sampling.reps-1)+2,vars)]),
                   "O_TMRCA_95L"=as.numeric(all["95% HPD Interval lower",seq(2,vars*(sampling.reps-1)+2,vars)]),
                   "O_TMRCA_95H"=as.numeric(all["95% HPD Interval higher",seq(2,vars*(sampling.reps-1)+2,vars)]))


    ## collate results before returning
    if(i==1){
      Full.Table <- data.frame(Result)
    } else {
      Full.Table <- rbind(Full.Table, data.frame(Result))
    }

  }
  res <- RMSE_calculation(Full.Table = Full.Table)
  return(res)
}
