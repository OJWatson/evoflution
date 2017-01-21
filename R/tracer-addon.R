#' Tracer Addon Function
#'
#' This function links with Tracer_Accumulation by creating the summary.txt file
#' needed, by loading the log files as needed and first writing it.
#'
#'
#' Tree Statistics pulled are internally coded to be the posterior, clock rate,
#' growthRate, ePop and TMRCA.
#'
#' The summary statistics pulled are internally coded to be the mean, standard error,
#' median, HPD 95% upper and lower estimate and effective sample size.
#'
#' The output is written to file for further use
#'
#'
#' @param cmd shellscript path name for loganalyser.cmd. Tracer Addon uses the internal
#' function system(), so path must be written accordingly.
#' @param file.dir directory where the sumamry.txt file will be written to
#' @param log.dir directory from which the log files are to be read from
#' @param reps number of reps within the log.dir
#'
#' @return writes summary.txt file
#'
#'
#' @export
#'
#'
Tracer_Addon <- function(cmd = paste("\"C:\\Users\\Oliver\\GoogleDrive\\AcademicWork\\Imperial\\O15-12\\BEAST\\BEAST v1.8.2\\bin\\loganalyser.cmd\""),
                         file.dir,log.dir,reps, replicate.name = "replicate"){

  ## create log directory where log files are
  log2.dir <- gsub(pattern = "/",replacement = "\\",fixed=TRUE,log.dir)
  tbcols <- c("mean","stdErr","median","hpdLower","hpdUpper","ESS","X50hpdLower","X50hpdUpper","*")

  for (i in 1:reps){

    sum.dir <- paste(log2.dir,replicate.name,"_summary",i,".txt",sep="")
    log3.dir <- paste(log2.dir,replicate.name,i,".log",sep="")
    sum.dir <- paste("\"",sum.dir,"\"",sep="")
    log3.dir <- paste("\"",log3.dir,"\"",sep="")

    cmd2 <- paste(cmd,log3.dir,sum.dir)
    system(cmd2)

    res <- data.frame(list(replicate.log..posterior=rep(0,6),  replicate.log..TreeHeight=0, replicate.log..clockRate=0,
                           replicate.log..ePopSize=0,   replicate.log..growthRate=0))
    rownames(res) <- c("mean","stdErr","median",
                       "95% HPD Interval lower","95% HPD Interval higher","effective sample size (ESS)")

    tb <- read.table(paste(log.dir,replicate.name,"_summary",i,".txt",sep=""),row.names = 1,skip = 1,col.names = tbcols,nrows = 10,sep = "\t")
    res["mean",] <-  tb$mean[c(1,5,6,9,10)]
    res["stdErr",] <- tb$stdErr[c(1,5,6,9,10)]
    res["median",] <- tb$median[c(1,5,6,9,10)]
    res["95% HPD Interval lower",] <- tb$hpdLower[c(1,5,6,9,10)]
    res["95% HPD Interval higher",] <- tb$hpdUpper[c(1,5,6,9,10)]
    res["effective sample size (ESS)",] <- tb$ESS[c(1,5,6,9,10)]

    if (i==1){
      all <- res
    } else {
      all <- cbind(all,res)
    }


  }

  write.table(as.data.frame(all),paste(file.dir,replicate.name,"_summary.txt",sep=""),sep="\t",quote=FALSE,col.names = TRUE,row.names = TRUE)

}
