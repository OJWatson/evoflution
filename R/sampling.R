#' Forwards genetic sequence simulator
#'
#' \code{New_sampling} uses the output of \code{Infection_Sim} and \code{Sequencing}
#' to sample a selection of sequences and write the annotated output top file as a fasta
#'
#' @param sim.out Output of \code{Infection_Sim}
#' @param seq.out Output of \code{Sequencing}
#' @param file Directory path where the nexus file of samples are to be saved to.
#' @param n.seq.samples Number of sequences to be sampled in total
#' @param sigma The interval parameter, where the interval of sequences is equal to 1/sigma. Default = 1, i.e. all sequences
#' @param bs.label String used to be added on to the created nex files. Default = "bs"
#' @param n.bs Number of bootstrap samples to create. Default = 50
#' @param oldestboolean Boolean if the oldest sample is included (TRUE). Default=TRUE
#' @param outgroup Boolean if the outgroup  is included (TRUE). Default=FALSE
#' @param discretised Boolean if discretised sampling occurs. Default = FALSE
#' @param proportional Boolean if proportional sampling occurs. Default = FALSE
#' @param t.cont Boolean if the simulation occurs in discrete (FALSE) or continuous time (TRUE), in turn directing
#' whether decimilisation is needed in uniform sampling strategies.
#' @param missed.generations Number of intiial generations not considered for sampling. Default = 0
#'
#' @note If the combination of prop.times and extant==FALSE creates a time window that does not include all clades, and all.clades
#' == TRUE, only all the possible clades within the time window will be sampled.
#'
#' @return save a nexus file of dated sequences
#'
#' @importFrom ape write.nexus.data
#' @export
#'
#'

New_sampling <- function(sim.out, seq.out, file, n.seq.samples=50, bs.label = "rep",
                         sigma=1, n.bs=50, random=TRUE, oldestboolean=TRUE, outgroup=TRUE, discretised=FALSE,
                         proportional=FALSE, t.cont=TRUE, missed.generations=0){

  ## HANDLE VARIABLES ##
  ##########################################################

  ## Check all arguments given are valid
  if(n.seq.samples>(sim.out$total.infections-1)) stop("Number of sequences chosen is greater than total cases")
  if(sigma < 1) stop("sigma must be greater than 1")

  ## catch directory ending
  if(tail(unlist(strsplit(file,"")),1)!="/") file <- paste(file,"/",sep="")

  ## handle arguments and create variables
  if(missed.generations){
    missed.generation.ids <- 0
    for(m in 1:missed.generations){
      missed.generation.ids <- missed.generation.ids + sum(sim.out$gens[[m]])
    }
  } else {
    missed.generation.ids <- 2
  }


  all.parents <- unlist(sim.out$parents)

  ## everryone including root
  all.individuals <- 1:length(seq.out$individuals.sequences)

  ## any missed generation, if not from 2 to end
  all.possible.seqs <- missed.generation.ids:length(seq.out$individuals.sequences)
  all.times <- unlist(sim.out$times)

  ## root time, all times, outgroup
  if(outgroup){
  world.times <- c(0,all.times,max(all.times))
  } else {
    world.times <- c(0,all.times)
  }

  all.sorted.times <- sort(world.times)
  all.possible.times <- world.times[all.possible.seqs]
  total.gens <- length(sim.out$gens)

  ## EXTRA FUNCTIONS ##
  ##########################################################

  ## function to evolve sequence in time explicit manner
  seq_time_evolve <- function(sequence,model,time){
    # first create the transition matrix P
    if (model == "JC69"){
      P <- JC69(mu,Tg)
    } else if (model == "HKY85"){
      P <- HKY85(mu = mu, t = time, kappa = kappa, pi = base.freq(sequence))
    }

    res <- t(nuc[(which(apply(P[(s2n(as.character(sequence),base4 = FALSE)),], 1, rmultinom, n=1, size=1)==1)-(0:(sequence.length-1))*4)])

    return(res)
  }

  ## uniform in time for discrete time tree bootstrap selection of sequences
  uniform_discrete_times <- function(all.possible.seqs,all.possible.times,n.bs){

    get_times <- function(a,b) {

      out <- cbind(a, bval = NA)

      for (i in seq_along(a)) {
        #which value of B is closest?
        whichB <- which.min(abs(b - a[i]))
        #Assign that value to the bval column
        out[i, "bval"] <- b[whichB]
        #Remove that value of B from being chosen again
        b <- b[-whichB]
      }

      return(out[,2])

    }

    get_positions <- function(c,d,all.possible.seqs){

      res <- vector(length=n.seq.samples-2)

      for (i in 1:length(c)){
        posibs <- which(d==c[i])
        if(length(posibs)==1){
          res[i] <- (posibs)
        } else {
          res[i] <- (sample(x=posibs,size=1))
        }
        d[res[i]] <- Inf
      }

      return(all.possible.seqs[res])
    }


    b <- seq(from=min(all.possible.times),to=max(all.possible.times),
             by=((max(all.possible.times) - min(all.possible.times))/(n.seq.samples-3)))
    c <- get_times(b,sort(all.possible.times))
    d <- all.possible.times

    res <- replicate(n.bs,get_positions(c,d,all.possible.seqs))


    return(res)

  }

  ## uniform in time for continuous time tree bootstrap selection of sequences
  uniform_continuous_times <- function(all.possible.seqs,all.possible.times,n.bs){

    get_times <- function(a,b) {

      out <- cbind(a, bval = NA)

      for (i in seq_along(a)) {
        #which value of B is closest?
        whichB <- which.min(abs(b - a[i]))
        #Assign that value to the bval column
        out[i, "bval"] <- b[whichB]
        #Remove that value of B from being chosen again
        b <- b[-whichB]
      }

      return(out[,2])

    }

    d <- all.possible.times
    sorted.possible.times <- sort(all.possible.times)

    for (i in 1:n.bs){

      tip.check <- TRUE
      while(tip.check){

        b <- sort(runif(n=n.seq.samples-2,min=min(all.possible.times),max=max(all.possible.times)))
        c <- get_times(b,sorted.possible.times)

        if(length(unique(c))==n.seq.samples-2) {
          tip.check <- FALSE
        }

      }

      if (i==1){
        res <- all.possible.seqs[match(c,d)]
      } else {
        res <- cbind(res,all.possible.seqs[match(c,d)])
      }
    }
    return(res)

  }

  ## MAIN FUNCTION ##
  ##########################################################

  ## initialise results
  Sample.IDs <- list(ID=rep(0,n.seq.samples))

  ## first ensure that the youngest and oldest are sampled
  youngests <- which(all.possible.times==max(all.possible.times))
  oldests <- which(all.possible.times==min(all.possible.times))

  ## youngest at position 2 i.e most recent
  if(length(youngests)==1){
    Sample.IDs$ID[2] <- all.possible.seqs[youngests]
  } else {
    Sample.IDs$ID[2] <- sample(all.possible.seqs[youngests],1)
  }

  ## oldest at position 1, i.e. first infection from root
  if(oldestboolean==TRUE){
    if(length(oldests)==1){
      Sample.IDs$ID[1] <- all.possible.seqs[oldests]
    } else {
      Sample.IDs$ID[1] <- sample(all.possible.seqs[oldests],1)
    }
  } else {
    Sample.IDs$ID[1] <- sample(all.possible.seqs[-youngests],1)
  }


  if (outgroup){
    Sample.IDs$ID[2] <- length(seq.out$individuals.sequences)
  }
  ## work out the available cases from which to draw the sample from, with the root obviously not being possible to sample

  already.chosen <- match(Sample.IDs$ID[1:2],all.possible.seqs)
  ## all possible sequences therefore equal to the first case:interval
  all.possible.seqs <- all.possible.seqs[-already.chosen[!is.na(already.chosen)]]
  all.possible.times <- all.possible.times[-already.chosen[!is.na(already.chosen)]]

  ## if it's discretised restrict the available sequences which will then be randomly sampled
  if(discretised==TRUE){
    time.cut.off <- all.sorted.times[floor(length(all.possible.seqs)/sigma) + 1]
    too.high <- all.possible.times < time.cut.off
    all.possible.seqs <- all.possible.seqs[too.high]
    all.possible.times <- all.possible.times[too.high]
  }

  ## Now do the actual sampling for the remaining sequences

  if (random==TRUE){
    if(proportional==FALSE){
      bs.sample.sequences <-  replicate(n.bs,sample(all.possible.seqs,size=(n.seq.samples-2)))
    } else {
      bs.sample.sequences <-  replicate(n.bs,sample(all.possible.seqs,size=(n.seq.samples-2),prob=(1/(1+(sigma-1)*(all.possible.times/max(all.possible.times))))))
    }
  } else {
    if(t.cont){
      bs.sample.sequences <-  uniform_continuous_times(all.possible.seqs = all.possible.seqs,all.possible.times = all.possible.times,n.bs = n.bs)
    } else {
      bs.sample.sequences <-  uniform_discrete_times(all.possible.seqs = all.possible.seqs,all.possible.times = all.possible.times,n.bs = n.bs)
    }
  }

  for(i in 1:n.bs) {

    file2 <- paste(file,bs.label,i,".nex",sep="")
    bs.Sample.IDs <- Sample.IDs

    bs.Sample.IDs$ID[3:n.seq.samples] <- bs.sample.sequences[,i]

    bs.Sample.IDs$ID <- unique(bs.Sample.IDs$ID)

    extended <- list(labels=0)

    ## create tip and node labels in form n/t#_sequenceID_Time

    t <- paste("t",seq(1:((length(all.individuals))-length(all.parents))),sep="")
    tip.label <- paste(t,all.individuals[-all.parents],world.times[-all.parents],sep="_")
    n <- paste("n",seq(1:length(all.parents)),sep="")
    node.label <- paste(n,all.individuals[all.parents],world.times[all.parents],sep="_")

    extended$labels[all.individuals[-all.parents]] <- tip.label
    extended$labels[all.individuals[all.parents]] <- node.label

    Sampled_Sequences <- seq.out$sequences[[seq.out$individuals.sequences[bs.Sample.IDs$ID][1]]]
    for(s in 2:n.seq.samples){
      Sampled_Sequences <- rbind(Sampled_Sequences,seq.out$sequences[[seq.out$individuals.sequences[bs.Sample.IDs$ID][s]]])
    }


    row.names(Sampled_Sequences) <- extended$labels[bs.Sample.IDs$ID]
    ## use list of tip ids to select relevant sequences
    ape::write.nexus.data(Sampled_Sequences, file = file2, datablock = FALSE, interleaved = FALSE)

  }

}

