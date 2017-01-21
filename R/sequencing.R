#' Forwards genetic sequence simulator
#'
#' \code{Sequencing} uses the output of \code{Infection_Sim} to simualate the evolutionary
#' process of the seqeunce from the root forwards
#'
#'
#' @param sim.out Output of \code{Infection_Sim}
#' @param seq.length Sequence length. Will be overwritten accordingly if root.path given. Default = 1701
#' @param mu Daily mutation rate. Default = 1.57e-05 sub / site / day
#' @param Tg mean generation time. Default = 2.6 days
#' @param model Nucleic subsitution model. Can be one of "JC69", "K80", "HKY85" (Default)
#' @param kappa Tv/Ts ration. Default = 2
#' @param root.path File path of saveRDS(root). Default = NULL
#' @param evolution.model String variable detailing evolutionary model. Can be any of "generational", "within.host" and "parallel"
#'
#' @return Returns a list of 4 lists
#'
#' @export
#'
#' @importFrom ape as.DNAbin base.freq makeLabel
#' @importFrom seqinr s2n
#'


Sequencing <- function(sim.out,seq.length=1701,mu=1.57e-05,Tg=2.6,model="HKY85",kappa=2,
                       root.path=NULL,evolution.model="generational",outgroup=TRUE){


  ## create reference bin
  nuc <- ape::as.DNAbin(c("a","c","g","t"))

  ## assign evoluttionary model
  generational <- within.host <- parallel <- FALSE
  model.choice <- match(evolution.model,c("generational","within.host","parallel"))
  if(model.choice==1){
    generational <- TRUE
  } else if(model.choice==2){
    within.host <- TRUE
  } else {
    parallel <- TRUE
  }


  ## EXTRA FUNCTIONS ##
  ##########################################################

  ## function to create sequence from scratch
  seq_generate <- function(nuc, seq.length){
    res <- t(sample(nuc, size=seq.length, replace=TRUE))
    #class(res) <- "DNAbin"
    return(res)
  }

  ## function that takes sequence, and mutates seq.change positions
  seq_evolve <- function(sequence, mu, time, model, seq.changes){
    # first create the transition matrix P
    if (model == "JC69"){
      P <- JC69(mu = mu,t = time)
    } else if (model == "HKY85"){
      P <- as.matrix(HKY85(mu = mu, t = time, kappa = kappa, pi = ape::base.freq(sequence)))
    }

    diag(P) <- 0

    random.nuc.sites <- sample(seq.length, size = seq.changes, replace = F)
    if(seq.changes == 1){
      sequence[random.nuc.sites] <- nuc[sample(4,1,prob = P[(s2n(as.character(sequence[random.nuc.sites]),base4 = FALSE)),])]
    } else {
      sequence[random.nuc.sites] <- nuc[(which(apply(P[(s2n(as.character(sequence[random.nuc.sites]),base4 = FALSE)),], 1, rmultinom, n=1, size=1)==1)-
                                           (0:(seq.changes-1))*4)]
    }

    return(sequence)
  }

  ## function to evaluate probability of mutation event
  mutation_prob <- function(sequence, mu, model, kappa, Tg){
    # first assign the class of DNAbin to sequence (funny)
    #class(sequence)<-"DNAbin"
    # second create the transition matrix P
    if (model == "JC69"){
      P <- JC69(t = Tg,mu = mu)
    } else if (model == "HKY85"){
      P <- HKY85(t = Tg, mu = mu, kappa=kappa, pi=ape::base.freq(sequence))
    } else if (model == "K80") {
      P <- K80(t = Tg, mu = mu, kappa=kappa)
    }
    # using the transition matrix and the absolute numbers of nucelotides caluclate the probability
    # that no mutations occur
    abs.freq <- ape::base.freq(sequence,freq=TRUE)
    res <- 1 - ( (P[1,1])^abs.freq[1]*(P[2,2])^abs.freq[2]*(P[3,3])^abs.freq[3]*(P[4,4])^abs.freq[4] )
    return(res)

  }

  ## uses vector of times from one parent to calcualte mutation probabilities for each time cumulatively
  within_host_prob <- function(time.vector,sequence,model="HKY85",kappa=2,sequence.length=seq.length){

    sorted.times <- sort(time.vector)

    # Create differences of times, i.e. first time, second time - first time, third time - second etc
    if(length(sorted.times)>1){
      relative.sorted.times <- c(sorted.times[1],sorted.times[2:length(sorted.times)] - sorted.times[1])
    } else {
      relative.sorted.times <- sorted.times
    }

    # Provide time vectors to seqeuence models which will yield equal length vector of probs of mutation.
    if (model == "JC69"){
      res <- sapply(relative.sorted.times,JC69,mu=mu,length=sequence.length)
    } else if (model == "HKY85"){
      res <- sapply(relative.sorted.times,HKY85,mu=mu,kappa=kappa, pi=ape::base.freq(sequence),length=sequence.length)
    } else if (model == "K80") {
      res <- sapply(relative.sorted.times,K80,mu=mu,kappa=kappa, length=sequence.length)
    }

    return(res)

  }


  ## function to evolve sequence in time explicit manner
  seq_time_evolve <- function(sequence,model="HKY85",kappa=2,time,sequence.length=seq.length){

    # first create the transition matrix P
    if (model == "JC69"){
      P <- JC69(mu = mu,t = time)
    } else if (model == "HKY85"){
      P <- HKY85(mu = mu, t = time, kappa = kappa, pi = ape::base.freq(sequence))
    }

    # Apply transition matrix across whole sequence
    res <- t(nuc[(which(apply(P[(s2n(as.character(sequence),base4 = FALSE)),], 1, rmultinom, n=1, size=1)==1)-(0:(sequence.length-1))*4)])

    return(res)
  }



  ## START MAIN FUNCTION ##
  ############################################################

  ## Create root
  if(!is.null(root.path)){
    root <- readRDS(root.path)
    seq.length <- dim(root)[2]
  } else {
    root <- seq_generate(nuc,seq.length)
  }

  ## sequence list
  sequences <- list()
  sequences[[1]] <- root

  ## what sequence id each individual has
  individuals.sequences <- vector(mode = "numeric",length=sim.out$total.infections)
  individuals.sequences[1] <- 1

  ## counters
  individual.counter <- 2
  parent.counter <- 1
  sequence.counter <- 1

  ## remove non-infecting individuals
  gens <- unlist(sim.out$gens)
  if(!sum(gens==0)){
    gens.no.zeros <- gens
  } else {
  gens.no.zeros <- gens[-which(gens==0)]
  }

  ## create vector of times and parents
  times <- c(0,unlist(sim.out$times))
  parents <- unlist(sim.out$parents)

  ## log counter
  print.counter <- seq(1,length(parents),1000)

  ## generational model
  if(generational){

    # creat list of mutation probailities, where each new sequence will have its probablity of mutation generated.
    # Will ultimately look very similar as the generational model assumes the same length of evoltuion time of Tg
    mutation.probs <- list()
    mutation.probs[[1]] <- mutation_prob(sequence = root, mu = mu, Tg = Tg, model = model, kappa = kappa)

    for(i in 1:length(parents)){

      # messsage loggers
      if(is.element(i,print.counter)){
        message(i)
      }

      # draw whether there are any mutations
      n.seq.changes <- rpois(1,mutation.probs[[individuals.sequences[parents[i]]]])

      # if changes then allocate and update sequence and individual counters
      if(n.seq.changes){

        sequence.counter <- sequence.counter + 1
        individuals.sequences[individual.counter:(individual.counter + gens.no.zeros[i]-1)] <- sequence.counter
        sequences[[sequence.counter]] <- seq_evolve(sequence = sequences[[individuals.sequences[parents[i]]]],
                                                    mu = mu, time = Tg, model = model, seq.changes = n.seq.changes)
        mutation.probs[[sequence.counter]] <- mutation_prob(sequences[[sequence.counter]],
                                                            mu = mu, Tg = Tg, model = model, kappa = kappa)
        individual.counter <- individual.counter + gens.no.zeros[i]

        # otherwise just update the individual counter
      } else {

        individuals.sequences[individual.counter:(individual.counter + gens.no.zeros[i]-1)] <- sequence.counter
        individual.counter <- individual.counter + gens.no.zeros[i]

      }

    }

  }

  ## within.host model
  if(within.host){

    for(i in 1:length(parents)){

      # message logger
      if(is.element(i,print.counter)){
        message(i)
      }

      # take offspring times and sequence of parent and produce the probability of mutation in each time window
      offspring.times <- times[individual.counter:(individual.counter + gens.no.zeros[i]-1)]

      starting.sequence <- sequences[[individuals.sequences[parents[i]]]]

      seq.change.probs <- within_host_prob(time.vector = offspring.times,sequence=starting.sequence,
                                           model=model,kappa=kappa,sequence.length=seq.length)

      # draw number of seq changes
      n.seq.changes <- as.numeric( sapply(seq.change.probs,FUN = function(x){return(rpois(n = 1,lambda = x))}) )

      # check for no mutations
      if(sum(n.seq.changes)){

        # cycle through the mutations
        for(j in 1:length(n.seq.changes)){

          # catch no mutations and allocate
          if(!n.seq.changes[j]){

            individuals.sequences[individual.counter] <- sequence.counter
            individual.counter <- individual.counter + 1

            # if there is a mutation then find out where and update as in generational model
          } else {

            sequence.counter <- sequence.counter + 1
            individuals.sequences[individual.counter] <- sequence.counter
            sequences[[sequence.counter]] <- seq_evolve(sequence = starting.sequence,seq.changes = n.seq.changes[j],
                                                        mu = mu, time = Tg, model = model)
            starting.sequence <- sequences[[sequence.counter]]
            individual.counter <- individual.counter + 1

          }

        }

      } else {

        individuals.sequences[individual.counter:(individual.counter + gens.no.zeros[i]-1)] <- sequence.counter
        individual.counter <- individual.counter + gens.no.zeros[i]

      }

    }

  }

  ## create labels
  seq.labels <- ape::makeLabel(rep("seq",sequence.counter+1),make.unique = TRUE)
  for(l in 1:sequence.counter){
    row.names(sequences[[l]]) <- seq.labels[[l]]
  }

  ## if you want an outgroup then create sequence evolved over length time equal to greatest root tip length
  if(outgroup){
    sequence.counter <- sequence.counter + 1
    individuals.sequences <- c(individuals.sequences, sequence.counter)
    sequences[[sequence.counter]] <- seq_time_evolve(sequence = sequences[[1]],model=model,kappa=kappa,time=max(times))
    row.names(sequences[[sequence.counter]]) <- paste("OUTGROUP_",seq.labels[[sequence.counter]],sep="")
  }

  # return list of results
  return(list("sequences"=sequences,"individuals.sequences"=individuals.sequences,
              "seq.labels"=seq.labels,"sequence.counter"=sequence.counter))


}



