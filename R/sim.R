#' Forwards Infection Simulator
#'
#'
#' \code{Infection_sim} produces a list of generations detailing an infection process
#'
#'
#' @param N Target infected population size.  Default = 1e5
#' @param R0 basic reproduction number. Default = 2
#' @param Tg mean generation time. Default = 2.6 days
#' @param fixed.offspring Boolean dictating whether R0 offspring are always produced. Default = FALSE
#' @param scale rgamma scale parameter. Default = 1
#' @param size rnbinom size parameter Default = 1
#' @param real.time Boolean dictating. whether an extra generation is simulated producing
#' extra for real.time consideration.Default = FALSE
#' @param continous.time Boolean dictating whether gamma process is used for generation time distribution. Default = TRUE
#'
#' @return Returns a list of 5 lists
#'
#' @export
#'
#'
#'

Infection_sim <- function(N=1e5,R0=2,Tg=2.6,fixed.offspring=FALSE,
                          scale=1,size=1,real.time=FALSE, continous.time=TRUE){


  infections <- 0

  generations <- list()
  times <- list()
  parents <- list()
  infections <- list()

  generation <- 1

  if(fixed.offspring){
    first.generation.offspring <- 2
  } else {
  first.generation.offspring <- 0

  while(first.generation.offspring == 0){
    first.generation.offspring <- rnbinom(n = 1,size = size,mu = R0)
  }

  }

  generations[[1]] <- first.generation.offspring

  if(continous.time){
  times[[1]] <- rgamma(n=sum(first.generation.offspring),shape=Tg,scale=scale)
  } else {
    times[[1]] <- rpois(n=sum(first.generation.offspring),lambda = Tg)
  }

  infections[[1]] <- total.infections <- sum(first.generation.offspring)

  parents[[1]] <- 1

  next.parent <- 2

  if(real.time){
    time.sufficient <- TRUE
  } else {
    time.sufficient <- FALSE
  }

  while(total.infections < N | time.sufficient){

    if(total.infections > N){
      time.sufficient <- FALSE
    }

    generation <- generation + 1

    if(fixed.offspring){
      next_generation <- rep(2,sum(generations[[generation-1]]))
    } else {


    next_generation <- 0

    while(sum(next_generation) == 0){
      next_generation <- rnbinom(n = sum(generations[[generation-1]]),size = size,mu = R0)
    }

    }

    generations[[generation]] <- next_generation

    infections[[generation]] <- sum(next_generation)

    parents[[generation]] <- next.parent : (next.parent + (sum(generations[[generation]]>0)-1))


    next.parent <- total.infections + 2

    total.infections <- total.infections + infections[[generation]]


    ## time steps

    ## create a vector containing how many of each of the previous generations times are needed
    reps = generations[[generation]]
    x = 1:length(reps)
    t = reps > 0
    a = cumsum(reps[t])
    b = rep(0,max(a))
    b[a - reps[t] + 1] = 1
    x1 = x[t]
    out = x1[cumsum(b)]

    if(continous.time){
      times[[generation]] <- rgamma(n=infections[[generation]],shape=Tg,scale=scale) + times[[generation-1]][out]
    } else {
      times[[generation]] <- rpois(n=infections[[generation]],lambda = Tg) + times[[generation-1]][out]
    }

  }

  return(list("gens"=generations,"times"=times,"parents"=parents,
              "infections"=infections,"total.infections"=total.infections))

}
