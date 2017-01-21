#' JC69 Substiution Model
#'
#' \code{JC69} creates the rate matrix of the JC69 model before calculating the
#' related transition matrix, P, using the mutation and branch length/time from parent
#' variables.
#'
#'
#' @param mu - mutation rate
#' @param t - branch length, i.e. specific generation time
#' @param length Length of sequence that was used in creating the base frequencies.
#' If provided the output will instead be the probability of mutation occuring
#'
#' @export
#'
#' @return This function returns a \code{dgeMatrix} containing the transition matrix P, except if length is provided
#'
#' @importFrom Matrix expm
#'
JC69 <- function (t,mu,length=NULL){

  ## check supplied variables
  if(!is.numeric(mu)) stop("Mutation rate not valid")
  if(!is.numeric(t)) stop("Branch length not valid")

  ## create rate matrix
  a = 0
  b = 1
  Q = matrix(data = rbind(c(a,b,b,b),c(b,a,b,b),c(b,b,a,b),c(b,b,b,a)), ncol = 4)

  ## convert rate matrix to transition matrix
  x = 1/4 + 3/4 * exp(-4*mu*t)
  y = 1/4 - 1/4 * exp(-4*mu*t)
  a = x
  b = y
  P = matrix(data = rbind(c(a,b,b,b),c(b,a,b,b),c(b,b,a,b),c(b,b,b,a)), ncol = 4)

  if(!is.null(length)){
    P <- 1 - ( (P[1,1])^length )
  }


  ##return transition matrix
  return(P)
}

