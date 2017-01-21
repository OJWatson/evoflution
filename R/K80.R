#' K80 Substiution Model
#'
#' \code{K80} This function creates the rate matrix of the K80 model before calculating the
#' related transition matrix, P, using the mutation (mu), branch length/time (t) and transition
#' transversion ratio (kappah)
#'
#'
#' @param mu Mutation rate
#' @param t Branch length, i.e. specific generation time
#' @param kappa ransition/transversion rate ratio
#' @param pi Base Frequencies
#' @param length Length of sequence that was used in creating the base frequencies.
#' If provided the output will instead be the probability of mutation occuring
#'
#' @export
#'
#' @return This function returns a \code{dgeMatrix} containing the transition matrix P, except if length is provided
#'
#' @importFrom Matrix expm
#'

K80 <- function (t, mu, kappa, pi=c(0.25,0.25,0.25,0.25),length=NULL){

  ## check supplied variables
  if(!is.numeric(mu)) stop("Mutation rate not valid")
  if(!is.numeric(t)) stop("Branch length not valid")
  if(!is.numeric(kappa)) stop("Transition/transversion rate ratio not valid")
  if(length(pi)==1) stop("Epidemic died out")


  ## create the generator matrix R incorporating the transition/transversion rate
  a = 0
  b = 1
  R = matrix(data = rbind(c(a,b,kappa,b),c(b,a,b,kappa),c(kappa,b,a,b),c(b,kappa,b,a)), ncol = 4)

  ## diagonalise the frequency vector
  freq = diag(pi);

  ## create the scaled Q matrix so that qii = ???Sum(j=0, j!=i)(4) qij
  Q <-(R%*%freq)-diag(apply(R%*%freq,1,sum));
  scaleQ <-sum(freq%*%R%*%freq);
  Q <- Q/scaleQ;

  ## Create the transition matrix P, by calculating the matrix exponential of the
  ## branch length * mu * Q.

  P <- expm(Q*mu*t);

  if(!is.null(length)){
    P <- 1 - ( (P[1,1])^length )
  }


  ##return transition matrix
  return(P)
}

