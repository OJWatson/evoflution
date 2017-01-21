Mu_Bias_Error_Calculator <- function(data,tree.reps,dirs,sample.reps){

  Bias_Calculator <- function(data,tree.reps,dirs,sample.reps){

    bias <- matrix(nrow = tree.reps,ncol = dirs)
    count <- 0
    for(j in 1:dirs){
      for(i in 1:tree.reps){
        bias[i,j] <- mean(data$E_Mu[((count*sample.reps)+1):((count*sample.reps)+sample.reps)]- data$O_Mu_Mean[((count*sample.reps)+1):((count*sample.reps)+sample.reps)])/ data$E_Mu[(count+1)*sample.reps]
        count <- count + 1
      }
    }

    return(bias)

  }

  Error_Calculator <- function(data,tree.reps,dirs,sample.reps){

    error <- matrix(nrow = tree.reps,ncol = dirs)
    count <- 0
    for(j in 1:dirs){
      for(i in 1:tree.reps){
        error[i,j] <- sqrt(mean((data$E_Mu[((count*sample.reps)+1):((count*sample.reps)+sample.reps)]- data$O_Mu_Mean[((count*sample.reps)+1):((count*sample.reps)+sample.reps)])^2))/data$E_Mu[(count+1)*sample.reps]
        count <- count + 1
      }
    }

    return(error)

  }

  bias <- colMeans(Bias_Calculator(data,tree.reps=tree.reps,dirs=dirs,sample.reps=sample.reps))
  error <- colMeans(Error_Calculator(data,tree.reps=tree.reps,dirs=dirs,sample.reps=sample.reps))

  res <- list("Error"=error,"Bias"=bias)

  return(res)

}
