
directories <- c("Within_Host_No_Outgroup","Within_Host_Outgroup","Generational_No_Outgroup","Generational_Outgroup")
sampling_schemes <- c("Uniform","Oldest_Bar_One","Random")
outgroups <- c(FALSE,TRUE,FALSE,TRUE)
evolutions <- c(rep("within.host",2), rep("generational",2))
N <- 1e5
mu=1.57e-05
size.reps <- c(0,10,1,0.1)
sample.reps=20
tree.reps=10
R0=2
Tg=2.6

for(i in 1:4){

  dir <- paste("C:/Users/Oliver/GoogleDrive/AcademicWork/Imperial/git/evoflution/Results/",directories[i],sep="")
  print(i)

  BEAST_R0_Full(N=N,
                R0=R0,
                Tg=Tg,
                t.cont=TRUE,
                real.time=FALSE,
                scale=1,
                size.reps=size.reps,
                ## sequencing parameters
                mu=mu,
                root.path="C:/Users/Oliver/GoogleDrive/AcademicWork/Imperial/git/evoflution/inst/extdata/A_California_07_2009(H1N1)_HA.rds",
                seq.length=1600,
                model="HKY85",
                evolution.model=evolutions[i],
                kappa=2,
                outgroup=outgroups[i],
                ## sampling parameters
                n.seq.samples=50,
                sample.reps=sample.reps,
                tree.reps=tree.reps,
                sigma.reps=c(1,5,9,13,18), # won't do anything yet
                missed.generations=0,
                oldestboolean=TRUE,
                template.xml = "C:/Users/Oliver/GoogleDrive/AcademicWork/Imperial/git/evoflution/inst/extdata/template.xml",
                dir=dir,
                network.path="\\\\fi--san02.dide.ic.ac.uk\\homes\\olw13\\Data\\")

}


statistics_results <- list()
length(statistics_results) <- length(directories)
names(statistics_results) <- directories


for(directory.no in 1:4){


  param.dir <- paste("mu_",mu,"_R0_",R0,"_Tg_",Tg,sep="")
  local.dir <- paste("C:/Users/Oliver/GoogleDrive/AcademicWork/Imperial/git/evoflution/Results/",directories[directory.no],"/",sep="")

  if(directory.no>2){
  network.dir <- paste("\\\\fi--san02.dide.ic.ac.uk\\homes\\olw13\\Data\\",directories[directory.no],"\\",sep="")
  } else {
    network.dir <- paste("C:/Users/Oliver/Desktop/Data\\",directories[directory.no],"\\",sep="")
  }





    for(sampling in 1:3){

      for(size in 1:4){

      full.local.dir <- paste(local.dir,param.dir,"/Size",size.reps[size],"/",sep="")
      full.network.dir <- paste(network.dir,param.dir,"\\Size",size.reps[size],"\\",sep="")

      tracer_pull <- BEAST_R0_Network_Pull_Tracer(sampling.reps=sample.reps,
                                   rep.name=sampling_schemes[sampling],
                                   outgroup=outgroups[directory.no],
                                   tree.reps=tree.reps,
                                   vars=5,
                                   var.R0=size.reps[size],
                                   sigma=1,
                                   N.v=c(rep(N,tree.reps)),
                                   mu.v=c(rep(mu,tree.reps)),
                                   Tg.v=c(rep(2.6,tree.reps)),
                                   R0.v=c(rep(2,tree.reps)),
                                   network.root.dir=full.network.dir,
                                   local.root.dir=full.local.dir)

      statistics_results[[directory.no]] <- rbind(statistics_results[[directory.no]],tracer_pull)

    }

  }

}

res <- statistics_results
for(i in 1:length(res)){
res[[i]]$name <- names(res)[i]
}

 TMRCA_Bias_Matrix <- lapply(res,Bias_Error_Calculator,10,12,20) %>%
   lapply(function(x){return(x$Bias)}) %>%
   unlist() %>%
  matrix(ncol = 4,byrow=T)

 Clock_Rate_Bias_Matrix <- lapply(res,Mu_Bias_Error_Calculator,10,12,20) %>%
   lapply(function(x){return(x$Bias)}) %>%
   unlist() %>%
   matrix(ncol = 4,byrow=T)
