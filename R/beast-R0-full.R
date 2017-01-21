#' BEAST_R0_Full
#'
#' Main simulation function to investigate the effects of sampling strategy against TMRCA/Clock Rate Bias. As opposed
#' to \code{BEAST_R0_Check}, this function considers all sampling schemes.
#'
#' \code{BEAST_R0_Full} proceeds by simulating a tree in a local directory, sampling from it and placing the nexus samples
#' within the local directory, before then creating the BEAST2.3xml file in a parallel network directory, and
#' then running the BEAST analysis using the batch collation file.
#'
#' Highly prescriptive function declaration, which considers all the simulation parameters needed, the range of
#' sigma parameters, the BEAST template required, directories etc.
#'
# SIMULATION PARAMETERS
#' @param N Target infected population size.  Default = 1e5
#' @param R0 basic reproduction number. Default = 2
#' @param Tg mean generation time. Default = 2.6 days
#' @param size.reps Vector of different size parameters in negative binomial offspring distribution. Default=c(0,10,1,0.1)
#' @param scale rgamma scale parameter. Default = 1
#' @param real.time Boolean dictating. whether an extra generation is simulated producing extra for real.time consideration.Default = FALSE
#' @param t.cont Boolean dictating whether gamma process is used for generation time distribution. Default = TRUE
#' @param size.reps - Vector of offspring distribution parameters. Default=c(0,10,1,0.1)
# SEQUENCING PARAMETERS
#' @param mu Clock Rate. Default=1e-4
#' @param seq.length Sequence length. Will be overwritten accordingly if root.path given. Default = 1701
#' @param model Nucleic subsitution model. Can be one of "JC69", "K80", "HKY85" (Default)
#' @param kappa Tv/Ts ration. Default = 2
#' @param root.path File path of saveRDS(root). Default = NULL
#' @param evolution.model String variable detailing evolutionary model. Can be any of "generational", "within.host" and "parallel"
#' @param outgoup Boolean concerning whether outgroup is simulated. Default = FALSE
# SAMPLING/REPETITION PARAMETERS
#' @param sample.reps number of tree subsamples under each scheme. Default = 10
#' @param tree.reps number of unique tree simulations conducted under each scheme. Default = 20
#' @param sigma.reps variables total. Default = c(1,5,9,13,18)
#' @param n.seq.samples number of sequences to be sampled. Default = 50
#' @param oldestboolean boolean if oldest sequence not root is included. Default=TRUE
#' @param missed.generations Number of intiial generations not considered for sampling. Default = 0
#' @param template.xml path for template xml. Default is a strict molecular clock, exponential coalescent, HKY85+G4 xml
#' found within the data folder called test.xml
#' @param dir Local Directory where simulations tree and nexus files are stored. Defaults specific to author.
#' @param network.path Network Directory which a cluster can see to access the BEAST xml files. Defaults specific to author.
#'
#' @return Returns a data frame of collated statstics
#'
#' @export

BEAST_R0_Full <- function(N=1e5,
                          R0=2,
                          Tg=2,
                          t.cont=TRUE,
                          real.time=FALSE,
                          scale=1,
                          size.reps=c(0,10,1,0.1),
                          ## sequencing parameters
                          mu=1.57e-05,
                          root.path=NULL,
                          seq.length=1600,
                          model="HKY85",
                          evolution.model="within.host",
                          kappa=2,
                          outgroup=FALSE,
                          ## sampling parameters
                          missed.generations=0,
                          n.seq.samples=50,
                          sample.reps=10,
                          tree.reps=20,
                          sigma.reps=c(1,5,9,13,18),
                          oldestboolean=FALSE,
                          template.xml = "C:/Users/Oliver/GoogleDrive/AcademicWork/Imperial/git/evoflution/inst/extdata/template.xml",
                          dir="C:/Users/Oliver/GoogleDrive/AcademicWork/Imperial/git/evoflution/Results/Within_Host_No_Outgroup",
                          network.path="\\\\fi--san02.dide.ic.ac.uk\\homes\\olw13\\Data\\"){

  ## HANDLE VARIABLES ##
  ##########################################################
  ptm <- proc.time()
  ## catch directory ending
  if(tail(unlist(strsplit(dir,"")),1)!="/") dir <- paste(dir,"/",sep="")

  ## CREATE RESULTS STORAGE ##
  ##########################################################
  dir.create(path=dir)
  network <- paste(network.path,tail(strsplit(dir,"/")[[1]],1),"\\",sep="")
  dir.create(path=network)
  network.dir <- paste(network.path,tail(strsplit(dir,"/")[[1]],1),"\\","mu_",mu,"_R0_",R0,"_Tg_",Tg,"\\",sep="")
  dir.create(path=network.dir)
  param.dir <- paste(dir,"mu_",mu,"_R0_",R0,"_Tg_",Tg,"/",sep="")
  dir.create(path=param.dir)

  ## MAIN LOOP ##
  ##########################################################

  ## set up for loop for changes in variation

  for (a in 1:length(size.reps)){

    ## create progress bar
    pb <- winProgressBar(title=paste(a,": Tree Reps",sep=""), label="0% done", min=0, max=tree.reps, initial=0)


    size <- size.reps[a]
    if(a==1){
      var.dir <- paste(param.dir,"Size",size,"/",sep="")
      network.var.dir <- paste(network.dir,"Size",size,"\\",sep="")
    } else {
      var.dir <- paste(param.dir,"Size",size,"/",sep="")
      network.var.dir <- paste(network.dir,"Size",size,"\\",sep="")
    }


    dir.create(var.dir)
    dir.create(network.var.dir)

    # R0 broken, size = 10, size = 0.1

    for (b in 1:tree.reps){

      sub.var.dir <- paste(var.dir,"rep",b,"/",sep="")
      network.sub.var.dir <- paste(network.var.dir,"rep",b,"\\",sep="")
      dir.create(sub.var.dir)
      dir.create(network.sub.var.dir)

      ## for loop for repetitions
      # simulate 10 trees


      ## simulate tree
      if(a==1){

        sim.out <- Infection_sim(N=N, R0=R0, Tg=Tg, fixed.offspring=TRUE, scale=scale,
                                 size=size, real.time=FALSE, continous.time=t.cont)
        seq.out <- Sequencing(sim.out, seq.length=seq.length, mu=1.57e-05, Tg=2.6, model="HKY85",
                              kappa=2, root.path=NULL, evolution.model="within.host", outgroup=outgroup)

      } else {

        sim.out <- Infection_sim(N=N, R0=R0, Tg=Tg, fixed.offspring=FALSE, scale=scale,
                                 size=size, real.time=FALSE, continous.time=t.cont)
        seq.out <- Sequencing(sim.out,seq.length=seq.length,mu=1.57e-05,Tg=2.6,model="HKY85",kappa=2,
                              root.path=NULL, evolution.model="within.host",outgroup=outgroup)

      }

      pseudo.case.tree <- list("sim.out" = sim.out, "seq.out" = seq.out)
      saveRDS(pseudo.case.tree,file=paste(sub.var.dir,"pseudo_case_tree.rds",sep=""))

      ###############################################
      ## UNIFORM
      ###############################################
      # sample sample.rep times in a uniform fashion

      bs.label = "Uniform"

      New_sampling(sim.out, seq.out, file=sub.var.dir, n.seq.samples=n.seq.samples, n.bs = sample.reps, outgroup=outgroup,
                   sigma=1, random=FALSE, oldestboolean=oldestboolean, bs.label=bs.label, t.cont=t.cont, proportional=FALSE,
                   missed.generations = missed.generations)


      ## create beast xml files

      for(j in 1:sample.reps){

        ## create BEAST xml in the parent data directory
        Nex_2_BEASTxml(template.xml = template.xml,
                       nex = paste(sub.var.dir,"Uniform",j,".nex",sep=""),file = paste(network.sub.var.dir,bs.label,j,".beast.xml",sep=""),
                       rep.str = paste("Uniform",j,sep=""),cont = t.cont)
      }

      ## create necessary submit bat in the network share
      fileConn<-file(paste(network.sub.var.dir,"treerep_",b,bs.label,".cluster.bat",sep=""))
      writeLines(paste("job submit /scheduler:fi--dideclusthn.dide.ic.ac.uk /stdout:",
                       network.sub.var.dir,bs.label,".out",b,".txt /stderr:",network.sub.var.dir,bs.label,".err",b,".txt ",
                       "/numcores:1-1 /jobtemplate:GeneralNodes ", paste(network.sub.var.dir,"treerep_",b,bs.label,".beast.bat",sep=""),sep=""),
                 fileConn)
      close(fileConn)

      ## create necessary beast bat in the network share
      fileConn2<-file(paste(network.sub.var.dir,"treerep_",b,bs.label,".beast.bat",sep=""))
      writeLines(paste("call setjava64 \n",
                       "FOR /L %%i IN (1 1 ",sample.reps,") do ( \n",
                       "java -Djava.library.path=\"\\\\fi--san02.dide.ic.ac.uk\\homes\\olw13\\libhmsbeagle\" -jar",
                       " \"\\\\fi--san02.dide.ic.ac.uk\\homes\\olw13\\BEAST\\lib\\beast.jar\"",
                       " -working -overwrite -beagle -beagle_CPU -beagle_SSE -beagle_instances 8 -threads 8 ",
                       paste("\"",network.sub.var.dir,bs.label,"%%i%.beast.xml\"\n",
                             ")",sep=""),
                       sep=""),fileConn2)
      close(fileConn2)

      ## call cluster BEAST
      system(paste(network.sub.var.dir,"treerep_",b,bs.label,".cluster.bat",sep=""))

      ##################################
      ## OLDEST BAR ONE
      ##################################
      # sample sample.rep times in an oldest bar.one fashion

      oldest.sigma <- length(seq.out$individuals.sequences)/n.seq.samples

      bs.label <- "Oldest_Bar_One"
      New_sampling(sim.out, seq.out, file=sub.var.dir, n.seq.samples=n.seq.samples, n.bs = sample.reps, outgroup=outgroup,
                   sigma=oldest.sigma, random=TRUE, oldestboolean=oldestboolean, bs.label=bs.label, t.cont=t.cont, proportional=FALSE,
                   missed.generations = missed.generations)


      ## create beast xml files

      for(j in 1:sample.reps){

        ## create BEAST xml in the parent data directory
        Nex_2_BEASTxml(template.xml = template.xml,
                       nex = paste(sub.var.dir,bs.label,j,".nex",sep=""),file = paste(network.sub.var.dir,bs.label,j,".beast.xml",sep=""),
                       rep.str = paste(bs.label,j,sep=""),cont = t.cont)
      }

      ## create necessary submit bat in the network share
      fileConn<-file(paste(network.sub.var.dir,"treerep_",b,bs.label,".cluster.bat",sep=""))
      writeLines(paste("job submit /scheduler:fi--dideclusthn.dide.ic.ac.uk /stdout:",
                       network.sub.var.dir,bs.label,".out",b,".txt /stderr:",network.sub.var.dir,bs.label,".err",b,".txt ",
                       "/numcores:1-1 /jobtemplate:GeneralNodes ", paste(network.sub.var.dir,"treerep_",b,bs.label,".beast.bat",sep=""),sep=""),
                 fileConn)
      close(fileConn)

      ## create necessary beast bat in the network share
      fileConn2<-file(paste(network.sub.var.dir,"treerep_",b,bs.label,".beast.bat",sep=""))
      writeLines(paste("call setjava64 \n",
                       "FOR /L %%i IN (1 1 ",sample.reps,") do ( \n",
                       "java -Djava.library.path=\"\\\\fi--san02.dide.ic.ac.uk\\homes\\olw13\\libhmsbeagle\" -jar",
                       " \"\\\\fi--san02.dide.ic.ac.uk\\homes\\olw13\\BEAST\\lib\\beast.jar\"",
                       " -working -overwrite -beagle -beagle_CPU -beagle_SSE -beagle_instances 8 -threads 8 ",
                       paste("\"",network.sub.var.dir,bs.label,"%%i%.beast.xml\"\n",
                             ")",sep=""),
                       sep=""),fileConn2)
      close(fileConn2)

      ## call cluster BEAST
      system(paste(network.sub.var.dir,"treerep_",b,bs.label,".cluster.bat",sep=""))

      ##################################
      ## RANDOM
      ##################################
      # sample sample.rep times in a random fashion


      bs.label <- "Random"

      New_sampling(sim.out, seq.out, file=sub.var.dir, n.seq.samples=n.seq.samples, n.bs = sample.reps, outgroup=outgroup,
                   sigma=1, random=TRUE, oldestboolean=oldestboolean, bs.label=bs.label, t.cont=t.cont, proportional=FALSE,
                   missed.generations = missed.generations)

      ## create beast xml files

      for(j in 1:sample.reps){

        ## create BEAST xml in the parent data directory
        Nex_2_BEASTxml(template.xml = template.xml,
                       nex = paste(sub.var.dir,bs.label,j,".nex",sep=""),file = paste(network.sub.var.dir,bs.label,j,".beast.xml",sep=""),
                       rep.str = paste(bs.label,j,sep=""),cont = t.cont)
      }

      ## create necessary submit bat in the network share
      fileConn<-file(paste(network.sub.var.dir,"treerep_",b,bs.label,".cluster.bat",sep=""))
      writeLines(paste("job submit /scheduler:fi--dideclusthn.dide.ic.ac.uk /stdout:",
                       network.sub.var.dir,bs.label,".out",b,".txt /stderr:",network.sub.var.dir,bs.label,".err",b,".txt ",
                       "/numcores:1-1 /jobtemplate:GeneralNodes ", paste(network.sub.var.dir,"treerep_",b,bs.label,".beast.bat",sep=""),sep=""),
                 fileConn)
      close(fileConn)

      ## create necessary beast bat in the network share
      fileConn2<-file(paste(network.sub.var.dir,"treerep_",b,bs.label,".beast.bat",sep=""))
      writeLines(paste("call setjava64 \n",
                       "FOR /L %%i IN (1 1 ",sample.reps,") do ( \n",
                       "java -Djava.library.path=\"\\\\fi--san02.dide.ic.ac.uk\\homes\\olw13\\libhmsbeagle\" -jar",
                       " \"\\\\fi--san02.dide.ic.ac.uk\\homes\\olw13\\BEAST\\lib\\beast.jar\"",
                       " -working -overwrite -beagle -beagle_CPU -beagle_SSE -beagle_instances 8 -threads 8 ",
                       paste("\"",network.sub.var.dir,bs.label,"%%i%.beast.xml\"\n",
                             ")",sep=""),
                       sep=""),fileConn2)
      close(fileConn2)

      ## call cluster BEAST
      system(paste(network.sub.var.dir,"treerep_",b,bs.label,".cluster.bat",sep=""))

    }


    ###############################################
    ## update progress bar
    info <- sprintf("%d%% done", round(((((b)/tree.reps)*100))))
    setWinProgressBar(pb, round((((b)/tree.reps)*tree.reps)), label=info)




  }

  ## print times
  close(pb)
  print(paste(tree.reps,"trees analysed, concerning ",sample.reps," samples",((proc.time() - ptm)[3]),"secs"))
}
