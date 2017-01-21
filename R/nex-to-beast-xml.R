
#' Take nex and turn into BEAST 2.3 xml file
#'
#' This function takes a generated nex file, and creates a related BEAST 2.3 xml
#' that uses another BEAST xml file as a template.
#'
#' @param template.xml - the template xml already parsed as an xml object via xmlParse > xmlToList
#' @param nex - the nex file path to be used to create the new beast.xml
#' @param file - the full file path where the created beast xml is to be written to
#' @param rep.str - the id string that is to replace the id string from the template xml, usually replicate#
#' @param cont - boolean concerning whether the simulation is time explicit or not. Default = FALSE
#' @note Function has been written to be used ideally within a loop of using sample_build directly afterwards,
#' and in this instance it makes more sense to predetermine the template.xml and template.txt. Hence the
#' template.path is default null, and function written intending user to allocated the xml and txt templates.
#'
#' @return written BEAST.xml
#'
#'
#'


Nex_2_BEASTxml <- function(template.xml, nex, file, rep.str, cont = FALSE,relaxed=FALSE){

  ## HANDLE VARIABLES ##
  ##########################################################
  template.xml=XML::xmlToList(XML::xmlParse(template.xml))
  ## Change dates

  nex2 <- ape::read.nexus.data(nex)
  names <- names(nex2)
  pos <- max(sapply(strsplit(names, split = "_"), length))
  spl <- sapply(strsplit(names,"_"), `[`, pos)

  for (i in 1:length(names)){

    if(cont==FALSE){
    if (i == 1) {
      new.dates <- paste("\n                ",names[1],"=",spl[1],".0",sep="")
    } else if (i == length(names)) {
      new.dates <- paste(new.dates,",\n",names[i],"=",spl[i],".0                ",sep="")
    }
    else {
      new.dates <- paste(new.dates,",\n",names[i],"=",spl[i],".0",sep="")
    }
    } else {
      if (i == 1) {
        new.dates <- paste("\n                ",names[1],"=",spl[1],"",sep="")
      } else if (i == length(names)) {
        new.dates <- paste(new.dates,",\n",names[i],"=",spl[i],"                ",sep="")
      }
      else {
        new.dates <- paste(new.dates,",\n",names[i],"=",spl[i],"",sep="")
      }
    }
  }

  ## set up new xml and replace dates

  new.xml <- template.xml
  new.xml$run$state$tree$trait$text <- new.dates

  ## extract sequences

  new.data <- rep(new.xml$data[1],length(names))
  new.data$.attrs <- new.xml$data$.attrs
  for (i in 1:length(names)){
    new.data[i]$sequence[1] <- paste("seq",names[i],sep="_")
      new.data[i]$sequence[2] <- names[i]
      new.data[i]$sequence[4] <- paste(unlist(nex2[i]),collapse="")
  }

  new.xml$data <-  new.data
  XML::saveXML(List_To_Xml(new.xml,"beast"),file = file)

  ## read Lines and do fine changes
  txt <- readLines(file)

    ## pull str to be replaced
    str2rep <- template.xml$data$.attrs[1]

    ## Brute txt replace of replicate number
    txt <- gsub(str2rep,rep.str,txt)
    txt <- gsub("</text>","",txt)
    txt <- gsub("<text>","",txt)
    txt <- gsub("<sequence>","",txt)
    txt <- gsub("</sequence>","",txt)
    txt <- gsub(paste("<alignment>",rep.str,"</alignment>",sep=""),paste("<alignment idref=\"",rep.str,'"/>',sep=""),txt)
    txt <- gsub(paste(" <taxonset>TaxonSet.",rep.str,"</taxonset>",sep=""),paste("<taxonset idref=\"","TaxonSet.",rep.str,'"/>',sep=""),txt)
    txt <- gsub(paste("<up>clockRate.c:",rep.str,"</up>",sep=""),paste("<up idref=\"","clockRate.c:",rep.str,'"/>',sep=""),txt)
    txt <- gsub(paste("<down>Tree.t:",rep.str,"</down>",sep=""),paste("<down idref=\"","Tree.t:",rep.str,'"/>',sep=""),txt)

    non.t <- c("posterior","likelihood","prior")
    for (i in 1:length(non.t)){
      txt <- gsub(paste("<log>",non.t[i],"</log>",sep=""),paste("<log idref=\"",non.t[i],'"/>',sep=""),txt)
    }

    non.t2 <- c("treeLikelihood.","clockRate.c:","gammaShape.s:","kappa.s:","CoalescentExponential.t:","ePopSize.t:","growthRate.t:")
    for (i in 1:length(non.t2)){
      txt <- gsub(paste("<log>",non.t2[i],rep.str,"</log>",sep=""),paste("<log idref=\"",non.t2[i],rep.str,'"/>',sep=""),txt)
    }

if (relaxed == TRUE){
  txt <- gsub(paste("<up>ucldMean.c:",rep.str,"</up>",sep=""),paste("<up idref=\"ucldMean.c:",rep.str,'"/>',sep=""),txt)
  txt <- gsub(paste("<log>ucldMean.c:",rep.str,"</log>",sep=""),paste("<log idref=\"ucldMean.c:",rep.str,'"/>',sep=""),txt)
  txt <- gsub(paste("<log>ucldStdev.c:",rep.str,"</log>",sep=""),paste("<log idref=\"ucldStdev.c:",rep.str,'"/>',sep=""),txt)
}

  ## replace actual data
  ## write the finished xml
  writeLines(txt,file(file))

}


