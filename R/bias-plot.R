#' Bias plotting
#'
#' This function takes the output of a pipeline and plots posterior parameter estimate biases
#'
#'
#' @param dataset Output of pipeline
#'
#' @return plots 4 x 1 boxplots and returns the ggplots
#'
#'
#' @export
#'
#'

Bias_Plot <- function(dataset){

  ## EXTRA FUNCTIONS AND PACKAGES
  require(ggplot2)
  blankPlot <- ggplot()+geom_blank(aes(1,1)) +
    cowplot::theme_nothing()

  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  windows()
  ## start
ACWUP_1 <- dataset

ACWUP_1$Var[which(ACWUP_1$Var==0)] <- "Perfect"
levels(ACWUP_1$Var) <- levels(c("Perfect",10,1,0.1))

acwup1 <- ggplot(ACWUP_1, aes(x=factor(Var,levels = c("Perfect",10,1,0.1),ordered = TRUE),
                              y=((O_TMRCA_Med-E_TMRCA)/E_TMRCA)*-1,
                              colour=(factor(Sampling,levels = c("Uniform","Random","Oldest_Bar_One"),ordered=TRUE))))
acwup1_mu <- ggplot(ACWUP_1, aes(x=factor(Var,levels = c("Perfect",10,1,0.1),ordered = TRUE),
                                 y=((O_Mu_Med-E_Mu)/E_Mu)*-1,
                                 colour=(factor(Sampling,levels = c("Uniform","Random","Oldest_Bar_One"),ordered = TRUE))))
acwup1_Ne <- ggplot(ACWUP_1, aes(x=factor(Var,levels = c("Perfect",10,1,0.1),ordered = TRUE),
                                 y=((O_ePop_Med-E_ePop)/E_ePop)*-1,
                                 colour=(factor(Sampling,levels = c("Uniform","Random","Oldest_Bar_One"),ordered = TRUE))))
acwup1_gr <- ggplot(ACWUP_1, aes(x=factor(Var,levels = c("Perfect",10,1,0.1),ordered = TRUE),
                                 y=((O_Gr_Med-E_Gr)/E_Gr)*-1,
                                 colour=(factor(Sampling,levels = c("Uniform","Random","Oldest_Bar_One"),ordered = TRUE))))

acwup2 <- acwup1 + geom_hline(yintercept = 0) +
  geom_boxplot(aes(colour=(factor(Sampling,levels = c("Uniform","Random","Oldest_Bar_One"),ordered = TRUE)))) +
  xlab("Negative binomial Dispersion Parameter") + ylab("TMRCA Bias") +
  theme_bw(base_size = 14) +
  scale_color_discrete(name="Sampling", labels=c("Uniform", "Random", "Oldest Bar One")) +
  scale_fill_grey(start = 0.2, end = 0.8, na.value = "grey50",name="Sigma") +
  ylim(c(-5,0.2))

acwup2_mu <- acwup1_mu + geom_hline(yintercept = 0) +
  geom_boxplot(aes(colour=(factor(Sampling,levels = c("Uniform","Random","Oldest_Bar_One"),ordered = TRUE)))) +
  xlab("Negative binomial Dispersion Parameter") + ylab("Clock Rate Bias") +
  theme_bw(base_size = 14) +
  scale_color_discrete(name="Sampling", labels=c("Uniform", "Random", "Oldest Bar One")) +
  scale_fill_grey(start = 0.2, end = 0.8, na.value = "grey50",name="Sigma")

acwup2_Ne <- acwup1_Ne + geom_hline(yintercept = 0) +
  geom_boxplot(aes(colour=(factor(Sampling,levels = c("Uniform","Random","Oldest_Bar_One"),ordered = TRUE))),outlier.shape=NA) +
  xlab("Negative binomial Dispersion Parameter") + ylab("Ne Bias") +
  theme_bw(base_size = 14) +
  scale_color_discrete(name="Sampling", labels=c("Uniform", "Random", "Oldest Bar One")) +
  scale_fill_grey(start = 0.2, end = 0.8, na.value = "grey50",name="Sigma") + coord_cartesian(ylim=c(-100, 100))

acwup2_gr <- acwup1_gr + geom_hline(yintercept = 0) +
  geom_boxplot(aes(colour=(factor(Sampling,levels = c("Uniform","Random","Oldest_Bar_One"),ordered = TRUE)))) +
  xlab("Negative binomial Dispersion Parameter") + ylab("Growth Rate Bias") +
  theme_bw(base_size = 14) +
  scale_color_discrete(name="Sampling", labels=c("Uniform", "Random", "Oldest Bar One")) +
  scale_fill_grey(start = 0.2, end = 0.8, na.value = "grey50",name="Sigma")

gridExtra::grid.arrange(acwup2,acwup2_mu,acwup2_Ne,acwup2_gr, top=grid::textGrob(dataset$name,gp=gpar(fontsize=20,font=3)))

return(list("TMRCA"=acwup2,"Mu"=acwup2_mu,"Ne"=acwup2_Ne,"Gr"=acwup2_gr))

}

