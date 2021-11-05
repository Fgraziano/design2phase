#'@importFrom reshape2 melt
#'@importFrom grDevices dev.off pdf
#'@importFrom graphics abline boxplot
#'@importFrom stats interaction.plot
#'@importFrom ggplot2 ggplot aes geom_line geom_point ylab theme_bw ylim
NULL
#'Plot power based on estimate obtained from PowerIIphase function
#'@description This function gives a plot power based on estimate obtained from objects of class PowerIIphase
#'
#'@param x a PowerIIphase object from the PowerIIphase function.
#'@return plot of powers based on different sampling designs and sample size of II phase
#'@seealso \url{PowerIIphase} function for \code{x} term
#'@export

PlotPower<-function(x){
  powP=NULL
  n<-NULL; Power<-NULL; design<-NULL;
  for (i in 1:length(x$n)) {
  powP<-cbind(powP,data.frame(x$PhaseII_Performance[i])[,4])
  }
  powP<-data.frame(powP)
  colnames(powP)<-paste(x$n, sep="")
  powP$id<-x$designs
  mdata <- melt(powP, id.vars="id")
  colnames(mdata)<-c("design","n","Power")

  ggplot(mdata, aes(x=n, y=100*Power,
                    group=design, color=design)) +
    #ylim(0,100)+
    geom_line() + geom_point() +ylab("Power (%)")+
    theme_bw()
  }

