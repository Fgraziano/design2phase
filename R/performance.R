
performanceDesign<-function(d, BETA,B){

  Blength<-length(d$sample)
  nm<-round(sum(d$n)/B,0)

  nevent<-round(sum(d$n.event)/B,0)

  pvalue<-d$pvalue
  d$sign<- ifelse( pvalue<0.05,1,0)
  Power<-sum(d$sign)/B

  ###DEFF
  Deff<-d$DEFF
  deff<-sum(Deff)/B
  RelEff<-d$RelEff
  RelEff<-sum(RelEff)/B

  Perform_measure<-data.frame (Blength=Blength , sample= nm, nevent=nevent,
                                 Power=Power, deff=deff)

  rownames(Perform_measure)<-"perf"
  return(Perform_measure)
}
