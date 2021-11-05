Comp_design<-function(BB, BETA,B, pstrata,  acc.aux, design2p){
  perf_design<-NULL
  for (i in 1: (length(BB)))
  {
    d<-data.frame(BB[i])
    colnames(d)<-c("sample", "n", "n.event",  "se.coef","pvalue", "DEFF")
    DD<-round(performanceDesign(d,BETA,B),4)
    perf_design<- rbind(perf_design,DD)
  }
  ###NO STRATA NO AUX
  if(is.null(pstrata) & is.null(acc.aux) & "NCC" %in% design2p & "PPS" %in% design2p){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "PPS event", "NCC")}

    ###NO STRATA NO AUX
  if(is.null(pstrata) & is.null(acc.aux) & !("NCC" %in% design2p) & "PPS" %in% design2p){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "PPS event")}

  ###NO STRATA SI AUX
  if(is.null(pstrata) & !is.null(acc.aux) & !("NCC" %in% design2p) & !("PPS" %in% design2p) &!("CM" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat aux")}

  ###NO STRATA SI AUX
  if(is.null(pstrata) & !is.null(acc.aux) & !("NCC" %in% design2p) & ("PPS" %in% design2p) &!("CM" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat aux", "PPS strat event", "PPS strat aux")}

  ###NO STRATA SI AUX
  if(is.null(pstrata) & !is.null(acc.aux) & ("NCC" %in% design2p) & !("PPS" %in% design2p) &!("CM" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat aux", "NCC")}
  ###NO STRATA SI AUX
  if(is.null(pstrata) & !is.null(acc.aux) & ("NCC" %in% design2p) & ("PPS" %in% design2p) &!("CM" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat aux", "PPS strat event", "PPS strat aux", "NCC")}

  ###NO STRATA SI AUX
  if(is.null(pstrata) & !is.null(acc.aux) & ("NCC" %in% design2p) & ("PPS" %in% design2p) &("CM" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat aux", "PPS strat event", "PPS strat aux", "NCC","CM")}

  ###NO STRATA SI AUX
  if(is.null(pstrata) & !is.null(acc.aux) & !("NCC" %in% design2p) & ("PPS" %in% design2p) & ("CM" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat aux", "PPS strat event", "PPS strat aux", "CM")}

     ###NO STRATA NO AUX
  if(is.null(pstrata) & is.null(acc.aux) & "NCC" %in% design2p & !("PPS" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "NCC")}
  ###NO STRATA NO AUX
  if(is.null(pstrata) & is.null(acc.aux) & !("NCC" %in% design2p) & !("PPS" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event")}

  ##########STRATA 1 & NO SURROGATE
  if( length(pstrata)==1 & is.null(acc.aux) & ("NCC" %in% design2p) & ("PPS" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum",
                             "PPS event", "PPS event and stratum", "NCC")}


  if(length(pstrata)==1 & is.null(acc.aux) & !("NCC" %in% design2p) & ("PPS" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum",
                             "PPS event", "PPS event and stratum")}

    if(length(pstrata)==1 & is.null(acc.aux) & ("NCC" %in% design2p) & !("PPS" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum",
                             "NCC")}
      if(length(pstrata)==1 & is.null(acc.aux) & !("NCC" %in% design2p) & !("PPS" %in% design2p)){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum")}


  ######NO STRATA e surrogate

  if(is.null(pstrata) & !is.null(acc.aux) & !("NCC" %in% design2p) & !("PPS" %in% design2p) & ("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and aux","CM")}

  if(is.null(pstrata) & !is.null(acc.aux) & ("NCC" %in% design2p) & !("PPS" %in% design2p) & ("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and aux","NCC","CM")}

  ##########strata 1 e surrogate ########################

  if(length(pstrata)==1 & !is.null(acc.aux) & !("NCC" %in% design2p) & !("PPS" %in% design2p) & ("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum","CC strat event and aux","CM")}

  if(length(pstrata)==1 & !is.null(acc.aux) & !("NCC" %in% design2p) & ("PPS" %in% design2p) & !("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum","CC strat event and aux",
                             "PPS event", "PPS event and stratum","PPS event and aux")}

  if(length(pstrata)==1 & !is.null(acc.aux) & ("NCC" %in% design2p) & !("PPS" %in% design2p) & !("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum","CC strat event and aux","NCC")}

  if(length(pstrata)==1 & !is.null(acc.aux) & !("NCC" %in% design2p) & !("PPS" %in% design2p) & !("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum","CC strat event and aux")}

  if(length(pstrata)==1 & !is.null(acc.aux) & "NCC" %in% design2p & ("PPS" %in% design2p) & !("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum","CC strat event and aux",
                             "PPS event", "PPS event and stratum","PPS event and aux","NCC")}

  if(length(pstrata)==1 & !is.null(acc.aux) & !("NCC" %in% design2p) & ("PPS" %in% design2p) & ("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum","CC strat event and aux",  "PPS event", "PPS event and stratum","PPS event and aux","CM")}

  if(length(pstrata)==1 & !is.null(acc.aux) & "NCC" %in% design2p & !("PPS" %in% design2p) & ("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event",
                             "CC strat event and stratum","CC strat event and aux", "NCC","CM")}

  if(length(pstrata)==1 & !is.null(acc.aux) & "NCC" %in% design2p & ("PPS" %in% design2p) & ("CM" %in% design2p) ){
    perf_design<-na.omit(perf_design)
    rownames(perf_design)<-c("SRS","CC strat event","CC strat event and stratum","CC strat event and aux",
                             "PPS event", "PPS event and stratum","PPS event and aux","NCC","CM")}

return(list(perf_design=perf_design,names=rownames(perf_design)))
}


