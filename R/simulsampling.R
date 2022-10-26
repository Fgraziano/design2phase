simulsampling<-function(betaBM, pBM, pstrata,betastrata,acc.aux, cens,tau,N,n,lambda, k, B,seed, design=c("PPS","NCC","CM"))
  {

  i<-1
  #Full_Cohort<-NULL
  Rsimple<-NULL
  CC_ev<- NULL
  CC_stra_elfin<-NULL
  CC_stra_surr<-NULL
  PPS_EVe<-NULL
  PPS_EVe_elf<-NULL
  PPS_EVsurr<-NULL
  NCC_1<-NULL
  CM<-NULL
  censor<-NULL
  eventfull<-NULL

  if(missing(n) ) {
    return("Missing n value!")
  }
  if(missing(B) ) {
    return("Missing n of simulation (B)!")
  }
  if(missing(betaBM) ) {
    return("Missing beta biomarker!")
  }
  if(missing(pBM) ) {
    return("Missing prevalence of biomarker!")
  }
  if(missing(tau) ) {
    return("Missing tau!")
  }

  X<-dataX(N=N, pBM=pBM, pstrata=pstrata,surrogate=acc.aux)

  for (i in 1:B) {

    dati<-genSAMPLEexpWeiC(N, lambda, k,betaBM, betastrata,X, cens, tau)
    eventfull[[i]]<-round(table(dati$censor)[2],1)

    ev=0
    idsampled<-sample(dati$id, n, replace = FALSE, prob = NULL)
    campione<-dati[dati$id %in% idsampled,]
    ev=sum(campione$censor==1)
    dati$incl<-ifelse(dati$id %in% idsampled,1,0)
    DATse<-dati[dati$incl==1,]
    Fit<-suppressWarnings(survival::coxph(Surv(time,censor)~x1, DATse, robust=TRUE))

   invisible(capture.output(sample.est<- summary(Fit)$coefficients[,c(4,6)]))

    py<-round(table(DATse$censor)[2],1)
    VAR_srs<-sample.est[1]^2
    DEFF<- 1

    Rsimple<-rbind(Rsimple,c(sample=i,length(DATse$x1),py, sample.est,DEFF) )
    colnames(Rsimple)<-c("sample", "n", "n.event",  "se.coef", "pvalue","DEFF")

    ###CASE CONTROL

    cases<-NULL;cntl<-NULL
    if(length(dati$id[dati$censor==1])>=n/2){
      cases<-sample(dati$id[dati$censor==1], n/2, replace = FALSE, prob = NULL)
      cntl<-sample(dati$id[!dati$censor==1],n-length(cases), replace = FALSE, prob = NULL)
    }
    if (length(dati$id[dati$censor==1])<n/2) cases<-dati$id[dati$censor==1]
    cntl<-sample(dati$id[!dati$censor==1],n/2, replace = FALSE, prob = NULL)

    idsampled<-c(cases,cntl);

    campione<-dati[dati$id %in% idsampled,]
    dati$incl<-ifelse(dati$id %in% idsampled,1,0);

    DATse<-dati[dati$incl==1,]
    py<-round(table(DATse$censor)[2],1)

    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~censor),data=dati)

    Fit<-suppressWarnings(svycoxph(Surv(time,censor)~x1 , des))

    invisible(capture.output(sample.est<- summary(Fit)$coefficients[,c("robust se","Pr(>|z|)")]))

    VAR_d<-sample.est[1]^2
    DEFF<- VAR_srs/VAR_d

     CC_ev<-rbind(CC_ev,c(sample=i,length(DATse$x1),py, sample.est,DEFF))
    colnames(CC_ev)<-c("sample", "n", "n.event","se.coef", "pvalue","DEFF")

    if("x2" %in% colnames(dati))
    {

      cases0<-NULL;cntl0<-NULL; cases1<-NULL;cntl1<-NULL #length(dati$id[dati$censor==1])

      if(length(dati$id[dati$censor==1 & dati$x2==0])>=n/4){
        cases0<-sample(dati$id[dati$censor==1 & dati$x2==0], n/4 ,replace = FALSE, prob = NULL) }
      if(length(dati$id[dati$censor==1 & dati$x2==0])<n/4) cases0<-dati$id[dati$censor==1 & dati$x2==0]

      cntl0<-sample(dati$id[!dati$censor==1 & dati$x2==0], n/4 ,replace = FALSE , prob = NULL)

      if(length(dati$id[dati$censor==1 & dati$x2==1])>=n/4){
        cases1<-sample(dati$id[dati$censor==1 & dati$x2==1], n/4,replace = FALSE , prob = NULL)
      }
      if(length(dati$id[dati$censor==1 & dati$x2==1])<n/4) cases1<-dati$id[dati$censor==1 & dati$x2==1]

      cntl1<-sample(dati$id[!dati$censor==1 & dati$x2==1], n/4 , replace = FALSE, prob = NULL)

      idsampled<-c(cases0,cases1,cntl0,cntl1); length(idsampled)

      campione<-dati[dati$id %in% idsampled,]
      dati$incl<-ifelse(dati$id %in% idsampled,1,0);

      DATse<-dati[dati$incl==1,]

      py<-round(table(DATse$censor)[2],1)

      des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x2, censor)),data=dati)

      Fit<-suppressWarnings(svycoxph(Surv(time,censor)~x1 , des))

      invisible(capture.output(sample.est<- summary(Fit)$coefficients[,c("robust se","Pr(>|z|)")]))
      VAR_d<-sample.est[1]^2
      DEFF<- VAR_srs/VAR_d

      CC_stra_elfin<-rbind(CC_stra_elfin,c(sample=i,length(DATse$x1),py,  sample.est,DEFF))
      colnames(CC_stra_elfin)<-c("sample", "n", "n.event", "se.coef", "pvalue","DEFF")
    }
    else {
      CC_stra_elfin<-rbind(CC_stra_elfin,c(sample=i,NA,NA,NA,NA,NA))
      colnames(CC_stra_elfin)<-c("sample", "n", "n.event", "se.coef", "pvalue","DEFF")
    }


    if("x3" %in% colnames(dati))
    {
      cases0<-NULL;cntl0<-NULL; cases1<-NULL;cntl1<-NULL

      if(length(dati$id[dati$censor==1 & dati$x3==0])>=n/4){
        cases0<-sample(dati$id[dati$censor==1 & dati$x3==0], n/4 ,replace = FALSE, prob = NULL) }
      if(length(dati$id[dati$censor==1 & dati$x3==0])<n/4) cases0<-dati$id[dati$censor==1 & dati$x3==0]

      cntl0<-sample(dati$id[!dati$censor==1 & dati$x3==0], n/4 ,replace = FALSE , prob = NULL)

      if(length(dati$id[dati$censor==1 & dati$x3==1])>=n/4){
        cases1<-sample(dati$id[dati$censor==1 & dati$x3==1], n/4,replace = FALSE , prob = NULL)}
      if(length(dati$id[dati$censor==1 & dati$x3==1])<n/4) cases1<-dati$id[dati$censor==1 & dati$x3==1]

      if(length(dati$id[!dati$censor==1 & dati$x3==1])>=n/4){
        cntl1<-sample(dati$id[!dati$censor==1 & dati$x3==1], n/4 , replace = FALSE, prob = NULL)}
      if(length(dati$id[!dati$censor==1 & dati$x3==1])<n/4) cntl1<-dati$id[!dati$censor==1 & dati$x3==1]

      idsampled<-c(cases0,cases1,cntl0,cntl1)

      campione<-dati[dati$id %in% idsampled,]
      dati$incl<-ifelse(dati$id %in% idsampled,1,0)

      DATse<-dati[dati$incl==1,]

      py<-round(table(DATse$censor)[2],1)

      des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x3, censor)),data=dati)

      Fit<-suppressWarnings(svycoxph(Surv(time,censor)~x1, des))

      invisible(capture.output(sample.est<- summary(Fit)$coefficients[,c("robust se","Pr(>|z|)")]))
      VAR_d<-sample.est[1]^2
      DEFF<- VAR_srs/VAR_d

      CC_stra_surr<-rbind(CC_stra_surr,c(sample=i,length(DATse$x1),py, sample.est,DEFF))
      colnames(CC_stra_surr)<-c("sample", "n", "n.event",  "se.coef", "pvalue", "DEFF")
    }
    else {
      CC_stra_surr<-rbind(CC_stra_surr,c(sample=i,NA,NA,NA,NA,NA))
      colnames(CC_stra_surr)<-c("sample", "n", "n.event",  "se.coef", "pvalue", "DEFF")
    }


    if("PPS" %in% design )
    {
    T1<-table(dati$censor)
    cases<-sample(dati$id[dati$censor==1], size= T1[2]*n/N, replace = FALSE)
    cntl<-sample(dati$id[!dati$censor==1], size= T1[1]*n/N, replace = FALSE)
    idsampled<-c(cases,cntl)

    campione<-dati[dati$id %in% idsampled,]
    dati$incl<-ifelse(dati$id %in% idsampled,1,0)

    des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~censor), data=dati)

    #Fit<-svycoxph(Surv(time,cens)~x1, des)
    #sample.est<- c(summary(Fit)$coefficients[,c(1:2,4:6)],summary(Fit)$conf.int[1,3:4])

    Fit<-suppressWarnings(svycoxph(Surv(time,censor)~x1, des))

    invisible(capture.output(sample.est<- summary(Fit)$coefficients[,c("se(coef)","Pr(>|z|)")]))
    VAR_d<-sample.est[1]^2
    DEFF<- VAR_srs/VAR_d

    DATse<-dati[dati$incl==1,]

    py<-round(table(DATse$censor)[2],1)

    PPS_EVe<-rbind( PPS_EVe,c(sample=i,length(DATse$x1),py, sample.est, DEFF) )
    colnames(PPS_EVe)<-c("sample", "n", "n.event", "se.coef", "pvalue", "DEFF")
    }
    else {
      PPS_EVe<-rbind(PPS_EVe,c(sample=i,NA,NA,NA,NA,NA))
      colnames(PPS_EVe)<-c("sample", "n", "n.event", "se.coef","pvalue", "DEFF")
    }

    if("x2" %in% colnames(dati) & "PPS" %in% design )
    {
      cases0<-NULL; cases1<-NULL; cntl0<- NULL ; cntl1<-NULL

      T1<-table(dati$censor,dati$x2); #str(T1)
      cases0<-sample(dati$id[dati$censor==1 & dati$x2==0], T1[2,1]*n/N ,replace = FALSE, prob = NULL)
      cntl0<-sample(dati$id[!dati$censor==1 & dati$x2==0], T1[1,1]*n/N ,replace = FALSE , prob = NULL)
      cases1<-sample(dati$id[dati$censor==1 & dati$x2==1], T1[2,2]*n/N,replace = FALSE , prob = NULL)
      cntl1<-sample(dati$id[!dati$censor==1 & dati$x2==1],T1[1,2]*n/N , replace = FALSE, prob = NULL)

      idsampled<-c(cases0,cases1,cntl0,cntl1)

      dati$incl<-ifelse(dati$id %in% idsampled,1,0);

      des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x2, censor)), data=dati)

      Fit<-suppressWarnings(svycoxph(Surv(time,censor)~x1 , des))

      invisible(capture.output(sample.est<- summary(Fit)$coefficients[,c("se(coef)","Pr(>|z|)")]))
      VAR_d<-sample.est[1]^2
      DEFF<- VAR_srs/VAR_d

      DATse<-dati[dati$incl==1,]

      py<-round(table(DATse$censor)[2],1)

      PPS_EVe_elf<-rbind(PPS_EVe_elf,c(sample=i,length(DATse$x1),py,  sample.est, DEFF))
      colnames(PPS_EVe_elf)<-c("sample", "n", "n.event", "se.coef","pvalue", "DEFF")
    }
    else {
      PPS_EVe_elf<-rbind(PPS_EVe_elf,c(sample=i,NA,NA,NA,NA,NA))
      colnames(PPS_EVe_elf)<-c("sample", "n", "n.event", "se.coef","pvalue", "DEFF")
    }
    if("x3" %in% colnames(dati) & "PPS" %in% design)
    {
      T1<-table(dati$censor,dati$x3)
      cases0<-sample(dati$id[dati$censor==1 & dati$x3==0], T1[2,1]*n/N ,replace = FALSE, prob = NULL)
      cntl0<-sample(dati$id[!dati$censor==1 & dati$x3==0], T1[1,1]*n/N,replace = FALSE , prob = NULL)
      cases1<-sample(dati$id[dati$censor==1 & dati$x3==1], T1[2,2]*n/N,replace = FALSE , prob = NULL)
      cntl1<-sample(dati$id[!dati$censor==1 & dati$x3==1],T1[1,2]*n/N, replace = FALSE, prob = NULL)

      idsampled<-c(cases0,cases1,cntl0,cntl1)

      dati$incl<-ifelse(dati$id %in% idsampled,1,0)

      des<-twophase(id=list(~id,~id),subset=~incl==1,strata=list(NULL,~interaction(x3, censor)), data=dati)

      Fit<-suppressWarnings(svycoxph(Surv(time,censor)~x1 , des))

      invisible(capture.output(sample.est<- summary(Fit)$coefficients[,c("se(coef)","Pr(>|z|)")]))
      VAR_d<-sample.est[1]^2
      DEFF<- VAR_srs/VAR_d

      DATse<-dati[dati$incl==1,]
      py<-round(table(DATse$censor)[2],1)
      PPS_EVsurr<-rbind(PPS_EVsurr,c(sample=i,length(DATse$x1),py, sample.est, DEFF) )
      colnames( PPS_EVsurr)<-c("sample", "n", "n.event", "se.coef", "pvalue", "DEFF")
    }
    else {
      PPS_EVsurr<-rbind(PPS_EVsurr,c(sample=i,NA,NA,NA,NA,NA))
      colnames(PPS_EVsurr)<-c("sample", "n", "n.event", "se.coef", "pvalue", "DEFF")
    }

    if( "NCC" %in% design)
    {
    campione<-suppressWarnings(ccwc(entry=0,exit=time,fail=censor,controls=1,data=dati,include=id,silent=TRUE))
    if((length(campione$Set)/2)>=n/2){
      SetId<-sample(1:((length(campione$Set))/2), n/2, replace = FALSE, prob = NULL)
    }
    if ((length(campione$Set)/2)<n/2) SetId<-unique(campione$Set)

    IDCASE<-campione$id[campione$Set %in% SetId & campione$Fail=="1"]
    IDContr<-campione$id[campione$Set %in% SetId & campione$Fail=="0"]
    dati$sampleidstat<-ifelse(dati$id %in% IDCASE,2,ifelse(dati$id %in% IDContr,1,0))
    dati$pro<-KMprob(dati$time,dati$sampleidstat,m=1);

    DATse<-dati[dati$incl==1,]
    py<-round(table(DATse$censor)[2],1)
    pCas<-table(DATse$censor)[2]/(length(dati$id[dati$censor=="1"]));
    dati$pro<-ifelse(dati$censor=="1", pCas,dati$pro)
    dati$incl<-ifelse(dati$id %in% c(IDCASE,IDContr),1,0)

    des<-twophase(id=list(~id,~id),subset=~incl==1,probs=list(NULL,~pro),data=dati)

    Fit<-suppressWarnings(svycoxph(Surv(time,censor)~x1 , des))

    invisible(capture.output(sample.est<- summary(Fit)$coefficients[,c("robust se","Pr(>|z|)")]))
    VAR_d<-sample.est[1]^2
    DEFF<- VAR_srs/VAR_d

    NCC_1<-rbind(NCC_1,c(sample=i,round(length(DATse$x1),1),py,  sample.est,DEFF))
    colnames(NCC_1)<-c("sample", "n", "n.event", "se.coef",  "pvalue","DEFF")}
    else {
      NCC_1<-rbind( NCC_1,c(sample=i,NA,NA,NA,NA,NA))
      colnames(NCC_1)<-c("sample", "n", "n.event", "se.coef", "pvalue", "DEFF")
    }

    if("x3" %in% colnames(dati) & "CM" %in% design)
    {
      dati<-dati[order(dati$time),]
      dati$surr<-dati$x3+1
      colstrat<-which( colnames(dati)=="surr")
      coltime<-which( colnames(dati)=="time")
      colcensor<-which( colnames(dati)=="censor")

     ccount<-cmsample1(data=dati, m=c(1,1), colstrat)

      rownames(ccount)<-1:nrow(ccount);
      ccount<-data.frame(ccount);

      risize<-rep(0,sum(ccount$censor=="1"));
      id<-1:nrow(ccount)

      risize<-rep(2, nrow(ccount)/2)
      ccount$riskset<-rep((1:sum(ccount$x3)), (risize));

      if((length(ccount$riskset)/2)>=n/2){
        SetId<-sample(1:((length(ccount$riskset))/2), n/2, replace = FALSE, prob = NULL)
      }
      if ((length(ccount$riskset)/2)<n/2) SetId<-unique(ccount$riskset)

      IDsample<-ccount$id[ccount$riskset %in% SetId]; length(IDsample)
      dati$incl<-ifelse(dati$id %in% IDsample,1,0);

      dati$pro<-Prob(m=c(1,1), dati, coltime,colcensor, colstrat)$p
      pCas<-length(dati$id[dati$incl==1 & dati$censor=="1"])/(length(dati$id[dati$censor=="1"]))
      dati$pro<-ifelse(dati$censor=="0", dati$pro, pCas)

      des<-twophase(id=list(~id,~id),subset=~incl==1, probs = list (NULL, ~pro),data=dati)
      Fit<-suppressWarnings(svycoxph(Surv(time,censor)~x1 , des))

      invisible(capture.output(sample.est<- summary(Fit)$coefficients[,c("robust se","Pr(>|z|)")]))
      VAR_d<-sample.est[1]^2
      DEFF<- VAR_srs/VAR_d

      DATse<-dati[dati$incl==1,]

      py<-round(table(DATse$censor)[2],1)

      CM<-rbind(CM,c(sample=i,round(length(DATse$x1),1),py, sample.est, DEFF) )
      colnames(CM)<-c("sample", "n", "n.event",  "se.coef",  "pvalue","DEFF")
    }
    else {
      CM<-rbind(CM,c(sample=i,NA,NA,NA,NA,NA))
      colnames(CM)<-c("sample", "n", "n.event",  "se.coef",  "pvalue","DEFF")
    }
    ###################################
    print(i)
  }
  return(c(list(eventfull=eventfull),list(Rsimple=data.frame(Rsimple),CC_ev=data.frame(CC_ev),
   CC_stra_elfin=data.frame(CC_stra_elfin),CC_stra_surr=data.frame(CC_stra_surr),
    PPS_EVe=data.frame(PPS_EVe), PPS_EVe_elf= data.frame(PPS_EVe_elf), PPS_EVsurr= data.frame(PPS_EVsurr),NCC_1=data.frame(NCC_1), CM=data.frame(CM))))
  }


