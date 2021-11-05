Prob<-function(m, data, coltime,colcensor, colstrat, t0=NULL,...){
  data$strat<-data[, colstrat]
  iter<-dim(data)[1]; colnames(data)[coltime]<-'t'; colnames(data)[colcensor]<-'d';dat<-data[order(data$t),];
  strat<-data$strat;k<-1:iter;t<-1:iter;
  if(is.null(t0)){t0=rep(min(data$t), iter)}
  output = apply(as.matrix(1:iter,iter,1), 1,function(x) Pkt_f(x,dat,iter, t0, m,  colstrat))
  Pkt<-matrix(output,iter,iter); p<-numeric(iter)
  daux<-dat$d
  for(j  in 1:iter){
    if(dat$d[j]==1){  p[j]<-1}
    if(dat$d[j]==0){  taux<- min(data$t[t0[j]<= data$t]);bj<-(1:iter)[data$t==taux]
    pl<-((1-Pkt[j,])[bj:(j)]);  pl[j]<-1 ; p[j]<- 1-prod(pl)}}
  list(p=p )}


Pkt_f<-function(t, dat,iter, t0=NULL,m, colstrat){
  dat$strat<-dat[, colstrat]
  if(is.null(t0)){t0=rep(min(data$t), iter)} # If t0 is null, then it assumes all subjects enter the study at the same time
  p_t<-numeric(iter)
  if(dat$d[t]==1) {
    c_t<- dat$strat[t]; c_k<- dat$strat;  L<- length(unique(dat$strat));  n_t<-numeric(L);  len<- numeric(L)
    for(r in 1:L)   {n_t[r]<-sum((rep(1, iter-t+1))[dat$strat[t:iter]==r &  t0[t:iter]<=dat$t[t]  ]  )
    if( n_t[r]>0) len[r]<-0    else  len[r]<-1    } ;   m_t<-m

    if( n_t[c_t]>=m_t[c_t]) {p_t[c_t==c_k ]<- (m_t[c_t]-1)/(n_t[c_t]-1)}
    if(n_t[c_t]>=m_t[c_t]){ p_t[c_t==c_k]<- (m_t[c_t]-1)/(n_t[c_t]-1) }

    p_t[c_t==c_k & n_t[c_k]<m_t[c_k]& n_t[c_k]>1]<-1
    p_t[c_t==c_k & n_t[c_k]<=1]<-0

    a<-as.numeric(c_t!=c_k & n_t[c_k]>=m_t[c_k])
    b<-c_k[a==1]
    p_t[c_t!=c_k & n_t[c_k]>=m_t[c_k]]<- (m_t[b])/(n_t[b])
    p_t[c_t!=c_k & n_t[c_k]<m_t[c_k]& n_t[c_k]>=1]<-1
    p_t[c_t!=c_k & n_t[c_k]==0]<-0
    p_t[1:t]<- 0
  }
  p_t}

cmsample1<- function(data, m, colstrat)
{
  data$strat<-data[, colstrat]
  iter<-nrow(data)
  mat<-numeric()
  for(j  in (1:iter)[data$cens==1]){
    wj<-numeric()
    st_j<-as.vector(table(data$strat[j:iter]));
    L<-length(st_j);
    stat<-1:L
    c_j<-data$strat[j];
    dl<-numeric();
    m_j<-m ;
    m_ja<-m_j;
    m_j[c_j]<- m_j[c_j]-1;
    w<- st_j/m_ja;
    sj<-numeric(); staj<-numeric()
    for(l in 1:L){
      saux<-(((j+1):iter)[data$strat[(j+1):iter]==stat[l]  ]); si<- as.numeric(sample(as.character(saux), m_j[l])); sj<-c(sj,si); w_jl<-rep(st_j[l]/m_ja[l], nrow(data.frame(si)))
      staj<-c(staj,data$strat[si]); sl<-data[si,]; dl<-rbind(dl,sl); wj<-c(wj, w_jl)
    }
    mat<-rbind(mat, cbind(t(cbind(t(data[j,]), t(dl ) ) ) ,  c(st_j[c_j]/(m_ja[c_j] ), wj )))}
  mat
}

Pkt_f<-function(t, dat,iter, t0=NULL,m, colstrat){
  dat$strat<-dat[, colstrat]
  if(is.null(t0)){t0=rep(min(data$t), iter)} # If t0 is null, then it assumes all subjects enter the study at the same time
  p_t<-numeric(iter)
  if(dat$d[t]==1) {
    c_t<- dat$strat[t]; c_k<- dat$strat;  L<- length(unique(dat$strat));  n_t<-numeric(L);  len<- numeric(L)
    for(r in 1:L)   {n_t[r]<-sum((rep(1, iter-t+1))[dat$strat[t:iter]==r &  t0[t:iter]<=dat$t[t]  ]  )
    if( n_t[r]>0) len[r]<-0    else  len[r]<-1    } ;   m_t<-m

    if( n_t[c_t]>=m_t[c_t]) {p_t[c_t==c_k ]<- (m_t[c_t]-1)/(n_t[c_t]-1)}
    if(n_t[c_t]>=m_t[c_t]){ p_t[c_t==c_k]<- (m_t[c_t]-1)/(n_t[c_t]-1) }

    p_t[c_t==c_k & n_t[c_k]<m_t[c_k]& n_t[c_k]>1]<-1
    p_t[c_t==c_k & n_t[c_k]<=1]<-0

    a<-as.numeric(c_t!=c_k & n_t[c_k]>=m_t[c_k])
    b<-c_k[a==1]
    p_t[c_t!=c_k & n_t[c_k]>=m_t[c_k]]<- (m_t[b])/(n_t[b])
    p_t[c_t!=c_k & n_t[c_k]<m_t[c_k]& n_t[c_k]>=1]<-1
    p_t[c_t!=c_k & n_t[c_k]==0]<-0
    p_t[1:t]<- 0
  }
  p_t}

