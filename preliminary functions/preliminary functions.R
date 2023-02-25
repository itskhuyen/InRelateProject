#function for getting K from indinvq file
getSubpop <- function(indivqfile){
  subpop <- indivqfile[,6:ncol(indivqfile)]
  return(subpop)
}

#function for reading in indivq
readIndivq <- function(indivqfile){
  indivq <- getSubpop(indivqfile);
  return(indivq)
}

#function for reading in pkla
readPkla <- function(pklafile){
  pkla <- pklafile[,3:ncol(pklafile)];
  return(pkla)
}

#function for reading in str
readStr <- function(strfile){
  str <- strfile;
  return(str)
}

#get unique alleles
getAlleles <- function(pklafile){
  groupAlleles <- aggregate(x = pklafile$V2, by = list(pklafile$V1), FUN = function(x) length(unique(x)));
  return(groupAlleles) #get the unique alleles within each loci
}

#get loci
getLoci <- function(strfile){
  loci <- strfile[,3:ncol(strfile)];
  return(loci)
}

#get individual
getIndiv <- function(strfile){
  individual <- as.numeric(unique(strfile$V1));
  return(individual)
}

#compute the number of pairs
computenumpairs <-function(indivqfile){
  numpairs <-getIndiv(indivqfile*(indivqfile-1)/2);
  return(numpairs);
}

#get IBS mode
getIBS <- function(indiv1a, indiv1b, indiv2a, indiv2b){
  for(l in 3:ncol(hs)) #ncol(hs) gives total number of columns in datafram +2, does this work or need another formula	
    # Indentify IBS mode.
    if(indiv1a[l]==indiv1b[l] && indiv1b[l]==indiv2a[l] && indiv2a[l] ==indiv2b[l])
      IBS<-1
    if(indiv1a[l]==indiv1b[l] && indiv2a[l]==indiv2b[l] && indiv1a[l]!=indiv2a[l])
      IBS<-2
    if((indiv1a[l]==indiv1b[l] && indiv1a[l]==indiv2a[l] && indiv2a[l]!=indiv2b[l]) || (indiv1a[l]==indiv1b[l] && indiv1a[l]==indiv2b[l] && indiv1a[l]!=indiv2a[l]))
      IBS<-3
    if((indiv1a[l]==indiv1b[l] && indiv2a[l]!=indiv2b[l] && indiv1a[l]!=indiv2a[l] && indiv1a[l]!=indiv2b[l]))
      IBS<-4
    if((indiv1a[l]!=indiv1b[l] && indiv1a[l]==indiv2a[l] && indiv2a[l]==indiv2b[l]) || (indiv1a[l]!=indiv1b[l] && indiv1b[l]==indiv2a[l] && indiv2a[l]==indiv2b[l]))
      IBS<-5
    if((indiv2a[l]==indiv2b[l] && indiv1a[l]!=indiv1b[l] && indiv1a[l]!=indiv2a[l] && indiv1b[l]!=indiv2b[l]))
      IBS<-6
    if((indiv1a[l]!=indiv1b[l] && indiv2a[l]!=indiv2b[l] && indiv1a[l]==indiv2a[l] && indiv1b[l]==indiv2b[l]) || (indiv1a[l]!=indiv1b[l] && indiv2a[l]!=indiv2b[l] && indiv1a[l]==indiv2b[l] && indiv1b[l]==indiv2a[l]))
      IBS<-7
    if(indiv1a[l]!=indiv1b[l] && indiv2a[l]!=indiv2b[l] && ((indiv1a[l]==indiv2a[l] && indiv1b[l]!=indiv2b[l]) || (indiv1a[l]==indiv2b[l] && indiv1b[l]!=indiv2a[l]) || (indiv1b[l]==indiv2a[l] && indiv1a[l]!=indiv2b[l]) || (indiv1b[l]==indiv2b[l] && indiv1a[l]!=indiv2a[l])))
      IBS<-8
    if((indiv1a[l]!=indiv1b[l] && indiv2a[l]!=indiv2b[l] && indiv1a[l]!=indiv2a[l] && indiv1a[l]!=indiv2b[l] && indiv1b[l]!=indiv2a[l] && indiv1b[l]!=indiv2b[l]))
      IBS<-9
    return(IBS)
}
  
#compute relatedness over all pairs of individuals
computeRelatedness <-function(indiv1.etaik, indiv2.etaik){
  for(i in 1:50) {
    for (j in i:50) {
      for(x in 1:K)	{
        indiv1.etaik[x]<-etaik[i,x+5]
      }
      for(x in 1:K)	{
        indiv2.etaik[x]<-etaik[j,x+5]
      }
    }
  }
}

#calculating coefficients
calcCoeff <- function(IBS, l, indiv1a, indiv1.etaik, indiv2a, indiv2.etaik, locus1, K) {
    if(IBS == 1) {
      a <- indiv1a[l]
      x <- which(locus1$V2 == a[,])
      allelei <- locus1[x,]
      
      pklallelei <- array(dim=c(K))
      for(x in 1:K) {
        pklallelei[x] <- allelei[1, x+2]
      }
      
      zi1 <- 0
      zi2 <- 0
      for(x in 1:K) {
        zi1 <- zi1 + pklallelei[x]*indiv1.etaik[x]
        zi2 <- zi2 + pklallelei[x]*indiv2.etaik[x]
      }
      
      # Under S1:
      d1 <- (zi1+zi2)/2
      d2 <- zi1*zi2
      d3 <- zi1*zi2
      d4 <- (zi1*zi2*zi2+zi1*zi1*zi2)/2
      d5 <- zi1*zi2
      d6 <- (zi1*zi2*zi2+zi1*zi1*zi2)/2
      d7 <- zi1*zi2
      d8 <- (zi1*zi1*zi2+zi2*zi2*zi1)/2
      d9 <- (zi1*zi1*zi2*zi2)
    }
    
    if(IBS == 2) {
      a <- indiv1a[l]
      b <- indiv2a[l]
      x <- which(locus1$V2 == a[,])
      y <- which(locus1$V2 == b[,])
      allelei <- locus1[x,]
      allelej <- locus1[y,]
      
      pklallelei <- array(dim=c(K))
      pklallelej <- array(dim=c(K))
      for(x in 1:K) {
        pklallelei[x] <- allelei[1, x+2]
        pklallelej[x] <- allelej[1, x+2]
      }
      
      zi1 <- 0
      zi2 <- 0
      zj1 <- 0
      zj2 <- 0
      for(x in 1:K) {
        zi1 <- zi1 + pklallelei[x]*indiv1.etaik[x]
        zi2 <- zi2 + pklallelei[x]*indiv2.etaik[x]
        zj1 <- zj1 + pklallelej[x]*indiv1.etaik[x]
        zj2 <- zj2 + pklallelej[x]*indiv2.etaik[x]
      }
      
      # Under S2:
      d1 <- 0
      d2 <- (zi1*zj2+zj1*zi2)/2
      d3 <- 0
      d4 <- (zi1*zj2*zj2+zj1*zi2*zi2)/2
      d5 <- 04
      d6 <- (zi1*zj2*zj2+zj1*zi2*zi2)/2
      d7 <- 0
      d8 <- 0
      d9 <- (zi1*zi1*zj2*zj2+zj1*zj1*zi2*zi2)/2
    }}

#do optimization 
doOptimization <- function(ibds) {
  mle <- solnp(ibds, fun = loglikibd, eqfun = eqn1, eqB = c(1), 
               ineqfun = ineq1, ineqLB = c(0, 0, 0, 0, 0, 0, 0, 0, 0), 
               ineqUB = c(1, 1, 1, 1, 1, 1, 1, 1, 1))
  return(mle)}

#calculate relatedness
calcRelatedness1 <- function(mle) {
  2*(mle$pars[1]+0.5*(mle$pars[3]+mle$pars[5]+mle$pars[7])+0.25*(mle$pars[8]))}

calcrRlatedness2 <- function(mle) {
  2*(mle$pars[7]*0.5+0.25*mle$pars[8])}
    
    
