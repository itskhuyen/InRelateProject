#function for getting K from indinvq file
getSubpop <- function(indivqfile){
  subpop <- indivqfile[,6:ncol(indivqfile)]
  return(subpop)
}

#function for reading in indinvq
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
getAlles <- function(pklafile){
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
  individual <- unique(strfile$V1);
  return(individual)
}

#compute the number of pairs
computenumpairs <-function(indivqfile){
  return(getIndiv(indivqfile*(indivqfile-1)/2));
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

#calculate Fst
Fst <- function(x, pop = NULL, quiet = TRUE, na.alleles = ""){
  if (any(.checkPloidy(x) != 2))
    stop("Fst() requires diploid data")
  
  NAMESX <- names(x)
  if (is.null(pop)) {
    ipop <- which(NAMESX == "population")
    if (!length(ipop)) stop("no 'population' column in x")
  } else {
    if (is.numeric(pop) && length(pop) == 1) {
      ipop <- pop
    } else {
      x$populationforthisanalysis <- factor(pop)
      ipop <- length(x)
    }
  }
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