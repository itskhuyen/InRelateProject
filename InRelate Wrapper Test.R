# define readIndivq function
readIndivq <- function(indivqfile){
  indivq <- getSubpop(indivqfile);
  return(indivq)
}

getSubpop <- function(indivqfile){
  subpop <- indivqfile[,6:ncol(indivqfile)]
  return(subpop)
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

# define get_data function wrap all read in functions
upload_data <- function() {
  strfile <- read.table("example.str", header = TRUE)
  indivqfile <- read.table("example.indivq")
  pklafile <- read.table("example.pkla")
  
  return(list(strfile = strfile, indivqfile = indivqfile, pklafile = pklafile))
}

# define read_data function wrap all read in functions needed to extract certain values
#such as total pairs of individuals or K
read_data <- function(strfile, indivqfile, pklafile) {
  
  # read in the files
  str <- readStr(strfile)
  indiv <- readIndivq(indivqfile)
  pkla <- readPkla(pklafile)
  subpop <- getSubpop(indivqfile)
  
  # return a list with the data frames
  return(list(str = str, indiv = indiv, pkla = pkla, subpop=subpop))
}

# get the data
data_list <- upload_data()

# read in and extract the data
extracted_data <- read_data(data_list$strfile, data_list$indivqfile, data_list$pklafile)

# print the extracted data
print(extracted_data)

# Setting initial values for optimization.
# You can change this if you like, if you know some of the related pairs. Else leave it be.
ibds<-c(0.0,0.0,0.0,0.0,0.0,0.0,0.25,0.5,0.25)

#also takes in predictor??
calcRelatedness <- function(ibds, filename1) {
  library(Rsolnp)
  
  loglikibd<-function(ibds) {
    # array has to be of dim=c(1, number of loci)
    loglik<-array(0,dim=c(1,300))
    for(l in 1:300){
      sumibdloci<-0
      for(i in 1:9){
        sumibdloci<-sumibdloci+predictors[l,i]*ibds[i]
      }
      
      loglik[l]<-log(sumibdloci)
    }
    return(-sum(loglik))
  }
  
  eqn1<-function(ibds){
    sum=ibds[1]+ibds[2]+ibds[3]+ibds[4]+ibds[5]+ibds[6]+ibds[7]+ibds[8]+ibds[9]
    return(sum)
  }
  
  ineq1<-function(ibds){
    z1=ibds[1]
    z2=ibds[2]
    z3=ibds[3]
    z4=ibds[4]
    z5=ibds[5]
    z6=ibds[6]
    z7=ibds[7]
    z8=ibds[8]
    z9=ibds[9]
    return(c(z1,z2,z3,z4,z5,z6,z7,z8,z9))
  }
  
  deltafile <- paste(filename1, ".deltas", sep = "")
  relatfile <- paste(filename1, ".relat", sep = "")
  
  relat <- matrix(nrow = 1, ncol = 3)
  g <- 1
  
  for (i in 1:8) {
    for (j in (i+1):9) {
      decrementi <- 1
      decrementj <- 0
      
      ibds_sub <- ibds[-c(i, j)]
      
      #do optimization
      mle <- solnp(ibds_sub, fun = loglikibd, eqfun = eqn1, eqB = c(1), 
                   ineqfun = ineq1, ineqLB = c(0, 0, 0, 0, 0, 0, 0, 0), 
                   ineqUB = c(1, 1, 1, 1, 1, 1, 1, 1))
      
      cat(c(mle$pars[1], mle$pars[2], mle$pars[3], mle$pars[4], mle$pars[5], mle$pars[6], mle$pars[7], mle$pars[8], mle$pars[9]), "\n", file = deltafile, append = TRUE)
      
      #build relatedness table
      relat[g, 1] <- sprintf("%d_%d", i, j)
      relat[g, 2] <- calcRelatedness1(mle)
      relat[g, 3] <- calcRelatedness2(mle)
      
      #decrementj <- decrementj + 2
      #decrementi <- decrementi + 2
      g <- g + 1
    }
  }
  
  calcRelatedness1 <- function(mle) {
    2*(mle$pars[1]+0.5*(mle$pars[3]+mle$pars[5]+mle$pars[7])+0.25*(mle$pars[8]))
  }
  
  calcRelatedness2 <- function(mle) {
    2*(mle$pars[7]*0.5+0.25*mle$pars[8])
  }
  
  colnames(relat) <- c("Pair", "Relatedness1", "Relatedness2")
  write.table(relat, relatfile)
}

