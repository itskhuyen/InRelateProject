library(Rsolnp);
hs<-read.table("africapairs.str",header=TRUE);
etaik<-read.table("onlyrecs.indivq");
pkla<-read.table("onlyrecs_pkla");


#Initial values for IBD coefficients - deltas

#ibds<-c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2)

#Log Likelihood function that is to be maximized. Note the negative return
#for maximization
loglikibd<-function(ibds) {

loglik<-array(0,dim=c(1,783))
for(l in 1:783){
sumibdloci<-0
for(i in 1:9){
sumibdloci<-sumibdloci+predictors[l,i]*ibds[i]
}

loglik[l]<-log(sumibdloci)
}
return(-sum(loglik))
}


#equality constraint, sum of deltas = 1
eqn1<-function(ibds){
 sum=ibds[1]+ibds[2]+ibds[3]+ibds[4]+ibds[5]+ibds[6]+ibds[7]+ibds[8]+ibds[9]
 return(sum)
 }
 
 
#inequality constraints: all should be >0, <1 
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


K<-7



#sols<-array(0,dim=c(,9))
relat<-array(0,dim=c(24,3))
indiv1.etaik<-array(dim=c(K))
indiv2.etaik<-array(dim=c(K))
predictors<-array(0,dim=c(783,9))
pklallelei<-array(dim=c(K))
pklallelej<-array(dim=c(K))
pklallelek<-array(dim=c(K))
pklallelel<-array(dim=c(K))
#g<-1

#constraints<-matrix(c(1,1,1,1,1,1,1,1,1,
#1,0,0,0,0,0,0,0,0,
#0,1,0,0,0,0,0,0,0,
#0,0,1,0,0,0,0,0,0,
#0,0,0,1,0,0,0,0,0,
#0,0,0,0,1,0,0,0,0,
#0,0,0,0,0,1,0,0,0,
#0,0,0,0,0,0,1,0,0,
#0,0,0,0,0,0,0,1,0,
#0,0,0,0,0,0,0,0,1),nrow=10,ncol=9,byrow=TRUE)

#cons.dir<-c("=",">",">",">",">",">",">",">",">",">")

#rhs<-c(1,0,0,0,0,0,0,0,0,0)

	g<-1
	decrementi<-0
	decrementj<-1

i<-1
#compute over pairs of individuals
for(i in seq(1,96,by=4)) {
for(x in 1:K){
indiv1.etaik[x]<-etaik[i-decrementi,x+5]
}
for(x in 1:K){
indiv2.etaik[x]<-etaik[i+2-decrementj,x+5]
}

indiv1a<-hs[i,]
indiv1b<-hs[i+1,]
indiv2a<-hs[i+2,]
indiv2b<-hs[i+3,]

#compute across loci

for(l in 3:785) {

#here 3 indexes the locus 
#indiv1, indiv2
#a allele and b allele

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


#number of uniquealleles in locus 1
#here 0 indexes the locus
length(which(pkla$V1==l-3))

#save locus1 as a subset
#same here - 0 indexs the locus
locus1<-subset(pkla,pkla$V1==l-3)

#calculate coefficients under S1:
#let's say that k indexes the subpopulations

if(IBS==1) {
a<-indiv1a[l]
x<-which(locus1$V2==a[,])
allelei<-locus1[x,]

for(x in 1:K){
pklallelei[x]<-allelei[1,x+2]
}

zi1=0
zi2=0

for(x in 1:K){
zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
}

#Under S1:
d1<-(zi1+zi2)/2
d2<-zi1*zi2
d3<-zi1*zi2
d4<-(zi1*zi2*zi2+zi1*zi1*zi2)/2
d5<-zi1*zi2
d6<-(zi1*zi2*zi2+zi1*zi1*zi2)/2
d7<-zi1*zi2
d8<-(zi1*zi1*zi2+zi2*zi2*zi1)/2
d9<-(zi1*zi1*zi2*zi2)

}

if(IBS==2) {

allelei<-array(dim=c(K))
allelej<-array(dim=c(K))

a<-indiv1a[l]
b<-indiv2a[l]
x<-which(locus1$V2==a[,])
y<-which(locus1$V2==b[,])

allelei<-locus1[x,]
allelej<-locus1[y,]


for(x in 1:K){
pklallelei[x]<-allelei[1,x+2]
pklallelej[x]<-allelej[1,x+2]
}

zi1=0
zi2=0
zj1=0
zj2=0

for(x in 1:K){
zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
}

#Under S2
d1<-0
d2<-(zi1*zj2+zj1*zi2)/2
d3<-0
d4<-(zi1*zj2*zj2+zj1*zi2*zi2)/2
d5<-0
d6<-(zi1*zj2*zj2+zj1*zi2*zi2)/2
d7<-0
d8<-0
d9<-(zi1*zi1*zj2*zj2+zj1*zj1*zi2*zi2)/2

}

if(IBS==3) {


pklallelei<-array(dim=c(K))
pklallelej<-array(dim=c(K))


a<-indiv1a[l]
b<-indiv2b[l]
x<-which(locus1$V2==a[,])
y<-which(locus1$V2==b[,])

allelei<-locus1[x,]
allelej<-locus1[y,]

for(x in 1:K){
pklallelei[x]<-allelei[1,x+2]
pklallelej[x]<-allelej[1,x+2]
}


zi1=0
zi2=0
zj1=0
zj2=0

for(x in 1:K){
zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
}

#Under S3
d1<-0
d2<-0
d3<-(zi1*zj2+zj1*zi2)/2
d4<-(zi1*zi2*zj2+zj1*zi2*zj2)/2
d5<-0
d6<-0
d7<-0
d8<-(zi1*zi2*zj2+zj1*zi2*zj2)/2
d9<-(zi1*zi1*zi2*zj2+zj1*zj1*zj2*zi2)/2

}

if(IBS==4) {


a<-indiv1a[l]
b<-indiv2a[l]
c<-indiv2b[l]
x<-which(locus1$V2==a[,])
y<-which(locus1$V2==b[,])
z<-which(locus1$V2==c[,])

allelei<-locus1[x,]
allelej<-locus1[y,]
allelek<-locus1[z,]

pklallelei<-array(dim=c(K))
pklallelej<-array(dim=c(K))
pklallelek<-array(dim=c(K))

for(x in 1:K){
pklallelei[x]<-allelei[1,x+2]
pklallelej[x]<-allelej[1,x+2]
pklallelek[x]<-allelek[1,x+2]
}

zi1=0
zi2=0
zj1=0
zj2=0
zk1=0
zk2=0

for(x in 1:K){
zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
zk1<-zk1+pklallelek[x]*indiv1.etaik[x]
zk2<-zk2+pklallelek[x]*indiv2.etaik[x]
}

#Under S4

d1<-0
d2<-0
d3<-0
d4<-(zi1*zj2*zk2+zi2*zj1*zk1)/2
d5<-0
d6<-0
d7<-0
d8<-0
d9<-(zi1*zi1*zj2*zk2+zi2*zi2*zj1*zk1)/2
}


if(IBS==5) {
#figure out index of allele index
a<-indiv2b[l]
b<-indiv1a[l]
x<-which(locus1$V2==a[,])
#7
y<-which(locus1$V2==b[,])
#4

#pull out allele frequencies for that allele 
allelej<-locus1[x,]
allelei<-locus1[y,]


pklallelei<-array(dim=c(K))
pklallelej<-array(dim=c(K))

for(x in 1:K){
pklallelei[x]<-allelei[1,x+2]
pklallelej[x]<-allelej[1,x+2]
}


#need to define zi1,zi2,zj2,zj1

zi1=0
zi2=0
zj1=0
zj2=0

for(x in 1:K){
zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
}

#now to get conditional probabilities:
#Under S5
d1<-0
d2<-0
d3<-0
d4<-0
d5<-(zi1*zj2+zj1*zi2)/2
d6<-(zi1*zi2*zj2+zj1*zi2*zj2)/2
d7<-0
d8<-(zi1*zi2*zj2+zj1*zi2*zj2)/2
d9<-(zi1*zi1*zi2*zj2+zj1*zj1*zj2*zi2)/2
}

if(IBS==6) {


a<-indiv1a[l]
b<-indiv1b[l]
c<-indiv2b[l]
x<-which(locus1$V2==a[,])
y<-which(locus1$V2==b[,])
z<-which(locus1$V2==c[,])

allelei<-locus1[z,]
allelej<-locus1[x,]
allelek<-locus1[y,]


pklallelei<-array(dim=c(K))
pklallelej<-array(dim=c(K))
pklallelek<-array(dim=c(K))

for(x in 1:K){
pklallelei[x]<-allelei[1,x+2]
pklallelej[x]<-allelej[1,x+2]
pklallelek[x]<-allelek[1,x+2]
}

zi1=0
zi2=0
zj1=0
zj2=0
zk1=0
zk2=0

for(x in 1:K){
zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
zk1<-zk1+pklallelek[x]*indiv1.etaik[x]
zk2<-zk2+pklallelek[x]*indiv2.etaik[x]
}

#Under S6
d1<-0
d2<-0
d3<-0
d4<-0
d5<-0
d6<-(zi1*zj2*zk2+zi2*zj1*zk1)/2
d7<-0
d8<-0
d9<-(zi1*zi1*zj2*zk2+zi2*zi2*zj1*zk1)/2

}

if(IBS==7) {

a<-indiv1a[l]
b<-indiv1b[l]

x<-which(locus1$V2==a[,])
y<-which(locus1$V2==b[,])

allelei<-locus1[x,]
allelej<-locus1[y,]
pklallelei<-array(dim=c(K))
pklallelej<-array(dim=c(K))


for(x in 1:K){
pklallelei[x]<-allelei[1,x+2]
pklallelej[x]<-allelej[1,x+2]
}

zi1=0
zi2=0
zj1=0
zj2=0

for(x in 1:K){
zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
}

#Under S7
d1<-0
d2<-0
d3<-0
d4<-0
d5<-0
d6<-0
d7<-(zi1*zj2+zj1*zi2)/2
d8<-(zj1*zj2*(zi1+zi2*0.5)+(zi1*zi2*(zj1+zj2)*0.5))/2
d9<-zi1*zi2*zj1*zj2
}

if(IBS==8) {


a<-indiv1a[l]
b<-indiv1b[l]
c<-indiv2b[l]

x<-which(locus1$V2==a[,])
y<-which(locus1$V2==b[,])
z<-which(locus1$V2==c[,])

allelei<-locus1[x,]
allelej<-locus1[y,]
allelek<-locus1[z,]

pklallelei<-array(dim=c(K))
pklallelej<-array(dim=c(K))
pklallelek<-array(dim=c(K))

for(x in 1:K){
pklallelei[x]<-allelei[1,x+2]
pklallelej[x]<-allelej[1,x+2]
pklallelek[x]<-allelek[1,x+2]
}



zi1=0
zi2=0
zj1=0
zj2=0
zk1=0
zk2=0

for(x in 1:K){
zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
zk1<-zk1+pklallelek[x]*indiv1.etaik[x]
zk2<-zk2+pklallelek[x]*indiv2.etaik[x]
}

#Under S8
d1<-0
d2<-0
d3<-0
d4<-0
d5<-0
d6<-0
d7<-0
d8<-0.5*(((zi1+zi2)/2)*(zj1*zk2+zk1*zj2))
d9<-0.5*((zi1*zi2)*(zj1*zk2+zj2*zk1))

}

if(IBS==9) {



a<-indiv1a[l]
b<-indiv1b[l]
c<-indiv2a[l]
d<-indiv2b[l]

x<-which(locus1$V2==a[,])
y<-which(locus1$V2==b[,])
z<-which(locus1$V2==c[,])
z1<-which(locus1$V2==d[,])

allelei<-locus1[x,]
allelej<-locus1[y,]
allelek<-locus1[z,]
allelel<-locus1[z1,]

pklallelei<-array(dim=c(K))
pklallelej<-array(dim=c(K))
pklallelek<-array(dim=c(K))
pklallelel<-array(dim=c(K))

for(x in 1:K){
pklallelei[x]<-allelei[1,x+2]
pklallelej[x]<-allelej[1,x+2]
pklallelek[x]<-allelek[1,x+2]
pklallelel[x]<-allelel[1,x+2]
}


zi1=0
zi2=0
zj1=0
zj2=0
zk1=0
zk2=0
zl1=0
zl2=0

for(x in 1:K){
zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
zk1<-zk1+pklallelek[x]*indiv1.etaik[x]
zk2<-zk2+pklallelek[x]*indiv2.etaik[x]
zl1<-zl1+pklallelel[x]*indiv1.etaik[x]
zl2<-zl2+pklallelel[x]*indiv2.etaik[x]
}

#Under S9

d1<-0
d2<-0
d3<-0
d4<-0
d5<-0
d6<-0
d7<-0
d8<-0
d9<-(zi1*zj1*zk2*zl2+zk1*zl1*zi2*zj2)/2
}


#trying flexmix
predictors[l-2,1]<-d1
predictors[l-2,2]<-d2
predictors[l-2,3]<-d3
predictors[l-2,4]<-d4
predictors[l-2,5]<-d5
predictors[l-2,6]<-d6
predictors[l-2,7]<-d7
predictors[l-2,8]<-d8
predictors[l-2,9]<-d9

}
ibds<-c(0.0,0.0,0.0,0.0,0.0,0.0,0.25,0.5,0.25)
mle=solnp(ibds,fun=loglikibd,eqfun=eqn1,eqB=c(1),ineqfun=ineq1,ineqLB=c(0,0,0,0,0,0,0,0,0),ineqUB=c(1,1,1,1,1,1,1,1,1))
#mle=solnp(ibds,fun=loglikibd,eqfun=eqn1,eqB=c(1,0,0,0,0,0,0),ineqfun=ineq1,ineqLB=c(0,0,0),ineqUB=c(1,1,1))
#em<-multmixEM(predictors,lambda=NULL,theta=NULL,k=9,epsilon=1e-06)

#relatedness<-0
#for (h in 1:30){
#relatedness<-relatedness+(em$posterior[h,1]+0.5*(em$posterior[h,3]+em$posterior[h,5]+em$posterior[h,7])+0.25*em$posterior[h,8])*2
#}
#relatedness=relatedness/30
relatedness1<-2*(mle$pars[1]+0.5*(mle$pars[3]+mle$pars[5]+mle$pars[7])+0.25*(mle$pars[8]))
cat(c(mle$pars[1],mle$pars[2],mle$pars[3],mle$pars[4],mle$pars[5],mle$pars[6],mle$pars[7],mle$pars[8],mle$pars[9]),"\n",file="africak7.deltas",append=TRUE)
relatedness2<-2*(mle$pars[7]*0.5+0.25*mle$pars[8])
r<-i-decrementi
s<-i+2-decrementj
relat[g,1]<-sprintf("%d_%d",r,s)
relat[g,2]<-relatedness1
relat[g,3]<-relatedness2
decrementj<-decrementj+2
decrementi<-decrementi+2
g<-g+1
}

write.table(relat,"africak7.relat")

