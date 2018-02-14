#' Sample Size Determination for Longitudinal Designs with Binary Outcome
#'
#' Provides the necessary sample size for a longitudinal study with binary outcome in 
#' order to attain a pre-specified power while strictly maintaining the Type I error rate.
#' The sample size computation requires the user to define a column of design matrix 
#' relating to the slope of time as a monotonic function of time, such as linear, log, 
#' sqrt etc., along with the respective beta parameters. The underlying model is assumed 
#' to be a two-level logistic mixed-effects regression model with random intercept and/or 
#' slope of time to account for within-subject correlations and between-subject variability.
#' Gaussian quadrature is used to compute the marginal likelihood integrals and to evaluate 
#' Fisher Information matrix. 
#'
#' Attrition vector: This package allows for the specification of different attrition
#' vectors for the control and treatment group. The element of attrition vector should
#' sum to 1.
#'
#' @param nt number of time-points.
#' @param Xd design column for the slope of time (monotonic function of time).
#' @param betap vector of beta parameters (b0=Intercept, b1=slope of time for control,
#' b3=group difference at time 0 between treatment and control groups, 
#' b4=main parameter of interest which captures difference between the slope parameters of 
#' treatment and control groups).
#' @param var.ri variance of random intercept.
#' @param var.rs variance of random slope.
#' @param cov.is  covariance of intercept and slope.
#' @param ratio proportion of subjects in the control group out of the total sample.
#' @param xi1 attrition vector of the control group. The elements of attrition vector should sum to 1.
#' @param xi2 attrition vector of the treatment group. The elements of attrition vector should sum to 1.
#' @param ... optional arguments \code{alpha}, \code{power}, \code{tail}, \code{num.quad}.
#' @export ssrm.logmer
#' @return results
#' @examples
#' ssrm.logmer(nt=4,Xd=c(0,1,2,3),betap=c(1,0,0.1,0.3),var.ri=0.5,
#'             ratio=0.5,xi1=c(0,0,0,1),xi2=c(0.1,0.1,0.2,0.6))
#' ssrm.logmer(nt=4,Xd=c(0,1,2,3),betap=c(1,0,0.1,0.3),var.ri=0.5,
#'             var.rs=0.25,cov.is=0.1,power=0.90,tail=1,alpha=0.025)
#' @references Kapur K, Bhaumik R, Charlene Tang X, Hur K, Reda DJ, Bhaumik DK (2014) <doi:10.1002/sim.6203>.
#' Sample size determination for longitudinal designs with binary response. Stat Med 33(22):3781-3800.

#Imports: 
#     statmod
#     stats
#     sfsmisc
 
ssrm.logmer<-function(nt=NULL,Xd=NULL,betap=NULL,var.ri=NULL,var.rs=NULL,cov.is=NULL,ratio=NULL,xi1=NULL,xi2=NULL,...){
num_re<-2

argList<-list(...);
	  if(is.null(argList$num.quad)==TRUE){argList$num.quad<-10};
	  if(is.null(argList$num.rand)==TRUE){argList$num.rand<-2};
	  if(is.null(argList$alpha)==TRUE){argList$alpha<-0.05};
        if(is.null(argList$power)==TRUE){argList$power<-0.80};
 	  if(is.null(argList$tail)==TRUE){argList$tail<-2};

logfunc<-function(xbeta,zbeta){
y<-exp(xbeta+zbeta)/(1+exp(xbeta+zbeta))
return(y)
}

logfunc2<-function(xbeta,zbeta){
y<-(exp(xbeta+zbeta))/((1+exp(xbeta+zbeta)))^2
return(y)
}

alpha<-argList$alpha
tail<-argList$tail
cprob<-1-alpha
if(tail==1){ 
cprob<-1-alpha/2
}
cf<-stats::qchisq(cprob,df=1,ncp=0,lower.tail = TRUE,log.p=FALSE)

num_quad<-argList$num.quad
quadherm<-statmod::gauss.quad.prob(num_quad,"normal")
q.w1<-kronecker(quadherm$weights,quadherm$weights)
q.n11<-rep(quadherm$nodes,each = num_quad)
q.n12<-rep(quadherm$nodes,times = num_quad)
q.n1<-matrix(cbind(q.n11,q.n12),ncol=num_re)

eps<-1e-12
#default values
if(is.null(var.ri)==TRUE){var.ri<-eps}
if(is.null(var.rs)==TRUE){var.rs<-eps}
if(is.null(cov.is)==TRUE){cov.is<-0}

if(is.null(betap)==TRUE){betap<-c(1,0,0.1,0.3)}

if(is.null(Xd)==TRUE){Xd<-c(0,1,2)}
if(is.null(nt)==TRUE){nt<-length(Xd)}

if(is.null(xi1)==TRUE){
xi1<-numeric(length(Xd))
xi1[length(Xd)]<-1
}
if(is.null(xi2)==TRUE){
xi2<-numeric(length(Xd))
xi2[length(Xd)]<-1
}
if(is.null(ratio)==TRUE){ratio<-0.5}

if(length(Xd)!=nt){stop("Column of Design Matrix not equal to Number of Time-Points")}
if(length(betap)!=4){stop("Beta vector length != 4 (beta0, beta1, beta2, beta3)")}
if(sum(Xd)==0){stop("Not a valid Design Matrix")}
#if(Xd[1]!=0){warning("Warning: In Design Matrix Column, the start time point has been shifted to 0")}
if(length(Xd)>=20){stop("Currently, the sample size calculation is restricted <=20 time points")}
erXd<-0
for(i in 1:(nt-1)){
if(Xd[i]>Xd[i+1]){erXd<-1}
}
if(erXd==1){stop("Column of Design matrix is not monotonically increasing")}
if(ratio<=0){stop("Not valid - Ratio")}
if(ratio>=1){stop("Not valid - Ratio")}
if(sum(xi1)!=1){stop("Missing vector should sum up to 1")}
if(sum(xi2)!=1){stop("Missing vector should sum up to 1")}
if(length(Xd)!=length(xi1)){stop("Mis-match in specification of Design and Missing vectors")}
if(length(Xd)!=length(xi2)){stop("Mis-match in specification of Design and Missing vectors")}
if(length(xi1)!=length(xi2)){stop("Mis-match in specification of Missing vectors")}
if(var.ri<=eps*eps){stop("Not valid - variance")}
if(var.rs<=eps*eps){stop("Not valid - variance")}

p.n1<-ratio
p.n2<-1-p.n1


re_vcov<-matrix(c(var.ri,cov.is,cov.is,var.rs),nrow=num_re,ncol=num_re)
if(eigen(re_vcov)$values[1]<=0|eigen(re_vcov)$values[2]<=0){stop("Sigma(random-effects) not positive definite")}
T_chol<-t(chol(re_vcov))
rd<-c(0,0)
nparm<-length(betap)
ntim<-length(Xd)
xd1<-rep(1,ntim)
xd2<-Xd
#if(Xd[1]!=0){xd2<-Xd-Xd[1]}
xd3<-rep(0,ntim)
Xdesg1<-matrix(c(xd1,xd2,xd3,xd3),nrow=ntim,ncol=length(betap))
Xdesg2<-matrix(c(xd1,xd2,xd1,xd2),nrow=ntim,ncol=length(betap))
Zdesg<-matrix(c(xd1,xd2),nrow=ntim)
XmuiJg1<-(Xdesg1%*%betap)
XmuiJg2<-(Xdesg2%*%betap)
ZmuiJg<-Zdesg%*%rd
ZmuiJ<-(Zdesg%*%T_chol%*%t(q.n1))

#print(Xd)
#print(betap)
#print(num_quad)
#print(cprob)
#print(eigen(re_vcov)$values)
#print(T_chol)
#print(length(betap))
#print(ntim)

xivec1<-xi1
xivec2<-xi2

y.out<-array(0, dim=c(ntim,2^ntim,ntim))
for (it in 2:ntim){
for(i in 0:(2^it-1)){
ty<-sfsmisc::digitsBase(i, base=2, it)
for(j in 1:it){
y.out[it,(i+1),j]<-ty[j,1]
}
}
}
y.out[1,1,1]<-1

Ib<-array(0, dim = c(nparm,nparm))
Ig1<-array(0, dim = c(ntim,nparm,nparm))
Ig2<-array(0, dim = c(ntim,nparm,nparm))

for(it in 1:ntim){

hg1<-numeric(2^it)
hg2<-numeric(2^it)
dhg1<-array(0, dim = c(nparm,2^it))
dhg2<-array(0, dim = c(nparm,2^it))
d2hg1<-array(0, dim = c(nparm,nparm,2^it))
d2hg2<-array(0, dim = c(nparm,nparm,2^it))

for (p in 1:(2^it)){

lg1<-numeric(num_quad^num_re)
lg2<-numeric(num_quad^num_re)
dlg1<-array(0,dim = c(nparm,num_quad^num_re))
dlg2<-array(0,dim = c(nparm,num_quad^num_re))
for (i in 1:it){

mg1<-(logfunc(XmuiJg1[i],ZmuiJ[i,]))
mg2<-(logfunc(XmuiJg2[i],ZmuiJ[i,]))

templmg1<-log(mg1)
templmg2<-log(mg2)
templmg1[which(templmg1==-Inf)]<- -1e-30
templmg2[which(templmg2==-Inf)]<- -1e-30

templogmg1<-log(1-mg1)
templogmg2<-log(1-mg2)
templogmg1[which(templogmg1==-Inf)]<- -1e-30
templogmg2[which(templogmg2==-Inf)]<- -1e-30

temp1<-y.out[it,p,i]*templmg1 + (1-y.out[it,p,i])*templogmg1
temp2<-y.out[it,p,i]*templmg2 + (1-y.out[it,p,i])*templogmg2
temp1[which(temp1==-Inf)]<- -1e-30
temp2[which(temp2==-Inf)]<- -1e-30

lg1 <- lg1 + temp1 
lg2 <- lg2 + temp2
lg1[which(lg1==-Inf)]<- -1e-30
lg2[which(lg2==-Inf)]<- -1e-30

for(r in 1:nparm){
dlg1[r,1:num_quad^num_re] <- dlg1[r,1:num_quad^num_re] +  (y.out[it,p,i]-mg1)*Xdesg1[i,r]
dlg2[r,1:num_quad^num_re] <- dlg2[r,1:num_quad^num_re] +  (y.out[it,p,i]-mg2)*Xdesg2[i,r]
}
}

lg1 = exp(lg1)
lg2 = exp(lg2)

hg1[p]<-lg1%*%(q.w1)
hg2[p]<-lg2%*%(q.w1)
if(sum(as.numeric(is.nan(lg1)))>0){stop("Not a feasible Design (likelihood possibly does not exist). Try rescaling the Column of the Design Matrix or Beta vector.")}
if(sum(as.numeric(is.nan(lg2)))>0){stop("Not a feasible Design (likelihood possibly does not exist). Try rescaling the Column of the Design Matrix or Beta vector.")}

for (r1 in 1:nparm){
dhg1[r1,p]<- (dlg1[r1,1:num_quad^num_re]*lg1)%*%(q.w1)
dhg2[r1,p]<- (dlg2[r1,1:num_quad^num_re]*lg2)%*%(q.w1)
}

d2hg1[1:nparm,1:nparm,p]<- dhg1[1:nparm,p]%*%t(dhg1[1:nparm,p])
d2hg2[1:nparm,1:nparm,p]<- dhg2[1:nparm,p]%*%t(dhg2[1:nparm,p])
d2hg1[1:nparm,1:nparm,p]<-d2hg1[1:nparm,1:nparm,p]/hg1[p]
d2hg2[1:nparm,1:nparm,p]<-d2hg2[1:nparm,1:nparm,p]/hg2[p]
}

for (p in 1:(2^it)){
Ig1[it,1:nparm,1:nparm]<- Ig1[it,1:nparm,1:nparm] + d2hg1[1:nparm,1:nparm,p]
Ig2[it,1:nparm,1:nparm]<- Ig2[it,1:nparm,1:nparm] + d2hg2[1:nparm,1:nparm,p]
}
}

for (it in 1:ntim){
Ib<- Ib + p.n1*(Ig1[it,1:nparm,1:nparm])*xivec1[it] + p.n2*(Ig2[it,1:nparm,1:nparm])*xivec2[it]
}

if(det(Ib)==0){stop("Fisher Information Matrix is singular")}

Ibinv<-solve(Ib)

ns1<-0
ns2<-1e12
eps<-1e-6
pvmid<-0
while((pvmid-argList$power)^2>eps){ 

nsm<-(ns1+ns2)/2
pv1<-stats::pchisq(cf, 1, ncp = ns1*betap[4]^2/(Ibinv[4,4]), lower.tail = FALSE, log = FALSE)
pv2<-stats::pchisq(cf, 1, ncp = ns2*betap[4]^2/(Ibinv[4,4]), lower.tail = FALSE,log = FALSE)
pvmid<-stats::pchisq(cf, 1, ncp = nsm*betap[4]^2/(Ibinv[4,4]), lower.tail = FALSE, log = FALSE)

if(pvmid<argList$power){
ns1<-nsm
}
if(pvmid>=argList$power){
ns2<-nsm
}
}

N1<-ceiling(ratio*nsm)
N2<-ceiling((1-ratio)*nsm)
TotalN<-N1+N2
#Output

Sigma<-re_vcov

if(abs(re_vcov[1,2])<=eps*eps){
Sigma[1,2]<-0
Sigma[2,1]<-Sigma[1,2]
}

if(re_vcov[1,1]<=eps){
Sigma[1,1]<-0
}

if(re_vcov[2,2]<=eps){
Sigma[2,2]<-0
}

pow<-argList$power

results<-list(betap,Sigma,nt,xi1,xi2,Xd,num_quad,alpha,pow,tail,ratio,N1,N2,TotalN)
names(results)<-c("Beta","VarCovRE","NumTime","AttritionGrp1","AttritionGrp2","ColumnDesign","Quadpoints","alpha","power","tail","Ratio","N.group1","N.group2","Total.N") 
return(results)

#End function
}


