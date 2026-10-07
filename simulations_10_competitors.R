library(logitFD)
library(fda.usc)
library(pROC)
library(nsROC)
# source("gROC.R") # Actualizar

AUC<- function(M, GR) {return(mean(ecdf(M[GR==0])(M[GR==1])))}
pts = 10000
t<- -1+ 2*(1:pts)/pts

#####
#####
##### Modelos I
simulIa<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
pts<- length(t)
f1<- 1.4*sin(pi*t)
f0<- sin(pi*t)
## A
# Negativos
for (j in 1:n0) {
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT0[j,]<- f0 + dis}
# Positivos
for (j in 1:n1) {
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT1[j,]<- f0 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#########

simulIb<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
pts<- length(t)
f1<- 1.4*sin(pi*t)
f0<- sin(pi*t)
## A
# Negativos
for (j in 1:n0) {
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT0[j,]<- f0 + dis}
# Positivos
for (j in 1:n1) {
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT1[j,]<- f1 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#########

simulIc<- function(n0,n1,t) {
pts<- length(t)
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
pts<- length(t)
f1<- 1.4*sin(pi*t)
f0<- sin(pi*t)
## A
# Negativos
for (j in 1:n0) {
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT0[j,]<- f0 + dis}
# Positivos
for (j in 1:n1) {
dis = rnorm(pts, 0, 0.01);
dis = cumsum(dis);
DT1[j,]<- f1 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#########

simulId<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
pts<- length(t)
f1<- 1.4*sin(pi*t)
f0<- sin(pi*t)
## A
# Negativos
for (j in 1:n0) {
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT0[j,]<- f0 + dis}
# Positivos
for (j in 1:n1) {
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT1[j,]<- f1 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#######
#######
####### Modelo II
f11<- -2*t*2
f12<- -f11
f01<- -2*t^2*(t<=0) + 2*t^2*(t>0)
f02<- -f01

#########

simulIIa<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
for (j in 1:n0) {
a<- ifelse(runif(1)<0.5, rnorm(1,-2,0.25), rnorm(1,2,0.25) )
DT0[j,]<- a*t^2}
# Positivos
for (j in 1:n1) {
a<- ifelse(runif(1)<0.5, rnorm(1,-2,0.25), rnorm(1,2,0.25) )
DT1[j,]<- a*t^2*(t<=0) - a*t^2*(t>0)}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#########

simulIIb<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
for (j in 1:n0) {
a<- ifelse(runif(1)<0.5, rnorm(1,-2,0.25), rnorm(1,2,0.25) )
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT0[j,]<- a*t^2 + dis}
# Positivos
for (j in 1:n1) {
a<- ifelse(runif(1)<0.5, rnorm(1,-2,0.25), rnorm(1,2,0.25) )
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT1[j,]<- a*t^2*(t<=0) - a*t^2*(t>0) + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#########

simulIIc<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
for (j in 1:n0) {
a<- ifelse(runif(1)<0.5, rnorm(1,-2,0.25), rnorm(1,2,0.25) )
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT0[j,]<- a*t^2 + dis}
# Positivos
for (j in 1:n1) {
a<- ifelse(runif(1)<0.5, rnorm(1,-2,0.25), rnorm(1,2,0.25) )
dis = rnorm(pts, 0, 0.01);
dis = cumsum(dis);
DT1[j,]<- a*t^2*(t<=0) - a*t^2*(t>0) + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#########

simulIId<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
for (j in 1:n0) {
a<- ifelse(runif(1)<0.5, rnorm(1,-2,0.25), rnorm(1,2,0.25) )
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT0[j,]<- a*t^2 + dis}
# Positivos
for (j in 1:n1) {
a<- ifelse(runif(1)<0.5, rnorm(1,-2,0.25), rnorm(1,2,0.25) )
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT1[j,]<- a*t^2*(t<=0) - a*t^2*(t>0) + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#######
#######
####### Modelo III

simulIIIa<- function(n0,n1,t) {
f1<- dnorm(t, 0.25, 0.65)
f0<- dnorm(t, -0.15, 0.65)
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
for (j in 1:n0) {
a<- rnorm(1, -0.15, 0.1)
b<- 0.5 + abs(rnorm(1,0,0.2))
DT0[j,]<- dnorm(t,a,b)}
# Positivos
for (j in 1:n1) {
a<- rnorm(1, 0.15, 0.1)
b<- 0.5 + abs(rnorm(1,0,0.2))
DT1[j,]<- dnorm(t,a,b)}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#####

simulIIIb<- function(n0,n1,t) {
f1<- dnorm(t, 0.15, 0.65)
f0<- dnorm(t, -0.15, 0.65)
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
for (j in 1:n0) {
a<- rnorm(1, -0.15, 0.1)
b<- 0.5 + abs(rnorm(1,0,0.2))
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT0[j,]<- dnorm(t,a,b) + dis}
# Positivos
for (j in 1:n1) {
a<- rnorm(1, 0.15, 0.1)
b<- 0.5 + abs(rnorm(1,0,0.2))
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT1[j,]<- dnorm(t,a,b) + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

##################

simulIIIc<- function(n0,n1,t) {
f1<- dnorm(t, 0.15, 0.65)
f0<- dnorm(t, -0.15, 0.65)
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
for (j in 1:n0) {
a<- rnorm(1, -0.15, 0.1)
b<- 0.5 + abs(rnorm(1,0,0.2))
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT0[j,]<- dnorm(t,a,b) + dis}
# Positivos
for (j in 1:n1) {
a<- rnorm(1, 0.15, 0.1)
b<- 0.5 + abs(rnorm(1,0,0.2))
dis = rnorm(pts, 0, 0.01);
dis = cumsum(dis);
DT1[j,]<- dnorm(t,a,b) + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

############################

simulIIId<- function(n0,n1,t) {
f1<- dnorm(t, 0.15, 0.65)
f0<- dnorm(t, -0.15, 0.65)
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
# Negativos
for (j in 1:n0) {
a<- rnorm(1, -0.15, 0.1)
b<- 0.5 + abs(rnorm(1,0,0.2))
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT0[j,]<- dnorm(t,a,b) + dis}
# Positivos
for (j in 1:n1) {
a<- rnorm(1, 0.15, 0.1)
b<- 0.5 + abs(rnorm(1,0,0.2))
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT1[j,]<- dnorm(t,a,b) + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR, real= 77)}

#####
#####
##### Modelos IV
simulIVa<- function(n0,n1,t) {
f1<- 0.5*t^3/(t-2)^2
f0<- -0.5*t^3/(t-2)^2
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
pts<- length(t)
## A
# Negativos
for (j in 1:n0) {
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT0[j,]<- f0 + dis}
# Positivos
for (j in 1:n1) {
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT1[j,]<- f0 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#########

simulIVb<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
pts<- length(t)
f1<- 0.5*t^3/(t-2)^2
f0<- -0.5*t^3/(t-2)^2
## A
# Negativos
for (j in 1:n0) {
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT0[j,]<- f0 + dis}
# Positivos
for (j in 1:n1) {
dis = rnorm(pts, 0, 0.01);
dis = cumsum(dis);
DT1[j,]<- f1 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

simulIVc<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
pts<- length(t)
f1<- 0.5*t^3/(t-2)^2
f0<- -0.5*t^3/(t-2)^2
## A
# Negativos
for (j in 1:n0) {
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT0[j,]<- f0 + dis}
# Positivos
for (j in 1:n1) {
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT1[j,]<- f1 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

simulIVd<- function(n0,n1,t) {
DT0<- matrix(NA, ncol=pts, nrow= n0)
DT1<- matrix(NA, ncol=pts, nrow= n1)
GR<- c(rep(0,n0),rep(1,n1))
pts<- length(t)
f1<- 0.5*t^3/(t-2)^2
f0<- -0.5*t^3/(t-2)^2
## A
# Negativos
for (j in 1:n0) {
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT0[j,]<- f0 + dis}
# Positivos
for (j in 1:n1) {
dis = rexp(pts, 100) - 0.01;
dis = cumsum(dis);
DT1[j,]<- f1 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR)}

#
# Diferentes metodos: logit (lineal, quadratico). farma (maximo, minimo, integral). Propuesto. -CALCULANDO CURVA ROC-
#

logit<- function(DT, GR, t) {
FourierBasis <- create.fourier.basis(rangeval = c(min(t),max(t)), nbasis=7)
curves.fd <- Data2fd(argvals = t, y=t(DT),basisobj = FourierBasis)
curves2.fd <- Data2fd(argvals = t, y=t(DT^2),basisobj = FourierBasis)
lin <- logitFD.pc(Response= GR, FDobj=list(curves.fd), ncomp = 2)
quad <- logitFD.pc(Response= GR, FDobj=list(curves.fd, curves2.fd), ncomp = c(2,2))
tr<- seq(0,1, 0.001)
a<- as.numeric(lin$ROC[[9]])/100
b<- as.numeric(quad$ROC[[9]])/100
Rlin<- approxfun(1 - lin$ROC$specificities/100, lin$ROC$sensitivities/100)(tr)
Rqua<- approxfun(1 - quad$ROC$specificities/100, quad$ROC$sensitivities/100)(tr)
list( t= tr, Rlin = Rlin, Rqua= Rqua, Alin = a, Aqua = b)
}

farma<- function(DT, GR, t) {
fmax<- rep(NA,nrow(DT))
fmin<- rep(NA,nrow(DT))
int<- rep(NA,nrow(DT))

for (i in 1:n) {
 fmax[i]<- max(DT[i,])
 fmin[i]<- min(DT[i,])
 int[i]<- abs(t[2]-t[1])*sum(DT[i,]) }
 
 tr<- seq(0,1, 0.001)
 a<- AUC(fmax, GR)
 ROCmax<- ecdf(1-ecdf(fmax[GR==0])(fmax[GR==1]))(tr)
 b<- AUC(fmin, GR)
 ROCmin<- ecdf(1-ecdf(fmin[GR==0])(fmin[GR==1]))(tr)
 c<- AUC(int, GR)
 ROCint<- ecdf(1-ecdf(int[GR==0])(int[GR==1]))(tr)
 list( t= tr, Rmin = ROCmin, Rmax=ROCmax, Rint=ROCint, Amin = b, Amax = a, Aint = c)
}

PBC<- function(DTn, GRn, t, DTo, GRo) {
nn<- length(GRn)
no<- length(GRo)
AX<- matrix(0, ncol= no, nrow= nn)
for (j in 1:nn) {for (i in 1:no) {AX[j,i]<- abs(t[2]-t[1])*sum( (DTo[i,] - DTn[j,])^2)}}
P1<- matrix(NA,ncol=2, nrow=nn)
a<- rep(NA,2)
for (i in 1:nn)  {
P1[i,1]<- as.numeric(gROC(AX[i,], GRo, side="both")$auc)
P1[i,2]<- as.numeric(gROC(AX[i,], GRo, side="left")$auc)
}

a[1]<- AUC(P1[,1], GRn)
a[2]<- AUC(P1[,2], GRn)
cr<- ifelse( (max(a) - min(a))<0.15, which.min(a),  which.max(a))
P<- P1[,cr]
area<- a[cr]
tr<- seq(0,1, 0.001)
ROC<- ecdf(1-ecdf(P[GRn==0])(P[GRn==1]))(tr)
crit<- ifelse(cr==1, "Both", "Left")
list(distance= AX, probs= P, auc= area, criterio=crit, t= tr, ROC=ROC)}

# Ranking procedures [Estevez and Vieu] - Doing using AI
library(pROC)

rank_fd <- function(functional_sample, GR, t) {
  curves_A  <- functional_sample[GR == 1, ]
  curves_NA <- functional_sample[GR == 0, ]
  
  mean_A  <- colMeans(curves_A)
  mean_NA <- colMeans(curves_NA)
  pilot_direction <- mean_A - mean_NA
  
  rank_scores <- apply(functional_sample, 1, function(curve) {
    sum(curve * pilot_direction) / sum(pilot_direction^2)
  })

  final_ranks <- rank(rank_scores)
  p<- seq(0,1,0.001)
  roc_curve <- roc(GR, rank_scores)
  p<- seq(0,1,0.001)
  ROC= approxfun(1- roc_curve$spec, roc_curve$sens)(p)
  return(list(score= rank_scores, gr= GR, t= p, roc= ROC, auc= roc_curve$auc ))}

# A kind of KNN

KNNf<- function(DT, GR, k=10, prob= 0.333) {
       ix1<- which(GR== 1); n1<- length(ix1); n1m<- round(n1*prob); n1t<- n1 - n1m
       ix0<- which(GR== 0); n0<- length(ix0); n0m<- round(n0*prob); n0t<- n0 - n0m
       sel<- rep(1, (n1+n0))
       sel[c(sample(ix1, n1m), sample(ix0, n0m))]= 0
       DTo<- DT[sel==0,]; GRo<- GR[sel==0]
       DTn<- DT[sel==1,]; GRn<- GR[sel==1]

       funct <- class::knn(DTo, DTn, GRo, k, prob=TRUE)
       punct<- attr(funct, "prob")
       
       roc_curve <- roc(GRn, punct)
       p<- seq(0,1,0.001)
       ROC= approxfun(1- roc_curve$spec, roc_curve$sens)(p)
       return(list(score= punct, gr= GR, t= p, roc= ROC, auc= roc_curve$auc ))}

## Depth FM

depFM<- function(DT, GR, t){
 fdata <- fdata(DT, argvals = t)
 fdata0<- fdata(DT[GR==0, ], argvals = t)
 pFM0 <- depth.FM(fdata, fdata0)$dep
 GMscores <- 1 - pFM0
 tr<- seq(0,1, 0.001)
 ROC<- ecdf(1-ecdf(GMscores[GR==0])(GMscores[GR==1]))(tr)
 return(list(score= GMscores, t=tr, ROC=ROC, auc= tr[2]*sum(ROC)))}

### Depth mode

depmode<- function(DT, GR, t){
 fdata <- fdata(DT, argvals = t)
 fdata0<- fdata(DT[GR==0, ], argvals = t)
 pFM0 <- depth.mode(fdata, fdata0)$dep
 GMscores <- 1 - pFM0
 tr<- seq(0,1, 0.001)
 ROC<- ecdf(1-ecdf(GMscores[GR==0])(GMscores[GR==1]))(tr)
 return(list(score= GMscores, t=tr, ROC=ROC, auc= tr[2]*sum(ROC)))}

### Simulaciones
##### Parametros

n0<- 200
n1<- 100
n<- n0+n1
pts = 10000
t<- -1+ 2*(1:pts)/pts

B= 1000
Areas<- matrix(-1, ncol= 5, nrow=B)

for (b in 1:B) {
dtO<- simulIVa(n0,n1,t)
Areas[b,1]<- rank_fd(dtO$DT, dtO$GR)$auc
Areas[b,2]<- KNNf(dtO$DT, dtO$GR)$auc
Areas[b,3]<- depFM(dtO$DT, dtO$GR, t)$auc
Areas[b,4]<- depmode(dtO$DT, dtO$GR, t)$auc
print(c(b, colMeans(Areas[1:max(2,b),1:4])) )
}

aux<- read.csv("C:/Users/Pablo Martinez-Cambl/Desktop/DHMC/Otros/Articulos/functionalROC_E/simulaciones_V2/Estimation/mdIVa.sizes4.csv")
Areas[,5]<- aux$Prob
summary(Areas)

Areas<- as.data.frame(Areas)
names(Areas)<- c("RK", "kNN", "D-FM", "D-MD", "PBC")
Areas$Model<- "Model IV-a"
Areas$size<- "(200,100)"
summary(Areas)
write.csv(Areas, "mdIVa.sizes4b.csv", row.names=FALSE)



























