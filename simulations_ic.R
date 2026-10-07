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
list(DT=DT, GR=GR, real=0.5)}

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
list(DT=DT, GR=GR, real=0.96)}

#########

simulIc<- function(n0,n1,t) {
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
list(DT=DT, GR=GR, real=0.96)}

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
list(DT=DT, GR=GR, real=0.96)}

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
list(DT=DT, GR=GR, real= 1.00)}

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
list(DT=DT, GR=GR, real=0.991)}

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
list(DT=DT, GR=GR, real=0.991)}

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
list(DT=DT, GR=GR, real= 0.991)}

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
list(DT=DT, GR=GR, real=0.98)}

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
list(DT=DT, GR=GR, real=0.76)}

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
list(DT=DT, GR=GR, real=0.79)}

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
list(DT=DT, GR=GR, real= 0.5)}

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
dis = rnorm(pts, 0, 0.005);
dis = cumsum(dis);
DT1[j,]<- f1 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR, real=0.86)}

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
list(DT=DT, GR=GR, real=0.90)}

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
dis = rexp(pts, 200) - 0.005;
dis = cumsum(dis);
DT1[j,]<- f1 + dis}
# Todos
DT<- rbind(DT0, DT1)
list(DT=DT, GR=GR, real=0.87)}


#
# Diferentes metodos: logit (lineal, quadratico). farma (maximo, minimo, integral). Propuesto
# Tambien la cosa del IC (training/testing approach)
#

logit<- function(DT, GR, t) {
FourierBasis <- create.fourier.basis(rangeval = c(min(t),max(t)), nbasis=7)
curves.fd <- Data2fd(argvals = t, y=t(DT),basisobj = FourierBasis)
curves2.fd <- Data2fd(argvals = t, y=t(DT^2),basisobj = FourierBasis)
lin <- logitFD.pc(Response= GR, FDobj=list(curves.fd), ncomp = 2)
quad <- logitFD.pc(Response= GR, FDobj=list(curves.fd, curves2.fd), ncomp = c(2,2))
a<- as.numeric(lin$ROC[[9]])/100
b<- as.numeric(quad$ROC[[9]])/100
return(c(a,b))}

farma<- function(DT, GR, t) {
fmax<- rep(NA,nrow(DT))
fmin<- rep(NA,nrow(DT))
int<- rep(NA,nrow(DT))
for (i in 1:n) {
 fmax[i]<- max(DT[i,])
 fmin[i]<- min(DT[i,])
 int[i]<- abs(t[2]-t[1])*sum(DT[i,]) }
a<- AUC(fmax, GR)
b<- AUC(fmin, GR)
c<- AUC(int, GR)
return ( c(a,b,c) )}

PBC<- function(DT, GR, t) {
n<- length(GR)
AX<- matrix(0, ncol= n, nrow= n)
for (j in 1:(n-1)) {for (i in (j+1):n) {AX[j,i]<- abs(t[2]-t[1])*sum( (DT[i,] - DT[j,])^2); AX[i,j]=AX[j,i]}}
P1<- matrix(NA,ncol=2, nrow=n)
a<- rep(NA,2)
for (i in 1:n)  {
mk<- AX[i,-i]
gr<- GR[-i]
P1[i,1]<- as.numeric(gROC(mk,gr, side="both")$auc)
P1[i,2]<- as.numeric(gROC(mk,gr, side="left")$auc)
}
a[1]<- AUC(P1[,1], GR)
a[2]<- AUC(P1[,2], GR)

cr<- ifelse( (max(a) - min(a))<0.15, which.min(a),  which.max(a))
P<- P1[,cr]
area<- a[cr]
crit<- ifelse(cr==1, "both", "left")
list(distance= AX, probs= P, auc= area, criterio=crit)}

#
# Checo la calidad del CI
#

PBC.IC<- function(DT, GR, t, prob= 0.333)
{
 ix1<- which(GR== 1); n1<- length(ix1); n1m<- round(n1*prob); n1t<- n1 - n1m
 ix0<- which(GR== 0); n0<- length(ix0); n0m<- round(n0*prob); n0t<- n0 - n0m
 sel<- rep(1, (n1+n0))
 sel[c(sample(ix1, n1m), sample(ix0, n0m))]= 0
 DTo<- DT[sel==0,]; GRo<- GR[sel==0]
 DTn<- DT[sel==1,]; GRn<- GR[sel==1]
 md<- PBC(DTo, GRo, t)
 no<- length(GRo); nn<- length(GRn)
 AX<- matrix(0, ncol= no, nrow= nn)
 for (i in 1:no) {for (j in 1:nn) {AX[j,i]<- abs(t[2]-t[1])*sum( (DTo[i,] - DTn[j,])^2)}}
 P<- rep(NA, nn)
 for (i in 1:nn)   P[i]<- as.numeric(gROC(AX[i,], GRo, side= md$criterio)$auc)
 area<- AUC(P, GRn)
 sd<- ((1/n1t)*var(ecdf(P[GRn==0])(P[GRn==1])) + (1/n0t)*var(ecdf(P[GRn==1])(P[GRn==0])))^0.5
 sd<- max(sd, 1/(n1t+n0t))
 ci<- c(area - 1.96*sd, area + 1.96*sd)
 list(DT= DTo, GR= GRo, ic=ci, criterio= md$criterio, auc=area)
}

PBCreal<- function(DTn, GRn, DTo, GRo, t, CI, direction, real)
{
 nn<- length(GRn)
 no<- length(GRo)
 AX<- matrix(0, ncol= no, nrow= nn)
 for (i in 1:no) {for (j in 1:nn) {AX[j,i]<- abs(t[2]-t[1])*sum( (DTo[i,] - DTn[j,])^2)}}
 P<- rep(NA, nn)
 for (i in 1:nn)   P[i]<- as.numeric(gROC(AX[i,], GRo, side= direction)$auc)
 realS<- AUC(P, GRn)
 result1<- ifelse(realS<CI[1] | realS> CI[2], 0, 1)
 result2<- ifelse(real<CI[1] | real> CI[2], 0, 1)
 return(c(realS, real, result1, result2))
}

### Simulaciones
##### Parametros

n0<- 200
n1<- 100
n<- n0+n1
B<- 1000
Results<- matrix(-1, ncol= 7, nrow=B)
Areas<- matrix(-1, ncol= 6, nrow=B)
dtN<- simulIVd(2500,2500,t)

for (b in 1:B) {
dtO<- simulIVd(n0,n1,t)
Areas[b,1:3]<- farma(dtO$DT, dtO$GR, t)
Areas[b,4:5]<- logit(dtO$DT, dtO$GR, t)
pbc<- PBC.IC(dtO$DT,dtO$GR, t)
Areas[b,6]<- pbc$auc
Results[b,]<- c(PBCreal(dtN$DT, dtN$GR, pbc$DT, pbc$GR, t, pbc$ic, pbc$criterio, dtO$real), pbc$auc, pbc$ic[1:2])
print(round(c(b, Results[b,], colMeans(Results[1:max(2,b),3:4]),colMeans(Areas[1:max(2,b),1:6])),3))
}

Areas<- as.data.frame(Areas)
names(Areas)<- c("Max", "Min", "Integral", "Linear", "Quadratic", "Prob")
Areas$Model<- "Model IV-d"
Areas$size<- "(200,100)"
write.csv(Areas, "mdIVd.sizes4.csv", row.names=FALSE)

Results<- as.data.frame(Results)
names(Results)<- c("Real.sample", "Real", "coverage1", "coverage2", "AUC", "IC.L", "IC.U")
Results$Model<- Areas$Model
Results$size<- Areas$size
write.csv(Results, "icIVd.sizes4.csv", row.names=FALSE)


















