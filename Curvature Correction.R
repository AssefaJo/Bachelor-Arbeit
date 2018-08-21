library(matlib)
library(shapes)
library(Morpho)

norm_vec <- function(x) sqrt(sum(x^2))

#Gibt die euklidische Norm eines 2d-arrays zurückkk
norm_arr <- function(x) norm_vec(as.vector(x))

#Exponentialabbilung
expo <- function(x,v){

  #Norm von v und x
  nv <- sqrt(sum(v^2))
  nx <- sqrt(sum(x^2))
  
  #Berechnung des Pre-Shape auf S^k_m
  pshape <- cos(nv)*x+((sin(nv)*nx)/nv)*v
  
  return(pshape)
}  
  

#Riemannscher Logarithmus
loga <- function(x,y){
  
  #Norm von y und x
  ny <- sqrt(sum(y^2))
  nx <- sqrt(sum(x^2))

  #theta
  t <- acos(abs(sum(x*y))/(nx*ny))
  
  #Projektion von y auf x
  pi <- (x*sum(x*y))/(nx^2)
  
  #Berechnung des Pre-Shape Vektors
  pshape_v <- (t*(y-pi))/sqrt(sum((y-pi)^2))
  
  return(pshape_v)
}


#Funktion zur Schätzung der Kovarainzmatrix von 
#normalverteilten pre-Shapes mit vorgegebener Kovarianz im Tangentialraum
#einmal exakt und einmal orthogonal projiziert
pca <- function(N,dimL=2,Landmarks=3,r=0.2){
  
  if(Landmarks<=dimL) {
    return(message("Die Anzahl der Landmarks muss groesser sein als die Dimension der Landmarks"))
  }
  
  #Dimension des umliegenden Raumes
  d = dimL * (Landmarks-1)
  
  #Dimension des Pre-Shape Space/Tangentialraums
  c = d - 1
  
  #Parameter Sigma, den es zu schaetzen gilt
  Sigma <- diag(rep(r,times=c))
  
  #Generiere array e mit N standardnormalverteilten c-dimensionalen Vektoren
  e<-array(0,c(c,N))
  
  for(i in 1:N) {
    e[,i] <- as.vector(rnorm(c,mean=0,sd=1))
  }
  
  #E ist N(0,Sigma%*%t(Sigma))-verteilt
  E<-array(0,c(c,N))
  
  for(i in 1:N) {
    E[,i] <- Sigma %*% e[,i]
  }
  
  #Mean Shape x als "Nordpol" der Einheitssphäre
  x <- array(c(rep(0,times=c),1))
  
  #In Koordinaten dargestellte orthonormale Basis des Tangentialraums an x
  b <- diag(rep(1,times=d))[,-d]
  
  #Vektoren erzeugen als Linearkombination
  v <- array(0,c(d,N))
  
  for(j in 1:N) {
    for(i in 1:c) {
      v[,j] <- v[,j] + (b[,i] * E[i,j])
    }
  }
  
  #Falls die Laenge der Vektoren laenger als pi/2=1.57 wird, macht eine orthogonale Projektion keinen Sinn.
  for(i in 1:N) {
    if(sqrt(sum(v[,i]^2))>pi/2) {
      return(message("Standardabweichung zu gross"))
    }
  }
  
  #Anwendung der Exponentialabbildung liefert Pre-Shapes ps auf M
  ps <- array(0,c(d,N))
  
  for(i in 1:N) {
    ps[,i] <- expo(x,v[,i])
  }
  
  #Bestimmung des Fehlers der orthogonalen Projektion mittels simulierten Daten:
  
  
  
  #Der mit log bestimmte exakte Wert von ps im Tangentialraum. 
  ps_exakt <- array(0,c(d,N))
  
  for(i in 1:N) {
    ps_exakt[,i] <- loga(x,ps[,i])+x
  }
  
  #Projiziere ps orthogonal in den Tangentialraum von x.
  ps_orp <- ps

  for (i in 1:N) {
    ps_orp[d,i] <- 1
  }

  
  #Annihilieren der letzten Koordinate liefert den Vektor in T_xM
  
  
  
  #Orthogonal projizierter Vektor
  v_orp <- array(0,c(d-1,N))
  
  for (i in 1:N) {
    v_orp[,i] <- ps_orp[,i][-d]
  }
  
  #Exakt projizierter Vektor
  v_exakt <- array(0,c(d-1,N))
  
  for (i in 1:N) {
    v_exakt[,i] <- ps_exakt[,i][-d]
  }
  
  
  
  
  #Berechnen der Kovarianzmatrix der exakten Daten
  cov_exakt <- matrix(0,nrow = c,ncol = c)
  
  for(i in 1:N) {
    cov_exakt <- cov_exakt + (1/N) * (v_exakt[,i] %*% t(v_exakt[,i]))
  }
  
  
  #Berechnen der Kovarianzmatrix der approximierten Daten
  cov_orp <- matrix(0,nrow = c,ncol = c)
  
  for(i in 1:N) {
    cov_orp <- cov_orp + (1/N) * (v_orp[,i] %*% t(v_orp[,i]))
  }
  
  
  
  
  
  
  #Singulaerwertzerlegung
  EDe <- eigen(cov_exakt)
  
  EDo <- eigen(cov_orp)
  
  cov_ex <- EDe$vectors %*% diag(sqrt(EDe$values)) %*% t(EDe$vectors)
  cov_ap <- EDo$vectors %*% diag(sqrt(EDo$values)) %*% t(EDo$vectors)
  
  A <- cov_ex - cov_ap
  
  diff <- sqrt(max(eigen(t(A)%*%A)$values))
  
#  meandiff <- 0
#  for (i in 1:c) {
#    for (j in 1:c) {
#      meandiff <- meandiff+((1/(c*c))*diff[i,j])
#    }
#  }
    
  out = (list(cov_ex = cov_ex,cov_ap = cov_ap,cov = Sigma,meandiff=diff))
   
  return(out)
}

pca(1000,r=0.2)
deriva<-seq(0.01,0.3,0.002)
x<-rep(0,times=length(deriva))
for (k in 1:length(deriva)) {
  x[k]<-pca(1000,r = deriva[k])$meandiff
}

plot(deriva,x,xlab="Standardabweichung",ylab="")



N=10
d=4
c=d-1
e<-array(0,c(c,N))

for(i in 1:N){e[,i]<-as.vector(rnorm(c,mean=0,sd=1))}

#c\times c Kovarianzmatrix
cov<-diag(rep(0.01,times=c))
#E ist zu N(0,cov*t(cov)) verteilt
E<-array(0,c(c,N)) 
for(i in 1:N){E[,i]<-cov%*%e[,i]}


#Basis b von T_xM mit PCs:
#Dimension der Basisvektoren
b<-diag(rep(1,times=d))[,-d]

#b<-proc$PCs
#Erhalte nun einen pre shape X im Tangentialraum aus Linearkombination:
#Erstelle leeres 16 dim. array
t<-array(0,c(d/2,2,N))
#t<-(b_1 * E_11 + .. + b_m* E_m1)
x<-array(c(rep(0,times=d-1),1),c(d/2,2))

x
for(j in 1:N){for(i in 1:c){t[,,j]<-t[,,j]+(b[,i]*E[i,j])}}
t
#Unsere generierten Vektoren sind sogar orthogonal zum mean-shape.
#t(t[,5])%*%as.vector(x)#10^-15
#x <- cbind(c(1,0),c(0,0))
z <- array(0,c(d/2,2,N))
#10 shapes auf meinem pre shape space
for(i in 1:N){z[,,i]<-expo(x,t[,,i])}
z
#Bestimmung des Fehlers der orthogonalen Projektion mittels simulierten Daten:

#Der mit log bestimmte exakten Wert von z im Tangentialraum. 
z_l <- array(0,c(d,N))
for(i in 1:N){z_l[,i]<-loga(x,z[,i])}
#Falls die Länge meiner Vektoren länger als pi/2=1.57 wird, macht dies keinen Sinn mehr.
#for(i in 1:N){message(norm_vec(t[,i]))}
z_l
z
#Damit folgt für den shape im Tangentialraum:
z_exakt<-array(0,c(d,N))
for(i in 1:N){z_exakt[,i]<-z_l[,i]+x}
Morpho:::orp(z,mshape=x)[,,2]
z_orp<-z
#Projiziere z orthogonal in den Tangentialraum von x
for (i in 1:N) {z_orp[d,i]<-1}

#Transferiere die orthogonalen Shapes in Koordinaten wieder in den Tangentialraum
z_ovec<-array(0,c(d-1,N))
for (i in 1:N) {z_ovec[,i]<-z_orp[,i][-d]}

#Die Tangentialvektoren der exakten shapes sind in Koordinaten grade E
z_evec<-E

plotshapes(proc$rotated[,,1:10],color = "red")
points(x,col="green",pch=19)



