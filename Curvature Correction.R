library(matlib)
library(shapes)
library(Morpho)

norm_vec <- function(x) sqrt(sum(x^2))

#Gibt die euklidische Norm eines 2d-arrays zurückkk
norm_arr <- function(x) norm_vec(as.vector(x))

#Exponentialabbilung
expo <- function(x,v){
  
  k <- nrow(x)
  m <- ncol(x)
  
  #x und v zu Vektor umschreiben
  x_vec <- as.vector(x)
  v_vec <- as.vector(v)
  
  #Norm von v und x
  nv <- norm_vec(v_vec)
  nx <- norm_vec(x_vec)
  
  #Berechnung des k*m dimensionalen shapes als vektor
  e <- cos(nv)*x_vec+((sin(nv)*nx)/nv)*v_vec
  e
  #Umwandlung in ein 2D array    
  #matrix(e, nrow = k, ncol = m, byrow=FALSE)
}  
  

#Riemannscher Logarithmus
loga <- function(x,y){
  
  k <- nrow(x)
  m <- ncol(x)
  
  #x und y zu Vektor umschreiben
  x_vec <- as.vector(x)
  y_vec <- as.vector(y)
  
  #Norm von y und x
  ny <- norm_vec(y_vec)
  nx <- norm_vec(x_vec)
  
  #y normieren
  
  #y_vec <- y_ve/n
  #ny <- norm_vec(y_vec)
  
  #theta
  t <- acos(abs(sum(x_vec*y_vec))/(nx*ny))
  
  #Projektion von y auf x
  pi <- (x_vec*sum(x_vec*y_vec))/(nx^2)
  
  #Berechnung und Ausgabe des Vektors im Tangentialraum von x
  vec<-(t*(y_vec-pi))/norm_vec(y_vec-pi)
  vec
  #Umwandlung in ein 2D array
  #matrix(vec, nrow = k, ncol = m, byrow=FALSE)
}


#Funktion zur Schätzung der Kovarainzmatrix von 
#normalverteilten pre-Shapes mit vorgegebener Kovarianz im Tangentialraum
#einmal exakt und einmal orthogonal projiziert
pca <- function(N,dimL=2,Landmarks=3,derivation=0.3){
  
  #Dimension des Pre-Shape Space
  d=dimL*(Landmarks-1)
  
  #Dimension des Tangentialraums
  c=d-1
  
  #Generiere array e mit N standardnormalverteilten d-dimensionalen Vektoren
  e<-array(0,c(c,N))
  
  for(i in 1:N){e[,i]<-as.vector(rnorm(c,mean=0,sd=1))}
  
  #c\times c Kovarianzmatrix
  cov<-diag(rep(derivation,times=c))
  #E ist zu N(0,cov*t(cov)) verteilt
  E<-array(0,c(c,N))
  for(i in 1:N){E[,i]<-cov%*%e[,i]}
  
  
  #Basis b von T_xM mit PCs:
  #Dimension der Basisvektoren
  b<-diag(rep(1,times=d))[,-d]
  
  #b<-proc$PCs
  #Erhalte nun einen pre shape X im Tangentialraum aus Linearkombination:
  #Erstelle leeres 16 dim. array
  t<-array(0,c(d,N))
  #t<-(b_1 * E_11 + .. + b_m* E_m1)
  x<-array(c(rep(0,times=d-1),1))
  
  x
  for(j in 1:N){for(i in 1:c){t[,j]<-t[,j]+(b[,i]*E[i,j])}}
  
  #Unsere generierten Vektoren sind sogar orthogonal zum mean-shape.
  #t(t[,5])%*%as.vector(x)#10^-15
  #x <- cbind(c(1,0),c(0,0))
  z <- array(0,c(d,N))
  #10 shapes auf meinem pre shape space
  for(i in 1:N){z[,i]<-expo(x,t[,i])}
  
  #Bestimmung des Fehlers der orthogonalen Projektion mittels simulierten Daten:
  
  #Der mit log bestimmte exakten Wert von z im Tangentialraum. 
  z_l <- array(0,c(d,N))
  for(i in 1:N){z_l[,i]<-loga(x,z[,i])}
  #Falls die Länge meiner Vektoren länger als pi/2=1.57 wird, macht dies keinen Sinn mehr.
  #for(i in 1:N){message(norm_vec(t[,i]))}
  
  #Damit folgt für den shape im Tangentialraum:
  z_exakt<-array(0,c(d,N))
  for(i in 1:N){z_exakt[,i]<-z_l[,i]+x}
  
  z_orp<-z
  #Projiziere z orthogonal in den Tangentialraum von x
  for (i in 1:N) {z_orp[d,i]<-1}
  
  #Transferiere die orthogonalen Shapes in Koordinaten wieder in den Tangentialraum
  z_ovec<-array(0,c(d-1,N))
  for (i in 1:N) {z_ovec[,i]<-z_orp[,i][-d]}
  
  #Die Tangentialvektoren der exakten shapes sind in Koordinaten grade E
  z_evec<-E
  #In den zwei plots sieht man den Unterschied von z_orp und z
  
  cov_e<-matrix(0,nrow = c,ncol = c)
  cov_o<-matrix(0,nrow = c,ncol = c)
  
  
  #Kovarianzmatrix von den exakten und orthogonalen shapes
  for(i in 1:N){cov_e<-cov_e+(1/N)*(z_evec[,i]%*%t(z_evec[,i]))}
  for(i in 1:N){cov_o<-cov_o+(1/N)*(z_ovec[,i]%*%t(z_ovec[,i]))}
  
  EDe<-eigen(cov_e)
  EDo<-eigen(cov_o)
  A<-EDe$vectors%*%diag(sqrt(EDe$values))%*%t(EDe$vectors)
  B<-EDo$vectors%*%diag(sqrt(EDo$values))%*%t(EDo$vectors)
  
  out= (list(Cov_ex = A,Cov_ap = B,cov = cov))
  return(out)
}
pca(10000,derivation=0.01)
procSym(gorf.dat)






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
