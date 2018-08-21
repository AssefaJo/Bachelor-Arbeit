library(matlib)
library(Morpho)
norm_vec <- function(x) sqrt(sum(x^2))

#Gibt die euklidische Norm eines 2d-arrays zurückkk
norm_arr <- function(x) norm_vec(as.vector(x))

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
  
  matrix(e, nrow = k, ncol = m, byrow=FALSE)
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
  
  #Umwandlung in ein 2D array
  matrix(vec, nrow = k, ncol = m, byrow=FALSE)
}


data(gorf.dat)

proc <- procSym(gorf.dat,orp=TRUE,scale=F)

x <- proc$mshape

N<-10
#Anzahl der Koordinaten der shapes 
d<-length(x)
#Dimension des Shape-Space
c<-length(proc$PCs[1,])

r<-0.06

#Plots für den Vergleich von unterschiedlichen Standardabweichungen
e<-array(0,c(c,N))

for(i in 1:N){e[,i]<-as.vector(rnorm(c,mean=0,sd=1))}

#c\times c Kovarianzmatrix
cov <- diag(rep(r,times=c))
#E ist zu N(0,cov*t(cov)) verteilt
E<-array(0,c(c,N)) 
for(i in 1:N){E[,i]<-cov%*%e[,i]}

#Basis von der tangentialen Hypereben an der Mean Shape
b <- proc$PCs

t<-array(0,c(d/2,2,N))
for(j in 1:N){for(i in 1:c){t[,,j]<-t[,,j]+(b[,i]*E[i,j])}}

z <- array(0,c(d/2,2,N))
#N shapes auf Shape Space
for(i in 1:N){z[,,i]<-expo(x,t[,,i])}


#Plots:

#4.2
plotshapes(z)
points(x,col="green",pch=19)
plotshapes(proc$rotated[,,1:N])
points(x,col="green",pch=19)

#2.1
plotshapes(gorf.dat[,,1:N])
plotshapes(proc$rotated[,,1:N])
points(x,col="green",pch=19)



