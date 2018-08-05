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
  
  #Umwandlung in ein 2D array    
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
plotshapes(gorf.dat[,,1:30])
plotshapes
procSym
#Procrustes Superimposition
proc <- procSym(gorf.dat,orp=TRUE,scale=F)
plotshapes(proc$rotated[,,1:30])
#Mean shape x
x <- proc$mshape

#Beispielshape
y <- proc$rotated[,,8]
y_orp <- proc$orpdata[,,8]
norm_arr(y)

#Beispielrechnungen:

#Tangentialvektor von x
v_orp <- y_orp-x

v_orp - loga(x,expo(x,as.vector(v_orp)))
y - expo(x,loga(x,y))
#Identität Abweichung 10^-16

(x+loga(x,y))-proc$orpdata[,,8]#10^-5
#Abweichung des mit log berechneten, exakten shapes im Tangentialraum von der othogonalen Projektion im Tangentialraum.


plotshapes(gorf.dat[,,3:15])
plotshapes(proc$rotated[,,3:15], color = 3)
deformGrid2d(proc$orpdata[,,3],x+loga(x,proc$rotated[,,3]),wireframe = c(1,6:2,8:6))
deformGrid2d(x,proc$rotated[,,1],wireframe = c(1,6:2,8:6))


#Folgende plots gilt es zu vergleichen:
plotshapes(proc$orpdata[,,3])
plotshapes(x+loga(x,proc$rotated[,,3]))