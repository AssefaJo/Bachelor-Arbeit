library(geomorph)
library(shapes)
library(Morpho)

data(gorf.dat)

norm_vec <- function(x) sqrt(sum(x^2))

#Exponentialabbilung
expo <- function(x,v){
  
  k <- nrow(x)
  m <- ncol(x)
  
  #x zu Vektor umschreiben
  x_vec <- as.vector(t(x))
  
  #Bedingung: v muss orthogonal zu x sein.
  if (sum(x_vec*v)==0){
    
    #Norm von v und x
    nv <- norm_vec(v)
    nx <- norm_vec(x_vec)
    
    #Berechnung des k*m dimensionalen shapes als vektor
    e <- cos(nv)*x_vec+sin(nv)*(nx*v)/nv
    
    #Umwandlung in eine Configuration Matrix
    matrix(e, nrow = k, ncol = m, byrow=TRUE)
    
  }
  else{message('v nicht im Tangentialraum von x')}
  
}

#Riemannscher Logarithmus
loga <- function(x,y){
  
  #x und y zu Vektor umschreiben
  x_vec <- as.vector(t(x))
  y_vec <- as.vector(t(y))
  
  #Norm von y und x
  ny <- norm_vec(y_vec)
  nx <- norm_vec(x_vec)
  
  #theta
  t <- acos(abs(sum(x_vec*y_vec))/(nx*ny))
  
  #Projektion von y auf x
  pi <- x_vec*sum(x_vec*y_vec)/nx^2
  
  #Berechnung und Ausgabe des Vektors im Tangentialraum von x
  (t*(y_vec-pi))/norm_vec(y_vec-pi)
}


data(boneData)
proc <- procSym(boneLM)
pop_sex <-name2factor(boneLM, which=3:4)

pcaplot3d(proc,pcshow=1:3,mag=-3)

proc$PCs
proc$PCscores

mean <- proc$mshape
#mean ist Punkt auf dem exp definiert ist, 
#benötige noch implementierung normalverteilter
#Tangentialvektoren

