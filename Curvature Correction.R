library(geomorph)
library(shapes)
library(Morpho)

norm_vec <- function(x) sqrt(sum(x^2))

#Exponentialabbilung
expo <- function(x,v){
  
  k <- nrow(x)
  m <- ncol(x)
  
  #x und v zu Vektor umschreiben
  x_vec <- as.vector(t(x))
  v_vec <- as.vector(t(v))
  
  #Bedingung: v muss orthogonal zu x sein.
  if (sum(x_vec*v_vec)<0.00000001){
    
    #Norm von v und x
    nv <- norm_vec(v_vec)
    nx <- norm_vec(x_vec)
    
    #Berechnung des k*m dimensionalen shapes als vektor
    e <- cos(nv)*x_vec+sin(nv)*(nx*v_vec)/nv
    
    #Umwandlung in ein 3D array
    matrix(e, nrow = k, ncol = m, byrow=TRUE)
    
  }
  else{message('v nicht im Tangentialraum von x')}
  
}

#Riemannscher Logarithmus
loga <- function(x,y){
  
  k <- nrow(x)
  m <- ncol(x)
  
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
  vec<-(t*(y_vec-pi))/norm_vec(y_vec-pi)
  
  #Umwandlung in ein 3D array
  matrix(vec, nrow = k, ncol = m, byrow=TRUE)
}


data(boneData)
proc <- procSym(boneLM,orp=FALSE)

#Wäre das nun ein Vektor im Tangentialraum vom Meanshape?
#der Mean shape selbst subtrahiert 
#von der ersten orthogonalen Projektion 
#(quasi der verbindungsvektor von Projektion und mean shape)
v<-Morpho:::orp(proc$rotated)[,,2]-proc$mshape
sum(diag(v%*%t(x)))#prüfe orthogonalistät (nahezu orthogonal)


#Beispielrechnung für expo und loga:
x<-proc$mshape
#Beispiel shape y
y<-proc$rotated[,,3]

#Erhalte v als Vektor im Tangentialraum von x
v<-loga(x,y)
v
#Vergleiche expo(x,v) mit dem tatsächlichen Wert von y
y-expo(x,v)#Wieso ist dies nicht null? Fehler Im code oder Ungenauigkeit des Algorithmus?


