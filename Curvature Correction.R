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
  v_ve <- as.vector(t(v))
  
  #v_ve orthogonalisieren zu x
  v_vec <- v_ve-(sum(v_ve*x_vec)*x_vec)
  
  #Norm von v und x
  nv <- norm_vec(v_vec)
  nx <- norm_vec(x_vec)
  #Berechnung des k*m dimensionalen shapes als vektor
  e <- cos(nv)*x_vec+sin(nv)*(nx*v_vec)/nv
    
  #Umwandlung in ein 3D array    
  matrix(e, nrow = k, ncol = m, byrow=TRUE)
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


data(gorf.dat)
proc <- procSym(gorf.dat,orp=FALSE)
x<-proc$mshape
#Wäre das nun ein Vektor im Tangentialraum vom Meanshape?
#der Mean shape selbst subtrahiert 
#von der ersten orthogonalen Projektion 
#(quasi der verbindungsvektor von Projektion und mean shape)
v_orp<-Morpho:::orp(proc$rotated)[,,3]-proc$mshape
sum(diag(v%*%t(x)))#prüfe orthogonalistät (nahezu orthogonal)


#Beispielrechnung für expo und loga:
x<-proc$mshape
#Beispiel shape y
y<-proc$rotated[,,3]

#Erhalte v als 3D array im Tangentialraum von x
v<-loga(x,y)

#Vergleiche expo(x,v)(also der Umkehrfunktion von loga) mit dem tatsächlichen Wert von y
y-expo(x,v)#Das müsste null sein, da expo(loga()) Identität
#Wieso ist dies nicht null? Fehler Im code oder Ungenauigkeit des Algorithmus

#Vergleiche v mit v_orp
v-v_orp#Abweichung im Bereich 10^-6

#Vergleiche y mit Projektion von y
y_orp <- Morpho:::orp(proc$rotated)[,,3]
y-y_orp#Abweichung im Bereich von 10^-4

#Vergleiche den durch exp und loga berechneten Wert von y mit Projektion von y
y_f <- expo(x,v)
y_f-y_orp#Abweichung im Bereich von 10^-4

#wenn die obigen Implementierungen korrekt sind, kann man folgern, dass die Anwendung von loga und exp, also der Identität auf dem Tangentialraum
#zu keiner singifikanten Abweichung führt (im Bereich 10^-6)

m <- proc$orpdata[,,2]
n <- Morpho:::orp(proc$rotated)[,,2]
#Wieso sind m und n nicht identisch?
