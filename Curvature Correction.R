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
  
  #v_ve nochmal zusätzlich orthogonalisieren zu x, um Fehler klein zu halten
  #v_vec <- v_ve-(sum(v_ve*x_vec)*x_vec)
  
  #Norm von v und x
  nv <- norm_vec(v_vec)
  nx <- norm_vec(x_vec)
  #Berechnung des k*m dimensionalen shapes als vektor
  e <- cos(nv)*x_vec+sin(nv)*(nx*v_vec)/nv
    
  #Umwandlung in ein 2D array    
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
  
  #Umwandlung in ein 2D array
  matrix(vec, nrow = k, ncol = m, byrow=TRUE)
}


data(gorf.dat)
#Performing Procrustes Superimposition on gorf.dat
proc <- procSym(gorf.dat,orp=TRUE)
#Mean shape x
x<-proc$mshape

#Wäre das nun ein Vektor im Tangentialraum vom Meanshape?
#der Mean shape selbst subtrahiert 
#von der ersten orthogonalen Projektion 
#(quasi der verbindungsvektor von Projektion und mean shape)
v_orp<-proc$orpdata[,,3]-proc$mshape
sum(diag(v%*%t(x)))#prüfe orthogonalistät von v_orp und x 


#Beispielrechnung für expo und loga:
x<-proc$mshape
#Beispiel shape y
y<-proc$rotated[,,3]

#Erhalte v als 2D array im Tangentialraum von x
v<-loga(x,y)


#Hauptproblem:
#Vergleiche expo(x,v)(also der Umkehrfunktion von loga) mit dem tatsächlichen Wert von y
y-expo(x,v)#Das müsste null sein, da expo(loga()) Identität
#Wieso ist dies nicht null? Fehler Im code oder Ungenauigkeit des Algorithmus
#Ich habe den Code für sehr leichte Zahlenbeispiele getestet
#Dabei erhalte ich durchaus das Ergebnis, dass expo(loga()) die Identität ist.


#Vergleiche v mit v_orp
v-v_orp#Abweichung im Bereich 10^-6
#Dies wäre also wenn alles richtig implementiert und definiert wurde,
#der Fehler der orthogonalen Projektion im Vergleich zum riemannschen Logarithmus.(?)

plotshapes(proc$rotated[,,3:7], color = 3)
plotshapes(proc$mshape,color = 4)
#Wie kann ich diese zwei Plots in einem Schaubild plotten?

#Gibt es einen Befehl, der mich direkt zu der Kovarianzmatrix führt, 
#oder muss man diese indirekt berechnen mit den EW und PCs?
