library(geomorph)
library(shapes)
library(Morpho)

norm_vec <- function(x) sqrt(sum(x^2))

#Norm_arr normiert ein 2D-array x auf Norm 1
norm_arr <- function(x) sqrt(sum(diag(t(x)%*%x)))

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
  e <- cos(nv)*x_vec+((sin(nv)*nx)/nv)*v_vec
    
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
  matrix(vec, nrow = k, ncol = m, byrow=TRUE)
}


data(gorf.dat)
#Performing Procrustes Superimposition on gorf.dat
proc <- procSym(gorf.dat,orp=TRUE)
#Mean shape x
x<-proc$mshape

#Wäre das nun ein Vektor im Tangentialraum vom Meanshape?
#der Mean shape selbst subtrahiert 
#von der dritten orthogonalen Projektion 
#(quasi der verbindungsvektor von Projektion und mean shape)
v_orp<-proc$orpdata[,,8]-proc$mshape

#Beispielrechnung für expo und loga:
x<-proc$mshape
#Beispiel shape y
y<-proc$rotated[,,8]
#y normieren, um Fehler gering zu halten
p<-y/norm_arr(y)

#Erhalte v als 2D array im Tangentialraum von x
#Wende Logarithmus einmal auf y und einmal auf y normiert (p) an.
v<-loga(x,p)
w<-loga(x,y)
v-w#Im Bereich 10^-16, also ist die Normierung vernachlässigbar.

#Hauptproblem:
#Wenn man jedoch erneut expo anwendet scheint die Normierung einen signifikanten Unterschied zu machen.
p-expo(x,v)#Im Bereich 10^-15
y-expo(x,w)#Im Bereich 10^-4


#Vergleiche v mit v_orp
v-v_orp
w-v_orp#Beides Abweichung im Bereich 10^-6
#Dies wäre also wenn alles richtig implementiert und definiert wurde,
#der Fehler der orthogonalen Projektion im Vergleich zum riemannschen Logarithmus.

(x+loga(x,p))-proc$orpdata[,,8]#10^-5
#Abweichung der orthogonalen Projektion des shape in den Tangentialraum.


plotshapes(gorf.dat[,,3:15])
plotshapes(proc$rotated[,,3:15], color = 3)


#Folgende plots gilt es zu vergleichen:
plotshapes(proc$orpdata[,,3])
plotshapes(x+loga(x,proc$rotated[,,3]))
#lineplot?
#Wie kann ich diese zwei Plots in einem Schaubild plotten?



#Gibt es einen Befehl, der mich direkt zu der Kovarianzmatrix führt, 
#oder muss man diese indirekt berechnen mit den EW und PCs?
