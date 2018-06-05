library(geomorph)
library(shapes)
library(Morpho)

norm_vec <- function(x) sqrt(sum(x^2))

#>Gibt die euklidische Norm eines 2d-arrays zurück
norm_arr <- function(x) sqrt(sum(diag(t(x)%*%x)))

#Exponentialabbilung
expo <- function(x,v){
  
  k <- nrow(x)
  m <- ncol(x)
  
  #x und v zu Vektor umschreiben
  x_vec <- as.vector(t(x))
  v_vec <- as.vector(t(v))
  
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
x <- proc$mshape

#Beispielshape
y <- proc$rotated[,,8]
y_orp <- proc$orpdata[,,8]



#Beispielrechnungen:

#Wäre das nun ein Vektor im Tangentialraum vom Meanshape?
#der Mean shape selbst subtrahiert 
#von der orthogonalen Projektion eines shape
#(quasi der verbindungsvektor von Projektion und mean shape)
v_orp <- y_orp-x

v_orp - loga(x,expo(x,v_orp))
#Die identität hat eine Abweichung im Bereich 10^-17. Scheint also richtig zu sein.

norm_arr(proc$rotated[,,8])
#Wieso ist die Norm von meinen Shapes nicht 1? Skalierungsinvarianz?
#Oder geht man hier von dem Procrustes distance aus?

(x+loga(x,y))-proc$orpdata[,,8]#10^-5
#Abweichung des mit log berechneten, exakten shapes im Tangentialraum von der othogonalen Projektion im Tangentialraum.



#Fragen zu den plots:

#Wie kann ich folgende zwei Plots in einem Schaubild plotten?
plotshapes(gorf.dat[,,3:15])
plotshapes(proc$rotated[,,3:15], color = 3)

#lines(c(-200,250),c(0,0))
#lines(c(0,0),c(-100,300))

#Folgende plots gilt es zu vergleichen:
plotshapes(proc$orpdata[,,3])
plotshapes(x+loga(x,proc$rotated[,,3]))
#Wie kann ich diese zwei Plots in einem Schaubild plotten? Mit linien.
#lineplot?



#Kovarianz und PCs:

proc$PCs
#Liefern mir die PCs eine Basis im Tangentialraum meines mean shape?
t(proc$PCs[,2])%*%as.vector(x)#skalarprodukt von mean shape und PC. 10^-16
#mean shape und PC stehen also orthogonal zueinander.
#Also PCs Basis von Tangentialraum(?).

#Gibt es einen Befehl, der mich direkt zu der Kovarianzmatrix von gorf.dat führt, 
#oder muss man diese indirekt berechnen mit den EW und PCs?




#Simulation von shapes im Tangentialraum (mit Annahme PCs bilden Basis vom Tangentialraum):

#Generiere array e mit 10 normalverteilten 12-dimensionalen Vektoren, die verteilt sind zu N(0,id).
e<-array(0,c(12,10))
for(i in 1:10){e[,i]<-as.vector(rnorm(12,mean=0,sd=1))}

#12x12 Kovarianzmatrix mit Bsp.: 0.5 auf Diagonalen
cov<-diag(rep(0.03,times=12))

#Multipliziere jeden generierten Zufallsvektor aus e mit cov und erhalte neue Vektoren 
#die verteilt sind zu N(0,cov*t(cov))
E<-array(0,c(12,10))
for(i in 1:10){E[,i]<-cov%*%e[,i]}

#Die wichtigste Frage bleibt: Bilden die PCs aus gorf.dat eine Basis des Tangentialraums?
#Basis b von T_xM mit PCs:
b<-proc$PCs

#Erhalte nun einen shape X im Tangentialraum aus folgender Linearkombination:
# X_1 = x+(b_1 * E_11 + .. + b_m* E_m1)
#x ist der mean shape.
#b_i ist ein Basisvektor, also die i-te Spalte von b.
#E_ij ist der i-j-te Eintrag der Matrix E. Also nehmen wir für einen shape einen Spaltenvektor aus E
#und jeder Eintrag ist ein Koeffizient für den zugehörigen Basisvektor (-> Harms).

#Erstelle leeres 16 dim. array
w<-array(0,c(16,10))
t<-array(0,c(16,10))

#t<-(b_1 * E_11 + .. + b_m* E_m1)
for(j in 1:10){for(i in 1:12){t[,j]<-t[,j]+(b[,i]*E[i,j])}}

#Unsere generierten Vektoren sind sogar orthogonal zum mean-shape.
t(t[,1])%*%as.vector(x)#10^-15

#Ich habe Angst hier durcheinander zu kommen mit der Vektornotation von shapes und der Matrixnotation. 
#Da ich oft die Matrix umschrieben muss als Vektor, aber ich mir nicht sicher bin ob zeilenweise oder spaltenweise. Tipp?

#Addiere mean shape
for(i in 1:10){w[,i]<-x+t[,i]}

#Erhalte 10 shapes im Tangentialraum. Tada.
w

z<-array(0,c(8,2,10))
#Und damit 10 shapes auf meinem shape space
for(i in 1:10){z[,,i]<-expo(x,t[,i])}

#Wenn alles richtig ist habe ich nun mittels der PCs aus gorf.dat und der neuen Kovarianzmatrix neue shapes simuliert.
#Sieht sogar so ähnlich aus wie gorf.dat. 
plotshapes(z)
plotshapes(proc$rotated[,,1:10])
#Umso geringer die Varianzen, desto genauer werden die shapes.
#Wieder die Frage: Wie plotte ich das in einem Fenster?

#Nun gilt es mit den simulierten Daten den Fehler der orthogonalen Projektion zu bestimmen.
#Wie bestimme ich nun die orthogonale Projektion von z? 
#....

v_orp



#Ignorieren:
#plot(c(1,1,2,1),c(1,1.5,1,1),pch=20,col="blue",xlim = c(-1,2),ylim = c(-1,2),xlab = "x",ylab = "y")
#lines(c(1,1,2,1),c(1,1.5,1,1))
#lines(c(0,0),c(-2,3))
#lines(c(-2,3),c(0,0))
#points(c(-1/3,-1/3,2/3,-1/3),c(1-(3.5/3),1.5-(3.5/3),1-(3.5/3),1-(3.5/3)),pch=20,col="red")
#lines(c(-1/3,-1/3,2/3,-1/3),c(1-(3.5/3),1.5-(3.5/3),1-(3.5/3),1-(3.5/3)))
#arrows(0.9,0.9,0.3,0.3,angle = 30,length=0.18)