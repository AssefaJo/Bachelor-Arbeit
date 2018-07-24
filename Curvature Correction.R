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

v_orp - loga(x,expo(x,v_orp))
y - expo(x,loga(x,y))
#Identität Abweichung 10^-16

(x+loga(x,y))-proc$orpdata[,,8]#10^-5
#Abweichung des mit log berechneten, exakten shapes im Tangentialraum von der othogonalen Projektion im Tangentialraum.


plotshapes(gorf.dat[,,3:15])
plotshapes(proc$rotated[,,3:15], color = 3)
deformGrid2d(proc$orpdata[,,3],x+loga(x,proc$rotated[,,3]),wireframe = c(1,6:2,8:6))
deformGrid2d(proc$mshape,proc$rotated[,,1],wireframe = c(1,6:2,8:6))

#Folgende plots gilt es zu vergleichen:
plotshapes(proc$orpdata[,,3])
plotshapes(x+loga(x,proc$rotated[,,3]))



#Simulation von shapes im Tangentialraum:
N=11
#Generiere array e mit N standardnormalverteilten 12-dimensionalen Vektoren
e<-array(0,c(12,N))

for(i in 1:N){e[,i]<-as.vector(rnorm(12,mean=0,sd=1))}

#12x12 Kovarianzmatrix
cov<-diag(rep(0.04,times=12))

#E ist zu N(0,cov*t(cov)) verteilt
E<-array(0,c(12,N))
for(i in 1:N){E[,i]<-cov%*%e[,i]}

E
#Basis b von T_xM mit PCs:
b<-proc$PCs

#Erhalte nun einen pre shape X im Tangentialraum aus Linearkombination:
#Erstelle leeres 16 dim. array
w<-array(0,c(16,N))
t<-array(0,c(16,N))

#t<-(b_1 * E_11 + .. + b_m* E_m1)
for(j in 1:N){for(i in 1:12){t[,j]<-t[,j]+(b[,i]*E[i,j])}}

#Unsere generierten Vektoren sind sogar orthogonal zum mean-shape.
#t(t[,5])%*%as.vector(x)#10^-15
z
z <- array(0,c(8,2,N))
#10 shapes auf meinem pre shape space
for(i in 1:N){z[,,i]<-expo(x,t[,i])}

#Vergleiche die plots der generierten shapes mit den Beispielashapes
plotshapes(z)
plotshapes(proc$rotated[,,1:N])
#Wie plotte ich das in einem Fenster?


#Bestimmung des Fehlers der orthogonalen Projektion mittels simulierten Daten:

#Der mit log bestimmte exakten Wert von z im Tangentialraum. 
z_t <- array(0,c(8,2,N))
for(i in 1:N){z_t[,,i]<-loga(x,z[,,i])}

#Falls die Länge meiner Vektoren länger als pi/2=1.57 wird, macht dies keinen Sinn mehr.
for(i in 1:N){message(norm_vec(t[,i]))}

#Damit folgt für den shape im Tangentialraum:
z_exakt<-array(0,c(8,2,N))
for(i in 1:N){z_exakt[,,i]<-z_t[,,i]+x}

#Projiziere z orthogonal in den Tangentialraum von x
z_orp <- Morpho:::orp(z, mshape =x)

#In den zwei plots sieht man den Unterschied von z_orp und z
plotshapes(z_orp)
plotshapes(z_exakt)

#Der Abstand von z_orp und dem exakten Wert des shapes im Tangentialraum.
for(i in 1:N){message(norm_arr(z_orp[,,i]-z_exakt[,,i]))}


#Damit wäre die Abweichung der orthogonalen Projektion gezeigt. 
#Dies kann man für unterschiedlichstes cov ausführen und somit 
#stärkere Abweichungen für größere Varianz und schwäche Abweichungen für geringere Varianz feststellen.

cov_e<-matrix(0,nrow = 16,ncol = 16)
cov_o<-matrix(0,nrow = 16,ncol = 16)

# mit den simulierten Daten erhalten wir eine 16x16 kovarianzmatrix.
#Unsere ursprüngliche zu schätzende kovarianzmatrix ist aber 12x12. 
#?

#Kovarianzmatrix von den exakten und orthogonalen shapes
for(i in 1:10){cov_e<-cov_e+(1/10)*(as.vector(z_exakt[,,i]-x)%*%t(as.vector(z_exakt[,,i]-x)))}
for(i in 1:10){cov_o<-cov_o+(1/10)*(as.vector(z_orp[,,i]-x)%*%t(as.vector(z_orp[,,i]-x)))}
cov%*%t(cov)
#Spectral Composition
EDe<-eigen(cov_e)
EDo<-eigen(cov_o)

#Untersuchung der Eigenwerte
EDe$values
EDo$values
