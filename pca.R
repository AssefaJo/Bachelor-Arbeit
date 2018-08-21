#Exponentialabbilung
expo <- function(x,v){
  
  #Norm von v und x
  nv <- sqrt(sum(v^2))
  nx <- sqrt(sum(x^2))
  
  #Berechnung des Pre-Shape auf S^k_m
  pshape <- cos(nv)*x+((sin(nv)*nx)/nv)*v
  
  return(pshape)
}  


#Riemannscher Logarithmus
loga <- function(x,y){
  
  #Norm von y und x
  ny <- sqrt(sum(y^2))
  nx <- sqrt(sum(x^2))
  
  #theta
  t <- acos(abs(sum(x*y))/(nx*ny))
  
  #Projektion von y auf x
  pi <- (x*sum(x*y))/(nx^2)
  
  #Berechnung des Pre-Shape Vektors
  pshape_v <- (t*(y-pi))/sqrt(sum((y-pi)^2))
  
  return(pshape_v)
}

pca <- function(N,dimL=2,Landmarks=3,r=0.2){
  
  if(Landmarks<=dimL) {
    return(message("Die Anzahl der Landmarks muss groesser sein als die Dimension der Landmarks"))
  }
  
  #Dimension des umliegenden Raumes
  d = dimL * (Landmarks-1)
  
  #Dimension des Pre-Shape Space/Tangentialraums
  c = d - 1
  
  #Parameter Sigma, den es zu schaetzen gilt
  Sigma <- diag(rep(r,times=c))
  #Generiere array e mit N standardnormalverteilten c-dimensionalen Vektoren
  e<-array(0,c(c,N))
  
  for(i in 1:N) {
    e[,i] <- as.vector(rnorm(c,mean=0,sd=1))
  }
  
  #Multiplikation mit Sigma
  E<-array(0,c(c,N))
  
  for(i in 1:N) {
    E[,i] <- Sigma %*% e[,i]
  }
  #Mean Shape x als "Nordpol" der Einheitssphaere
  x <- array(c(rep(0,times=c),1))
  
  #In Koordinaten dargestellte orthonormale Basis des Tangentialraums an x
  b <- diag(rep(1,times=d))[,-d]
  
  #Vektoren erzeugen als Linearkombination der Basis b und den Eintraegen von E als Koeffizienten
  v <- array(0,c(d,N))
  
  for(j in 1:N) {
    for(i in 1:c) {
      v[,j] <- v[,j] + (b[,i] * E[i,j])
    }
  }
  
  #Falls die Laenge der Vektoren laenger als pi/2=1.57 wird, macht eine orthogonale Projektion keinen Sinn.
  for(i in 1:N) {
    if(sqrt(sum(v[,i]^2))>pi/2) {
      return(message("Standardabweichung zu gross"))
    }
  }
  
  #Anwendung der Exponentialabbildung liefert Pre-Shapes ps auf M
  ps <- array(0,c(d,N))
  
  for(i in 1:N) {
    ps[,i] <- expo(x,v[,i])
  }
  
  
  #Der mit log bestimmte exakte Wert von ps im Tangentialraum. 
  ps_exakt <- array(0,c(d,N))
  
  for(i in 1:N) {
    ps_exakt[,i] <- loga(x,ps[,i])+x
  }
  
  #Projiziere ps orthogonal in den Tangentialraum von x.
  ps_orp <- ps
  
  for (i in 1:N) {
    ps_orp[d,i] <- 1
  }
  
  #Orthogonal projizierter Vektor
  v_orp <- array(0,c(d-1,N))
  
  for (i in 1:N) {
    v_orp[,i] <- ps_orp[,i][-d]
  }
  
  #Exakt projizierter Vektor
  v_exakt <- array(0,c(d-1,N))
  
  for (i in 1:N) {
    v_exakt[,i] <- ps_exakt[,i][-d]
  }
  
  #Berechnen der Kovarianzmatrix der exakten Daten
  cov_exakt <- matrix(0,nrow = c,ncol = c)
  
  for(i in 1:N) {
    cov_exakt <- cov_exakt + (1/N) * (v_exakt[,i] %*% t(v_exakt[,i]))
  }
  
  
  #Berechnen der Kovarianzmatrix der approximierten Daten
  cov_orp <- matrix(0,nrow = c,ncol = c)
  
  for(i in 1:N) {
    cov_orp <- cov_orp + (1/N) * (v_orp[,i] %*% t(v_orp[,i]))
  }
  
  #Singulaerwertzerlegung
  EDe <- eigen(cov_exakt)
  
  EDo <- eigen(cov_orp)
  
  Sigma_ex <- EDe$vectors %*% diag(sqrt(EDe$values)) %*% t(EDe$vectors)
  Sigma_ap <- EDo$vectors %*% diag(sqrt(EDo$values)) %*% t(EDo$vectors)
  
  D <- Sigma_ex - Sigma_ap
  
  #Abweichung der Schaetzungen
  Delta <- sqrt(max(eigen(t(D) %*% D)$values))
  
  out = (list(Sigma = Sigma , Sigma_ex = Sigma_ex , Sigma_ap = Sigma_ap , Differenz = Delta))
  
  return(out)
}


#Plots:

#4.1
deriv<-seq(0.01,0.30,0.002)
x<-rep(0,times=length(deriv))
for (k in 1:length(deriv)) {
  x[k]<-pca(1000,r = deriv[k])$Differenz
}
plot(deriv,x,xlab="r",ylab=expression(Delta^r))

