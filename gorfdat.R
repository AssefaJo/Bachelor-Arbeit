library(matlib)
library(shapes)
library(Morpho)


data(gorf.dat)

proc <- procSym(gorf.dat,orp=TRUE,scale=F)

x <- proc$mshape

#Folgende zwei plots möchte ich in einem Koordinatensystem darstellen
#deformgrid2d liefert mir die Darstellung von lediglich 2 shapes
plotshapes(proc$rotated[,,1:10],color = "red")
plotshapes(x,color = "green")






#deformGrid2d(proc$orpdata[,,3],x+loga(x,proc$rotated[,,3]),wireframe = c(1,6:2,8:6))
#deformGrid2d(x,proc$rotated[,,1],wireframe = c(1,6:2,8:6))
