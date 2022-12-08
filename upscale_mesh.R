# Author name: Yi-Ying Chen 
# Email: yiyingchen@gate.sinica.edu.tw
# 
# load R library 
library("sp")
library("rgdal") #readOGR
library("rgeos") #gCentroid
library("viridis") 
# load finish net at 500m by 500m spacing 

mesh.500m = readOGR(verbose = FALSE, 
            "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/fishnet/mesh/500m/taiwan_raster_t97.shp")
mesh.12km = readOGR(verbose = FALSE,
             "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/fishnet/mesh/12km/MESH_Taiwan.shp")
#convert the projection into twd97
mesh.12km =  spTransform(mesh.12km, sp::CRS("+init=epsg:3826"))  

# get gelocation of center point  of the grid.

cent.xy.500m = gCentroid(mesh.500m, byid=TRUE)
# create dataframe for xy coordinate
df.500m <- as.data.frame(cent.xy.500m)

cent.xy.12km = gCentroid(mesh.12km, byid=TRUE)


# test to show the points in the mesh
plot(mesh.12km, axes=TRUE, xlim=c(27.0e+4,32.0e+4), ylim=c(27.5e+5,28.0e+5)) 

#my.color<- colorRampPalette(c("brown","yellow","gray","gray","green","forestgreen"))(16)
  my.col <- viridis(n=259, alpha=1, direction=1, option="H") # D for viridus  H for turbo/rainbow


for (i in 1:259) {
  # the grid box xy
  #   pt4    pt3
  #  
  ##  pt1/5    pt2   
  xmax=max(mesh.12km@polygons[[i]]@Polygons[[1]]@coords[,1])
  xmin=min(mesh.12km@polygons[[i]]@Polygons[[1]]@coords[,1])
  #
  ymax=max(mesh.12km@polygons[[i]]@Polygons[[1]]@coords[,2])
  ymin=min(mesh.12km@polygons[[i]]@Polygons[[1]]@coords[,2])
  #
  print(paste("xmin:",xmin, "xmax:",xmax, sep=" "))
  print(paste("ymin:",ymin, "ymax;",ymax, sep=" ")) 
 
  gd.pts <- subset (df.500m, (df.500m$x > xmin & df.500m$x < xmax & df.500m$y > ymin & df.500m$y < ymax )) 
  print(paste("Total grid @500m: ", length(gd.pts$x),sep=""))
  points( gd.pts, bg=my.col[i],pch=22, col=NA,cex=1.0)

}
#plot 500 points
points (x=df.500m$x, y=df.500m$y, cex=0.3, pch=16)
print(paste("Total points:", length(df.500m$y), sep="")) 



#read geo tif 
lulcc_236_209 <- raster(x = "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/vivian_code/results/2016/all2016_236_209_taiwanclassification.tif")
lulcc_236_210 <- raster(x = "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/vivian_code/results/2016/all2016_236_210_taiwanclassification.tif")
lulcc_235_209 <- raster(x = "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/vivian_code/results/2016/all2016_235_209_taiwanclassification.tif")
lulcc_235_210 <- raster(x = "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/vivian_code/results/2016/all2016_235_210_taiwanclassification.tif")


#overlap the plot 


plot( lulcc_235_209, add=T)
plot( lulcc_235_210, add=T)
plot( lulcc_236_209, add=T)
plot( lulcc_236_210, add=T)



