# Author name: Yi-Ying Chen 
# Email: yiyingchen@gate.sinica.edu.tw
# 
# load R library 
library("sp")
library("raster")
library("rgdal") #readOGR
library("rgeos") #gCentroid
library("viridis") 
#library("tidyverse")

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

# Using add_column()
df.500m.share <- data.frame(x=df.500m$x, y=df.500m$y,forest = 100, agri=0, water=0, built=0, other=0  )


cent.xy.12km = gCentroid(mesh.12km, byid=TRUE)


# worlk on the mesh tables

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
  print(paste("Total 500m-grid mech for each SPOT image : ", length(gd.pts$x),sep=""))
  #points( gd.pts, bg=my.col[i],pch=22, col=NA,cex=1.0)

}

#plot 500 points
#points (x=df.500m$x, y=df.500m$y, cex=0.3, pch=16)
print(paste("Total points:", length(df.500m$y), sep="")) 

# set color palette
  my.col.5c <- c("#439c6e","#e86d5f","#8bd2e8","#f0d86e","#999999") #dark green for forest, red for builtup, blue for water, orange for agri, grey for unknown   
  my.col.forest<- colorRampPalette(c("gray","lightgreen","#439c6e"))(101)
  #my.col <- viridis(n=259, alpha=1, direction=1, option="H") # D for viridus  H for turbo/rainbow

#set xy id for the wrking image 
xid <- formatC("236",format="s",width=3)
yid <- formatC("210",format="s",width=3) 

wrk_img <- raster(x = paste("/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/vivian_code/results/2016/all2016","_",xid,"_",yid,"_taiwanclassification.tif",sep=""))
# test to show the points in the mesh
brd=200
#plot(mesh.12km, axes=TRUE, xlim=c(30.0e+4,31.20e+4), ylim=c(27.55e+5,27.67e+5)) 
plot(mesh.12km, axes=TRUE, xlim=c(wrk_img@extent@xmin-brd,wrk_img@extent@xmax+brd), 
                           ylim=c(wrk_img@extent@ymin-brd,wrk_img@extent@ymax+brd))
plot(wrk_img, col=my.col.5c, add=T , zlim=c(1,5))

# find id within the image
g_xmax <- wrk_img@extent@xmax
g_xmin <- wrk_img@extent@xmin
g_ymax <- wrk_img@extent@ymax
g_ymin <- wrk_img@extent@ymin

df.500m.share.gd <- subset( df.500m.share, (df.500m.share$x < g_xmax & df.500m.share$x > g_xmin & df.500m.share$y < g_ymax & df.500m.share$y > g_ymin))

#create spatial point obj
ref.pt <- SpatialPoints(coords=cbind(df.500m.share.gd$x,df.500m.share.gd$y), proj4string= sp::CRS("+init=epsg:3826")  ) 
 
tot.pt <- length(ref.pt) 
for ( i in 1:length(ref.pt) ) {
#  print(paste("progresssing: ", formatC((i/tot.pt)*100, digits=1, width = 4, format = "fg"),
#               "%",sep=""))
  #set buffer distance as 250m from the center 
   pp <- extract( wrk_img,  ref.pt[i] , buffer=200)
  #classify as forest=1, builtup=2, water=3, agri=4, unknown=5
  ##forest type
  tot.n <- length(pp[[1]])
  df.500m.share.gd$forest[i]= length(which(pp[[1]]==1))/ tot.n
    if (is.na(df.500m.share.gd$forest[i])) df.500m.share.gd$forest[i]=0
  df.500m.share.gd$agri[i]  = length(which(pp[[1]]==4))/ tot.n
    if (is.na(df.500m.share.gd$agri[i])) df.500m.share.gd$agri[i]=0
  df.500m.share.gd$water[i] = length(which(pp[[1]]==3))/ tot.n
    if (is.na(df.500m.share.gd$water[i])) df.500m.share.gd$water[i]=0
  df.500m.share.gd$built[i] = length(which(pp[[1]]==2))/ tot.n
    if (is.na(df.500m.share.gd$built[i])) df.500m.share.gd$built[i]=0
  df.500m.share.gd$other[i] = length(which(pp[[1]]>=5))/ tot.n
    if (is.na(df.500m.share.gd$other[i])) df.500m.share.gd$other[i]=0

  #combine other and forest 
   #df.500m.share.gd$forest[i] = df.500m.share.gd$other[i] + df.500m.share.gd$forest[i]+df.500m.share.gd$agri[i]
   #print(pp)

   lulcc_gd <- c(df.500m.share.gd$forest[i], df.500m.share.gd$built[i], df.500m.share.gd$water[i], df.500m.share.gd$agri[i], df.500m.share.gd$other[i])
   lulcc.col <- which.max(lulcc_gd) 
   print( paste("Dominate land type:", lulcc.col, " share percentage", formatC(max(lulcc_gd)*100, digits=1, width = 4, format = "fg"),"%" ,sep="")) 

#output some statistics 
   points(x= df.500m.share.gd$x[i],  df.500m.share.gd$y[i],
     pch=21, bg=my.col.5c[lulcc.col],
     #bg=my.col.forest[ as.integer(ceiling(df.500m.share.gd$forest[i]*100))+1 ],
     cex=2.0, col=NA)
   print(paste("Forest coverage: ", formatC(df.500m.share.gd$forest[i]*100, digits=1, width = 4, format = "fg"),
              "%",sep=""))
   tot.frac = (df.500m.share.gd$forest[i]+ df.500m.share.gd$agri[i]+
               df.500m.share.gd$water[i]+df.500m.share.gd$built[i])*100
   if (tot.frac < 100) { 
       print(paste("Total coverage: ", formatC( tot.frac, digits=1, width = 4, format = "fg"),
              "%",sep="")) }
}




# combine table to the whole island

#overlap the plot 
#my.col<- colorRampPalette(c("white","forestgreen"))(100)
#my.col <- viridis(n=100, alpha=1, direction=1, option="H") # D for viridus  H for turbo/rainbow


#points(x= df.500m.share.gd$x,  df.500m.share.gd$y,pch=22, bg=my.col[round(df.500m.share.gd$built*100)],cex=0.5, col=NA)

#plot(lulcc_236_209, col=my.col.5c, add=T )



#plot( lulcc_235_209, add=T)
#plot( lulcc_235_210, add=T)
#plot( lulcc_236_209, add=T)
#plot( lulcc_236_210, add=T)



