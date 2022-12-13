# Author name: Yi-Ying Chen 
# Email: yiyingchen@gate.sinica.edu.tw
# 
# load R library 
library("sp")
library("raster")
library("rgdal") #readOGR
library("rgeos") #gCentroid
library("viridis") 
library("proj4")
library("snow")
#library("tidyverse")

ii = 0




# load finish net at 500m by 500m spacing 

mesh.500m = readOGR(verbose = FALSE, 
            "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/fishnet/mesh/500m/taiwan_raster_t97.shp")
#mesh.500m = readOGR(verbose = FALSE, 
#            "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/fishnet/mesh/500m/taiwan_raster_WGS84.shp")
#convert the projection into twd97
#mesh.500m =  spTransform(mesh.500m, sp::CRS("+init=epsg:3826"))  

mesh.12km = readOGR(verbose = FALSE,
             "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/fishnet/mesh/12km/MESH_Taiwan.shp")
#convert the projection into twd97
mesh.12km =  spTransform(mesh.12km, sp::CRS("+init=epsg:3826"))  

# get gelocation of center point  of the grid.
cent.xy.500m = gCentroid(mesh.500m, byid=TRUE)
# create dataframe for xy coordinate
df.500m <- as.data.frame(cent.xy.500m)
# Using add_column()
df.500m.share <- data.frame(x=df.500m$x, y=df.500m$y,forest = 0, agri=0, water=0, built=0, other=0  )
cent.xy.12km = gCentroid(mesh.12km, byid=TRUE)

#set image path 
#image_path = c("/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/vivian_code/results/2016/all2016")
image_path = c("/work/vivianlin0921/PCA/classificationresults/taiwan/2021/all/2021")
image_subname = c("_taiwanclassification.tif")

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
  my.col.5c <- c("#439c6e","#e86d5f","#8bd2e8","#f0d86e","#999999") 
      #dark green for forest, red for builtup, blue for water, orange for agri, grey for unknown   
  my.col.forest<- colorRampPalette(c("gray","lightgreen","#439c6e"))(101)
  #my.col <- viridis(n=259, alpha=1, direction=1, option="H") # D for viridus  H for turbo/rainbow



df.500m <- data.frame()
wrk_img <- list()
ref.pt <- list()
df.500m.share.wrk <- list()

#for (imesh in 1:259) {
for (imesh in 1:20)  {

print(paste("Working on imesh:",imesh,".", sep=" "))
xid=mesh.12km@data$XID[imesh]
yid=mesh.12km@data$YID[imesh]

print( paste("XID:", xid, " YID:", yid, sep="") )


#set xy id for the wrking image 
xid <- formatC(xid,format="s",width=3)
yid <- formatC(yid,format="s",width=3) 
img_path <- c( paste(image_path,"_",xid,"_",yid, image_subname,sep="") )

  print(paste("Working on the image file:",img_path,sep=""))

# check classificed image file
    if( file.exists(img_path) ) {
 
    #initiated the ii index based on image 
    ii = ii + 1 
    wrk_img[[ii]] <- raster(x = img_path)
 
   } else{
     print(paste("can't find the image file:",img_path,sep=""))
     ##---exist the loop--- for next iteration 
     next
   }

# test to show the points in the mesh
brd=200
set_buf = 1
#plot(mesh.12km, axes=TRUE, xlim=c(30.0e+4,31.20e+4), ylim=c(27.55e+5,27.67e+5)) 
plot(mesh.12km, axes=TRUE, xlim=c(wrk_img[[ii]]@extent@xmin-brd,wrk_img[[ii]]@extent@xmax+brd), 
                           ylim=c(wrk_img[[ii]]@extent@ymin-brd,wrk_img[[ii]]@extent@ymax+brd))
plot(wrk_img[[ii]], col=my.col.5c, add=T , zlim=c(1,5))

# find id within the image
g_xmax <- wrk_img[[ii]]@extent@xmax
g_xmin <- wrk_img[[ii]]@extent@xmin
g_ymax <- wrk_img[[ii]]@extent@ymax
g_ymin <- wrk_img[[ii]]@extent@ymin

#     #do nothing and go to next interation
     if (length(df.500m.share$x) == 0 ){
     #reset ii =ii -1
     ii = ii-1 
     print("no referenc grid @500m, reset index")
     next
    }

df.500m.share.wrk[[ii]] <- subset( df.500m.share, (df.500m.share$x < g_xmax & df.500m.share$x > g_xmin & df.500m.share$y < g_ymax & df.500m.share$y > g_ymin))

     if (length(df.500m.share.wrk[[ii]]$x) == 0 ) {
      #reset ii =ii -1
      ii = ii-1 
      print("no reference points @6m, reset index")
      next
     }

print(paste("wrking images:",ii, sep="") )
#print(df.500m.share.wrk[[ii]])

  # x.wrk <- df.500m.share.wrk[[img]]$x[igd]
    # y.wrk <- df.500m.share.wrk[[img]]$y[igd]  
   
#create spatial point obj
ref.pt[[ii]] <- SpatialPoints(coords=cbind(df.500m.share.wrk[[ii]]$x,df.500m.share.wrk[[ii]]$y), proj4string= sp::CRS("+init=epsg:3826")  ) 
#
tot.pt <- length(ref.pt[[ii]])
#

ld.do  <- FALSE

if (ld.do) {
print(paste("singal cpu start: ", Sys.time(),sep="")) 

  for ( j in 1: length(ref.pt[[ii]]) ) {
  #  print(paste("progresssing: ", formatC((i/tot.pt)*100, digits=1, width = 4, format = "fg"),
  #               "%",sep=""))
  #set buffer distance as 250m from the center :wq!
  #      pp <- extract( wrk_img[[ii]],  ref.pt[[ii]][j] , buffer=250)
  #      pp <- extract( wrk_img[[ii]],  ref.pt[[ii]][j] , buffer=set_buf)
         pp <- extract( wrk_img[[ii]],  ref.pt[[ii]][j] )
  #classify as forest=1, builtup=2, water=3, agri=4, unknown=5
  ##forest type
  tot.n <- length(pp[[1]])
  #forest 
  df.500m.share.wrk[[ii]]$forest[j]= length(which(pp[[1]]==1))/ tot.n
    if (is.na(df.500m.share.wrk[[ii]]$forest[j])) df.500m.share.wrk[[ii]]$forest[j]=0
  #builtup
  df.500m.share.wrk[[ii]]$built[j] = length(which(pp[[1]]==2))/ tot.n
    if (is.na(df.500m.share.wrk[[ii]]$built[j])) df.500m.share.wrk[[ii]]$built[j]=0
  #water 
   df.500m.share.wrk[[ii]]$water[j] = length(which(pp[[1]]==3))/ tot.n
    if (is.na(df.500m.share.wrk[[ii]]$water[j])) df.500m.share.wrk[[ii]]$water[j]=0
  #agri
  df.500m.share.wrk[[ii]]$agri[j]  = length(which(pp[[1]]==4))/ tot.n
    if (is.na(df.500m.share.wrk[[ii]]$agri[j])) df.500m.share.wrk[[ii]]$agri[j]=0
  #unknown
  df.500m.share.wrk[[ii]]$other[j] = length(which(pp[[1]]>=5))/ tot.n
    if (is.na(df.500m.share.wrk[[ii]]$other[j])) df.500m.share.wrk[[ii]]$other[j]=0
 
  #combine other and forest 
   #df.500m.share.gd$forest[i] = df.500m.share.gd$other[i] + df.500m.share.gd$forest[i]+df.500m.share.gd$agri[i]
   print(paste("sampling number for each grid:", tot.n, sep=" "))

  } #end for j 

#### combine table to the whole island ###
#update the share percentage in the original dataframe "df.500m.share" by the subset table "df.500m.share.gd"
#
#    for (igd in 1: length(df.500m.share.wrk[[ii]]$x) ) {
#   x.wrk <- df.500m.share.wrk[[ii]]$x[igd]
#   y.wrk <- df.500m.share.wrk[[ii]]$y[igd]  
   #replace the value for the grid  
#   df.500m.share[(df.500m.share$x== x.wrk & df.500m.share$y== y.wrk), ] <- df.500m.share.wrk[[ii]][igd,]
   #print(df.500m.share.gd[igd,])
#   }

   #ouput result of working image  
   print(df.500m.share.wrk[[ii]])
  #


} # end if ld.do 
} # end imesh



#save/ouput the df.500m.share table 
save(df.500m.share, file = "data.frame.500m.share.rda")


#dev.new()
#dev.off()
print(paste("singal cpu end: ", Sys.time(),sep=" ")) 


ld.do <-  TRUE

if (ld.do) {  
#library("snowfall")
# Now, create a R cluster using all the machine cores minus one
#sfInit(parallel=TRUE, cpus=parallel:::detectCores()-1)
#sfInit(parallel=TRUE, cpus=2)


# Load the required packages inside the cluster
#sfLibrary(raster)
#sfLibrary(sp)

# Run parallelized 'extract' function and stop cluster
#e.df <- sfSapply(s.list, extract, y=pts)
#sfStop()
#pp_c <- sfRapply(wrk_img, fun=raster::extract, MoreArgs=list(ref.pt, buffer=set_buf) )
  
#sfStop()
 
#assign cpu numbers regarding to (ii) numbers of images 
cluster_n = ii
print(paste("cluster cpus", cluster_n,"start:", Sys.time(),sep=" "))
# start to stack the wrking images & reference points
# do R cluster base on the number of assign cpus max set to 10
 beginCluster(n=cluster_n)
    #using mapply to do the extract funtion with  image/raster list and point list, and buffer condition 
 #   pp_c <- mapply(wrk_img, FUN=raster::extract, ref.pt, buffer=set_buf )
     pp_c <- mapply(wrk_img, FUN=raster::extract, ref.pt )
 endCluster()

   #working on calculating the share fraction
    for (img in 1:cluster_n) {
    #classify as forest=1, builtup=2, water=3, agri=4, unknown=5
    
    # point id within a img get point information from ref.pt list    
    for (ipt in 1: length(ref.pt[[img]])) {
    
    ## total sample n from the reference point 
     tot.n <- length(!is.na((pp_c[[img]][ipt])))
     #print(paste("total sampling pixels for each grid:",tot.n,sep=""))
    #forest 
     df.500m.share.wrk[[img]]$forest[ipt]= length(which(unlist(pp_c[[img]][ipt])==1))/ tot.n
     if (is.na(df.500m.share.wrk[[img]]$forest[ipt])) df.500m.share.wrk[[img]]$forest[ipt]=0
    #builtup
     df.500m.share.wrk[[img]]$built[ipt] = length(which(unlist(pp_c[[img]][ipt])==2))/ tot.n
     if (is.na(df.500m.share.wrk[[img]]$built[ipt])) df.500m.share.wrk[[img]]$built[ipt]=0
    #water 
     df.500m.share.wrk[[img]]$water[ipt] = length(which(unlist(pp_c[[img]][ipt])==3))/ tot.n
     if (is.na(df.500m.share.wrk[[img]]$water[ipt])) df.500m.share.wrk[[img]]$water[ipt]=0
    #agri
     df.500m.share.wrk[[img]]$agri[ipt]  = length(which(unlist(pp_c[[img]][ipt])==4))/ tot.n
     if (is.na(df.500m.share.wrk[[img]]$agri[ipt])) df.500m.share.wrk[[img]]$agri[ipt]=0
    #unknown
     df.500m.share.wrk[[img]]$other[ipt] = length(which(unlist(pp_c[[img]][ipt])>=5))/ tot.n
     if (is.na(df.500m.share.wrk[[img]]$other[ipt])) df.500m.share.wrk[[img]]$other[ipt]=0
  
        x.wrk <- df.500m.share.wrk[[img]]$x[ipt]
        y.wrk <- df.500m.share.wrk[[img]]$y[ipt]  
   
    # print(paste("working on grid:", "x:", x.wrk, "y:", y.wrk,sep=" "))
 
     #replace the value for the grid  
     df.500m.share[(df.500m.share$x== x.wrk & df.500m.share$y== y.wrk), ] <- unlist(df.500m.share.wrk[[img]][ipt,])
    
      } #end for ipt

     # update vshare percentage for each grid 
    for (igd in 1: length(df.500m.share.wrk[[img]]$x) ) {
     x.wrk <- df.500m.share.wrk[[img]]$x[igd]
     y.wrk <- df.500m.share.wrk[[img]]$y[igd]  
     #replace the value for the grid  
     df.500m.share[(df.500m.share$x== x.wrk & df.500m.share$y== y.wrk), ] <- unlist(df.500m.share.wrk[[img]][igd,])
    }
    print(df.500m.share.wrk[[img]])
 
    #update the data.frame 
    df.500m <- rbind(df.500m, df.500m.share.wrk[[img]])

    } #end for img     
} #end if ld.do 

#print(paste("cluster cpus", cluster_n,"end:", Sys.time(),sep=" "))
 
#conver dataframe to raster object 
raster.lulcc <- rasterFromXYZ(df.500m, crs=sp::CRS("+init=epsg:3826"))
#save the ratser as netCDf file
writeRaster(raster.lulcc, "tt_lulcc_500m_twd97.nc", overwrite=TRUE, format="CDF", varname="LU", varunit="fraction", 
        longname="Landuse/cover type derived from SPOT images 6m and upscale to 500m grid, Forest=1, Builtup=2, Water=3, Agri=4, Unkn.=5 ",
        xname="x", yname="y", zname="Coverage")



#create raster object based on the data.frame
#spg <- df.500m.share
#coordinates(spg) <- ~ x + y
#gridded(spg) <- TRUE
#raster.lulcc.500m_1 <- raster(spg)
#crs(raster.lulcc.500m_1) <- CRS("+init=epsg:3826")
ld.do <- TRUE

if( ld.do ) {

#conver dataframe to raster object 
raster.lulcc.500m <- rasterFromXYZ(df.500m.share, crs=sp::CRS("+init=epsg:3826"))
#save the ratser as netCDf file
writeRaster(raster.lulcc.500m, "test_lulcc_500m_twd97.nc", overwrite=TRUE, format="CDF", varname="LU", varunit="fraction", 
        longname="Landuse/cover type derived from SPOT images 6m and upscale to 500m grid, Forest=1, Builtup=2, Water=3, Agri=4, Unkn.=5 ",
        xname="x", yname="y", zname="Coverage")

#reproject to WGS84 but the grid spacing will be exactly the same
#raster.lulcc.500m.wgs84 <- projectExtent(raster.lulcc.500m, crs=as.character(CRS("+init=epsg:4326")) ) 
#adjust resolution to have the same dimension for (747 x 900) as the previous one 
#res(raster.lulcc.500m.wgs84) <- c(0.00499,0.00454)
#reproject the raster 
#raster.lulcc.500m.wgs84 <- projectRaster(raster.lulcc.500m,  res=c(0.00499,0.00454),  crs=as.character(CRS("+init=epsg:4326")))
raster.lulcc.500m.wgs84 <- projectRaster(raster.lulcc.500m,  res=c(0.00506,0.004598),  crs=as.character(CRS("+init=epsg:4326")))
#write raster file
writeRaster(raster.lulcc.500m.wgs84, "test_lulcc_500m_wgs84.nc", overwrite=TRUE, format="CDF", varname="LU", varunit="fraction",
        longname="Landuse/cover type derived from SPOT images 6m and upscale to 500m grid, Forest=1, Builtup=2, Water=3, Agri=4, Unkn.=5 ",
        xname="Longitude", yname="Latitude", zname="Coverage")
}

#convert the projection to WGS84: EPSG:4326

# Source data
#xy <- data.frame(x=df.500m.share$x, y=df.500m.share$y)

# Transform coordinate from TM2 xy data back to WGS84
#proj4text <- as.character(CRS("+init=epsg:3826"))
#pj <- proj4::project(xy, proj=proj4text, inverse=TRUE)
#df.500m.share$x=pj$x
#df.500m.share$y=pj$y



#create raster object based on the data.frame
#spg <- df.500m.share
#coordinates(spg) <- ~ x + y
#gridded(spg) <- TRUE
#raster.lulcc.500m.wgs84 <- raster(spg)
#crs(raster.lulcc.500m.wgs84) <- CRS("+init=epsg:4326")

#mesh.500m.wgs84 = readOGR(verbose = FALSE, 
#            "/lfs/home/ychen/scripts/R/Rscripts/SPOT_CLASS/fishnet/mesh/500m/taiwan_raster_WGS84.shp")

#convert the projection into twd97
#mesh.500m.wgs84 =  spTransform(mesh.500m.wgs84, sp::CRS("+init=epsg:4326"))  

#cent.lalo.500m = gCentroid(mesh.500m.wgs84, byid=TRUE)
# create dataframe for xy coordinate
#df.500m.lalo <- as.data.frame(cent.lalo.500m)

# Using add_column()
#df.500m.share <- data.frame(x=df.500m.lalo$x, y=df.500m.lalo$y,forest = 0, agri=0, water=0, built=0, other=0  )

# replace the x y to lalo 
#df.500m.share$x=df.500.lalo$x
#df.500m.share$y=df.500.lalo$y

#raster.lulcc.500m.wgs84 <- rasterFromXYZ(df.500m.share, crs=sp::CRS("+init=epsg:4326")) 

#save the ratser as netCDf file
#writeRaster(raster.lulcc.500m.wgs84, "lulcc_500m_wgs84.nc", overwrite=TRUE, format="CDF", varname="LU", varunit="fraction", 
#        longname="Landuse/cover type derived from SPOT images 6m and upscale to 500m grid, Forest=1, Builtup=2, Water=3, Agri=4, Unkn.=5 ",
#        xname="Longitude", yname="Latitude", zname="Coverage")


#overlap the plot 
#my.col<- colorRampPalette(c("white","forestgreen"))(100)
#my.col <- viridis(n=100, alpha=1, direction=1, option="H") # D for viridus  H for turbo/rainbow


#points(x= df.500m.share.gd$x,  df.500m.share.gd$y,pch=22, bg=my.col[round(df.500m.share.gd$built*100)],cex=0.5, col=NA)

#plot(lulcc_236_209, col=my.col.5c, add=T )



#plot( lulcc_235_209, add=T)
#plot( lulcc_235_210, add=T)
#plot( lulcc_236_209, add=T)
#plot( lulcc_236_210, add=T)



