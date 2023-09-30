######################################################
#
# Elmira Berjisian,  elmiraberjisian@gmail.com
# Estimate Elevation and Grade from Crowd-sourced GPS data
#
# 
#
########################################################
#########################################Required Packages To Load#########################################
library(rgdal)
library(rgeos)
library(geosphere)
library(maptools)
library(sf)
library(rgeos)
library(sp)
library(tidyr)
library(SimilarityMeasures)
library(TSdist)
library(dbscan)
library(cluster)
library(kmlShape)
library(imputeTS)
library(dplyr)
library(zoo)
library(dtwclust)
library(move)
library(geosphere)
library(signal)
library(matrixStats) 
library(Metrics)
library(pracma)
library(creditmodel)
#########################################Load Crowdsourced Trips and Map-matched Routes Associated with Them#########################################
samplestr<-read.csv("Input/GT_Burrard Bridge.csv") #includes ground truth sampling points, for every point at one meter equidistant along the location we use linear interpolation (line 84)
sampletrajcom<-read.csv("Input/rawtrip_Burrard Bridge.csv") #includes all  Endomondo trips crossing a bridge
highRes<-read.csv("Input/highRes_Burrard Bridge.csv") #includes high resolution elevation estimates (here from Lidar data and using Boyko method to process lidar data)
strtrwgps<-read.csv("Input/RWGPS_Burrard Bridge.csv") #includes RWGPS elevation estimates
nme<-"Burrard Bridge" #name your location for graphs and output generation

###################################################METHOD 1: Eeffort>Filter-Smooth-Clean>prototype>Smooth>smoothed prototype>Grade#########################################
 napts<-samplestr
#step 1 identify all trips that their associated route has the sample location in and detect which records fall into a 50 meters buffer around our sample location

sampletrajcom$inside<-FALSE
s_sampletrajcom<-SpatialPointsDataFrame(coords = sampletrajcom[,c(which(colnames(sampletrajcom)=="longitude"),
                                                                    which(colnames(sampletrajcom)=="latitude"))], data = sampletrajcom,
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
s_sampletrajcom_projected<-spTransform(s_sampletrajcom,CRSobj = CRS("+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  
samplestr<-samplestr[,c("long","lat")]
samplestrpoints <- sp::SpatialPoints(samplestr)
samplestr_line <- as(samplestrpoints,"SpatialLines")
proj4string(samplestr_line)<-"+proj=longlat +datum=WGS84 +no_defs"
samplestr_line_projected<-spTransform(samplestr_line,CRSobj = CRS("+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
gLength(samplestr_line_projected)
samplestrbuffer<-gBuffer(samplestr_line_projected,width = 50,capStyle = "SQUAR")
sampletrajcom$inside[as.numeric(rownames(over(samplestrbuffer, s_sampletrajcom_projected, T)[[1]]))]<-TRUE
#step 5 create efforts from records in the buffer around our sample location
sample_effort<-sampletrajcom[which(sampletrajcom$inside==TRUE),]
sample_effort$x<-proj4::project(sample_effort[,c(which(colnames(sample_effort)=="longitude"),which(colnames(sample_effort)=="latitude"))],"+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=WGS84",degrees = T)$x
sample_effort$y<- proj4::project(sample_effort[,c(which(colnames(sample_effort)=="longitude"),which(colnames(sample_effort)=="latitude"))],"+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=WGS84",degrees = T)$y
sample_effort_single<-split(sample_effort,f=as.factor(sample_effort$workoutID))
length(sample_effort_single)
#step 6 filter out efforts with low number of records
# Find indices of elements with NA altitudes
na_indices <- sapply(sample_effort_single, function(x) any(is.na(x$altitude)))
# Filter out elements with NA altitudes
sample_effort_single[na_indices] <- lapply(sample_effort_single[na_indices], function(x) x[!is.na(x$altitude),])
# Calculate the filter value 
nrows <- sapply(sample_effort_single, function(x) nrow(x))
fltr <- max(round(median(nrows) - 2 * sd(nrows)), 1)
# Find indices of elements with fewer rows than fltr
low_row_indices <- which(nrows <= fltr)
# Remove elements with fewer rows than fltr
if (length(low_row_indices) > 0) {
  sample_effort_single <- sample_effort_single[-low_row_indices]
}

length(sample_effort_single)

#step 7 generates points every meter along our  location and estimate z based on linear interpolation in between where you have GT elevations
gt_str_onem<-gInterpolate(samplestr_line_projected,seq(0,gLength(samplestr_line_projected),1))
sample_gt_elev<-as.data.frame(gt_str_onem@coords)
samplestrpoints <- sp::SpatialPoints(samplestr)
proj4string(samplestrpoints)<-"+proj=longlat +datum=WGS84 +no_defs"
samplestrpoints_projected<-spTransform(samplestrpoints,CRSobj = CRS("+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
sample_gt_elev$zGT<-NA
sample_gt_elev$zGT[napts$along]<-napts$z
sample_gt_elev$zGT<-na.approx(sample_gt_elev$zGT)
sample_gt_elev$highRes<-highRes$zhighRes
if (round(gLength(samplestr_line_projected)/5)%%2 ==0) {sg<-round(gLength(samplestr_line_projected)/5)-1}else{sg<-round(gLength(samplestr_line_projected)/5)}
sample_gt_elev$dist<-gProject(samplestr_line_projected,gt_str_onem)
#smooth our best estimate for elevation from sources above (lidar/dem)
sample_gt_elev$ShighRes<-sgolayfilt(sample_gt_elev$highRes, n= sg)
#step 8 add elevation from each effort assigned to the sampling point along the location
sample_gt_elevpts <-gt_str_onem
for (m in 1:length(sample_effort_single)) {
  # Extract longitude and latitude columns
  lon_lat_cols <- c("longitude", "latitude")
  effort_subset <- sample_effort_single[[m]][, lon_lat_cols]
  
  # Convert to SpatialLines
  effort_line <- as(sp::SpatialPoints(effort_subset),"SpatialLines")
  
  # Create the ID dataframe
  p.df <- data.frame(ID = 1:length(effort_line), row.names = sapply(slot(effort_line, "lines"), function(x) slot(x, "ID")))
  
  # Create SpatialLinesDataFrame
  l <- sp::SpatialLinesDataFrame(effort_line, p.df)
  sp::proj4string(l) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  # Transform to UTM coordinates
  main_l_p <- sp::spTransform(l, CRSobj = sp::CRS("+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  
  # Initialize sample_gt_elev$tmp
  sample_gt_elev$tmp <- NA
  
  # Loop through sample_gt_elev
  for (k in 1:nrow(sample_gt_elev)) {
    sp_project <- snapPointsToLines(sample_gt_elevpts[k, ], main_l_p, maxDist = NA, withAttrs = FALSE, idField = NA)
    
    # Create a dataframe of coordinates
    points_df <- data.frame(coordinates(main_l_p))
    points_df$ID <- seq(1, nrow(sample_effort_single[[m]]))
    
    # Add the snapped point to the dataframe
    points_df <- rbind(points_df, c(sp_project@coords[1, 1], sp_project@coords[1, 2], 999999))
    
    # Create a SpatialPointsDataFrame
    s_pointsdf <- sp::SpatialPointsDataFrame(coords = points_df[, c(1, 2)], data = points_df,
                                             proj4string = sp::CRS("+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
    
    # Calculate distances
    points_df$dist <- gProject(main_l_p, s_pointsdf)
    
    # Order by distance
    points_df <- points_df[order(points_df$dist), ]
    indx <- which(points_df$ID == 999999)
    
    if (indx - 1 != 0 && !is.na(points_df$ID[indx + 1])) {
      prior_id <- points_df$ID[indx - 1]
      post_id <- points_df$ID[indx + 1]
      sample_gt_elev$tmp[k] <- (((sample_effort_single[[m]]$altitude[post_id] - sample_effort_single[[m]]$altitude[prior_id]) /
                                   (points_df$dist[indx + 1] - points_df$dist[indx - 1])) * (points_df$dist[indx] - points_df$dist[indx - 1])) + sample_effort_single[[m]]$altitude[prior_id]
    } else {
      if (indx - 1 == 0) {
        sample_gt_elev$tmp[k] <- sample_effort_single[[m]]$altitude[1]
      } else {
        sample_gt_elev$tmp[k] <- sample_effort_single[[m]]$altitude[nrow(sample_effort_single[[m]])]
      }
    }
  }
  # Set column name for sample_gt_elev
  colnames(sample_gt_elev)[length(colnames(sample_gt_elev))] <- paste("z_eff", m, sep = "")
  print(m)
}
#step 9 adding what MJ paper did with respect to estimate mean and median from elevation obtained above
colnames(sample_gt_elev)[c(which(colnames(sample_gt_elev)=="z_eff1"):which(colnames(sample_gt_elev)==paste("z_eff",length(sample_effort_single),sep = "")))]#make sure this refers to the effort elevation columns
sample_gt_elev$meanZ<-rowMeans(sample_gt_elev[c(which(colnames(sample_gt_elev)=="z_eff1"):which(colnames(sample_gt_elev)==paste("z_eff",length(sample_effort_single),sep = "")))])
sample_gt_elev$meadianZ<-rowMedians(as.matrix(sample_gt_elev[c(which(colnames(sample_gt_elev)=="z_eff1"):which(colnames(sample_gt_elev)==paste("z_eff",length(sample_effort_single),sep = "")))]))
#step 10 where our contribution starts first we detect outliers using DBSCAN
sample_ls <- lapply(sample_gt_elev[, c(which(colnames(sample_gt_elev)=="z_eff1"):which(colnames(sample_gt_elev)==paste("z_eff",length(sample_effort_single),sep = "")))], sgolayfilt, n = sg)

#now we identify and remove outliers (i.e, cleaning step)
  sm<-do.call(cbind,sample_ls)
  # Calculate the Euclidean distance matrix using the dist function
  distm <- dist(t(sm),method = "euclidean")
 
#now let's do dbscan to see if it identifies noise
  
  db<-dbscan::dbscan(distm, eps=45, minPts = 2)
  if(length(which(db$cluster==0))!=0){
  sm<-sm[,-which(db$cluster==0)]
  sample_ls<-sample_ls[-which(db$cluster==0)]}
  length(sample_ls) #shows number of cleaned efforts
#step 11 we add Clustering methods: pam and kml respectively length(sample_ls) #this shows number of filtered and cleaned efforts
smw<-t(sm)
traj_long<-cldsWide(smw,times=seq(1:dim(smw)[2]),id=as.factor(as.character(seq(1:nrow(smw)))))
smitcls2<-kmlShape(traj_long, nbClusters = 1, timeScale =0.1, FrechetSumOrMax = "max", toPlot="both", parAlgo=parKmlShape())
distances_not_in_smitcls2 <- sample_gt_elev$dist[!(round(sample_gt_elev$dist) %in% smitcls2@trajMeans[, 2])]

# Create a data frame with times and traj columns
kmldf <- data.frame(times = smitcls2@trajMeans[, 2], traj = smitcls2@trajMeans[, 3])

# Add rows for distances not in smitcls2@trajMeans[, 2] with NA traj values
new_rows <- data.frame(times = distances_not_in_smitcls2, traj = NA)
kmldf <- rbind(kmldf, new_rows)

# Sort the data frame by times
kmldf <- kmldf[order(kmldf$times),]

# Interpolate NA values in the traj column using spline interpolation
kmldf$traj <- na.interpolation(kmldf$traj, option = "spline")
#step 12 now we create a data frame of all of our estimates for each point along sample loaction
  plot_sampleedf<-sample_gt_elev[,c(which(colnames(sample_gt_elev)=="dist"),
                                    which(colnames(sample_gt_elev)=="zGT"),
                                    which(colnames(sample_gt_elev)=="highRes"),
                                    which(colnames(sample_gt_elev)=="ShighRes"),
                                    which(colnames(sample_gt_elev)=="meanZ"),
                                    which(colnames(sample_gt_elev)=="meadianZ"))
  ] 
  colnames(plot_sampleedf)<-c("dist","GT","highRes","ShighRes","Fmean","Fmedian")
#step 13 add our estimates from Filtered and smoothed efforts
  plot_sampleedf$pam_dbscan<-as.numeric(pam(t(sm),1)$medoids)
  plot_sampleedf$mean_dbscan<-rowMeans(sm)
  plot_sampleedf$median_dbscan<-rowMedians(sm)  

  if(length(kmldf$traj[which(kmldf$times %in%  round(sample_gt_elev$dist))])!=nrow(plot_sampleedf)){
    plot_sampleedf$kml_dbscan<-NA
    for (i in (1:nrow(plot_sampleedf)))
    {
      plot_sampleedf$kml_dbscan[i]<-kmldf$traj[which(paste0(round(kmldf$times,2),".00",sep="")==(paste0(round(sample_gt_elev$dist[i],0),".00",sep="")))]
      
    }
  }else{ plot_sampleedf$kml_dbscan<-kmldf$traj[which(kmldf$times %in%  round(sample_gt_elev$dist))]}
  
colnames(plot_sampleedf)<-c("dist","GT","highRes","ShighRes","Fmean","Fmedian","FSCpam","FSCmean","FSCmedian","FSCkml")
 plot_sampleedf<-plot_sampleedf[,c(1,2,3,4,5,6,8,9,7,10)]
#step 14 rwgps profile
plot_sampleedf$RWGPS<-strtrwgps$altitude
#step 15 smoothing prototypes
  plot_sampleedf[,c(12:17)]<-NA
  colnames(plot_sampleedf)[12:17]<-as.vector(paste("S",colnames(plot_sampleedf)[5:10],sep=""))
  # Specify the columns to be smoothed
  columns_to_smooth <- 5:10
  
  # Apply the Savitzky-Golay filter to the selected columns
  plot_sampleedf[, (columns_to_smooth + 7)] <- lapply(plot_sampleedf[, columns_to_smooth], sgolayfilt, n = sg)
write.csv(plot_sampleedf, file = paste("Output//method1",nme,"elev.csv",sep=""))
#step 16 evaluation measures
# Columns to consider
columns_to_compare <- c(3:4, 12:17, 11)

  # RMSE
  rmsegrandvlbrdg <- sapply(columns_to_compare, function(i) rmse(plot_sampleedf$GT, plot_sampleedf[, i]))
  # Hausdorff distance
  hgrandvlbrdg <- sapply(columns_to_compare, function(i) {
    hausdorff_dist(
      matrix(c(plot_sampleedf$dist, plot_sampleedf$GT), ncol = 2),
      matrix(c(plot_sampleedf$dist, plot_sampleedf[, i]), ncol = 2)
    )
  })
  # EMD
  emdgrandvlbrdg <- sapply(columns_to_compare, function(i) {
    emd(
      SpatialPoints(matrix(c(plot_sampleedf$dist, min_max_norm(plot_sampleedf$GT)), ncol = 2)),
      SpatialPoints(matrix(c(plot_sampleedf$dist, min_max_norm(plot_sampleedf[, i])), ncol = 2))
    )
  })
  # Frechet Distance
  fregrandvlbrdg <- sapply(columns_to_compare, function(i) {
    FrechetDistance(plot_sampleedf[, i], plot_sampleedf$GT, plot_sampleedf$dist, plot_sampleedf$dist)
  })
  
  #build the evaluation table
  tbl<-data.frame(matrix(NA,nrow = 9,ncol = 4))
  rownames(tbl)<-colnames(plot_sampleedf)[c(3:4,12:17,11)]
  colnames(tbl)<-c("RMSE","H","EM","Frechet")
  tbl[,1]<-round(rmsegrandvlbrdg,2)
  tbl[,2]<-round(hgrandvlbrdg,2)
  tbl[,3]<-round(emdgrandvlbrdg,2)
  tbl[,4]<-round(fregrandvlbrdg,2)

  write.csv(tbl, file = paste("Output//EvaluationTables_method1",nme,"elev.csv",sep=""))
  #step 17 extend to grade
  # Columns to compute gradients for
  columns_to_compute_gradients <- c(
    "GT", "highRes", "ShighRes", "SFmean", "SFmedian",
    "SFSCmean", "SFSCmedian", "SFSCpam", "SFSCkml", "RWGPS"
  )
  for (col_name in columns_to_compute_gradients) {
    gradient_col_name <- paste0(col_name, "_g")
    plot_sampleedf[, gradient_col_name] <- c(NA,diff(plot_sampleedf[, col_name]) / diff(plot_sampleedf$dist))
  }

  
  plot_sampleegf<-plot_sampleedf[,c(1,18:27)]
  colnames(plot_sampleegf)<-c("dist","GT",colnames(plot_sampleedf)[c(3:4,12:17,11)])
 
  
  plot_sampleegf<-plot_sampleegf[-1,]
  write.csv(plot_sampleegf, file = paste("Output//method1",nme,"grade.csv",sep=""))
  #step 18 evaluation measures 
  # Columns to consider
  columns_to_compare <- c(3:11)
  
  # RMSE
  grmsegrandvlbrdg <- sapply(columns_to_compare, function(i) rmse(plot_sampleegf$GT, plot_sampleegf[, i]))
  # Hausdorff distance
  ghgrandvlbrdg <- sapply(columns_to_compare, function(i) {
    hausdorff_dist(
      matrix(c(plot_sampleegf$dist, plot_sampleegf$GT), ncol = 2),
      matrix(c(plot_sampleegf$dist, plot_sampleegf[, i]), ncol = 2)
    )
  })
  # EMD
  gemdgrandvlbrdg <- sapply(columns_to_compare, function(i) {
    emd(
      SpatialPoints(matrix(c(plot_sampleegf$dist, min_max_norm(plot_sampleegf$GT)), ncol = 2)),
      SpatialPoints(matrix(c(plot_sampleegf$dist, min_max_norm(plot_sampleegf[, i])), ncol = 2))
    )
  })
  # Frechet Distance
  gfregrandvlbrdg <- sapply(columns_to_compare, function(i) {
    FrechetDistance(plot_sampleegf[, i], plot_sampleegf$GT, plot_sampleegf$dist, plot_sampleegf$dist)
  })
  #build the table
  tbl<-data.frame(matrix(NA,nrow = 9,ncol = 4))
  rownames(tbl)<-colnames(plot_sampleedf)[3:11]
  colnames(tbl)<-c("RMSE","H","EM","Frechet")
  tbl[,1]<-round(grmsegrandvlbrdg,2)
  tbl[,2]<-round(ghgrandvlbrdg,2)
  tbl[,3]<-round(gemdgrandvlbrdg,2)
  tbl[,4]<-round(gfregrandvlbrdg,2)

  write.csv(tbl, file = paste("Output//EvaluationTables_method1",nme,"grade.csv",sep=""))


#########################################METHOD 2:Efforts>Filter-Smooth>GradeEfforts>Smooth-Clean>Gradeprototype#########################################
#step1 create grade efforts
  for (m in (1:length(sample_effort_single)))
  {

    q<-which(colnames(sample_gt_elev)==paste("z_eff",m,sep=""))
    if(length(q)!=0){
      
      stmp<-sgolayfilt(sample_gt_elev[,q], n= sg)
      gtmp<-c(NA,diff(stmp) / diff(sample_gt_elev$dist))
      sample_gt_elev <- cbind(sample_gt_elev, gtmp)
      colnames(sample_gt_elev)[length(colnames(sample_gt_elev))]<-paste("g_eff",m,sep="")
    }   
    
  }
  
  sample_gt_grade<-sample_gt_elev[,c(1:6,which(colnames(sample_gt_elev)=="g_eff1"):which(colnames(sample_gt_elev)==paste("g_eff",length(sample_effort_single),sep = "")))]
#step 2 create Ground truth grade effort

  sample_gt_grade$GGT <- c(NA, diff(sample_gt_elev$zGT) / diff(sample_gt_elev$dist))

#step 3 add our estimate from high res data source or any source you want to compare to
  sample_gt_grade$GhighRes<-c(NA, diff(sample_gt_elev$ShighRes) / diff(sample_gt_elev$dist))

#step 4 adding what MJ paper did with respect to estimate mean and median from elevation obtained above
  sample_gt_grade$meanG<-rowMeans(sample_gt_grade[,c(which(colnames(sample_gt_grade)=="g_eff1"):which(colnames(sample_gt_grade)==paste("g_eff",length(sample_effort_single),sep = "")))])
  sample_gt_grade$meadianG<-rowMedians(as.matrix(sample_gt_grade[,c(which(colnames(sample_gt_grade)=="g_eff1"):which(colnames(sample_gt_grade)==paste("g_eff",length(sample_effort_single),sep = "")))]))
#step 5 add RWGPS
  sample_gt_grade$GRWGPS<-c(NA, diff(strtrwgps$altitude) / diff(sample_gt_elev$dist))
#step 6 where our contribution starts first we detect outliers using DBSCAN

  sample_gt_grade<-sample_gt_grade[-1,]
  sample_ls <- lapply(sample_gt_grade[, c(which(colnames(sample_gt_grade)=="g_eff1"):which(colnames(sample_gt_grade)==paste("g_eff",length(sample_effort_single),sep = "")))], sgolayfilt, n = sg)
  
  #now we identify and remove outliers (i.e, cleaning step)
  sm<-do.call(cbind,sample_ls)
  # Calculate the Euclidean distance matrix using the dist function
  distm <- dist(t(sm),method = "euclidean")
  
  db<-dbscan::dbscan(distm, eps=0.35, minPts = 2)
  if(length(which(db$cluster==0))!=0){sm<-sm[,-which(db$cluster==0)]
  sample_ls<-sample_ls[-which(db$cluster==0)]}
  length(sample_ls)
  #step 7 we add Clustering methods: pam and kml respectively length(sample_ls) #this shows number of filtered and cleaned efforts
  smw<-t(sm)
  traj_long<-cldsWide(smw,times=seq(1:dim(smw)[2]),id=as.factor(as.character(seq(1:nrow(smw)))))
  smitcls2<-kmlShape(traj_long, nbClusters = 1, timeScale =0.1, FrechetSumOrMax ="max", toPlot="both", parAlgo=parKmlShape())
  distances_not_in_smitcls2 <- sample_gt_grade$dist[!(round(sample_gt_grade$dist) %in% smitcls2@trajMeans[, 2])]
  
  # Create a data frame with times and traj columns
  kmldf <- data.frame(times = smitcls2@trajMeans[, 2], traj = smitcls2@trajMeans[, 3])
  
  # Add rows for distances not in smitcls2@trajMeans[, 2] with NA traj values
  new_rows <- data.frame(times = distances_not_in_smitcls2, traj = NA)
  kmldf <- rbind(kmldf, new_rows)
  
  # Sort the data frame by times
  kmldf <- kmldf[order(kmldf$times),]
  
  # Interpolate NA values in the traj column using spline interpolation
  kmldf$traj <- na.interpolation(kmldf$traj, option = "spline")
  
  sample_gt_grade$Gtmp<-sgolayfilt( sample_gt_grade$GhighRes, n= sg)
  
  sample_gt_grade$GhighRes<-sample_gt_grade$Gtmp
  #step 8 now we create a data frame of all of our estimates for each point along sample loaction
  plot_sampleedf<-sample_gt_grade[,c(which(colnames(sample_gt_grade)=="dist"),
                                     which(colnames(sample_gt_grade)=="GGT"),
                                     which(colnames(sample_gt_grade)=="GhighRes"),
                                     which(colnames(sample_gt_grade)=="meanG"),
                                     which(colnames(sample_gt_grade)=="meadianG"))
  ] 
 
  colnames(plot_sampleedf)<-c("dist","GT","ShighRes","Fmean","Fmedian")
  #step 9 add our estimates from Filtered and smoothed efforts
  plot_sampleedf$pam_dbscan<-as.numeric(pam(t(sm),1)$medoids)
  plot_sampleedf$mean_dbscan<-rowMeans(sm)
  plot_sampleedf$median_dbscan<-rowMedians(sm)  

  if(length(kmldf$traj[which(kmldf$times %in%  round(sample_gt_grade$dist))])!=nrow(plot_sampleedf)){
    plot_sampleedf$kml_dbscan<-NA
    for (i in (1:nrow(plot_sampleedf)))
    {
      plot_sampleedf$kml_dbscan[i]<-kmldf$traj[which(paste0(round(kmldf$times,2),".00",sep="")==(paste0(round(sample_gt_grade$dist[i],0),".00",sep="")))]
      
    }
  }else{ plot_sampleedf$kml_dbscan<-kmldf$traj[which(kmldf$times %in%  round(sample_gt_grade$dist))]}
  
  
  
  colnames(plot_sampleedf)<-c("dist","GT","ShighRes","Fmean","Fmedian","FSCpam","FSCmean","FSCmedian","FSCkml")
  plot_sampleedf<-plot_sampleedf[,c(1,2,3,4,5,7,8,6,9)]
  plot_sampleedf$RWGPS<-sample_gt_grade$GRWGPS
  write.csv(plot_sampleedf, file = paste("Output//method2",nme,"grade.csv",sep=""))
  #step 16 evaluation measures
  # Columns to consider
  columns_to_compare <- 3:10
  
  # RMSE
  rmsegrandvlbrdg <- sapply(columns_to_compare, function(i) rmse(plot_sampleedf$GT, plot_sampleedf[, i]))
  
  # Hausdorff distance
  hgrandvlbrdg <- sapply(columns_to_compare, function(i) {
    hausdorff_dist(
      matrix(c(plot_sampleedf$dist, plot_sampleedf$GT), ncol = 2),
      matrix(c(plot_sampleedf$dist, plot_sampleedf[, i]), ncol = 2)
    )
  })
  # EMD
  library(move)
  emdgrandvlbrdg<-vector()
  emdgrandvlbrdg <- sapply(columns_to_compare, function(i) {
    if (any(is.nan(min_max_norm(plot_sampleedf[, i])))) {
      emd(SpatialPoints(matrix(c(plot_sampleedf$dist, min_max_norm(plot_sampleedf$GT)), ncol = 2)),
          SpatialPoints(matrix(c(plot_sampleedf$dist, rep(0, nrow(plot_sampleedf))), ncol = 2))
      )
    } else {
      emd(SpatialPoints(matrix(c(plot_sampleedf$dist, min_max_norm(plot_sampleedf$GT)), ncol = 2)),
          SpatialPoints(matrix(c(plot_sampleedf$dist, min_max_norm(plot_sampleedf[, i])), ncol = 2))
      )
    }
  })
  
  # Frechet Distance
  fregrandvlbrdg <- sapply(columns_to_compare, function(i) {
    FrechetDistance(plot_sampleedf[, i], plot_sampleedf$GT, plot_sampleedf$dist, plot_sampleedf$dist)
  })
  #build the table
  tbl<-data.frame(matrix(NA,nrow = 8,ncol = 4))
  rownames(tbl)<-colnames(plot_sampleedf)[3:10]
  colnames(tbl)<-c("RMSE","H","EM","Frechet")
  tbl[,1]<-round(rmsegrandvlbrdg,2)
  tbl[,2]<-round(hgrandvlbrdg,2)
  tbl[,3]<-round(emdgrandvlbrdg,2)
  tbl[,4]<-round(fregrandvlbrdg,2)
  
  write.csv(tbl, file = paste("Output//EvaluationTables_method2",nme,"grade.csv",sep=""))
  
 


