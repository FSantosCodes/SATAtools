#' @title Detect breakpoints from a GEE ts
#' @description This function detects breakpoints from GEE ts, detecting them and giving as output a shapefile with the correspondent dates of change. To calibrate the algorithm, it is possible to create a series of plots which aims to help to define its parameters (ie. change_threshold and segment_size). OSGeo4W is required to be installed and 'OSGeo4W.bat' declared as enviroment variable for its proper working.
#' @param GEE_ts GEE time series filename. The GEE time series to analyze and detect breakpoints. It should be projected in UTM WGS84 coordinates system.
#' @param change_threshold Numerical. A breakpoint detection change threshold based in the standard deviation from mean of GEE ts. Values can vary from 0.001 (more sensibility to changes) to > 1 (less sensibility)
#' @param segment_size Numerical. The minimun segment size to consider in brakepoint detection. Depends of the number of observations in the GEE ts.
#' @param min_area Numerical. The minimun area of change areas detected by the algorithm
#' @param calibration_mode Logical. Activates the calibration mode of the algorithm. This generates a series of plots of breakpoints based in samples points (requires a points based shapefile of samples).
#' @param samples_shp Shapefile filename. A set of samples for evaluate calibration settings of the algorithm.
#' @param column_id_name String. The name of the column of the unique identifier of each sample.
#' @param size_jpg Numeric. A value which defines the size of plots (in pixels). Plots are based in an square, so one dimension defines its size.
#' @param ncores Numeric. A number defining how many cores will be used in processing.
#' @param output_folder Foldername. A folder to use for store plots.
#' @return Plots outputs are stored in: *output_folder*. If calibration_mode is TRUE, a series of plots are create, otherwise a shapefile of detected breakpoints is generated.
#' @examples
#' #Defining variables for breakpoint detection
#' GEE_ts <- "C:/DATA/GAF/demo/sentinel1/gee_outputs/coca_S1_ts_VH_P50_I14.tif"
#' change_threshold <- 10
#' segment_size <- 3
#' min_area <- 4
#'
#' #For calibration mode
#' samples_shp <- "C:/DATA/GAF/demo/sentinel1/samples/muestras.shp"
#' column_id_name <- "clase"
#' size_jpg <- 3000
#'
#' #Variables for processing and outputs
#' ncores <- 3
#' output_folder <- "C:/DATA/GAF/demo/sentinel1/breakpoints"
#'
#' #running the function in CALIBRATION mode
#' detect_breakpoints(GEE_ts = GEE_ts,
#'                    change_threshold = change_threshold,
#'                    segment_size = segment_size,
#'                    min_area = min_area,
#'                    calibration_mode = T,
#'                    samples_shp = samples_shp,
#'                    column_id_name = column_id_name,
#'                    size_jpg = size_jpg,
#'                    ncores = ncores,
#'                    output_folder = output_folder)
#'
#' #running the function in PROCESSING mode
#' detect_breakpoints(GEE_ts = GEE_ts,
#'                    change_threshold = change_threshold,
#'                    segment_size = segment_size,
#'                    min_area = min_area,
#'                    calibration_mode = F,
#'                    column_id_name = column_id_name,
#'                    ncores = ncores,
#'                    output_folder = output_folder)
#' @export

##### MAIN FUNCTION #####

detect_breakpoints <- function(
  GEE_ts,
  change_threshold,
  segment_size,
  min_area = 4,
  calibration_mode = F,
  samples_shp,
  column_id_name,
  size_jpg = 1000,
  ncores = 1,
  output_folder
){

  #(0)#### TEST #####

  test <- F
  if(test){
    GEE_ts <- "C:/DATA/GAF/demo/sentinel1/gee_outputs/coca_S1_ts_VH_P50_I14.tif"
    change_threshold <- 10
    segment_size <- 3
    min_area <- 4
    calibration_mode <- T
    samples_shp <- "C:/DATA/GAF/demo/sentinel1/samples/muestras.shp"
    column_id_name <- "clase"
    size_jpg <- 3000
    ncores <- 3
    output_folder <- "C:/DATA/GAF/demo/sentinel1/breakpoints"
  }

  #(1)#### CONFIG #####

  check.libs <- c("rgdal","raster","foreach","doParallel",
                  "ggplot2","changepoint","gtools","plyr",
                  "lubridate","rgeos")
  get.libraries <- function(libraries){
    for (i in libraries){
      if(!require(i,character.only=T)){
        install.packages(i)
        library(i,character.only=T)}
    }
  }
  get.libraries(check.libs)

  #(2)#### READ DATA #####

  #read GEE ts
  if(!file.exists(GEE_ts)){
    stop("GEE_ts do not exists")
  }
  ts.band <- stack(GEE_ts)
  #read dates
  ts.dates <- gsub(".tif",".csv",basename(GEE_ts))
  ts.dates <- paste0(dirname(GEE_ts),"/",gsub("ts","dates",ts.dates))
  if(file.exists(ts.dates)){
    ts.dates <- read.csv(ts.dates,header=T,stringsAsFactors=F)$dates
    ts.dates <- unlist(strsplit(gsub("\\[|\\]| ","",ts.dates),"[,]"))
  }else{
    stop("GEE_ts .csv file is missing, can\'t get time series dates")
  }
  if(length(ts.dates)!=nlayers(ts.band)){
    stop("the number of observation in GEE_ts differ from the number of dates in .csv file")
  }
  #create ouput folder
  dir.create(output_folder,showWarnings = F, recursive = T)

  #(3)#### CALIBRATION MODE #####

  if(calibration_mode){
    #get & prepare shapefile
    samples.shp <- shapefile(samples_shp)
    samples.nam <- grep(column_id_name,names(samples.shp))
    if(length(samples.nam)==0){
      stop("column_id_name not found in the shapefile")
    }
    samples.nam <- samples.shp@data[,samples.nam]
    #project
    if(samples.shp@proj4string@projargs != ts.band@crs@projargs){
      samples.shp <- spTransform(samples.shp, ts.band@crs)
    }
    #extract data + standarize
    #samples.mat <- t(as.matrix(extract(ts.band,samples.shp,df=T)))
    #samples.mat <- t(scale(samples.mat[2:nrow(samples.mat),]))
    xy.coords <- coordinates(samples.shp)
    samples.mat <- apply(xy.coords,1,function(x){
      x <- system(paste0("osgeo4w && gdallocationinfo -valonly -geoloc ",GEE_ts," ",paste(x,collapse=" ")),intern=T)
      x <- x[grep("OSGEO",x):length(x)]
      x[1] <- unlist(strsplit(x[1],">"))[2]
      return(as.numeric(x))
    })
    samples.mat <- t(scale(samples.mat))
    samples.mat[is.na(samples.mat)] <- 0
    #run changepoint
    samples.chg <- cpt.mean(samples.mat,
                            penalty="Manual",
                            pen.value=change_threshold,
                            minseglen=segment_size,
                            class=F)[,1]
    dummy.val <- samples.chg == ncol(samples.mat)
    samples.chg[dummy.val] <- NA; rm(dummy.val)
    #save plots
    samples.num <- 1:nrow(samples.mat)
    samples.num <- split(samples.num, ceiling(seq_along(samples.num)/12))
    for(i in 1:length(samples.num)){
      #plot data
      plot.x <- samples.mat[samples.num[[i]],]
      plot.chg <- samples.chg[samples.num[[i]]]
      plot.names <- samples.nam[samples.num[[i]]]
      if(length(plot.names)!=12){
        plot.names <- c(plot.names,rep(NA,12-nrow(plot.x)))
      }
      #plot config
      out.name <- paste0("PLOT_",i,"_",unlist(strsplit(basename(GEE_ts),"[.]"))[1],".jpg")
      out.name <- paste0(output_folder,"/",out.name)
      jpeg(out.name,width=size_jpg,height=size_jpg,res=300)
      par(mfrow=c(4,3))
      for(j in 1:12){
        if(!is.na(plot.names[j])){
          plot(1:ncol(plot.x),
               plot.x[j,],
               main=plot.names[j],
               xlab="observations",
               ylab="z-scores")
          if(!is.na(plot.chg[j])){
            abline(v=plot.chg[j],col="red",lty=2)
            mtext(ts.dates[plot.chg[j]],3,at=plot.chg[j],col="red",cex=0.75)
          }
        }
      }
      p.text <- paste0("change_threshold = ",change_threshold," segment_size = ",segment_size)
      mtext(p.text, outer = T,cex=1, line=-1.5)
      dev.off()
    }
  }

  #(4)#### PROCESSING MODE #####

  if(!calibration_mode){
    #start time
    ptm <- proc.time()
    #change function
    chg.function <- function(data.input){
      data.input <- t(scale(data.input))
      data.input[is.na(data.input)] <- 0
      data.input <- cpt.mean(data.input,
                             penalty="Manual",
                             pen.value=change_threshold,
                             minseglen=segment_size,
                             class=F)[,1]
      return(data.input)
    }
    #prepare
    csv.folder <- paste0(output_folder,"/csv");dir.create(csv.folder,showWarnings = F, recursive = T)
    unlink(list.files(csv.folder,full.names=T),recursive = T, force = T)
    img.blocks <- raster::blockSize(ts.band)
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    #start parallel
    foreach(i=1:img.blocks$n,.packages=c("raster","changepoint")) %dopar% {
      #get data
      data.input <- t(raster::as.matrix(getValues(ts.band,row=img.blocks$row[i],nrows=img.blocks$nrows[i])))
      #function
      data.input <- chg.function(data.input)
      #save
      out.name <- paste0(csv.folder,'/breaks_',i,'.csv')
      write.csv(data.input,out.name,row.names = F)
    }
    stopCluster(cl)
    #get csv
    csv.files <- mixedsort(list.files(csv.folder,full.names=T,pattern=".csv"))
    data.input <- rbind.fill(lapply(csv.files, read.csv, header=TRUE))
    #remove dummy values
    data.input <- data.input[,1]
    dummy.val <- data.input == nlayers(ts.band)
    data.input[dummy.val] <- 0; rm(dummy.val)
    data.input[is.na(data.input)] <- 0
    #put in raster
    template <- raster(nrows=ts.band[[1]]@nrows,ncols=ts.band[[1]]@ncols,ext=extent(ts.band[[1]]),crs=ts.band[[1]]@crs)
    values(template) <- data.input
    ras.name <- unlist(strsplit(basename(GEE_ts),"[.]"))[1]
    ras.name <- paste0(output_folder,"/BREAKS_",ras.name,".tif")
    writeRaster(template,ras.name,overwrite=T,datatype="INT2U")
    #filter
    system(paste0("osgeo4w && gdal_sieve -st ",min_area," ",ras.name))
    template <- raster(ras.name)
    values(template)[values(template)==0] <- NA
    writeRaster(template,ras.name,overwrite=T,datatype="INT2U")
    #put in shapefile
    shp.name <- unlist(strsplit(basename(GEE_ts),"[.]"))[1]
    shp.name <- paste0(output_folder,"/BREAKS_",shp.name,".shp")
    if(file.exists(shp.name)){
      unlink(shp.name,recursive = T, force = T)
    }
    system(paste0("osgeo4w && gdal_polygonize ",ras.name," -f \"ESRI Shapefile\" ",shp.name))
    #prepare
    ts.start <- sapply(strsplit(ts.dates,"_"),head,1)
    ts.end <- sapply(strsplit(ts.dates,"_"),tail,1)
    shp.out <- shapefile(shp.name)
    shp.out@data$start <- ts.start[shp.out@data$DN]
    shp.out@data$end <- ts.end[shp.out@data$DN]
    shp.out@data$DN <- NULL
    shp.out@data$area <- gArea(shp.out,byid=T)
    shapefile(shp.out,shp.name,overwrite=T)
    unlink(ras.name,recursive = T, force = T)
    #end timer
    print(paste0(as.character(round((proc.time() - ptm)[3],2))," seconds"))
  }
  print("*** script finished succesfully ***")
}
