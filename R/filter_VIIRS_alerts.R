#' @title Filter VIIRS fire alerts
#' @description This function filters VIIRS fire alerts. It uses a shapefile of spatial units to generate plots and visualize some statistics of alerts by spatial units.
#' @param spatial_units_shp Shapefile filename. This is the shapefile of spatial units to analyze with SATA alerts. It should be projected in UTM WGS84 coordinates system.
#' @param column_units_name String. The name of the column indicating the names of the spatial units to analize. Case-sensitive applies.
#' @param VIIRS_shp Shapefile filename. The shapefile of VIIRS alerts from SATA (as points).
#' @param column_date_name String. The name of the column indicating the date of alerts. Case-sensitive applies.
#' @param column_frp_name String. The name of the column indicating the radiative power of alerts (normally as frp). Case-sensitive applies.
#' @param column_confidence_name String. The name of the column indicating the confidence of alerts (normally as confidence). Case-sensitive applies.
#' @param min_frp Numeric. A value defining the minimun radiative power in filtering.
#' @param min_confidence String. A category defining the minimun confidence in filtering (normally three: "high", "nominal" and "low").
#' @param size_jpg Numeric. Two values which define the width and height of boxplot (in cms).
#' @param output_folder Foldername. A folder to use for store alerts filtered.
#' @return Filtered alerts and plots are stored in: *output_folder*.
#' @examples
#' #Defining variables from spatial units shapefile
#' spatial_units_shp <- "C:/DATA/GAF/demo/shp/unidades_espaciales/provincias.shp"
#' column_units_name <- "dpa_despro" #case sensitive
#'
#' #Defining variables from VIIRS shapefile
#' VIIRS_shp <- "C:/DATA/GAF/demo/shp/alertas/pts_calor_VIIRS_30d.shp"
#' column_date_name <- "acq_date" #case sensitive
#' column_frp_name <- "frp" #case sensitive
#' column_confidence_name <- "confidence" #case sensitive
#' min_frp <- 10 #same units as area column
#' min_confidence <- "nominal" #all categories below this one are filtered
#'
#' #Defining variables of ouputs
#' size_jpg <- c(40,20) #first value defines width and second height
#' output_folder <- "C:/DATA/GAF/demo/estadisticas"
#'
#' #running the function
#' filter_VIIRS_alerts(spatial_units_shp = spatial_units_shp,
#'                             column_units_name = column_units_name,
#'                             VIIRS_shp = VIIRS_shp,
#'                             column_date_name = column_date_name,
#'                             column_frp_name = column_frp_name,
#'                             column_confidence_name = column_confidence_name,
#'                             min_frp = min_frp,
#'                             min_confidence = min_confidence,
#'                             size_jpg = size_jpg,
#'                             output_folder = output_folder)
#' @export

##### MAIN FUNCTION #####

filter_VIIRS_alerts <- function(
  spatial_units_shp,
  column_units_name,
  VIIRS_shp,
  column_date_name,
  column_frp_name,
  column_confidence_name,
  min_frp,
  min_confidence,
  size_jpg = c(10,10),
  output_folder
){

  #(0)#### TEST #####

  test <- F
  if(test){
    spatial_units_shp <- "C:/DATA/GAF/demo/shp/unidades_espaciales/provincias.shp"
    column_units_name <- "dpa_despro"
    VIIRS_shp <- "C:/DATA/GAF/demo/shp/alertas/pts_calor_VIIRS_30d.shp"
    column_date_name <- "acq_date"
    column_frp_name <- "frp"
    column_confidence_name <- "confidence"
    min_frp <- 10
    min_confidence <- "nominal"
    size_jpg <-  c(10,10)
    output_folder <- "C:/DATA/GAF/demo/estadisticas"
  }

  #(1)#### CONFIG #####

  check.libs <- c("rgdal","raster","sp","ggplot2",
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

  #check spatial units
  if(!file.exists(spatial_units_shp)){
    stop("spatial_units_shp do not exists")
  }
  spa.shp <- shapefile(spatial_units_shp)
  spa.prj <- spa.shp@proj4string@projargs
  spa.nam <- grep(column_units_name,names(spa.shp))
  if(length(spa.nam)==0){
    stop("column_units_name do not found in spatial_units_shp")
  }
  #create ouput folder
  dir.create(output_folder,showWarnings = F, recursive = T)

  #(2)#### GET VIIRS #####

  #check VIIRS fires alerts
  if(!file.exists(VIIRS_shp)){
    stop("VIIRS_shp not found")
  }
  #read shapefile
  viirs.shp <- shapefile(VIIRS_shp)
  viirs.shp <- viirs.shp[!duplicated(viirs.shp@data),]
  #date name
  viirs.date <- grep(column_date_name,names(viirs.shp))
  if(length(viirs.date)==0){
    stop(paste0("do not exists column_date_name: ", column_date_name))
  }
  names(viirs.shp)[viirs.date] <- "acq_date"
  #frp name
  viirs.frp <- grep(column_frp_name,names(viirs.shp))
  if(length(viirs.frp)==0){
    stop(paste0("do not exists column_frp_name: ", column_frp_name))
  }
  names(viirs.shp)[viirs.frp] <- "frp"
  #confidence name
  viirs.conf <- grep(column_confidence_name,names(viirs.shp))
  if(length(viirs.conf)==0){
    stop(paste0("do not exists column_confidence_name: ", column_confidence_name))
  }
  names(viirs.shp)[viirs.conf] <- "confidence"

  #(3)#### FUNCTIONS #####

  summary.stats <- function(z=x.int[[1]]){
    z.count <- nrow(z)
    z.start <- min(z$acq_date)
    z.end <- max(z$acq_date)
    z.stats <- data.frame(alerts=z.count,
                          start=z.start,
                          end=z.end,
                          frp_avg=round(mean(z$frp),2),
                          conf_mode=names(which.max(table(z$confidence))),
                          stringsAsFactors=F)
    return(z.stats)
  }

  get.viirs <- function(x=spa.shp,y=viirs.shp){
    #filter
    y <- y[y@data$frp >= min_frp,]
    if(min_confidence ==  "nominal"){
      y <- y[y@data$confidence != "low",]
    }else if(min_confidence ==  "high"){
      y <- y[y@data$confidence == "high",]
    }
    #project
    if(y@proj4string@projargs != spa.prj){
      y <- spTransform(y, CRS(spa.prj))
    }
    #get coordinates
    y@data$x_coord <- coordinates(y)[,1]
    y@data$y_coord <- coordinates(y)[,2]
    #intersect & process
    x.int <- sp::over(x,y,returnList=T)
    if(length(unlist(x.int))==0){
      x.ret <- "No VIIRS alerts found for these spatial units"
    }else{
      #remove null
      names(x.int) <- x@data[,spa.nam]
      x.int <- x.int[sapply(x.int,nrow)!=0]
      x.nam <- names(x.int)
      #SUMMARY ALERTS
      x.summary <- lapply(x.int,summary.stats)
      x.summary <- do.call("rbind",x.summary)
      x.summary <- cbind(data.frame(names=rownames(x.summary)),x.summary)
      rownames(x.summary) <- NULL
      x.summary <- x.summary[order(x.summary$alerts,decreasing=T),]
      x.summary$start <- gsub("Z","",x.summary$start)
      x.summary$end <- gsub("Z","",x.summary$end)
      #ALL ALERTS
      x.all <- do.call("rbind",x.int)
      rownames(x.all) <- NULL
      x.all <- SpatialPointsDataFrame(data.frame(x.all$x_coord,x.all$y_coord),
                                      proj4string=CRS(spa.prj),
                                      data.frame(names=rep(names(x.int),sapply(x.int,nrow)),
                                                 date=gsub("Z","",x.all$acq_date),
                                                 frp=x.all$frp,
                                                 conf=x.all$confidence,
                                                 stringsAsFactors=F))
      #return
      x.ret <- list(summary=x.summary,all=x.all)
    }
    return(x.ret)
  }

  #(4)#### RUN #####

  viirs.shp <- get.viirs(spa.shp,viirs.shp)

  #(5)#### SAVE #####

  #save
  if(length(viirs.shp)==2){
    #outputs
    out.name <- paste0(output_folder,"/VIIRS_stats.csv")
    write.csv(viirs.shp$summary,out.name,row.names=F)
    out.name <- paste0(output_folder,"/VIIRS_alerts.shp")
    shapefile(viirs.shp$all,out.name,overwrite=T)
    #plots1
    out.name <- paste0(output_folder,"/VIIRS_plot.jpg")
    area.p <- ggplot(viirs.shp$all@data,aes(x=names,y=frp)) +
      geom_boxplot() +
      coord_flip() +
      xlab("Spatial units") +
      ylab("Radiative power [megawatts]") +
      labs(title="VIIRS: radiative power by spatial unit")
    ggsave(out.name,area.p,width = size_jpg[1], height = size_jpg[2], units="cm")
  }else{
    warning(viirs.shp)
  }

  ###### close ######
  print("****script ends*****")

}
