#' @title Filter deforestation alerts
#' @description This function filters SATA deforestation alerts. It uses a shapefile of spatial units to generate plots and visualize some statistics of alerts by spatial units.
#' @param spatial_units_shp Shapefile filename. This is the shapefile of spatial units to analyze with SATA alerts. It should be projected in UTM WGS84 coordinates system.
#' @param column_units_name String. The name of the column indicating the names of the spatial units to analize. Case-sensitive applies.
#' @param deforestation_shp Shapefile filename. The shapefile of deforestation alerts from SATA (as polygons).
#' @param column_date_name String. The name of the column indicating the date of alerts. Case-sensitive applies.
#' @param column_area_name String. The name of the column indicating the area of alerts. Case-sensitive applies.
#' @param min_area Numeric. A value defining the minimun area in filtering. NA ignore filtering.
#' @param size_jpg Numeric. Two values which define the width and height of boxplot (in cms).
#' @param output_folder Foldername. A folder to use for store alerts filtered.
#' @return Filtered alerts and plots are stored in: *output_folder*.
#' @examples
#' #Defining variables from spatial units shapefile
#' spatial_units_shp <- "C:/DATA/GAF/demo/shp/unidades_espaciales/provincias.shp"
#' column_units_name <- "dpa_despro" #case sensitive
#'
#' #Defining variables from deforestation shapefile
#' deforestation_shp <- "C:/DATA/GAF/demo/shp/alertas/deforestation_2016_2017.shp"
#' column_date_name <- "inicio" #case sensitive
#' column_area_name <- "área" #case sensitive
#' min_area <- 1 #with same units as area column but here ignored
#'
#' #Defining variables of ouputs
#' size_jpg <- c(40,20) #first value defines width and second height
#' output_folder <- "C:/DATA/GAF/demo/estadisticas"
#'
#' #running the function
#' filter_deforestation_alerts(spatial_units_shp = spatial_units_shp,
#'                             column_units_name = column_units_name,
#'                             deforestation_shp = deforestation_shp,
#'                             column_date_name = column_date_name,
#'                             column_area_name = column_area_name,
#'                             min_area = min_area,
#'                             size_jpg = size_jpg,
#'                             output_folder = output_folder)
#' @export

##### MAIN FUNCTION #####

filter_deforestation_alerts <- function(
  spatial_units_shp,
  column_units_name,
  deforestation_shp,
  column_date_name,
  column_area_name,
  min_area = NA,
  size_jpg = c(10,10),
  output_folder
){

  #(0)#### TEST #####

  test <- F
  if(test){
    spatial_units_shp <- "C:/DATA/GAF/demo/shp/unidades_espaciales/provincias.shp"
    column_units_name <- "dpa_despro"
    deforestation_shp <- "C:/DATA/GAF/demo/shp/alertas/deforestation_2016_2017.shp"
    column_date_name <- "inicio"
    column_area_name <- "área"
    min_area <- 1
    size_jpg <- c(10,10)
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

  #(3)#### GET DEFORESTATION #####

  #check deforestation alerts
  if(!file.exists(deforestation_shp)){
    stop("deforestation_shp not found")
  }
  #read shape
  def.shp <- shapefile(deforestation_shp)
  def.shp <- def.shp[!duplicated(def.shp@data),]
  #date name
  def.date <- grep(column_date_name,names(def.shp))
  if(length(def.date)==0){
    stop(paste0("do not exists column_date_name: ", column_start_name))
  }
  names(def.shp)[def.date] <- "acq_date"
  #area name
  def.area <- grep(column_area_name,names(def.shp))
  if(length(def.area)==0){
    stop(paste0("do not exists column_area_name: ", column_area_name))
  }
  names(def.shp)[def.area] <- "area"

  #(4)#### FUNCTIONS #####

  summary.stats <- function(z=x.int[[1]]){
    z.count <- nrow(z)
    z.start <- min(z$acq_date)
    z.end <- max(z$acq_date)
    z.stats <- data.frame(alerts=z.count,
                          start=z.start,
                          end=z.end,
                          area_avg=round(mean(z$area),2),
                          stringsAsFactors=F)
    return(z.stats)
  }

  get.deforestation  <- function(x=spa.shp,y=def.shp){
    #filter + put as points
    if(!is.na(min_area)){
      y <- y[y@data$area >= min_area,]
    }
    y <- SpatialPointsDataFrame(coordinates(gCentroid(y,byid=T)),
                                      y@data,proj4string=y@proj4string)
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
      x.ret <- "No DEFORESTATION alerts found for these spatial units"
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
                                                 area=x.all$area,
                                                 stringsAsFactors=F))
      #return
      x.ret <- list(summary=x.summary,all=x.all)
    }
    return(x.ret)
  }

  #(5)#### RUN #####

  def.shp <- get.deforestation(spa.shp,def.shp)

  #(5)#### SAVE #####

  #save
  if(length(def.shp)==2){
    #outputs
    out.name <- paste0(output_folder,"/DEFORESTATION_stats.csv")
    write.csv(def.shp$summary,out.name,row.names=F)
    out.name <- paste0(output_folder,"/DEFORESTATION_alerts.shp")
    shapefile(def.shp$all,out.name,overwrite=T)
    #plots1
    out.name <- paste0(output_folder,"/DEFORESTACION_plot.jpg")
    area.p <- ggplot(def.shp$all@data,aes(x=names,y=area)) +
      geom_boxplot() +
      coord_flip() +
      xlab("Spatial units") +
      ylab("Area [ha]") +
      labs(title="DEFORESTATION: area by spatial unit")
    ggsave(out.name,area.p,width = size_jpg[1], height = size_jpg[2], units="cm")
  }else{
    warning(stats.df)
  }

  ###### close ######
  print("****script ends*****")

}
