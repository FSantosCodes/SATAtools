% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_deforestation_alerts.R
\name{filter_deforestation_alerts}
\alias{filter_deforestation_alerts}
\title{Filter deforestation alerts}
\usage{
filter_deforestation_alerts(spatial_units_shp, column_units_name,
  deforestation_shp, column_date_name, column_area_name, min_area = NA,
  size_jpg = c(10, 10), output_folder)
}
\arguments{
\item{spatial_units_shp}{Shapefile filename. This is the shapefile of spatial units to analyze with SATA alerts. It should be projected in UTM WGS84 coordinates system.}

\item{column_units_name}{String. The name of the column indicating the names of the spatial units to analize. Case-sensitive applies.}

\item{deforestation_shp}{Shapefile filename. The shapefile of deforestation alerts from SATA (as polygons).}

\item{column_date_name}{String. The name of the column indicating the date of alerts. Case-sensitive applies.}

\item{column_area_name}{String. The name of the column indicating the area of alerts. Case-sensitive applies.}

\item{min_area}{Numeric. A value defining the minimun area in filtering. NA ignore filtering.}

\item{size_jpg}{Numeric. Two values which define the width and height of boxplot (in cms).}

\item{output_folder}{Foldername. A folder to use for store alerts filtered.}
}
\value{
Filtered alerts and plots are stored in: *output_folder*.
}
\description{
This function filters SATA deforestation alerts. It uses a shapefile of spatial units to generate plots and visualize some statistics of alerts by spatial units.
}
\examples{
#Defining variables from spatial units shapefile
spatial_units_shp <- "C:/DATA/GAF/demo/shp/unidades_espaciales/provincias.shp"
column_units_name <- "dpa_despro" #case sensitive

#Defining variables from deforestation shapefile
deforestation_shp <- "C:/DATA/GAF/demo/shp/alertas/deforestation_2016_2017.shp"
column_date_name <- "inicio" #case sensitive
column_area_name <- "área" #case sensitive
min_area <- 1 #with same units as area column but here ignored

#Defining variables of ouputs
size_jpg <- c(40,20) #first value defines width and second height
output_folder <- "C:/DATA/GAF/demo/estadisticas"

#running the function
filter_deforestation_alerts(spatial_units_shp = spatial_units_shp,
                            column_units_name = column_units_name,
                            deforestation_shp = deforestation_shp,
                            column_date_name = column_date_name,
                            column_area_name = column_area_name,
                            min_area = min_area,
                            size_jpg = size_jpg,
                            output_folder = output_folder)
}
