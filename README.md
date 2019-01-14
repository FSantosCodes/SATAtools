# SATAtools

Tools for SATA (Sistema de Alertas Tempranas del Sistema Nacional de Monitoreo de Bosques del Ecuador) 

## Instalación

Para instalar la librería, es necesario instalar otras dependencias y tener la última version del OSGEO4W. Con esto se ejecutan las primeras líneas para la instalación:
```
#install required libraries
check.libs <- c("rgdal","raster","sp","foreach","doParallel",
                "lubridate","rgeos","scales","ggplot2","changepoint",
                "gtools","plyr","devtools","rlang")
get.libraries <- function(libraries){
  for (i in libraries){
    if(!require(i,character.only=T)){
      install.packages(i)
      library(i,character.only=T)}
  }
}
get.libraries(check.libs)
```
y la siguiente instala el SATAtools:
```
#install SATAtools
install_github("FSantosCodes/SATAtools")
library(SATAtools)
```
Los ejemplos que se muestran a continuación se basan en los resultados de la descarga de los algoritmos corridos en el GEE. El pos-procesamiento de estos son una muestra de los potenciales usos de esta herramienta. Las rutinas de GEE se muestran a continuación:
- Landsat 4-8 compositing
https://code.earthengine.google.com/6a2c42cdbebe701dd8f3b31287da29c7
- Sentinel 1 compositing
https://code.earthengine.google.com/777ef0b8a4aaa638b4d923fc45f04efd
- Sentinel 2 compositing
https://code.earthengine.google.com/612a604f0baf41581d225f4e8b2b06c0

### Filtrado de alertas

Las coberturas usadas en estos ejemplos usan las alertas de deforestacion del SATA:

```
##### FILTER DEFORESTATION #####

#Defining variables from spatial units shapefile
spatial_units_shp <- "C:/DATA/GAF/demo/shp/unidades_espaciales/provincias.shp"
column_units_name <- "dpa_despro" #case sensitive

#Defining variables from deforestation shapefile (output from SATA or 'detect_breakpoints' function)
deforestation_shp <- "C:/DATA/GAF/demo/shp/alertas/deforestation_2016_2017.shp"
column_date_name <- "inicio" #case sensitive
column_area_name <- "Ã¡rea" #case sensitive
min_area <- NA #with same units as area column but here ignored

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
```

aquí se lo hace con las alertas de incendios del MODIS:

```
#Defining variables from spatial units shapefile
spatial_units_shp <- "C:/DATA/GAF/demo/shp/unidades_espaciales/provincias.shp"
column_units_name <- "dpa_despro" #case sensitive

#Defining variables from MODIS shapefile
MODIS_shp <- "C:/DATA/GAF/demo/shp/alertas/pts_calor_MODIS_30d.shp"
column_date_name <- "acq_date" #case sensitive
column_frp_name <- "frp" #case sensitive
column_confidence_name <- "confidence" #case sensitive
min_frp <- 10 #with same units as area column
min_confidence <- 50 #with same units as area column

#Defining variables of ouputs
size_jpg <- c(40,20) #first value defines width and second height
output_folder <- "C:/DATA/GAF/demo/estadisticas"

#running the function
filter_MODIS_alerts(spatial_units_shp = spatial_units_shp,
                    column_units_name = column_units_name,
                    MODIS_shp = MODIS_shp,
                    column_date_name = column_date_name,
                    column_frp_name = column_frp_name,
                    column_confidence_name = column_confidence_name,
                    min_frp = min_frp,
                    min_confidence = min_confidence,
                    size_jpg = size_jpg,
                    output_folder = output_folder)
```

## Validación de alertas

Usando un shapefile de las alertas y sus fechas, esta rutina toma cada alerta (coodenadas, fecha) y la grafica usando la serie temporal generada con las rutinas del GEE. Estos graficos nos permiten evaluar el comportamiento antes/después de la alerta y observar la evidencia con que el algoritmo detectó el cambio. 

```
#Defining variables for alerts shapefile
alerts_shp <- "C:/DATA/GAF/demo/estadisticas/MODIS_alerts.shp"
column_dates_name <- "date"

#Defining variables for GEE ts files (normally as a TIF). Take into account that should exist in the same folder its version as CSV (dates files)
red_channel_GEE_ts <- "C:/DATA/GAF/demo/sentinel1/gee_outputs/coca_S1_ts_VV_P50_I14.tif"
green_channel_GEE_ts <- "C:/DATA/GAF/demo/sentinel1/gee_outputs/coca_S1_ts_VH_P50_I14.tif"
blue_channel_GEE_ts <- "C:/DATA/GAF/demo/sentinel1/gee_outputs/coca_S1_ts_VV_P50_I14.tif"

#other variables required
observations_evaluate <- 5 #could be less if processing is slow
buffer_distance <- 500 #consider that is in meters and measured from the alert point
stretch_perc <- 5 #recommendable if plots do not have good contrast
size_jpg <- 1000 #plot is based in an square, so one dimension defines its size
ncores <- 3 #modify if your resources are less than this

#output folder to save plots
output_folder <- "C:/DATA/GAF/demo/test"

#running the function
validate_alerts(alerts_shp = alerts_shp,
                column_dates_name = column_dates_name,
                red_channel_GEE_ts = red_channel_GEE_ts,
                green_channel_GEE_ts = green_channel_GEE_ts,
                blue_channel_GEE_ts = blue_channel_GEE_ts,
                observations_evaluate = observations_evaluate,
                buffer_distance = buffer_distance,
                stretch_perc = stretch_perc,
                size_jpg = size_jpg,
                ncores = ncores,
                output_folder = output_folder)
```
## Detección de cambios
Ejecución de rutinas de deteccion de cambios a partir de series temporales del GEE. La función calcula cambios en el promedio de la serie temporal estandarizada. El output que genera es un shapefile filtrado por área y fecha de cambio. Requiere que el OSGEO4W este previamente instalado.

##### DETECT BREAKPOINTS #####

```
##Defining variables for breakpoint detection
GEE_ts <- "C:/DATA/GAF/demo/sentinel1/gee_outputs/coca_S1_ts_VH_P50_I14.tif"
change_threshold <- 10
segment_size <- 3
min_area <- 4

#For calibration mode
samples_shp <- "C:/DATA/GAF/demo/sentinel1/samples/muestras.shp"
column_id_name <- "clase"
size_jpg <- 3000

#Variables for processing and outputs
ncores <- 3
output_folder <- "C:/DATA/GAF/demo/sentinel1/breakpoints"

#running the function in CALIBRATION mode
detect_breakpoints(GEE_ts = GEE_ts,
                   change_threshold = change_threshold,
                   segment_size = segment_size,
                   min_area = min_area,
                   calibration_mode = T,
                   samples_shp = samples_shp,
                   column_id_name = column_id_name,
                   size_jpg = size_jpg,
                   ncores = ncores,
                   output_folder = output_folder)

#running the function in PROCESSING mode
detect_breakpoints(GEE_ts = GEE_ts,
                   change_threshold = change_threshold,
                   segment_size = segment_size,
                   min_area = min_area,
                   calibration_mode = F,
                   column_id_name = column_id_name,
                   ncores = ncores,
                   output_folder = output_folder)

```

## Contribuciones

Por favor revíse el  [DESCRIPTION](https://github.com/FSantosCodes/SATAtools/blob/master/DESCRIPTION) para detalles sobre la implementación del SATAtools.

## Autores

* **Fabián Santos** - *Programación y concepto* - (https://www.researchgate.net/profile/Fabian_Santos2)

## Licencia

Este proyecto tiene licencia GNU-3 (General Public License v3.0) - mire el link (https://www.gnu.org/licenses/gpl-3.0.en.html) para más detalles.

## Agradecimientos

* GAFAG, Ministerio del Ambiente, KFW, ZFL.
