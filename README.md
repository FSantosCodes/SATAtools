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
Los ejemplos que se muestran a continuación se basan en los resultados de la descarga de los algoritmos corridos en el GEE. El pos-procesamiento de estos son los que constituyen esta muestra de los potenciales usos de esta herramienta.

### Filtrado de alertas

Las coberturas usadas en estos ejemplos usan las bases de detección de deforestacion del SATA:

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

esta otras con la base del MODIS y VIRRS:

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



## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
