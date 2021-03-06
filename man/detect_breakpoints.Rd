% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detect_breakpoints.R
\name{detect_breakpoints}
\alias{detect_breakpoints}
\title{Detect breakpoints from a GEE ts}
\usage{
detect_breakpoints(GEE_ts, change_threshold, segment_size, min_area = 4,
  calibration_mode = F, samples_shp, column_id_name, size_jpg = 1000,
  ncores = 1, output_folder)
}
\arguments{
\item{GEE_ts}{GEE time series filename. The GEE time series to analyze and detect breakpoints. It should be projected in UTM WGS84 coordinates system.}

\item{change_threshold}{Numerical. A breakpoint detection change threshold based in the standard deviation from mean of GEE ts. Values can vary from 0.001 (more sensibility to changes) to > 1 (less sensibility)}

\item{segment_size}{Numerical. The minimun segment size to consider in brakepoint detection. Depends of the number of observations in the GEE ts.}

\item{min_area}{Numerical. The minimun area of change areas detected by the algorithm}

\item{calibration_mode}{Logical. Activates the calibration mode of the algorithm. This generates a series of plots of breakpoints based in samples points (requires a points based shapefile of samples).}

\item{samples_shp}{Shapefile filename. A set of samples for evaluate calibration settings of the algorithm.}

\item{column_id_name}{String. The name of the column of the unique identifier of each sample.}

\item{size_jpg}{Numeric. A value which defines the size of plots (in pixels). Plots are based in an square, so one dimension defines its size.}

\item{ncores}{Numeric. A number defining how many cores will be used in processing.}

\item{output_folder}{Foldername. A folder to use for store plots.}
}
\value{
Plots outputs are stored in: *output_folder*. If calibration_mode is TRUE, a series of plots are create, otherwise a shapefile of detected breakpoints is generated.
}
\description{
This function detects breakpoints from GEE ts, detecting them and giving as output a shapefile with the correspondent dates of change. To calibrate the algorithm, it is possible to create a series of plots which aims to help to define its parameters (ie. change_threshold and segment_size). OSGeo4W is required to be installed and 'OSGeo4W.bat' declared as enviroment variable for its proper working.
}
\examples{
#Defining variables for breakpoint detection
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
}
