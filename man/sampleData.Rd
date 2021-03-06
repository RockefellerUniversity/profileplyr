% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allClasses.R
\name{sampleData}
\alias{sampleData}
\alias{sampleData,profileplyr-method}
\alias{sampleData<-}
\alias{sampleData<-,profileplyr,DataFrame-method}
\title{Retrieve and set sample data in profileplyr object}
\usage{
sampleData(object = "profileplyr")

\S4method{sampleData}{profileplyr}(object = "profileplyr")

sampleData(object) <- value

\S4method{sampleData}{profileplyr,DataFrame}(object) <- value
}
\arguments{
\item{object}{A profileplyr object}

\item{value}{DataFrame of sample information}
}
\value{
A DataFrame containing sample data

A DataFrame containing sample data to replace current sample data
}
\description{
Retrieve and set sample data in profileplyr object
}
\examples{
example <- system.file("extdata", "example_deepTools_MAT", package = "profileplyr") 
object <- import_deepToolsMat(example) 
sampleData(object)
sampleData(object)$scale <- c(1,10,1)
}
