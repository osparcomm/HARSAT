% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/information_functions.R
\name{get_basis_default}
\alias{get_basis_default}
\title{The default function for generating a basis}
\usage{
get_basis_default(data, info)
}
\arguments{
\item{data}{the data}

\item{info}{the information object}
}
\description{
This default is very simplistic, but works in all cases for sediment and water
most_common was used by AMAP in their mercury assessment and takes the most
common basis reported in each station, species (biota), matrix and
determinand_group combination
}
