% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/propensityModelFunction.r
\name{plotPropensityModel}
\alias{plotPropensityModel}
\title{Plot Propensity Model}
\usage{
plotPropensityModel(
  CountryAnalysis,
  indep_vars,
  treatment_binary,
  outcome_var,
  outcome_label,
  sample_density,
  lower_dist_bound = 0.01,
  upper_dist_bound = 0.5,
  maximum_match_diff = 0.9,
  debug = FALSE
)
}
\arguments{
\item{CountryAnalysis}{Takes a dataset with a 'distance' column...}

\item{indep_vars}{List of all independent variables to be used to model propensity and outcome}

\item{outcome_var}{Outcome variable (i.e., under 5 child mortality)}

\item{sample_density}{Total number of runs}

\item{lower_dist_bound}{Smallest distance (in decimal degrees) to calculate results for. Defaults to 0.01.}

\item{upper_dist_bound}{Largest distance (in decimal degrees) to calculate results for. Defaults to 0.5.}

\item{maximum_match_diff}{From 0 to 1, the largest difference in match quality allowed. Defaults to 0.9.}

\item{debug}{Set to TRUE if testing on specific original data. Defaults to FALSE. (to be removed)}
}
\description{
Matches samples of groups of treated and control observations that have a spatial distribution (i.e., latitude and longitude coordinates), so that matched groups have (1) similar covariate distributions \link{via a propensity score matching}, and (2) do not result in paired matches being geographically proximate.  The package also provides a distance-decay estimate of the treatment effect (i.e., the impact a treatment had on a defined outcome) following the procedure outlined in Runfola and Batra et al. (2020) \url{https://doi.org/10.3390/su12083225}.
}
\examples{
plotPropensityModel(read.csv("../CountryAnalysis.csv"), independent_vars, x,x,x, 30*16)
}