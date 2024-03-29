% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_select.R
\name{sim_select}
\alias{sim_select}
\title{Selective Simulation and Analysis for Logistic Regression}
\usage{
sim_select(
  Q = 500,
  sample_size = 2000,
  nsis = 2000/(4 * log(2000)),
  mediator_size = 2,
  p = 1000,
  a = rep(0, 2),
  r = 1,
  b = rep(0, 2),
  alpha1 = 0,
  Rsmethod = "McFadden",
  nm2 = NULL,
  nm1 = NULL,
  trainrate = 0.75,
  cov = NULL,
  FDR = 0.2,
  select = TRUE,
  gamma = 3,
  iter = FALSE,
  g = NULL,
  l = NULL
)
}
\arguments{
\item{Q}{Number of simulations.}

\item{sample_size}{Number of observations per simulation.}

\item{nsis}{Number of selections in Stepwise Information Criterion.}

\item{mediator_size}{Number of mediators.}

\item{p}{Total number of predictors.}

\item{a}{Coefficients for mediator interaction.}

\item{r}{Coefficient for the primary predictor.}

\item{b}{Coefficients for the mediators.}

\item{alpha1}{Intercept for the logistic model.}

\item{Rsmethod}{Method for R-squared calculation.}

\item{nm2}{Number of noise mediators to include.}

\item{nm1}{Number of additional noise mediators.}

\item{trainrate}{Proportion of data used for training.}

\item{cov}{Optional covariance matrix for data simulation.}

\item{FDR}{False Discovery Rate for mediator selection.}

\item{select}{Boolean to activate selection process.}

\item{gamma}{Gamma parameter for penalization in mediator selection.}

\item{iter}{Boolean to iterate selection process.}

\item{g}{Additional parameter for `sim_logistic_cov`.}

\item{l}{Additional parameter for `sim_logistic_cov`.}
}
\value{
A list of simulation results including R-squared values, selection performance metrics, and coefficients.
}
\description{
This function conducts selective simulations to evaluate the impact of various mediators and predictors on logistic regression outcomes. It allows for mediator selection based on specified criteria and calculates performance metrics.
}
\examples{
results <- sim_select()
}
