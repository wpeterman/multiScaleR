% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_scale_fn.R
\name{kernel_scale_fn}
\alias{kernel_scale_fn}
\title{Kernel scaling function}
\usage{
kernel_scale_fn(par, kernel_inputs, fitted_mod, mod_return = NULL)
}
\arguments{
\item{par}{list of parameters}

\item{kernel_inputs}{object created from `kernel_prep`}

\item{fitted_mod}{fitted model object}

\item{mod_return}{Default: NULL}
}
\value{
Either estimated parameters or the fitted model using provided parameters
}
\description{
Function for internal use with optim
}
\details{
For internal use
}
\keyword{internal}
