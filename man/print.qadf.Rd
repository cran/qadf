\name{print.qadf}
\alias{print.qadf}
\alias{summary.qadf}
\title{Print and Summary Methods for qadf Objects}
\usage{
\method{print}{qadf}(x, digits = 4L, ...)

\method{summary}{qadf}(object, ...)
}
\arguments{
\item{x}{An object of class \code{"qadf"}.}

\item{digits}{Integer. Number of significant digits for display.}

\item{...}{Further arguments (ignored).}

\item{object}{An object of class \code{"qadf"}.}
}
\value{
Invisibly returns the input object.
}
\description{
Print and summary methods for objects of class \code{"qadf"} returned
by \code{\link{qadf}}.
}
\examples{
set.seed(1)
y <- cumsum(rnorm(80))
res <- qadf(y, tau = 0.5)
print(res)
summary(res)
}
