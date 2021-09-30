\name{TopGearMPG}
\alias{TopGearMPG}
\docType{data}
\title{
Top Gear fuel consumption data
}
\description{
The data set contains information on fuel consumption of cars featured on the
website of the popular BBC television show \emph{Top Gear}, together with
car information and performance measures.  Only cars with complete information
are included in the data set.
}
\usage{data("TopGearMPG")}
\format{
  A data frame with 255 observations on the following 11 variables.
  \describe{
    \item{\code{Maker}}{factor; the car maker.}
    \item{\code{Model}}{factor; the car model.}
    \item{\code{Type}}{factor; the exact model type.}
    \item{\code{MPG}}{numeric; the combined fuel consuption (urban + extra
    urban; in miles per gallon).}
    \item{\code{Cylinders}}{numeric; the number of cylinders in the engine.}
    \item{\code{Displacement}}{numeric; the displacement of the engine (in cc).}
    \item{\code{BHP}}{numeric; the power of the engine (in bhp).}
    \item{\code{Torque}}{numeric; the torque of the engine (in lb/ft).}
    \item{\code{Acceleration}}{numeric; the time it takes the car to get from 0
    to 62 mph (in seconds).}
    \item{\code{TopSpeed}}{numeric; the car's top speed (in mph).}
    \item{\code{Weight}}{numeric; the car's curb weight (in kg).}
  }
}
\source{
The data were scraped from \code{http://www.topgear.com/uk/} on 2014-02-24.
}
\examples{
data("TopGearMPG")
plot(TopGearMPG[, -(1:3)])
}
\keyword{datasets}
