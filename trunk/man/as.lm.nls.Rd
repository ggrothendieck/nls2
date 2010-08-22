\name{as.lm.nls}
\alias{as.lm.nls}
\alias{as.lm}
\title{
Linearize an nls object.
}
\description{
Given an \code{"nls"} object, this method produces an \code{"lm"} object which approximates it.  
}
\usage{
\method{as.lm}{nls}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     \code{"nls"} object.
}
  \item{\dots}{
%%     Currently not used.
}
}
\details{
Given an object of class \code{"nls"}, \code{as.lm.nls} produces an \code{"lm"} object tangent to the \code{"nls"} object which approximates it.
}
\value{
The returned value is an \code{"lm"} class object.
}
\references{
Bates, DM and DG Watts, Nonlinear Regression Analysis and its Applications, Wiley, 1988.
}

\seealso{
\code{\link{nls}}.
}
\examples{
# data is from ?nls
DNase1 <- subset(DNase, Run == 1)
fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)

# these give same result
vcov(fm1DNase1)
vcov(as.lm(fm1DNase1))

# nls confidence and prediction intervals based on asymptotic approximation
# are same as as.lm confidence intervals.
predict(as.lm(fm1DNase1), interval = "confidence")
predict(as.lm(fm1DNase1), interval = "prediction")

}
\keyword{nonlinear}
\keyword{regression}
\keyword{models}
