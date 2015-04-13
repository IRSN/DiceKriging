\name{km-class}
\docType{class}
\alias{km-class}
\alias{coef,km-method}
%%\alias{plot,km-method}
%%\alias{predict,km-method}
%%\alias{show,km-method}
%%\alias{simulate,km-method}

\title{ Kriging models class }

\description{	S4 class for kriging models.}

\section{Objects from the Class}{
To create a \code{km} object, use \code{\link{km}}. See also this function for more details.
}

\section{Slots}{
	 \describe{
    \item{\code{d}:}{Object of class \code{"integer"}. The spatial dimension. }
    \item{\code{n}:}{Object of class \code{"integer"}. The number of observations. }
    \item{\code{X}:}{Object of class \code{"matrix"}. The design of experiments. }
    \item{\code{y}:}{Object of class \code{"matrix"}. The vector of response values at design points. }
    \item{\code{p}:}{Object of class \code{"integer"}. The number of basis functions of the linear trend. }
    \item{\code{F}:}{Object of class \code{"matrix"}. The experimental matrix corresponding to the evaluation of the linear trend basis functions at the design of experiments.}
    \item{\code{trend.formula}:}{Object of class \code{"formula"}. A formula specifying the trend as a linear model (no response needed). }
    \item{\code{trend.coef}:}{Object of class \code{"numeric"}. Trend coefficients. }
    \item{\code{covariance}:}{Object of class \code{"covTensorProduct"}. See \code{\link{covTensorProduct-class}}. }
    \item{\code{noise.flag}:}{Object of class \code{"logical"}. Are the observations noisy? }
    \item{\code{noise.var}:}{Object of class \code{"numeric"}. If the observations are noisy, the vector of noise variances. }
    \item{\code{known.param}:}{Object of class \code{"character"}. Internal use. One of: \code{"None", "All", "CovAndVar"} or \code{"Trend"}. }
    \item{\code{case}:}{Object of class \code{"character"}. Indicates the likelihood to use in estimation (Internal use). One of: \code{"LLconcentration_beta", "LLconcentration_beta_sigma2", "LLconcentration_beta_v_alpha"}. }
    \item{\code{param.estim}:}{Object of class \code{"logical"}. \code{TRUE} if at least one parameter is estimated, \code{FALSE} otherwise. }
    \item{\code{method}:}{Object of class \code{"character"}. \code{"MLE"} or \code{"PMLE"} depending on \code{penalty}. }
    \item{\code{penalty}:}{Object of class \code{"list"}. For penalized ML estimation. }
    \item{\code{optim.method}:}{Object of class \code{"character"}. To be chosen between \code{"BFGS"} and \code{"gen"}.}
    \item{\code{lower}:}{Object of class \code{"numeric"}. Lower bounds for covariance parameters estimation. }
    \item{\code{upper}:}{Object of class \code{"numeric"}. Upper bounds for covariance parameters estimation. }
    \item{\code{control}:}{Object of class \code{"list"}. Additional control parameters for covariance parameters estimation. }
    \item{\code{gr}:}{Object of class \code{"logical"}. Do you want analytical gradient to be used ? }
    \item{\code{call}:}{Object of class \code{"language"}. User call reminder. }
    \item{\code{parinit}:}{Object of class \code{"numeric"}. Initial values for covariance parameters estimation. }
    \item{\code{logLik}:}{Object of class \code{"numeric"}. Value of the concentrated log-Likelihood at its optimum. }
    \item{\code{T}:}{Object of class \code{"matrix"}. Triangular matrix delivered by the Choleski decomposition of the covariance matrix. }
    \item{\code{z}:}{Object of class \code{"numeric"}. Auxiliary variable: see \code{\link{computeAuxVariables}}. }
    \item{\code{M}:}{Object of class \code{"matrix"}. Auxiliary variable: see \code{\link{computeAuxVariables}}. }
  }
}

\section{Methods}{
  \describe{
    \item{coef}{\code{signature(x = "km")} Get the coefficients of the \code{km} object. }
    \item{plot}{\code{signature(x = "km")}: see \code{\link{plot,km-method}}. }
    \item{predict}{\code{signature(object = "km")}: see \code{\link{predict,km-method}}. }
    \item{show}{\code{signature(object = "km")}: see \code{\link{show,km-method}}. }
    \item{simulate}{\code{signature(object = "km")}: see \code{\link{simulate,km-method}}. }
	 }
}


\author{ O. Roustant, D. Ginsbourger}

\seealso{ \code{\link{km}} for more details about slots and to create a \code{km} object, \code{\link{covStruct.create}} to construct a covariance structure, and \code{\link{covTensorProduct-class}} for the S4 covariance class defined in this package.}

\keyword{classes}
