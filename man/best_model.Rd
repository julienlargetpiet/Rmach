% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_fun.R
\name{best_model}
\alias{best_model}
\title{best_model}
\usage{
best_model(
  inpt_datf,
  Degree,
  Coeff_v = NA,
  Powers = NA,
  Mth_symb,
  Numrtr_v = NA
)
}
\arguments{
\item{inpt_datf}{is the input dataframe, first column for the x values and second column for the y values}

\item{Degree}{is a vector containing all the degrees. Each degree represents how many coefficients the model has.}

\item{Coeff_v}{is a list containing the vector containing the coefficients for each model. The first value of each coefficient vector is always the constant, so it is not linked to any math symbol}

\item{Powers}{is a list containing all the values associated with the math symbols of mth_symb list for each model. Because you can have multiple models in the function, so Powers is separated with the "-" separator between the different powers values for each model like in the examples}

\item{Mth_symb}{is a list containing the vector of the different math symbols linked to the coefficients from the second value}

\item{Numrtr_v}{is a list containing the different numerator values for each math symbol for each model, see examples}
}
\description{
Returns the best input models. The coefficient of the best model can be found with the poly_model function
}
\examples{

print(best_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), Degree=c(2, 2), Coeff_v=c("-", 32.5, -3, "-", 32.5, -5, "-"), Powers=c("-", 1, "-", 1, "-"), Mth_symb=c("-", "x", "-", "x", "-"), Numrtr_v=NA))

[1] 2

print(best_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), Degree=c(2, 2), Coeff_v=c("-", c(32.5, -5), "-", c(32.5, -3), "-"), Powers=c("-", c(1), "-", c(1), "-"), Mth_symb=c("-", c("x"), "-", c("1/x"), "-"), Numrtr_v=NA))

[1] 1

print(best_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), Degree=c(2, 2), Coeff_v=c("-", c(32.5, -5), "-", c(32.5, -3), "-"), Powers=c("-", c(1), "-", list(c(1:length(mtcars$wt))), "-"), Mth_symb=c("-", c("x"), "-", c("1/x"), "-"), Numrtr_v=NA))

[1] 1

print(best_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), Degree=c(2, 2), Coeff_v=c("-", c(32.5, -5), "-", c(32.5, -3, 2), "-"), Powers=c("-", c(1), "-", list(c(1:length(mtcars$wt)), 2), "-"), Mth_symb=c("-", c("x"), "-", c("list/x", "x"), "-"), Numrtr_v=c("-", list(c(1:length(mtcars$wt))), "-")))

#' [1] 1

}
