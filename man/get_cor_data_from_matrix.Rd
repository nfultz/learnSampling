% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/samplingFunctions_r.R
\name{get_cor_data_from_matrix}
\alias{get_cor_data_from_matrix}
\title{Calculate a data frame where columns correlate as specified by a correlation matrix}
\usage{
get_cor_data_from_matrix(r.matrix, n)
}
\arguments{
\item{r.matrix}{Desired correlation matrix}

\item{n}{Sample size}
}
\value{
Data frame with sample data
}
\description{
Calculate a data frame where columns correlate as specified by a correlation matrix
}
\examples{
N <- 600
M <- matrix(c(1.00, 0.60, 0.30, 0.30,
            0.60, 1.00, 0.00, 0.60,
            0.30, 0.00, 1.00, 0.00,
            0.30, 0.60, 0.00, 1.00), nrow=4, ncol=4)
my.data <- get_cor_data_from_matrix(M,N)
library(apaTables)
apa.cor.table(my.data)
}

