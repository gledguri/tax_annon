% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_annon.R
\name{seq_dev}
\alias{seq_dev}
\title{Sequence developer}
\usage{
seq_dev(k, i)
}
\arguments{
\item{k}{the length of the sequence of numbers}

\item{i}{the loop factor}
}
\value{
A vector of a sequence of numbers k-long
}
\description{
This function creates a sequence of numbers from 1:k when i=0 and k+1:2k 
This is very handy when you want to slice data into fragments
}
\examples{
seq_dev(100,0)
seq_dev(100,1)
for(i in 0:3) vector[i] <- seq_dev(100,i)
}
