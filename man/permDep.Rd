\name{permDep}
\alias{permDep}
\title{Permutation test for general dependent truncation}
\description{
  Perform permutatoin test based on conditional or unconditional approach.
}
\usage{
permDep(trun, obs, permSize, cens, sampling = "cond", 
         kendallOnly = FALSE, minp1Only = FALSE, minp2Only = FALSE,
         nc = ceiling(detectCores() / 2), seed = NULL) 
}
\arguments{
  \item{trun}{Left truncation time}
  \item{obs}{Observed failure time, might be right-censored}
  \item{permSize}{Number of permutations}
  \item{cens}{The status indicator, normally 0 = censored, 1 = event. }
  \item{sampling}{A character string specifying the sampling method used
    in permutation. The following are permitted:
    \code{cond}: conditional permutation,
    \code{ucond}: unconditional permutation,
    \code{iscond}: importance sampling conditional permutation,  
    \code{isucond}: importance sampling unconditional permutation. 
  }
  \item{kendallOnly, minp1Only, minp2Only}{Logical values. Which tests
    will be performed? }
  \item{nc}{Number of cores. If nc > 1, parallel computing technique
    will be performed.}
  \item{seed}{An optional vector containing random seeds to be used in
    generating permutations. Random seeds will be used when left unspecified.}
}
%% \details{
%% }
\value{
  a list containing output with the following components:
  \item{obsKen}{observed Kendall's tau value}
  \item{obsP1}{observed minp1 p-value}
  \item{obsP2}{observed minp2 p-value}
  \item{obsTest1}{observed minp1 test statistics}
  \item{obsTest2}{observed minp2 test statistics}
  \item{permKen}{Kendall's tau values from permutations}
  \item{permP1}{minp1 p-values from permutations}
  \item{permP2}{minp2 p-values from permutations}
  \item{permDep1}{minp1 test statistics from permutations}
  \item{permDep2}{minp2 test statistics from permutations}
  \item{perm.data.Ken}{The `most' null permutation set based on
    Kendall's tau test}
  \item{perm.data.P1}{The `most' null permutation set based on
    minp1 test}
  \item{perm.data.P2}{The `most' null permutation set based on
    minp2 test}
}

\references{
  Chiou, S.H., Qian, J., and Betensky, R.A. (2017).
  Permutation Test for General Dependent Truncation
  \emph{Techreport} \bold{1}{1}: 1--10.
}

\seealso{\code{\link{permDep}}}
\examples{
simDat <- function(n) {
    k <- s <- 1
    tt <- xx <- yy <- cc <- delta <- rep(-1, n)
    while(k <= n){
      tt[k] <- runif(1, 0, 3.5)
      xx[k] <- 1.95 + 0.65 * (tt[k] - 1.25)^2 + rnorm(1, sd = 0.1)
      cc[k] <- runif(1, 0, 10)
      delta[k] <- (xx[k] <= cc[k])
      yy[k] <- pmin(xx[k], cc[k])
      s <- s + 1
      if(tt[k] <= yy[k]){k = k+1}
    }
    data.frame(list(trun = tt, obs = yy, delta = delta))
  }

set.seed(123)
dat <- simDat(50)
B <- 30

## Perform conditional permutation with Kendall's tau, minp1 and minp2
set.seed(123)
system.time(fit <- with(dat, permDep(trun, obs, B, delta, nc = 1)))
fit
}
\keyword{permDep}

