\name{qgsim}
\alias{qgsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qgsim(npops = 2, mig = 0, Ne = 500, theta0 = 0, z0 = 0, h2 = 0.5, Gcor = 0.2, omega11 = 1, omega22 = 1, omegaCor = 0, model = "Brownian", ngens = 100, tsd = 0.02)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npops}{
%%     ~~Describe \code{npops} here~~
}
  \item{mig}{
%%     ~~Describe \code{mig} here~~
}
  \item{Ne}{
%%     ~~Describe \code{Ne} here~~
}
  \item{theta0}{
%%     ~~Describe \code{theta0} here~~
}
  \item{z0}{
%%     ~~Describe \code{z0} here~~
}
  \item{h2}{
%%     ~~Describe \code{h2} here~~
}
  \item{Gcor}{
%%     ~~Describe \code{Gcor} here~~
}
  \item{omega11}{
%%     ~~Describe \code{omega11} here~~
}
  \item{omega22}{
%%     ~~Describe \code{omega22} here~~
}
  \item{omegaCor}{
%%     ~~Describe \code{omegaCor} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{ngens}{
%%     ~~Describe \code{ngens} here~~
}
  \item{tsd}{
%%     ~~Describe \code{tsd} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (npops = 2, mig = 0, Ne = 500, theta0 = 0, z0 = 0, h2 = 0.5, 
    Gcor = 0.2, omega11 = 1, omega22 = 1, omegaCor = 0, model = "Brownian", 
    ngens = 100, tsd = 0.02) 
{
    ntraits <- 2
    G <- matrix(c(h2, Gcor, Gcor, h2), nrow = 2, byrow = TRUE)
    omega <- matrix(c(omega11, omegaCor, omegaCor, omega22), 
        nrow = 2, byrow = TRUE)
    theta0 <- matrix(theta0, nrow = npops, ncol = ntraits, byrow = TRUE)
    z0 <- matrix(z0, nrow = npops, ncol = ntraits, byrow = TRUE)
    z_mat <- vector("list", npops)
    theta_mat <- vector("list", npops)
    for (j in 1:npops) {
        z_mat[[j]] <- matrix(z0[j, ], nrow = ngens, ncol = ntraits)
        theta_mat[[j]] <- matrix(theta0[j, ], nrow = ngens, ncol = ntraits)
    }
    for (i in 2:ngens) {
        for (j in 1:npops) {
            z <- z_mat[[j]][i - 1, ]
            theta <- theta_mat[[j]][i - 1, ]
            if (model == "Brownian") {
                theta <- theta + rnorm(ntraits, 0, tsd)
            }
            z <- z + (theta - z) \%*\% solve(omega) \%*\% G + rmnorm(1, 
                varcov = G/Ne)
            z_mat[[j]][i, ] <- z
            theta_mat[[j]][i, ] <- theta
        }
    }
    o <- list(z = z_mat, theta = theta_mat)
    return(o)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
