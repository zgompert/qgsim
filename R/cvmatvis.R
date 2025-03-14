#' Visualize variance/covariance matrices
#'
#' Visualize the G and omega variance/covariance matrices generated in qgsim
#' @param h2 the trait heritabilities, assumed to be the same for both traits (must be between 0 and 1).
#' @param Gcor the genetic correlation between the pair of traits (must be between -1 and 1).
#' @param omega11 denotes the intensity of selection (curvature of the adaptive landscape) with respect to trait 1, larger values result in weaker selection [default = 1].
#' @param omega22 denotes the intensity of selection (curvature of the adaptive landscape) with respect to trait 2, larger values result in weaker selection [default = 1].
#' @param omegaCor denotes the strength of correlational selection, larger values denote stronger selection for combinations of trait 1 and trait 2 (must be between -1 and 1, 0 denotes independent selection on each trait) [default = 0].
#' @export
#'
#' @return The plot showing 95% CI for both generated matrices
var_covar_mat_vis<-
function(h2=0.5, Gcor=0.2, omega11=1, omega22=1, omegaCor=0) {
	## construct the Gmatrix from the trait heritabilities, assumed to be the same for both traits,
	## and genetic correlations
	## this assumes the phenotypic variance is 1
	G<-matrix(c(h2,Gcor,Gcor,h2),nrow=2,byrow=TRUE)
	## construct omega matrix, curvature of the adaptive landscape
	omega<-matrix(c(omega11,omegaCor,omegaCor,omega22),nrow=2,byrow=TRUE)

  gen_ellipse_vals<-
  function(sigma) {
    # set theta
    theta<-c(0, 0)

    # build out starting points (unit circle)
    n_points<-100
    xy <- cbind(sin(seq(0, 2 * pi, length.out = n_points)),
                cos(seq(0, 2 * pi, length.out = n_points)))

    # scale dimensions
    ev<-eigen(sigma)
    xy[, 1] <- xy[, 1] * 1
    xy[, 2] <- xy[, 2] * sqrt(min(ev$values) / max(ev$values))

    # rotate
    phi <- atan(ev$vectors[2, 1] / ev$vectors[1, 1])
    R <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2)
    xy <- tcrossprod(R, xy)

    # compute quantile
    chi_val<-qchisq(0.95, 2)*max(ev$values)

    # generate ellipse lines
    r<-sqrt(chi_val)
    # lines(r * xy[1, ] + theta[1], r * xy[2, ] + theta[2], lty = 1)
    list(x=(r*xy[1, ]+theta[1]), y=(r*xy[2, ]+theta[2]))
  }

  ev_G<-gen_ellipse_vals(G)
  ev_omega<-gen_ellipse_vals(omega)

  par(mar = c(4.5, 4, .5, .5))
  plot(c(-5, 5), c(-5, 5), type="n", xlab="x", ylab="y")

  legend(x="topleft",legend=c("G matrix", "Omega matrix"), fill=c("blue", "orange"))
  lines(ev_G$x, ev_G$y, lty=1, col="blue")
  lines(ev_omega$x, ev_omega$y, lty=1, col="orange")
}
