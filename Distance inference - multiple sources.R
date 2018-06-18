# Global settings

# In your root directory put the notebooks (*.ipynb) into a directory called "notebook/" and the R codes (*.R) into a directory called "Rcode/".

rootDir <- "/Users/coryn/local/astrostats_2017_ESAC/" # directory containing notebook, Rcode, and data directories
rootDir <- "~/Documents/edas/2-Astrostats-2017-ESAC_Baille-Jones/"

# install.packages("PolynomF")
# install.packages(c("fields", "RColorBrewer"))

library("PolynomF")

setwd(paste(rootDir, "Rcode", sep=""))
conv <- pi/180 
source("general_functions.R")
source("distance_functions.R")
setwd("../")
library(MASS) # for truehist
library(fields) # for image.plot
library(RColorBrewer) # for colorRampPalette
mypalette <- colorRampPalette(brewer.pal(9, "Greys"), space="rgb", interpolate="linear", bias=2.5)
mycols <- mypalette(64)

## Read and check GDR1 data

# Gaia archive query for gdr1set03.csv (164 Pleiades members from GDR1 release paper)
# 
# select gaia.ra, gaia.dec, gaia.parallax, gaia.parallax_error, gaia.pmra, gaia.pmdec
# from gaiadr1.tgas_source as gaia
# where contains(point('ICRS',gaia.ra,gaia.dec),circle('ICRS',56.75,24.12,5)) = 1
# and sqrt(power(gaia.pmra-20.5,2)+power(gaia.pmdec+45.5,2)) < 6.0

# parallaxes are in mas, angles in degrees
dat <- read.csv(paste(rootDir, "data/gdr1set03.csv", sep=""), sep=",")
# To prevent hitting numerical limits in likelihood calculations, retain parallaxes in mas. 
# Thus distances are now kpc.
# For testing use smaller data set
#dat <- dat[1:25,]
# To investigate impact of larger uncertainties, scale parallax SD:
#dat$parallax_error <- dat$parallax_error*20

par(mfrow=c(2,2), mar=c(5,5,0.5,1), oma=c(0.1,0.1,0.5,0.1), mgp=c(2.2,0.8,0), cex=1.0) 
# Show source positions on sky along with proper motions
plot(dat$ra, dat$dec, xlim=rev(range(dat$ra)), xlab="RA [deg]", ylab="Dec [deg]")
sf <- 0.03 # scale factor for visualizing PM arrows
arrows(dat$ra, dat$dec, dat$ra+sf*dat$pmra, dat$dec+sf*dat$pmdec, length=0.05, lw=0.5)
# Histogram of fractional parallax uncertainties
truehist(dat$parallax_error/dat$parallax, nbins=25)
# Histogram of inverse parallaxes, for orientation.
truehist(1/dat$parallax, nbins=25, prob=FALSE, ylab="number")
segments(x0=1/dat$parallax, y0=0, x1=1/dat$parallax, y1=10, lw=0.3)

## Naive parallax combination

# It is incorrect to estimate the cluster distance via the mean of the inverse parallaxes
# (image some parallaxes were negative), and is meaningless when the fractional uncertainties are large.
# But we do it anyway to get a feel for the data:
cat("Mean and SD of inverse parallaxes: ", mean(1/dat$parallax), sd(1/dat$parallax), "kpc\n")
# Note that this is an SD in the sense of a spread, not a standard uncertainty in the mean

Nstar <- nrow(dat)
# Take inverse of mean of parallaxes and estimate variance by first order Taylor expansion
mpar <- mean(dat$parallax)
rNaive1 <- 1/mpar
rNaive1SD <- (1/Nstar)*sqrt(sum(dat$parallax_error^2)) / mpar^2
cat("Inverse of mean parallaxes: ", rNaive1, "+/-", rNaive1SD, "kpc\n")
# This gives an estimate of the distance to the cluster, together with the standard deviation in this mean.
cat("SD of inverse mean parallax in the sense of spread:", sd(dat$parallax_error)/mpar^2, "kpc\n")

# Take inverse of variance weighted mean of parallaxes 
# and again estimate variance by first order Taylor expansion
wfac <- sum(1/dat$parallax_error^2)
wmpar <- sum(dat$parallax/dat$parallax_error^2)/wfac
rNaive2 <- 1/wmpar
rNaive2SD <- 1/(wmpar^2 * sqrt(wfac))
cat("Inverse of weighted mean of parallaxes: ", rNaive2, "+/-", rNaive2SD, "kpc\n")

## Estimate distance to cluster (no correlation or cluster spread)

# Now adopt multivariate likelihood, assuming that each parallax is an independent estimate of the parallax
# of the cluster, i.e. 
# - no correlations between parallax measurements, and
# - no intrinsic spread in stellar parallaxes.
# Posterior densities can exceed numerical limits, so check here before trying to plot.
rlen <- 1e3
r <- seq(from=0.130, to=0.140, length.out=1e3)
#r <- seq(from=0.100, to=0.200, length.out=1e3) # with x20 uncertainties
# Inspect density values to make sure we're not running into numerical limits
(dense <- ud.distpost3multi(r=r, w=dat$parallax, wsd=dat$parallax_error, rlen=rlen))[1:10]
range(dense)

mom <- pdfmom(dense, r)
cat("Posterior mean, SD for rc =", mom$mean, mom$sd, "kpc\n")
plot(r, dense/mom$Z, type="l", xlab="distance r [kpc]", ylab="P(r | {w},{wsd},rlen) * const")
abline(v=c(rNaive1, rNaive2), col="red")
# Note that the mode is very close to rNaive2, because data are highly 
# informative compared to prior, so max. posterior is max. likelihood

## Invent simple model for spatial parallax correlations

# Function parcor computes correlations between parallaxes for sources according to
# amp*exp(-sep/len), where sep is the angular separation on the sky of two sources,
# 0<=amp<=1, and len>=0 is the length scale.
parcor(dat$ra[1], dat$dec[1], dat$ra[2], dat$dec[2], amp=0.5, len=0.5)
# Function parcovmat compute parallax covariance matrix for given sources plus parameters for parcor
(mat <- parcovmat(dat=dat[1:3,], amp=0.5, len=0.5))
# d.likemulti computes the density of an N-dimensional likelihood (i.e. for N sources) with
# specified covariance matrix (parcovmat). For independent parallaxes, use parcovmat=NULL.
d.likemulti(w=dat$parallax[1:3], r=0.12, wsd=dat$parallax_error[1:3], parcovmat=mat)

## Naive parallax combination with correlations

# This would just be maximum likelihood, which is equal to maximum posterior with a uniform (non-truncated) prior. 
# So we'll just skip this and get to the more general posterior solution.

## Distance estimation with correlations

# ud.distpost3multi accommodates parallax covariances when its argument parcovmat!=NULL
rlen <- 1 # kpc!
# test the function: the following two numbers should be identical
ud.distpost3multi(r=0.120, w=dat$parallax[1:4], wsd=dat$parallax_error[1:4], rlen=rlen)
ud.distpost3multi(r=0.120, w=dat$parallax[1:4], wsd=dat$parallax_error[1:4], 
parcovmat=diag(dat$parallax_error[1:4]^2), rlen=rlen)

# compute covariance matrix 
mat <- parcovmat(dat=dat, amp=0.5, len=0.5)

# compute posterior (scaled to have maximum=1) on a grid and plot
r <- seq(from=0.125, to=0.140, length.out=1e3)
#r <- seq(from=0.050, to=0.250, length.out=1e3) # with x20 uncertainties
dense <- ud.distpost3multi(r=r, w=dat$parallax, wsd=dat$parallax_error, parcovmat=mat, rlen=rlen)
mom <- pdfmom(dense, r)
cat("Posterior mean, SD for r =", mom$mean, mom$sd, "kpc\n")
plot(r, dense/mom$Z, type="l", xlab="distance r [kpc]", ylab="P(r | {w},{wsd},{wCov},rlen)")
# overplot result (obtained above) without correlations
dense <- ud.distpost3multi(r=r, w=dat$parallax, wsd=dat$parallax_error, rlen=rlen)
mom <- pdfmom(dense, r)
cat("Posterior mean, SD for r =", mom$mean, mom$sd, "kpc\n")
lines(r, dense/mom$Z, col="red")

## Infer cluster distance and size (no correlations)

# We now want posterior P*(rc,sc|{w},{wsd}), where rc and sc are cluster distance and size respectively. The cluster center in (RA, Dec) is assumed known. I nominally assume the true stellar distances are drawn from a 3D isotropic Gaussian of mean rc and standard deviation sc. In general the likelihood is a marginalization
# over the N unknown distances, i.e. an N-dimensional integral. Here I neglect any correlations between the parallax measurements to reduce this to a product of N one-dimensional integrals. I then proceed in three different ways, which corresponds to three different simplifications to this likelihood. 
# 1. Nelgect the angular size of the cluster. I can therefore assume that the true distances of the stars from the cluster center (z) are along the l.o.s, and thus drawn from a 1D Gaussian with mean 0 and stdev sc. To simplify the integral (which is otherwise intractable), I further perform a binomial approximation of (w-1/r) with r=r_c+x. The integral is can then be expressed in terms of an error function, which is fast. This is done by d.likecluster1.
# 2. As 1, but now without the binomial approximation, and using numerical integration (Gaussian quadrature) instead. (Motivation: the binomial approxiation seems to be poor.) This is done by d.listcluster2, with argument costheta=NULL.
# 3. Retain the full 3D model. For this we must specify the angular separations of the stars from the cluster centre, costheta. The integrals are done numerically (Gaussian quadrature). Sometimes the integrals don't converge (presumably due to limitations of the R function "integrate"), so for those specific cases (given star and rc,sc values) I revert to case 2 (for which there is no justification other than that's robust). It currently doesn't work too well. This is also done by d.likecluster2
                                                                                                                                                                                                                                                                                                                                                                            
# I use separable priors on cluster distance and cluster size. For the former I use the exponential decreasing space density prior. I compute the posterior on a regular 2D grid of (rc,sc).

### The cluster size prior (and set parameter of distance prior)

# I use the following gamma distribution prior
sc <- seq(from=0, to=0.050, length.out=1e3)
scPriorShape <- 2
scPriorScale <- 0.005
plot(sc, dgamma(sc, shape=scPriorShape, scale=scPriorScale), type="l", xlab="sc [kpc]")
# Set parameter of distance prior
rlen <- 1 # kpc!

### Use this block to generate simulated data in order to test likelihood 1 & 2 (overwrites dat!)

# True cluster is 1D along the l.o.s, so strictly theta=0 for all of them.
# Thus cannot expect results with Likelihood 3 below to be good (as that assumes a spherical cluster).
set.seed(12345)
rcTrue <- 0.120
scTrue <- 0.01
Nstar <-  25
rTrue <- rnorm(n=Nstar,mean=rcTrue,sd=scTrue)
wsd <- rep.int(x=0.3,times=Nstar) 
w <- rnorm(n=Nstar,mean=1/rTrue,sd=wsd) # generates one random variable for each row of cbind(w,wsd)
dat <- data.frame(parallax=w, parallax_error=wsd)
dat[1:10,]

### Likelihood 1 (1D approximation, binomial approximation for integration)

# This can be skipped, as it does not give good results due to the poor binomial approximation for many cases

# Evaluate 2D posterior on the following grid
# 1D marginals will be normalized numerically (see PBI section 6.3 for details). 
# For this to be correct, the range of rc and sc must encompass essentially all density.
rc <- seq(from=0.110, to=0.130, length.out=100)
sc <- seq(from=0.0001, to=0.03, length.out=100)
dense <- matrix(0, nrow=length(rc), ncol=length(sc))
cat("Of", length(sc), "outer loop steps, no. remaining:\n")
for(j in 1:length(sc)) {
cat(j," ")
for(i in 1:length(rc)) {
dense[i,j] <- d.likecluster1(w=dat$parallax, wsd=dat$parallax_error, r=rc[i], s=sc[j]) *
d.distprior3(r=rc[i], rlen=rlen)*dgamma(x=sc[j], shape=scPriorShape, scale=scPriorScale)   
}
}
cat("\n")

# Plot 2D posterior as well as the two marginal PDFs (achieved simply by summing grid, as it's regular)
range(dense)
dense <- dense/max(dense)
par(mfrow=c(2,2), mar=c(5,5,0.5,1), oma=c(0.1,0.1,0.5,0.1), mgp=c(2.2,0.8,0), cex=1.0) 
image.plot(z=dense, x=rc, y=sc, nlevel=1024, col=mycols)
scDense <- apply(dense, 2, sum)
mom <- pdfmom(scDense, sc)
cat("Posterior mean, SD for sc =", mom$mean, mom$sd, "kpc\n")
plot(sc, scDense/mom$Z, type="l", ylab="P(sc | {w},{wsd},priors)")
lines(sc, dgamma(sc, shape=2, scale=0.003), col="red") # overplot prior on sc
rcDense <- apply(dense, 1, sum)
mom <- pdfmom(rcDense, rc)
cat("Posterior mean, SD for rc =", mom$mean, mom$sd, "kpc\n")
plot(rc, rcDense/mom$Z, type="l", ylab="P(rc | {w},{wsd},priors)")

### Likelihood 2 (1D approximation, numerical integrations)

# Evaluate 2D posterior on the following grid.
# This takes about 2.5 minutes on my laptop using gdr1set03.csv (164 stars), 
# so first test it using 10-20 simulated stars (see code for this above)
rc <- seq(from=0.120, to=0.140, length.out=100)
sc <- seq(from=0.0001, to=0.030, length.out=100)
dense <- matrix(0, nrow=length(rc), ncol=length(sc))
cat("Of", length(sc), "outer loop steps, no. remaining:\n")
for(j in 1:length(sc)) {
cat(j," ")
for(i in 1:length(rc)) {
dense[i,j] <- d.likecluster2(w=dat$parallax, wsd=dat$parallax_error, rc=rc[i], sc=sc[j]) *
d.distprior3(r=rc[i], rlen=rlen)*dgamma(x=sc[j], shape=scPriorShape, scale=scPriorScale)   
}
}
cat("\n")

# Plot 2D posterior as well as the two marginal PDFs (achieved simply by summing grid, as it's regular)
range(dense)
dense <- dense/max(dense)
par(mfrow=c(2,2), mar=c(5,5,0.5,1), oma=c(0.1,0.1,0.5,0.1), mgp=c(2.2,0.8,0), cex=1.0) 
image.plot(z=dense, x=rc, y=sc, nlevel=1024, col=mycols)
scDense <- apply(dense, 2, sum)
mom <- pdfmom(scDense, sc)
cat("Posterior mean, SD for sc =", mom$mean, mom$sd, "kpc\n")
plot(sc, scDense/mom$Z, type="l", ylab="P(sc | {w},{wsd},priors)")
lines(sc, dgamma(sc, shape=2, scale=0.003), col="red") # overplot prior on sc
rcDense <- apply(dense, 1, sum)
mom <- pdfmom(rcDense, rc)
cat("Posterior mean, SD for rc =", mom$mean, mom$sd, "kpc\n")
plot(rc, rcDense/mom$Z, type="l", ylab="P(rc | {w},{wsd},priors)")
# compare result with parallaxes (!) derived in GDR1 overview and cluster papers

# ### Use this block to generate simulated data in order to test likelihood 3 (overwrites dat!)
# 
# set.seed(12345)
# rcTrue <- 0.120
# scTrue <- 0.01
# Nstar <-  10
# rTrue <- rnorm(n=Nstar,mean=rcTrue,sd=scTrue)
# wsd <- rep.int(x=0.3,times=Nstar) 
# w <- rnorm(n=Nstar,mean=1/rTrue,sd=wsd) # generates one random variable for each row of cbind(w,wsd)
# theta <- rnorm(n=Nstar, mean=0, sd=0.11) # in degrees
# dat <- data.frame(parallax=w, parallax_error=wsd)
# cbind(dat, theta)[1:10,]

### Likelihood 3 (3D, numerical integrations)

# This should be skipped, as the numerical integration I am using doesn't converge. See below for notes in case you want to try a better integrator...

# Compute angular separations of stars from nominal cluster centre, clusterCen, which is assumed known.
# Comment out this block out if using test data (as theta is already assigned above)
clustCen <- data.frame(ra=56.75, dec=24.12) # RA, Dec of cluster center in degrees
theta <- sqrt(((dat$ra-clustCen$ra)*cos(conv*clustCen$dec))^2 +
(dat$dec-clustCen$dec)^2) # on-sky separations from clustCen in degrees

# Evaluate 2D posterior on the following grid.
# This takes about 18 minutes on my laptop using gdr1set03.csv (164 stars) and gives meaningless results.
# So first test it using 10-20 simulated stars (see code for this above).
rc <- seq(from=0.120, to=0.140, length.out=200)
sc <- seq(from=0.0005, to=0.030, length.out=100)
dense <- matrix(0, nrow=length(rc), ncol=length(sc))
cat("Of", length(sc), "outer loop steps, no. remaining:\n")
for(j in 1:length(sc)) {
cat(j," ")
for(i in 1:length(rc)) {
dense[i,j] <- d.likecluster2(w=dat$parallax, wsd=dat$parallax_error, costheta=cos(conv*theta), 
rc=rc[i], sc=sc[j]) *
d.distprior3(r=rc[i], rlen=rlen)*dgamma(x=sc[j], shape=scPriorShape, scale=scPriorScale)   
}
}
cat("\n")

# Plot 2D posterior as well as the two marginal PDFs (achieved simply by summing grid, as it's regular)
range(dense)
dense <- dense/max(dense)
par(mfrow=c(2,2), mar=c(5,5,0.5,1), oma=c(0.1,0.1,0.5,0.1), mgp=c(2.2,0.8,0), cex=1.0) 
image.plot(z=dense, x=rc, y=sc, nlevel=1024, col=mycols)
scDense <- apply(dense, 2, sum)
mom <- pdfmom(scDense, sc)
cat("Posterior mean, SD for sc =", mom$mean, mom$sd, "kpc\n")
plot(sc, scDense/mom$Z, type="l", ylab="P(sc | {w},{wsd},priors)")
lines(sc, dgamma(sc, shape=2, scale=0.003), col="red") # overplot prior on sc
rcDense <- apply(dense, 1, sum)
mom <- pdfmom(rcDense, rc)
cat("Posterior mean, SD for rc =", mom$mean, mom$sd, "kpc\n")
plot(rc, rcDense/mom$Z, type="l", ylab="P(rc | {w},{wsd},priors)")

#### Notes

# With the above test data we get strange results: the 2D PDF has strong discontinuities as a function of rc, as though it were being cut. They do not appear when doing this with costheta=0 (1D cluster approximation, the previous block). They remain even if I reduce the approximation threshold (for avoiding divide by zero) in the line 

"sel <- which(x/rc<1e-6)"

# in d.likecluster2() to zero (divide by zero avoided, presumably as no theta are very small). Inserting the proper integrastion limits (0,Inf), rather than the truncated ones (which are meant to be faster) makes no difference.
# Some integrastions in all cases do fail, and if I replace those failures with zero (likelihood), rather than reverting to the 1D numerical case (i.e. costheta=NULL), then every element of dense (posterior) is zero. Thus we need to prevent integrate.func() (which is just a wrapper to integrate{stats}) failing. Experience tells me this can't be improved so needs to be replaced. Could do MCMC, but that's very slow.

