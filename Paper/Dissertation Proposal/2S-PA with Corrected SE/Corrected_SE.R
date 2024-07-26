library(lavaan)
library(R2spa)
library(numDeriv)
library(boot)
source("/Users/jimmy_z/R Projects/R2spa/R/get_fscore.R")
source("/Users/jimmy_z/R Projects/R2spa/R/tspa_corrected_se.R")
source("/Users/jimmy_z/R Projects/R2spa/R/helper.R")

num_obs <- 5000
# Simulate Latent Data
# path: 0.3, 0.3, 0.3
# cor_xm = 0
# x, m, ey
cov_xm_ey <- matrix(c(1, 0, 0,
                      0, 1, 0,
                      0, 0, 0.73), nrow = 3)
eta <- MASS::mvrnorm(num_obs, mu = rep(0, 3), Sigma = cov_xm_ey,
                     empirical = FALSE)
# xm
eta <- cbind(eta, eta[, 1] * eta[, 2])
# y
etay <- eta[, -3] %*% rep(0.3, 3) + eta[, 3]
# x, m, xm, y
eta <- as.data.frame(cbind(eta[,1], eta[,2], eta[,4], etay))
names(eta) <- c("x_j", "m_j", "xm_j", "y")

# Latent continuous vars
lambdax <- c(.60, .70, .80)
lambdam <- c(.30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85)
xstar <- as.data.frame(eta[, 1] %*% t(lambdax) + rnorm(num_obs * length(lambdax)))
names(xstar) <- paste0("x", 1:3, "j")
mstar <- as.data.frame(eta[, 2] %*% t(lambdam) + rnorm(num_obs * length(lambdam)))
names(mstar) <- paste0("m", 1:12, "j")
df <- cbind(xstar, mstar, etay)
names(df)[ncol(df)] <- "y"

# SEM Model
model = "x_j =~ x1j + x2j + x3j
         m_j =~ m1j + m2j + m3j + m4j + m5j + m6j +
                m7j + m8j + m9j + m10j + m11j + m12j
         y ~ x_j + m_j"

cfa_x_j <- cfa("x_j =~ x1j + x2j + x3j", data = df)
fs_x_j <- get_fs_lavaan(cfa_x_j, vfsLT = TRUE)
vldev_x_j <- attr(fs_x_j, which = "vfsLT") #?
cfa_m_j <- cfa("m_j =~ m1j + m2j + m3j + m4j + m5j + m6j +
                m7j + m8j + m9j + m10j + m11j + m12j", data = df)
fs_m_j <- get_fs_lavaan(cfa_m_j, vfsLT = TRUE)
vldev_m_j <- attr(fs_m_j, which = "vfsLT") #?

# Factor Scores
fs_dat <- data.frame(
  fs_x_j = fs_x_j$fs_x_j,
  fs_m_j = fs_m_j$fs_m_j
)
# Combine sampling variance of loading and error variance
# Note: loadings first, then error variance
vldev <- block_diag(vldev_x_j, vldev_m_j)[c(1, 3, 2, 4), c(1, 3, 2, 4)]




vldev1 <- attr(fs_x_j, which = "vfsLT")
cfa_dem60 <- cfa("dem60 =~ y1 + y2 + y3 + y4",
                 data = PoliticalDemocracy)
# Regression factor scores
fs2 <- get_fs_lavaan(cfa_dem60, vfsLT = TRUE)
# Delta method variance of (loading, intercept)
vldev2 <- attr(fs2, which = "vfsLT")
fs_dat <- data.frame(
  fs_ind60 = fs1$fs_ind60,
  fs_dem60 = fs2$fs_dem60
)
# Combine sampling variance of loading and error variance
# Note: loadings first, then error variance
vldev <- block_diag(vldev1, vldev2)[c(1, 3, 2, 4), c(1, 3, 2, 4)]

# ld <- block_diag(attr(fs_x_j, which = "fsL"),
#                  attr(fs_m_j, which = "fsL"))
# ev <- block_diag(attr(fs_x_j, which = "fsT"),
#                  attr(fs_m_j, which = "fsT"))
# tspa_fit <- tspa(model = "x_j ~ m_j",
#                  data = fs_dat,
#                  fsL = ld,
#                  fsT = ev)




