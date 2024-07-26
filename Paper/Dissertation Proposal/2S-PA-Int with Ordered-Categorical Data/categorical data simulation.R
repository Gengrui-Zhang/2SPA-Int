library(semTools)
library(lavaan)
library(ggplot2)

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
# Logistic Distribution
lambdax <- c(.60, .70, .80)
lambdam <- c(.30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85)
xstar <- as.data.frame(eta[, 1] %*% t(lambdax) + rnorm(num_obs * length(lambdax)))
names(xstar) <- paste0("x", 1:3, "j")
mstar <- as.data.frame(eta[, 2] %*% t(lambdam) + rnorm(num_obs * length(lambdam)))
names(mstar) <- paste0("m", 1:12, "j")

# # Test factor loadings
# model = "x_j =~ x1j + x2j + x3j"
# x_j_cfa <- sem(model, xstar, std.lv = T)
# summary(x_j_cfa)
# model = "m_j =~ m1j + m2j + m3j + m4j + m5j + m6j +
#                 m7j + m8j + m9j + m10j + m11j + m12j"
# m_j_cfa <- sem(model, mstar, std.lv = T)
# summary(m_j_cfa)
# # Reliability
# reliability(x_j_cfa)
# reliability(m_j_cfa)

# Categorical vars
# Symmetric
thresx_symm <- matrix(c(0, 0, 0), nrow = 1)  # binary
thresm_symm <- t(matrix(c(-1.5, -0.5, 0.5, 1.5), nrow = 12, ncol = 4, byrow = TRUE))  # ordinal with 5 categories
x_obs_symm <- vapply(
  seq_along(lambdax),
  FUN = \(i) {
    findInterval(xstar[, i], thresx_symm[, i])
  },
  FUN.VALUE = numeric(num_obs))
x_obs_symm <- as.data.frame(x_obs_symm)
names(x_obs_symm) <- paste0("x", 1:3, "j")
m_obs_symm <- vapply(
  seq_along(lambdam),
  FUN = \(i) {
    findInterval(mstar[, i], thresm_symm[, i])
  },
  FUN.VALUE = numeric(num_obs))
m_obs_symm <- as.data.frame(m_obs_symm)
names(m_obs_symm) <- paste0("m", 1:12, "j")

# Skewed
thresx_skew <- matrix(c(-2.2, -2.2, -2.2), nrow = 1)  # binary
thresm_skew <- t(matrix(c(.05, .75, 1.55, 2.55), nrow = 12, ncol = 4, byrow = TRUE))  # ordinal with 5 categories
x_obs_skew <- vapply(
  seq_along(lambdax),
  FUN = \(i) {
    findInterval(xstar[, i], thresx_skew[, i])
  },
  FUN.VALUE = numeric(num_obs))
x_obs_skew <- as.data.frame(x_obs_skew)
names(x_obs_skew) <- paste0("x", 1:3, "j")
m_obs_skew <- vapply(
  seq_along(lambdam),
  FUN = \(i) {
    findInterval(mstar[, i], thresm_skew[, i])
  },
  FUN.VALUE = numeric(num_obs))
m_obs_skew <- as.data.frame(m_obs_skew)
names(m_obs_skew) <- paste0("m", 1:12, "j")

