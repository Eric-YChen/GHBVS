# Usage (from simulation/): Rscript sim_polarised.R
suppressMessages({ library(MASS); library(mvtnorm); library(survival); library(glmnet) })
BVS_CORE <- normalizePath("../core", winslash = "/")
suppressMessages(source(file.path(BVS_CORE, "ADS_sampler.R")))

set.seed(123)
n <- 30; k <- 3
Sigma <- matrix(0.2, k, k); diag(Sigma) <- 0.7
X_uncensored <- mvrnorm(n = 6, mu = c(0, 0, 0), Sigma = Sigma)
X_censored <- mvrnorm(n = 24, mu = c(-0.05, 5.6, -0.05), Sigma = Sigma)
X <- rbind(X_uncensored, X_censored); colnames(X) <- c("X1", "X2", "X3")
y_uncensored <- log(rexp(6, rate = exp(-1 + X_uncensored %*% c(-1.8, 3, -0.7))))
y <- exp(c(y_uncensored, rep(-48, 24)))
delta <- c(rep(1, 6), rep(0, 24))
y[1:6] <- exp(c(-1.6, 1.2, 0.0, 1.8, 1.3, 3.6))
X[1:6, ] <- matrix(c(-0.4,  0.4,  0.6,
                      0.6,  0.2,  1.9,
                     -0.1, -0.1, -0.6,
                      0.5, -0.0, -0.2,
                     -1.1, -1.0, -0.9,
                      0.4, -0.6,  0.6), nrow = 6, byrow = TRUE)
dat <- list(y = y, x = X, xt = X, status = delta)

mcmc <- list(n.iter = 20000, n.burnin = 10000, n.thin = 2, n.chain = 1, printevery = 1e9)
res <- do.call(rbind, lapply(c("LCM", "Product"), function(pr) {
  set.seed(1)
  mp <- list(inits = list(nu = rnorm(1), theta0 = rnorm(1), theta = rnorm(k, 0, 0.5), eta = rnorm(k, 0, 0.5),
                          g.theta = 1, g.eta = 1, g = 1, a_g = NULL, gam = c(4L, 4L, 0L), nu.alpha = 0.01, nu.beta = 0.01),
             gam.param = list(layer2.mulnom.par = c(0, 0, 0, 1)), data = dat)
  ch <- BVS.MCMC(mp, mcmc, verbose = FALSE, gprior.mat = pr, gamprior = "custom")
  pp <- ch$post.prob
  incl <- do.call(rbind, lapply(strsplit(pp$Vars, "", fixed = TRUE), function(s) s != "0"))
  cat(sprintf("\n%s prior: %d models evaluated; models and posterior probabilities:\n", pr, nrow(pp)))
  print(data.frame(model = pp$Vars, prob = round(pp$prob, 3), visits = pp$visits), row.names = FALSE)
  if (length(ch$opt.log)) { cat("optimisation problems:\n"); print(ch$opt.log) }
  data.frame(method = paste0("AFT-", pr), X1 = sum(pp$prob * incl[, 1]), X2 = sum(pp$prob * incl[, 2]),
             X3 = sum(pp$prob * incl[, 3]), n.eval = nrow(pp))
}))
cat("\nPIPs (Table tab:motive-ex):\n"); print(transform(res, X1 = round(X1, 3), X2 = round(X2, 3), X3 = round(X3, 3)), row.names = FALSE)
dir.create("results", showWarnings = FALSE)
write.csv(res, file.path("results", "sim_polarised.csv"), row.names = FALSE)
