# design: "orig" (Kappa, Lambda) or "trans" (their sum and ratio)
prep.flchain <- function(design = c("orig", "trans"), file = "data/flchain.csv"){
  design <- match.arg(design)
  fl <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  d <- fl[, c("age", "sex", "kappa", "lambda", "creatinine", "mgus", "futime", "death")]
  d <- d[complete.cases(d), ]
  d <- d[d$futime > 0, ]
  sex <- as.numeric(d$sex == "M")
  if (design == "orig"){
    X <- cbind(age = scale(d$age)[, 1], kappa = scale(d$kappa)[, 1],
               lambda = scale(d$lambda)[, 1], sex = sex,
               creatinine = scale(d$creatinine)[, 1], mgus = d$mgus)
  } else {
    X <- cbind(age = scale(d$age)[, 1], "kap+lam" = scale(d$kappa + d$lambda)[, 1],
               "kap/lam" = scale(d$kappa / d$lambda)[, 1], sex = sex,
               creatinine = scale(d$creatinine)[, 1], mgus = d$mgus)
  }
  list(x = X, xt = X, status = d$death, y = d$futime / 365.24)
}

prep.nki <- function(file = "data/nki70.csv"){
  nk <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  X <- nk[, -(1:2)]
  X$Diam  <- as.numeric(X$Diam == ">2cm")
  X$N     <- as.numeric(X$N == ">=4")
  X$ER    <- as.numeric(X$ER == "Positive")
  X$Grade <- match(X$Grade, c("Well diff", "Intermediate", "Poorly diff")) - 1
  stopifnot(!anyNA(X$Grade))
  num <- setdiff(colnames(X), c("Diam", "N", "ER", "Grade"))
  X[num] <- lapply(X[num], function(v) scale(v)[, 1])
  X <- as.matrix(X)
  list(y = nk$time, status = nk$event, x = X, xt = X)
}

# noise: add a standard normal noise covariate generated with seed noise.seed
prep.ipilimumab <- function(noise = FALSE, noise.seed = 1, file = "data/CA184043_2A.csv"){
  df0 <- read.csv(file, header = TRUE)
  df <- data.frame(time = df0$Time, status = df0$Event,
                   trt = as.numeric(df0$Arm == "ipilimumab"))
  if (noise){
    df <- df[order(df$time), ]
    set.seed(noise.seed)
    df$spur <- as.numeric(scale(rnorm(nrow(df))))
  }
  X <- as.matrix(df[, setdiff(colnames(df), c("time", "status")), drop = FALSE])
  list(df = df, dat = list(x = X, xt = X, status = as.logical(df$status), y = df$time))
}
