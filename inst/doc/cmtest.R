## ----echo = TRUE--------------------------------------------------------------
lnL_bc <- function(coef, X, y, sum = FALSE, gradient = FALSE, hessian = FALSE){
    y <- y
    N <- length(y)
    K <- length(coef) - 3
    beta <- coef[1:(K + 1)]
    sigma <- coef[K + 2]
    lambda <- coef[K + 3]
    bX <- as.numeric(X %*% beta)
    if (lambda == 0){
        bc <- log(y)
        bc_l  <- 1 / 2 * log(y) ^ 2
        bc_ll <- 1 / 3 * log(y) ^ 3
    }
    else{
        bc <- (y ^ lambda - 1) / lambda
        bc_l <- bc * log(y) - (bc - log(y)) / lambda
        bc_ll <- bc_l * log(y) - (bc_l * lambda - (bc - log(y))) / lambda ^ 2
    }
    e <- bc - bX
    lnL <- - 1 / 2 * log(2 * pi) - log(sigma) - (1 - lambda) * log(y) -
        1 / (2 * sigma ^ 2) * e ^ 2

    if (gradient){
        g_beta <- 1 / sigma ^ 2 * e * X
        g_sigma <- - 1 / sigma + 1 / (sigma ^ 3) * e ^ 2
        g_lambda <- log(y) - 1 / (sigma ^ 2) * e * bc_l
        G <- cbind(g_beta, g_sigma,  g_lambda)
    }
    if (hessian){
        h_bb <- - 1 / sigma ^ 2 * crossprod(X)
        h_bs <- apply(- 2 / sigma ^ 3 * X * e, 2, sum)
        h_bl <- apply( 1 / sigma ^ 2 * bc_l * X, 2, sum)
        h_ss <- sum(1 / sigma ^ 2 - 3 / sigma ^ 4 * e ^ 2)
        h_sl <- sum(2 / sigma ^ 3 * e * bc_l)
        h_ll <- sum(- 1 / sigma ^ 2 * (e * bc_ll + bc_l ^ 2))
        H <- rbind(cbind(h_bb, sigma = h_bs, lambda = h_bl),
                   sigma = c(h_bs, h_ss, h_sl),
                   lambda = c(h_bl, h_sl, h_ll))
    }
    if (sum){
        lnL <- sum(lnL)
        if (gradient) G <- apply(G, 2, sum)
    }
    if (gradient) attr(lnL, "gradient") <- G
    if (hessian) attr(lnL, "hessian") <- H
    lnL
}

## -----------------------------------------------------------------------------
lin_lm <- lm(dist ~ speed + I(speed ^ 2), cars)
log_lm <- update(lin_lm, log(.) ~ .)

## -----------------------------------------------------------------------------
coefs_lin <- coef(lin_lm)
coefs_lin[1] <- coefs_lin[1] - 1

## -----------------------------------------------------------------------------
X <- model.matrix(lin_lm)
y <- model.response(model.frame(lin_lm))

## ----collapse = TRUE----------------------------------------------------------
sigma2 <- function(x) sigma(x) * sqrt(df.residual(x) / nobs(x))
coefs_lin <- c(coefs_lin, sigma = sigma2(lin_lm), lambda = 1)
coefs_log <- c(coef(log_lm), sigma = sigma2(log_lm), lambda = 0)
lnL_lin <- lnL_bc(coefs_lin, X, y, sum = TRUE, gradient = TRUE, hessian = FALSE)
lnL_log <- lnL_bc(coefs_log, X, y, sum = TRUE, gradient = TRUE, hessian = FALSE)
lnL_lin
lnL_log

## ----warning = FALSE, message = FALSE-----------------------------------------
library("maxLik")
bc_lin <- maxLik(lnL_bc, start = coefs_lin, X = X, y = y,
                 sum = TRUE, gradient = TRUE, hessian = TRUE)
bc_log <- maxLik(lnL_bc, start = coefs_log, X = X, y = y,
                 sum = TRUE, gradient = TRUE, hessian = TRUE)

## -----------------------------------------------------------------------------
cbind(coef(bc_lin), coef(bc_log))

## -----------------------------------------------------------------------------
summary(bc_lin)

## -----------------------------------------------------------------------------
lnL_lin <- lnL_bc(coefs_lin, X, y, sum = FALSE, gradient = TRUE, hessian = TRUE)
lnL_log <- lnL_bc(coefs_log, X, y, sum = FALSE, gradient = TRUE, hessian = TRUE)
G_lin <- attr(lnL_lin, "gradient")
H_lin <- attr(lnL_lin, "hessian")
G_log <- attr(lnL_log, "gradient")
H_log <- attr(lnL_log, "hessian")
g_lin <- apply(G_lin, 2, sum)
g_log <- apply(G_log, 2, sum)

## ----collapse = TRUE----------------------------------------------------------
as.numeric(g_lin %*% solve(- H_lin) %*% g_lin)
as.numeric(g_lin %*% solve(crossprod(G_lin)) %*% g_lin)
as.numeric(g_log %*% solve(- H_log) %*% g_log)
as.numeric(g_log %*% solve(crossprod(G_log)) %*% g_log)

## ----collapse = TRUE----------------------------------------------------------
iota <- rep(1, length(y))
lm_ones_lin <- lm(iota ~ G_lin - 1)
lm_ones_log <- lm(iota ~ G_log - 1)
summary(lm_ones_lin)$r.squared * nobs(lm_ones_lin)
summary(lm_ones_log)$r.squared * nobs(lm_ones_log)

## -----------------------------------------------------------------------------
cm_norm <- function(x, type = c("analytical", "opg", "reg")){
    type <- match.arg(type)
    X <- model.matrix(x)
    y <- model.response(model.frame(x))
    N <- length(y)
    K <- length(coef(x))
    coefs <- c(coef(x), sigma = sigma2(x))
    beta <- coefs[1:K]
    sigma <- coefs[K + 1]
    epsilon <- as.numeric(y - X %*% beta)
    G <- cbind(1 / sigma ^ 2 * X * epsilon,
               sigma = - 1  / sigma + 1 / sigma ^ 3 * epsilon ^ 2)
    M <- cbind(asym = epsilon ^ 3, kurt = epsilon ^ 4 - 3 * sigma ^ 4)
    m <- apply(M, 2, sum)
    if (type == "analytical"){
        I <- - rbind(cbind(- 1 / sigma ^ 2 * crossprod(X),
                           sigma = - 2 / sigma ^ 3 * apply(X * epsilon, 2, sum)),
                     sigma = c(- 2 / sigma ^ 3 * apply(X * epsilon, 2, sum),
                               1 / sigma ^ 2 - 3 / sigma ^ 4 * sum(epsilon ^ 2)))
        W <- - rbind(cbind(asym = - 3 * apply(epsilon ^ 2 * X, 2, sum),
                           kurt = - 4 * apply(epsilon ^ 3 * X, 2, sum)),
                     sigma = c(0, - 12 * N / sigma ^ 3))
    }
    if (type == "opg"){
        I <- crossprod(G)
        W <- crossprod(G, M)
    }
    if (type != "reg"){
        Q <- crossprod(M - G %*% solve(I) %*% W)
        stat <- as.numeric(m %*% solve(Q) %*% m)
    }
    else stat <- summary(lm(rep(1, N) ~ G + M - 1))$r.squared * N
    stat
}

## ----collapse = TRUE----------------------------------------------------------
cm_norm(lin_lm, "analytical")
cm_norm(lin_lm, "opg")
cm_norm(lin_lm, "reg")

## ----collapse = TRUE----------------------------------------------------------
cm_norm(log_lm, "analytical")
cm_norm(log_lm, "opg")
cm_norm(log_lm, "reg")

## ----eval = FALSE, include = FALSE--------------------------------------------
#  JB <- function(x){
#      e <- resid(x)
#      mu <- mean(e)
#      sigma <- sqrt( mean( (e - mu) ^ 2))
#      m3 <- mean( (e - mu) ^ 3)
#      m4 <- mean( (e - mu) ^ 4)
#      Asym <- m3 / sigma ^ 3
#      Kurt <- m4 / sigma ^ 4
#      df.residual(x) / 6 * (Asym ^ 2 + (Kurt - 3) ^ 2 / 4)
#  }

