# =============================================================================
# ФУНКЦИИ ОПТИМИЗАЦИИ ПОРТФЕЛЯ
# =============================================================================

# 1. Equal Weight
optimize_equal <- function(returns, Sigma = NULL) {
  n <- ncol(returns)
  w <- rep(1/n, n)
  names(w) <- colnames(returns)
  
  mu <- colMeans(returns)
  if (is.null(Sigma)) Sigma <- cov(returns)
  
  port_return <- as.numeric(crossprod(mu, w))
  port_risk <- as.numeric(sqrt(t(w) %*% Sigma %*% w))
  
  return(list(weights = w, return = port_return, risk = port_risk,
              sharpe = ifelse(port_risk > 0, port_return / port_risk, NA)))
}

# 2. Inverse Volatility
optimize_invvol <- function(returns, Sigma = NULL) {
  n <- ncol(returns)
  
  if (is.null(Sigma)) {
    Sigma <- cov(returns)
  }
  
  vols <- sqrt(diag(Sigma))
  w <- 1 / vols
  w <- w / sum(w)
  names(w) <- colnames(returns)
  
  mu <- colMeans(returns)
  port_return <- as.numeric(crossprod(mu, w))
  port_risk <- as.numeric(sqrt(t(w) %*% Sigma %*% w))
  
  return(list(weights = w, return = port_return, risk = port_risk,
              sharpe = ifelse(port_risk > 0, port_return / port_risk, NA)))
}

# 3. Risk Parity
optimize_riskparity <- function(returns, Sigma = NULL) {
  n <- ncol(returns)
  
  if (is.null(Sigma)) {
    Sigma <- cov(returns)
  }
  
  # Функция для расчета вклада в риск
  risk_contribution <- function(w, Sigma) {
    port_var <- as.numeric(t(w) %*% Sigma %*% w)
    if (port_var <= 0) return(rep(0, length(w)))
    (w * (Sigma %*% w)) / sqrt(port_var)
  }
  
  # Целевая функция - минимизация дисперсии вкладов в риск
  objective_rp <- function(w, Sigma) {
    w_sum <- sum(w)
    if (abs(w_sum) < 1e-6) return(1e6)
    w <- w / w_sum
    rc <- risk_contribution(w, Sigma)
    var(rc)
  }
  
  # Поиск решения с ограничениями long only
  best_obj <- Inf
  best_w <- NULL
  
  for (attempt in 1:30) {
    w0 <- runif(n, 0, 1)
    w0 <- w0 / sum(w0)
    
    opt <- tryCatch({
      optim(w0, objective_rp, Sigma = Sigma, 
            method = "L-BFGS-B",
            lower = rep(0, n),
            upper = rep(1, n),
            control = list(maxit = 2000))
    }, error = function(e) NULL)
    
    if (!is.null(opt) && opt$value < best_obj) {
      best_obj <- opt$value
      best_w <- opt$par
    }
  }
  
  if (is.null(best_w)) {
    w <- rep(1/n, n)
  } else {
    w <- best_w / sum(best_w)
  }
  
  names(w) <- colnames(returns)
  
  mu <- colMeans(returns)
  port_return <- as.numeric(crossprod(mu, w))
  port_risk <- as.numeric(sqrt(t(w) %*% Sigma %*% w))
  
  return(list(weights = w, return = port_return, risk = port_risk,
              sharpe = ifelse(port_risk > 0, port_return / port_risk, NA)))
}

# 4. Minimum Variance
optimize_minvar <- function(returns, Sigma = NULL) {
  n <- ncol(returns)
  
  if (is.null(Sigma)) {
    Sigma <- cov(returns)
  }
  
  if (requireNamespace("quadprog", quietly = TRUE)) {
    Dmat <- 2 * Sigma
    dvec <- rep(0, n)
    Amat <- cbind(rep(1, n), diag(n))
    bvec <- c(1, rep(0, n))
    sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
    w <- sol$solution
  } else {
    ones <- rep(1, n)
    Sigma_inv <- solve(Sigma + diag(1e-8, n))
    w <- Sigma_inv %*% ones
    w <- as.vector(w / sum(w))
    w <- pmax(w, 0)
    w <- w / sum(w)
  }
  
  names(w) <- colnames(returns)
  
  mu <- colMeans(returns)
  port_return <- as.numeric(crossprod(mu, w))
  port_risk <- as.numeric(sqrt(t(w) %*% Sigma %*% w))
  
  return(list(weights = w, return = port_return, risk = port_risk,
              sharpe = ifelse(port_risk > 0, port_return / port_risk, NA)))
}

# 5. Tangency (Max Sharpe)
optimize_tangency <- function(returns, Sigma = NULL) {
  if (is.null(returns) || (is.matrix(returns) && nrow(returns) == 0)) {
    n <- if (!is.null(Sigma)) ncol(Sigma) else 6
    w <- rep(1/n, n)
    names(w) <- if (!is.null(colnames(Sigma))) colnames(Sigma) else paste0("Asset", 1:n)
    
    if (is.null(Sigma)) Sigma <- diag(n)
    mu <- rep(0, n)
    
    port_return <- 0
    port_risk <- as.numeric(sqrt(t(w) %*% Sigma %*% w))
    
    return(list(weights = w, return = port_return, risk = port_risk, sharpe = 0))
  }
  
  if (!is.matrix(returns)) {
    returns <- as.matrix(returns)
  }
  
  n <- ncol(returns)
  
  if (is.null(Sigma)) {
    Sigma <- cov(returns) + diag(1e-8, n)
  }
  
  mu <- colMeans(returns)
  
  if (requireNamespace("quadprog", quietly = TRUE)) {
    tryCatch({
      Dmat <- 2 * Sigma
      dvec <- rep(0, n)
      Amat <- cbind(rep(1, n), diag(n), -diag(n))
      bvec <- c(1, rep(0, n), rep(-1, n))
      meq <- 1
      
      sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
      w <- sol$solution
      w <- pmax(w, 0)
      w <- pmin(w, 1)
      w <- w / sum(w)
    }, error = function(e) {
      w <- rep(1/n, n)
    })
  } else {
    tryCatch({
      Sigma_inv <- solve(Sigma + diag(1e-8, n))
      w <- Sigma_inv %*% mu
      w <- as.vector(w)
      w <- pmax(w, 0)
      w <- pmin(w, 1)
      if (sum(w) > 0) {
        w <- w / sum(w)
      } else {
        w <- rep(1/n, n)
      }
    }, error = function(e) {
      w <- rep(1/n, n)
    })
  }
  
  names(w) <- colnames(returns)
  
  port_return <- as.numeric(crossprod(mu, w))
  port_risk <- as.numeric(sqrt(t(w) %*% Sigma %*% w))
  
  return(list(weights = w, return = port_return, risk = port_risk,
              sharpe = ifelse(port_risk > 0, port_return / port_risk, 0)))
}

# 6. Risk Aversion (с ограничениями)
optimize_riskaversion <- function(returns, Sigma = NULL, gamma = 5) {
  n <- ncol(returns)
  
  if (is.null(Sigma)) {
    Sigma <- cov(returns)
  }
  
  mu <- colMeans(returns)
  
  if (requireNamespace("quadprog", quietly = TRUE)) {
    Dmat <- gamma * Sigma
    dvec <- mu
    Amat <- cbind(rep(1, n), diag(n), -diag(n))
    bvec <- c(1, rep(0, n), rep(-1, n))
    meq <- 1
    
    sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
    w <- sol$solution
  } else {
    Sigma_inv <- solve(Sigma + diag(1e-8, n))
    w_unconstrained <- (1/gamma) * Sigma_inv %*% mu
    w <- as.vector(w_unconstrained)
  }
  
  w <- pmax(w, 0)
  w <- pmin(w, 1)
  
  if (sum(w) > 0) {
    w <- w / sum(w)
  } else {
    w <- rep(1/n, n)
  }
  
  names(w) <- colnames(returns)
  
  port_return <- as.numeric(crossprod(mu, w))
  port_risk <- as.numeric(sqrt(t(w) %*% Sigma %*% w))
  
  return(list(weights = w, return = port_return, risk = port_risk,
              sharpe = ifelse(port_risk > 0, port_return / port_risk, NA),
              utility = port_return - (gamma/2) * port_risk^2))
}

# 7. Tobin (с ограничениями)
optimize_tobin <- function(returns, Sigma = NULL, gamma = 5) {
  if (is.null(returns) || (is.matrix(returns) && nrow(returns) == 0)) {
    n <- if (!is.null(Sigma)) ncol(Sigma) else 6
    w <- rep(0, n)
    names(w) <- if (!is.null(colnames(Sigma))) colnames(Sigma) else paste0("Asset", 1:n)
    
    return(list(weights = w, return = 0, risk = 0, sharpe = 0, risky_share = 0))
  }
  
  if (!is.matrix(returns)) {
    returns <- as.matrix(returns)
  }
  
  n <- ncol(returns)
  
  if (is.null(Sigma)) {
    Sigma <- cov(returns) + diag(1e-8, n)
  }
  
  mu <- colMeans(returns)
  
  tangency_result <- optimize_tangency(returns, Sigma)
  w_tangency <- tangency_result$weights
  
  mu_t <- sum(w_tangency * mu)
  sigma_t <- as.numeric(sqrt(t(w_tangency) %*% Sigma %*% w_tangency))
  
  risky_share <- 0
  if (sigma_t > 1e-8 && mu_t > 0) {
    risky_share <- mu_t / (gamma * sigma_t^2)
    risky_share <- min(max(risky_share, 0), 1)
  }
  
  w <- w_tangency * risky_share
  w <- pmax(w, 0)
  w <- pmin(w, 1)
  
  names(w) <- colnames(returns)
  
  port_return <- as.numeric(crossprod(mu, w))
  port_risk <- as.numeric(sqrt(t(w) %*% Sigma %*% w))
  
  return(list(weights = w, return = port_return, risk = port_risk,
              sharpe = ifelse(port_risk > 0, port_return / port_risk, 0),
              risky_share = risky_share))
}

# Словарь всех стратегий
optimization_strategies <- list(
  "Equal" = optimize_equal,
  "InvVol" = optimize_invvol,
  "MinVariance" = optimize_minvar,
  "RiskParity" = optimize_riskparity,
  "Tangency" = optimize_tangency,
  "RiskAversion" = optimize_riskaversion,
  "Tobin" = optimize_tobin
)