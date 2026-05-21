# =============================================================================
# GARCH МОДЕЛИ ДЛЯ ПРОГНОЗИРОВАНИЯ ВОЛАТИЛЬНОСТИ
# =============================================================================

# GARCH(1,1)
forecast_garch <- function(returns_series) {
  tryCatch({
    spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 0)),
      distribution.model = "std"
    )
    fit <- ugarchfit(spec, data = returns_series, solver = "hybrid")
    
    if (convergence(fit) == 0) {
      forecast <- ugarchforecast(fit, n.ahead = 1)
      return(list(vol = as.numeric(sigma(forecast)[1]), 
                  sign = sign(as.numeric(fitted(forecast))[1])))
    } else {
      return(list(vol = sd(returns_series, na.rm = TRUE), 
                  sign = sign(mean(returns_series, na.rm = TRUE))))
    }
  }, error = function(e) {
    return(list(vol = sd(returns_series, na.rm = TRUE), 
                sign = sign(mean(returns_series, na.rm = TRUE))))
  })
}

# EGARCH
forecast_egarch <- function(returns_series) {
  tryCatch({
    spec <- ugarchspec(
      variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 0)),
      distribution.model = "std"
    )
    fit <- ugarchfit(spec, data = returns_series, solver = "hybrid")
    
    if (convergence(fit) == 0) {
      forecast <- ugarchforecast(fit, n.ahead = 1)
      return(list(vol = as.numeric(sigma(forecast)[1]), 
                  sign = sign(as.numeric(fitted(forecast))[1])))
    } else {
      return(list(vol = sd(returns_series, na.rm = TRUE), 
                  sign = sign(mean(returns_series, na.rm = TRUE))))
    }
  }, error = function(e) {
    return(list(vol = sd(returns_series, na.rm = TRUE), 
                sign = sign(mean(returns_series, na.rm = TRUE))))
  })
}

# GJR-GARCH
forecast_gjr <- function(returns_series) {
  tryCatch({
    spec <- ugarchspec(
      variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 0)),
      distribution.model = "std"
    )
    fit <- ugarchfit(spec, data = returns_series, solver = "hybrid")
    
    if (convergence(fit) == 0) {
      forecast <- ugarchforecast(fit, n.ahead = 1)
      return(list(vol = as.numeric(sigma(forecast)[1]), 
                  sign = sign(as.numeric(fitted(forecast))[1])))
    } else {
      return(list(vol = sd(returns_series, na.rm = TRUE), 
                  sign = sign(mean(returns_series, na.rm = TRUE))))
    }
  }, error = function(e) {
    return(list(vol = sd(returns_series, na.rm = TRUE), 
                sign = sign(mean(returns_series, na.rm = TRUE))))
  })
}

# FIGARCH
forecast_figarch <- function(returns_series) {
  tryCatch({
    spec <- ugarchspec(
      variance.model = list(model = "fiGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 0)),
      distribution.model = "std"
    )
    fit <- ugarchfit(spec, data = returns_series, solver = "hybrid")
    
    if (convergence(fit) == 0) {
      forecast <- ugarchforecast(fit, n.ahead = 1)
      return(list(vol = as.numeric(sigma(forecast)[1]), 
                  sign = sign(as.numeric(fitted(forecast))[1])))
    } else {
      return(list(vol = sd(returns_series, na.rm = TRUE), 
                  sign = sign(mean(returns_series, na.rm = TRUE))))
    }
  }, error = function(e) {
    return(list(vol = sd(returns_series, na.rm = TRUE), 
                sign = sign(mean(returns_series, na.rm = TRUE))))
  })
}

# DCC-GARCH
forecast_dcc <- function(returns_matrix) {
  tryCatch({
    n_assets <- ncol(returns_matrix)
    
    if (nrow(returns_matrix) < 30) {
      Sigma <- cov(returns_matrix)
      return(list(cov_matrix = Sigma, 
                  vols = sqrt(diag(Sigma)),
                  signs = sign(colMeans(returns_matrix))))
    }
    
    uspec <- multispec(replicate(n_assets, 
                                 ugarchspec(
                                   variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(1, 0)),
                                   distribution.model = "std"
                                 )))
    
    spec <- dccspec(uspec, dccOrder = c(1, 1), distribution = "mvnorm")
    fit <- dccfit(spec, data = returns_matrix, solver = "hybrid")
    
    forecast <- dccforecast(fit, n.ahead = 1)
    cov_forecast <- rcov(forecast)[[1]][,,1]
    mu_forecast <- fitted(forecast)[[1]][1, ]
    
    return(list(cov_matrix = cov_forecast,
                vols = sqrt(diag(cov_forecast)),
                signs = sign(mu_forecast)))
    
  }, error = function(e) {
    Sigma <- cov(returns_matrix)
    return(list(cov_matrix = Sigma,
                vols = sqrt(diag(Sigma)),
                signs = sign(colMeans(returns_matrix))))
  })
}

# LSTM-GARCH (эмуляция)
forecast_lstm_garch <- function(returns_series, lookback = 24) {
  # Эмуляция LSTM через комбинацию GARCH моделей
  garch1 <- forecast_garch(returns_series)
  garch2 <- forecast_egarch(returns_series)
  garch3 <- forecast_gjr(returns_series)
  
  vol <- mean(c(garch1$vol, garch2$vol, garch3$vol), na.rm = TRUE)
  sign_pred <- sign(sum(c(garch1$sign, garch2$sign, garch3$sign), na.rm = TRUE))
  
  return(list(vol = vol, sign = sign_pred, type = "ensemble"))
}

# Словарь методов прогнозирования
forecast_methods <- list(
  "none" = function(returns) {
    list(vol = sd(returns), sign = sign(mean(returns)), cov_matrix = NULL)
  },
  "garch" = function(returns) {
    if (is.matrix(returns) && ncol(returns) > 1) {
      vols <- apply(returns, 2, function(x) forecast_garch(x)$vol)
      signs <- apply(returns, 2, function(x) forecast_garch(x)$sign)
      list(vols = vols, signs = signs, cov_matrix = cor(returns))
    } else {
      forecast_garch(returns)
    }
  },
  "egarch" = function(returns) {
    if (is.matrix(returns) && ncol(returns) > 1) {
      vols <- apply(returns, 2, function(x) forecast_egarch(x)$vol)
      signs <- apply(returns, 2, function(x) forecast_egarch(x)$sign)
      list(vols = vols, signs = signs, cov_matrix = cor(returns))
    } else {
      forecast_egarch(returns)
    }
  },
  "gjrgarch" = function(returns) {
    if (is.matrix(returns) && ncol(returns) > 1) {
      vols <- apply(returns, 2, function(x) forecast_gjr(x)$vol)
      signs <- apply(returns, 2, function(x) forecast_gjr(x)$sign)
      list(vols = vols, signs = signs, cov_matrix = cor(returns))
    } else {
      forecast_gjr(returns)
    }
  },
  "figarch" = function(returns) {
    if (is.matrix(returns) && ncol(returns) > 1) {
      vols <- apply(returns, 2, function(x) forecast_figarch(x)$vol)
      signs <- apply(returns, 2, function(x) forecast_figarch(x)$sign)
      list(vols = vols, signs = signs, cov_matrix = cor(returns))
    } else {
      forecast_figarch(returns)
    }
  },
  "dcc" = function(returns) {
    forecast_dcc(returns)
  },
  "lstm" = function(returns) {
    if (is.matrix(returns) && ncol(returns) > 1) {
      vols <- numeric(ncol(returns))
      signs <- numeric(ncol(returns))
      for (i in 1:ncol(returns)) {
        lstm_pred <- forecast_lstm_garch(returns[, i])
        vols[i] <- lstm_pred$vol
        signs[i] <- lstm_pred$sign
      }
      list(vols = vols, signs = signs, cov_matrix = cor(returns))
    } else {
      forecast_lstm_garch(returns)
    }
  }
)