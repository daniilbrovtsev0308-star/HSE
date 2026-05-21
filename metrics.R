# =============================================================================
# ФУНКЦИИ РАСЧЕТА МЕТРИК ПОРТФЕЛЯ
# =============================================================================

# Функция для расчета метрик портфеля (с CAGR и Сортино)
calculate_metrics <- function(weights, returns_matrix) {
  # Доходность портфеля
  port_returns <- returns_matrix %*% weights
  
  # Основные метрики
  total_return <- prod(1 + port_returns) - 1
  n <- length(port_returns)
  
  # CAGR (Compound Annual Growth Rate) - для недельных данных
  cagr <- (1 + total_return)^(52/n) - 1
  
  # Годовая волатильность
  annual_vol <- sd(port_returns) * sqrt(52)
  
  # Коэффициент Шарпа (предполагаем безрисковую ставку 0%)
  sharpe <- cagr / annual_vol
  
  # Коэффициент Сортино (используем только отрицательные отклонения)
  downside <- port_returns[port_returns < 0]
  downside_dev <- if(length(downside) > 1) sd(downside) * sqrt(52) else 0
  sortino <- if(downside_dev > 0) cagr / downside_dev else NA
  
  # Максимальная просадка
  cum_returns <- cumprod(1 + port_returns)
  running_max <- cummax(cum_returns)
  drawdown <- (cum_returns - running_max) / running_max
  max_dd <- min(drawdown)
  
  # Win Rate
  win_rate <- mean(port_returns > 0)
  
  # VaR
  var_95 <- quantile(port_returns, 0.05) * sqrt(52)
  
  # Концентрация
  hhi <- sum(weights^2)
  eff_n <- 1 / hhi
  
  return(c(
    CAGR = cagr * 100,
    Волатильность = annual_vol * 100,
    Шарп = sharpe,
    Сортино = sortino,
    Просадка = max_dd * 100,
    VaR_95 = var_95 * 100,
    Win_Rate = win_rate * 100,
    HHI = hhi,
    Эфф_N = eff_n
  ))
}

# Функция для расчета метрик портфеля с учетом безрисковой ставки
calculate_portfolio_metrics <- function(weights, returns_matrix, initial_capital = 1000000, risk_free_rate = 0.05) {
  
  # Доходность портфеля
  port_returns <- returns_matrix %*% weights
  n <- length(port_returns)
  
  # Общая доходность
  total_return <- prod(1 + port_returns) - 1
  
  # Годовая доходность (CAGR)
  cagr <- (1 + total_return)^(52/n) - 1
  
  # Годовая волатильность
  annual_vol <- sd(port_returns) * sqrt(52)
  
  # Коэффициент Шарпа
  sharpe <- (cagr - risk_free_rate) / annual_vol
  
  # Максимальная просадка
  cum_returns <- cumprod(1 + port_returns)
  running_max <- cummax(cum_returns)
  drawdown <- (cum_returns - running_max) / running_max
  max_dd <- min(drawdown)
  
  # Коэффициент Сортино
  downside <- port_returns[port_returns < 0]
  downside_dev <- if(length(downside) > 1) sd(downside) * sqrt(52) else 0
  sortino <- if(downside_dev > 0) (cagr - risk_free_rate) / downside_dev else NA
  
  # VaR и CVaR
  var_95 <- quantile(port_returns, 0.05)
  cvar_95 <- mean(port_returns[port_returns <= var_95])
  
  # Win Rate
  win_rate <- mean(port_returns > 0)
  
  # Коэффициент Калмара
  calmar <- ifelse(abs(max_dd) > 0, cagr / abs(max_dd), NA)
  
  # Концентрация (HHI)
  hhi <- sum(weights^2)
  eff_n <- 1 / hhi
  
  return(list(
    total_return = total_return * 100,
    cagr = cagr * 100,
    annual_vol = annual_vol * 100,
    sharpe = sharpe,
    sortino = sortino,
    max_drawdown = max_dd * 100,
    var_95 = var_95 * 100,
    cvar_95 = cvar_95 * 100,
    win_rate = win_rate * 100,
    calmar = calmar,
    hhi = hhi,
    eff_n = eff_n
  ))
}

# Функция для расчета метрик бэктестинга
calculate_backtest_metrics <- function(returns_series, portfolio_values, initial_capital, step_size = 4) {
  
  n <- length(returns_series)
  n_years <- (n * step_size) / 52
  
  # Финальная стоимость
  final_value <- portfolio_values[length(portfolio_values)]
  total_return <- (final_value / initial_capital - 1) * 100
  
  # Годовая доходность (CAGR)
  if (n_years > 0) {
    cagr <- ((final_value / initial_capital)^(1/n_years) - 1) * 100
  } else {
    cagr <- NA
  }
  
  # Годовая волатильность
  annual_vol <- sd(returns_series, na.rm = TRUE) * sqrt(52/step_size) * 100
  
  # Коэффициент Шарпа
  sharpe <- ifelse(annual_vol > 0, cagr / annual_vol, NA)
  
  # Коэффициент Сортино
  downside <- returns_series[returns_series < 0]
  downside_dev <- if(length(downside) > 1) sd(downside) * sqrt(52/step_size) * 100 else 0
  sortino <- if(downside_dev > 0) cagr / downside_dev else NA
  
  # Максимальная просадка
  cum_wealth <- portfolio_values / initial_capital
  running_max <- cummax(cum_wealth)
  drawdown <- (cum_wealth - running_max) / running_max * 100
  max_dd <- min(drawdown, na.rm = TRUE)
  
  # Win Rate
  win_rate <- mean(returns_series > 0, na.rm = TRUE) * 100
  
  return(list(
    final_value = final_value,
    total_return = total_return,
    cagr = cagr,
    annual_vol = annual_vol,
    sharpe = sharpe,
    sortino = sortino,
    max_drawdown = max_dd,
    win_rate = win_rate,
    n_periods = n,
    n_years = n_years
  ))
}

# Функция для расчета транзакционных издержек
calculate_transaction_costs <- function(weights_history, transaction_cost, portfolio_values) {
  if (is.null(weights_history) || nrow(weights_history) < 2) {
    return(list(
      avg_turnover = NA,
      total_turnover = NA,
      avg_commission_pct = NA,
      total_commission_pct = NA,
      total_commission_abs = NA,
      n_rebalances = 0
    ))
  }
  
  n_periods <- nrow(weights_history)
  turnovers <- numeric(n_periods - 1)
  
  for (i in 2:n_periods) {
    weight_changes <- abs(weights_history[i, ] - weights_history[i-1, ])
    turnovers[i-1] <- sum(weight_changes) / 2
  }
  
  avg_turnover <- mean(turnovers, na.rm = TRUE)
  total_turnover <- sum(turnovers, na.rm = TRUE)
  
  avg_commission_pct <- avg_turnover * transaction_cost * 100
  total_commission_pct <- total_turnover * transaction_cost * 100
  
  avg_portfolio_value <- mean(portfolio_values, na.rm = TRUE)
  total_commission_abs <- total_turnover * transaction_cost * avg_portfolio_value
  
  return(list(
    avg_turnover = avg_turnover * 100,
    total_turnover = total_turnover * 100,
    avg_commission_pct = avg_commission_pct,
    total_commission_pct = total_commission_pct,
    total_commission_abs = total_commission_abs,
    n_rebalances = n_periods - 1
  ))
}