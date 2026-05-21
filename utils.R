# =============================================================================
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# =============================================================================

# Функция для форматирования дат на графиках
format_dates_for_plot <- function(dates) {
  min_date <- min(dates)
  max_date <- max(dates)
  
  n_weeks <- as.numeric(difftime(max_date, min_date, units = "weeks"))
  
  if (n_weeks > 520) {
    by_period <- "1 year"
    date_format <- "%Y"
  } else if (n_weeks > 260) {
    by_period <- "6 months"
    date_format <- "%b %Y"
  } else if (n_weeks > 104) {
    by_period <- "3 months"
    date_format <- "%b %Y"
  } else if (n_weeks > 52) {
    by_period <- "2 months"
    date_format <- "%b %Y"
  } else {
    by_period <- "1 month"
    date_format <- "%b %Y"
  }
  
  all_dates <- seq(min_date, max_date, by = by_period)
  formatted <- format(all_dates, date_format)
  
  return(list(breaks = all_dates, labels = formatted))
}

# Функция для получения цветов портфеля
get_portfolio_colors <- function(n) {
  if (n <= 8) {
    return(RColorBrewer::brewer.pal(max(3, n), "Set2"))
  } else {
    return(viridis::viridis(n))
  }
}

# Функция для генерации эффективной границы
generate_efficient_frontier <- function(returns_matrix, n_portfolios = 1000) {
  n_assets <- ncol(returns_matrix)
  mu <- colMeans(returns_matrix) * 52
  Sigma <- cov(returns_matrix) * 52
  
  weights <- matrix(runif(n_portfolios * n_assets), ncol = n_assets)
  weights <- weights / rowSums(weights)
  
  port_returns <- weights %*% mu
  port_risks <- sqrt(diag(weights %*% Sigma %*% t(weights)))
  
  return(data.frame(
    risk = port_risks * 100,
    return = port_returns * 100
  ))
}

# Функция для обработки выбросов
replace_outliers <- function(data, method = "IQR", threshold = 3) {
  if (!is.numeric(data)) {
    stop("Данные должны быть числовым вектором")
  }
  
  cleaned_data <- data
  
  if (method == "IQR") {
    Q <- quantile(data, probs = c(0.25, 0.75), na.rm = TRUE)
    IQR <- Q[2] - Q[1]
    lower_bound <- Q[1] - threshold * IQR
    upper_bound <- Q[2] + threshold * IQR
    outliers <- which(data < lower_bound | data > upper_bound)
  } else if (method == "sd") {
    mean_val <- mean(data, na.rm = TRUE)
    sd_val <- sd(data, na.rm = TRUE)
    lower_bound <- mean_val - threshold * sd_val
    upper_bound <- mean_val + threshold * sd_val
    outliers <- which(data < lower_bound | data > upper_bound)
  } else {
    stop("Неподдерживаемый метод. Используйте 'IQR' или 'sd'")
  }
  
  for (i in outliers) {
    prev_val <- ifelse(i > 1, data[i - 1], NA)
    next_val <- ifelse(i < length(data), data[i + 1], NA)
    
    if (!is.na(prev_val) & !is.na(next_val)) {
      cleaned_data[i] <- mean(c(prev_val, next_val))
    } else if (!is.na(prev_val)) {
      cleaned_data[i] <- prev_val
    } else if (!is.na(next_val)) {
      cleaned_data[i] <- next_val
    }
  }
  
  return(cleaned_data)
}