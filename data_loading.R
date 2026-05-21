# =============================================================================
# ФУНКЦИИ ЗАГРУЗКИ ДАННЫХ (ОПТИМИЗИРОВАННАЯ ВЕРСИЯ)
# Версия 4.0 - быстрая агрегация
# =============================================================================

library(httr)
library(jsonlite)
library(tidyverse)
library(lubridate)
library(quantmod)
library(tidyquant)

# =============================================================================
# ФУНКЦИИ ДЛЯ MOEX
# =============================================================================

# Загрузка данных индексов с MOEX (оптимизированная)
get_moex_index_data <- function(ticker, start_date, end_date) {
  cat(paste("Загрузка индекса", ticker, "..."))
  
  tryCatch({
    url <- paste0("https://iss.moex.com/iss/history/engines/stock/markets/index/securities/", 
                  ticker, ".json")
    
    params <- list(
      from = format(as.Date(start_date), "%Y-%m-%d"),
      till = format(as.Date(end_date), "%Y-%m-%d"),
      start = 0,
      limit = 100
    )
    
    all_data <- data.frame()
    
    while(TRUE) {
      response <- GET(url, query = params, timeout(30))
      
      if (status_code(response) != 200) {
        cat(" ✗ Ошибка HTTP\n")
        return(NULL)
      }
      
      json_data <- fromJSON(content(response, "text"), flatten = TRUE)
      
      if (is.null(json_data$history) || is.null(json_data$history$data)) {
        cat(" ✗ Нет данных\n")
        return(NULL)
      }
      
      if (length(json_data$history$data) == 0) {
        break
      }
      
      columns <- unlist(json_data$history$columns)
      page_data <- as.data.frame(json_data$history$data, stringsAsFactors = FALSE)
      
      if (length(columns) == ncol(page_data)) {
        colnames(page_data) <- columns
        
        if ("TRADEDATE" %in% colnames(page_data) && "CLOSE" %in% colnames(page_data)) {
          temp_df <- data.frame(
            Date = as.Date(page_data$TRADEDATE),
            Price = as.numeric(gsub(",", ".", as.character(page_data$CLOSE)))
          )
          temp_df <- temp_df[!is.na(temp_df$Price) & temp_df$Price > 0, ]
          
          if (nrow(temp_df) > 0) {
            all_data <- rbind(all_data, temp_df)
          }
        }
      }
      
      if (nrow(page_data) < params$limit) {
        break
      }
      
      params$start <- params$start + params$limit
    }
    
    if (nrow(all_data) == 0) {
      cat(" ✗ Нет данных\n")
      return(NULL)
    }
    
    all_data <- all_data %>%
      distinct(Date, .keep_all = TRUE) %>%
      arrange(Date)
    
    colnames(all_data)[2] <- ticker
    
    cat(paste(" ✓ Загружено", nrow(all_data), "записей\n"))
    return(all_data)
    
  }, error = function(e) {
    cat(paste(" ✗ Ошибка:", e$message, "\n"))
    return(NULL)
  })
}

# Загрузка данных акций с MOEX
get_moex_shares_data <- function(ticker, start_date, end_date) {
  cat(paste("Загрузка акции", ticker, "..."))
  
  tryCatch({
    url <- paste0("https://iss.moex.com/iss/history/engines/stock/markets/shares/boards/tqbr/securities/", 
                  ticker, ".json")
    
    params <- list(
      from = format(as.Date(start_date), "%Y-%m-%d"),
      till = format(as.Date(end_date), "%Y-%m-%d"),
      start = 0,
      limit = 100
    )
    
    all_data <- data.frame()
    
    while(TRUE) {
      response <- GET(url, query = params, timeout(30))
      
      if (status_code(response) != 200) {
        cat(" ✗ Ошибка HTTP\n")
        return(NULL)
      }
      
      json_data <- fromJSON(content(response, "text"), flatten = TRUE)
      
      if (is.null(json_data$history) || is.null(json_data$history$data)) {
        cat(" ✗ Нет данных\n")
        return(NULL)
      }
      
      if (length(json_data$history$data) == 0) {
        break
      }
      
      columns <- unlist(json_data$history$columns)
      page_data <- as.data.frame(json_data$history$data, stringsAsFactors = FALSE)
      
      if (length(columns) == ncol(page_data)) {
        colnames(page_data) <- columns
        
        if ("TRADEDATE" %in% colnames(page_data) && "CLOSE" %in% colnames(page_data)) {
          temp_df <- data.frame(
            Date = as.Date(page_data$TRADEDATE),
            Price = as.numeric(gsub(",", ".", as.character(page_data$CLOSE)))
          )
          temp_df <- temp_df[!is.na(temp_df$Price) & temp_df$Price > 0, ]
          
          if (nrow(temp_df) > 0) {
            all_data <- rbind(all_data, temp_df)
          }
        }
      }
      
      if (nrow(page_data) < params$limit) {
        break
      }
      
      params$start <- params$start + params$limit
    }
    
    if (nrow(all_data) == 0) {
      cat(" ✗ Нет данных\n")
      return(NULL)
    }
    
    all_data <- all_data %>%
      distinct(Date, .keep_all = TRUE) %>%
      arrange(Date)
    
    colnames(all_data)[2] <- ticker
    
    cat(paste(" ✓ Загружено", nrow(all_data), "записей\n"))
    return(all_data)
    
  }, error = function(e) {
    cat(paste(" ✗ Ошибка:", e$message, "\n"))
    return(NULL)
  })
}

# ОСНОВНАЯ ФУНКЦИЯ ЗАГРУЗКИ MOEX (С БЫСТРОЙ АГРЕГАЦИЕЙ)
load_moex_complete <- function(tickers, start_date, end_date, period = "weekly") {
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""))
  cat("\n")
  cat("ЗАГРУЗКА ДАННЫХ MOEX\n")
  cat(paste(rep("=", 60), collapse = ""))
  cat("\n")
  cat(paste("Период:", start_date, "-", end_date, "\n"))
  cat(paste("Таймфрейм:", period, "\n\n"))
  
  indices <- c("MCFTR", "RGBITR", "IMOEX", "RGBI", "RTSI")
  
  data_list <- list()
  
  for (ticker in tickers) {
    if (ticker %in% indices) {
      ticker_data <- get_moex_index_data(ticker, start_date, end_date)
    } else {
      ticker_data <- get_moex_shares_data(ticker, start_date, end_date)
    }
    
    if (!is.null(ticker_data) && nrow(ticker_data) > 0) {
      data_list[[ticker]] <- ticker_data
    } else {
      cat(paste("✗ Не удалось загрузить", ticker, "\n"))
      return(NULL)
    }
  }
  
  if (length(data_list) == 0) {
    cat("\n✗ Не удалось загрузить данные\n")
    return(NULL)
  }
  
  cat("\nОбъединение данных...\n")
  
  # Быстрое объединение
  result <- data_list[[1]]
  if (length(data_list) > 1) {
    for (i in 2:length(data_list)) {
      result <- merge(result, data_list[[i]], by = "Date", all = TRUE)
    }
  }
  
  result <- result[order(result$Date), ]
  
  # Заполнение пропусков (оптимизированное)
  cat("Заполнение пропусков...\n")
  for (col in 2:ncol(result)) {
    result[[col]] <- zoo::na.locf(result[[col]], na.rm = FALSE)
  }
  
  # БЫСТРАЯ АГРЕГАЦИЯ ПО ТАЙМФРЕЙМУ (без циклов)
  if (period != "daily") {
    cat(paste("Агрегация:", period, "...\n"))
    
    if (period == "weekly") {
      result$Week <- floor_date(result$Date, "week", week_start = 1)
      
      # Быстрая агрегация через dplyr
      result <- result %>%
        group_by(Week) %>%
        summarise(
          Date = max(Date),
          across(all_of(tickers), ~ dplyr::last(na.omit(.x))),
          .groups = "drop"
        ) %>%
        filter(!is.na(Date)) %>%
        arrange(Date)
      
      result$Date <- as.Date(result$Date)
      
    } else if (period == "monthly") {
      result$YearMonth <- format(result$Date, "%Y-%m")
      
      result <- result %>%
        group_by(YearMonth) %>%
        summarise(
          Date = max(Date),
          across(all_of(tickers), ~ last(na.omit(.x))),
          .groups = "drop"
        ) %>%
        filter(!is.na(Date)) %>%
        arrange(Date)
      
      result$Date <- as.Date(result$Date)
    }
  }
  
  # Удаляем временные колонки
  result <- result %>%
    select(-any_of(c("Week", "YearMonth")))
  
  # Фильтруем по датам
  result <- result %>%
    filter(Date >= as.Date(start_date) & Date <= as.Date(end_date))
  
  # Удаляем строки с полностью пустыми данными
  result <- result[rowSums(!is.na(result[, -1, drop = FALSE])) > 0, ]
  
  if (nrow(result) == 0) {
    cat("\n✗ Нет данных после агрегации\n")
    return(NULL)
  }
  
  cat("\n")
  cat(paste(rep("-", 60), collapse = ""))
  cat("\n")
  cat("✓ ЗАГРУЗКА ЗАВЕРШЕНА\n")
  cat(paste("  Наблюдений:", nrow(result), "\n"))
  cat(paste("  Активов:", ncol(result) - 1, "\n"))
  cat(paste("  Диапазон дат:", format(min(result$Date), "%d.%m.%Y"), 
            "-", format(max(result$Date), "%d.%m.%Y"), "\n"))
  
  return(result)
}

# =============================================================================
# ФУНКЦИИ ДЛЯ МИРОВЫХ АКТИВОВ (YAHOO FINANCE)
# =============================================================================

robust_data_download <- function(ticker, symbol, from, to) {
  
  from_date <- as.Date(from)
  to_date <- as.Date(to)
  
  cat(paste("  Загрузка", ticker, "...\n"))
  
  # Попытка через tq_get
  cat("    tq_get: ")
  tryCatch({
    data <- tq_get(symbol, from = from_date, to = to_date)
    if (!is.null(data) && nrow(data) > 10) {
      if ("adjusted" %in% colnames(data)) {
        results <- data.frame(Date = data$date, Price = data$adjusted)
      } else if ("close" %in% colnames(data)) {
        results <- data.frame(Date = data$date, Price = data$close)
      } else {
        results <- data.frame(Date = data$date, Price = data$close)
      }
      results <- results[!is.na(results$Price), ]
      colnames(results)[2] <- ticker
      cat(paste("✓ Успешно (", nrow(results), "записей)\n", sep=""))
      return(results)
    }
    cat("✗\n")
  }, error = function(e) cat("✗\n"))
  
  # Попытка через quantmod
  cat("    quantmod: ")
  tryCatch({
    data <- getSymbols(symbol, from = from_date, to = to_date, auto.assign = FALSE, warnings = FALSE)
    if (!is.null(data) && nrow(data) > 10) {
      dates <- index(data)
      prices <- as.numeric(Cl(data))
      valid_idx <- !is.na(prices)
      if (sum(valid_idx) > 10) {
        results <- data.frame(Date = dates[valid_idx], Price = prices[valid_idx])
        colnames(results)[2] <- ticker
        cat(paste("✓ Успешно (", nrow(results), "записей)\n", sep=""))
        return(results)
      }
    }
    cat("✗\n")
  }, error = function(e) cat("✗\n"))
  
  cat("    ✗ Не удалось загрузить", ticker, "\n")
  return(NULL)
}

# Загрузка мировых данных (оптимизированная)
load_world_data <- function(tickers, start_date, end_date, period = "weekly") {
  
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""))
  cat("\n")
  cat("ЗАГРУЗКА МИРОВЫХ ДАННЫХ\n")
  cat(paste(rep("=", 60), collapse = ""))
  cat("\n")
  
  world_data_list <- list()
  
  for (ticker in tickers) {
    symbol <- switch(ticker,
                     "GOLD" = "GC=F",
                     "USDRUB" = "RUB=X", 
                     "Bitcoin" = "BTC-USD",
                     "SPXTR" = "^GSPC",
                     ticker
    )
    
    data <- robust_data_download(ticker, symbol, start_date, end_date)
    
    if (!is.null(data) && nrow(data) > 0) {
      world_data_list[[ticker]] <- data
    } else {
      cat(paste("✗ Не удалось загрузить", ticker, "\n"))
      return(NULL)
    }
  }
  
  if (length(world_data_list) == 0) {
    return(NULL)
  }
  
  # Быстрое объединение
  all_dates <- sort(unique(unlist(lapply(world_data_list, function(x) x$Date))))
  data_world <- data.frame(Date = all_dates)
  
  for (ticker in names(world_data_list)) {
    data_world <- left_join(data_world, world_data_list[[ticker]], by = "Date")
  }
  
  # Заполнение пропусков
  for (col in 2:ncol(data_world)) {
    data_world[[col]] <- zoo::na.locf(data_world[[col]], na.rm = FALSE)
  }
  
  # Быстрая агрегация
  if (period != "daily") {
    if (period == "weekly") {
      data_world$Week <- floor_date(data_world$Date, "week")
      data_world <- data_world %>%
        group_by(Week) %>%
        summarise(
          Date = max(Date),
          across(all_of(tickers), ~ last(na.omit(.x))),
          .groups = "drop"
        ) %>%
        filter(!is.na(Date)) %>%
        arrange(Date)
      data_world$Date <- as.Date(data_world$Date)
      
    } else if (period == "monthly") {
      data_world$YearMonth <- format(data_world$Date, "%Y-%m")
      data_world <- data_world %>%
        group_by(YearMonth) %>%
        summarise(
          Date = max(Date),
          across(all_of(tickers), ~ last(na.omit(.x))),
          .groups = "drop"
        ) %>%
        filter(!is.na(Date)) %>%
        arrange(Date)
      data_world$Date <- as.Date(data_world$Date)
    }
    
    data_world <- data_world %>% select(-any_of(c("Week", "YearMonth")))
  }
  
  cat("\n✓ ЗАГРУЗКА МИРОВЫХ ДАННЫХ ЗАВЕРШЕНА\n")
  cat(paste("  Наблюдений:", nrow(data_world), "\n"))
  cat(paste("  Активов:", ncol(data_world) - 1, "\n"))
  
  return(data_world)
}

# =============================================================================
# ОСНОВНАЯ ФУНКЦИЯ ЗАГРУЗКИ ВСЕХ ДАННЫХ
# =============================================================================

load_all_data <- function(rus_tickers, world_tickers, start_date, end_date, period = "weekly") {
  
  data_list <- list()
  
  if (length(rus_tickers) > 0) {
    data_rus <- load_moex_complete(rus_tickers, start_date, end_date, period)
    if (!is.null(data_rus)) data_list[["rus"]] <- data_rus
  }
  
  if (length(world_tickers) > 0) {
    data_world <- load_world_data(world_tickers, start_date, end_date, period)
    if (!is.null(data_world)) data_list[["world"]] <- data_world
  }
  
  if (length(data_list) == 0) {
    return(NULL)
  }
  
  result <- data_list[[1]]
  for (i in 2:length(data_list)) {
    result <- full_join(result, data_list[[i]], by = "Date")
  }
  
  result <- result[order(result$Date), ]
  
  # Заполнение пропусков
  for (col in 2:ncol(result)) {
    result[[col]] <- zoo::na.locf(result[[col]], na.rm = FALSE)
  }
  
  # Удаляем строки с полностью пустыми данными
  result <- result[rowSums(!is.na(result[, -1, drop = FALSE])) > 0, ]
  
  return(result)
}