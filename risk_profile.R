# =============================================================================
# МОДУЛЬ ОПРЕДЕЛЕНИЯ РИСК-ПРОФИЛЯ ИНВЕСТОРА
# На основе результатов Тобит-регрессии из ВКР
# Формулировки вопросов соответствуют опросу "Исследование связи финансовых знаний
# и инвестиционного поведения"
# Версия: 3.0
# =============================================================================

# -----------------------------------------------------------------------------
# 1. КОЭФФИЦИЕНТЫ ТОБИТ-МОДЕЛИ (из таблицы 2.4.2 ВКР)
# -----------------------------------------------------------------------------
TOBIT_COEFFICIENTS <- list(
  const                = 37.063,
  gender               = 6.963,
  self_fin_est         = -7.532,
  risk_aversion_status = 3.870,
  age_sq               = -0.018,
  income_sq            = -3.734,
  self_compare_herd    = 1.319,
  age_income           = 0.656,
  risk_aversion_divers = -0.001
)

# Примечание: коэффициент risk_aversion_divers оставлен, но diversification
# теперь рассчитывается на основе квалификации инвестора (квалифицированные
# инвесторы считаются более диверсифицированными)

# -----------------------------------------------------------------------------
# 2. ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# -----------------------------------------------------------------------------

#' Ограничение значения диапазоном
clamp <- function(x, min_val, max_val) {
  max(min_val, min(x, max_val))
}

#' Оператор для значений по умолчанию
`%||%` <- function(a, b) {
  if (is.null(a) || is.na(a)) b else a
}

#' Преобразование ответов в числовые значения с валидацией
#' @param input_list список с ответами пользователя
#' @return список с проверенными числовыми значениями
validate_and_prepare <- function(input_list) {
  result <- list()
  
  # Пол (0 - женский, 1 - мужской)
  result$gender <- as.numeric(input_list$q_gender %||% 1)
  if (!result$gender %in% c(0, 1)) result$gender <- 1
  
  # Возраст
  result$age <- as.numeric(input_list$q_age %||% 30)
  if (is.na(result$age) || result$age < 18) result$age <- 30
  if (result$age > 100) result$age <- 70
  
  # Доход (код категории 1-4)
  result$income <- as.numeric(input_list$q_income %||% 2)
  result$income <- clamp(result$income, 1, 4)
  
  # Самооценка финансовых знаний (по шкале 0-4)
  # Вопрос: "Как вы оцениваете свои знания в области инвестиций и финансов?"
  result$self_fin_est <- as.numeric(input_list$q_self_fin_est %||% 2)
  result$self_fin_est <- clamp(result$self_fin_est, 0, 4)
  
  # Вопрос: "Я хорошо разбираюсь в финансовых вопросах по сравнению с другими инвесторами" (1-5)
  result$self_compare <- as.numeric(input_list$q_self_compare %||% 3)
  result$self_compare <- clamp(result$self_compare, 1, 5)
  
  # Вопрос: "Когда я вижу, что многие покупают какую-то акцию, я тоже склонен её купить" (1-5)
  result$herd_effect <- as.numeric(input_list$q_herd_effect %||% 3)
  result$herd_effect <- clamp(result$herd_effect, 1, 5)
  
  # Вопрос: "Какой тип портфеля больше всего соответствует вашему подходу к инвестированию?" (1-4)
  result$risk_aversion <- as.numeric(input_list$q_risk_aversion %||% 2)
  result$risk_aversion <- clamp(result$risk_aversion, 1, 4)
  
  # Вопрос: "Имеете ли вы статус квалифицированного инвестора?"
  # Квалифицированные инвесторы считаются более диверсифицированными
  result$qualified_investor <- as.numeric(input_list$q_qualified %||% 0)
  result$qualified_investor <- ifelse(result$qualified_investor > 0, 1, 0)
  
  # Диверсификация (прокси через статус квалифицированного инвестора)
  # Квалифицированные инвесторы с большей вероятностью имеют диверсифицированный портфель
  result$diversification <- result$qualified_investor
  
  return(result)
}

#' Расчет доли рисковых активов по Тобит-модели
#' @param data список с подготовленными переменными
#' @return доля рисковых активов (0-100)
calculate_risk_share <- function(data) {
  # Производные переменные
  age_sq <- data$age^2
  income_sq <- data$income^2
  age_income <- data$age * data$income
  self_compare_herd <- data$self_compare * data$herd_effect
  risk_aversion_divers <- data$risk_aversion * data$diversification
  
  # Расчёт по формуле
  share <- TOBIT_COEFFICIENTS$const +
    TOBIT_COEFFICIENTS$gender * data$gender +
    TOBIT_COEFFICIENTS$self_fin_est * data$self_fin_est +
    TOBIT_COEFFICIENTS$risk_aversion_status * data$risk_aversion +
    TOBIT_COEFFICIENTS$age_sq * age_sq +
    TOBIT_COEFFICIENTS$income_sq * income_sq +
    TOBIT_COEFFICIENTS$self_compare_herd * self_compare_herd +
    TOBIT_COEFFICIENTS$age_income * age_income +
    TOBIT_COEFFICIENTS$risk_aversion_divers * risk_aversion_divers
  
  # Цензурирование в диапазон [0, 100] и округление
  return(clamp(round(share, 1), 0, 100))
}

#' Определение риск-профиля на основе доли рисковых активов
#' @param risk_share доля рисковых активов (0-100)
#' @return список с профилем, гаммой и рекомендуемой стратегией
determine_profile <- function(risk_share) {
  if (risk_share <= 20) {
    list(profile = "Консервативный", gamma = 14, recommended = "MinVariance")
  } else if (risk_share <= 40) {
    list(profile = "Умеренно-консервативный", gamma = 11, recommended = "RiskParity")
  } else if (risk_share <= 60) {
    list(profile = "Сбалансированный", gamma = 8, recommended = "RiskParity")
  } else if (risk_share <= 80) {
    list(profile = "Умеренно-агрессивный", gamma = 5, recommended = "Tangency")
  } else {
    list(profile = "Агрессивный", gamma = 2, recommended = "RiskAversion")
  }
}

# -----------------------------------------------------------------------------
# 3. UI МОДУЛЯ
# -----------------------------------------------------------------------------
risk_profile_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("📊 Оценка инвестиционного профиля"),
    p("Ответьте на вопросы, чтобы определить ваш риск-профиль. 
      Результат будет использован для персонализации стратегии пенсионного портфеля."),
    hr(),
    
    wellPanel(
      h4("А. Ваш профиль"),
      
      radioButtons(ns("q_gender"), "А1. Ваш пол:",
                   choices = c("Женский" = 0, "Мужской" = 1), 
                   selected = 1, inline = TRUE),
      
      sliderInput(ns("q_age"), "А2. Ваш возраст (лет):",
                  min = 18, max = 70, value = 30, step = 1),
      
      radioButtons(ns("q_income"), "А3. Ваш среднемесячный доход:",
                   choices = c("до 50 000 руб." = 1, 
                               "50 000 – 100 000 руб." = 2,
                               "100 000 – 200 000 руб." = 3, 
                               "более 200 000 руб." = 4), 
                   selected = 2),
      
      radioButtons(ns("q_qualified"), "А4. Имеете ли вы статус квалифицированного инвестора?",
                   choices = c("Нет" = 0, "Да" = 1), 
                   selected = 0, inline = TRUE)
    ),
    
    wellPanel(
      h4("Б. Оценка знаний и поведения"),
      
      radioButtons(ns("q_self_fin_est"), "Б1. Как вы оцениваете свои знания в области инвестиций и финансов?",
                   choices = c("Очень низкие (почти ничего не понимаю)" = 0, 
                               "Скорее низкие (базовые понятия)" = 1,
                               "Средние (понимаю основные инструменты)" = 2,
                               "Высокие (активно инвестирую, разбираюсь в рынке)" = 3,
                               "Очень высокие (профессиональные знания)" = 4), 
                   selected = 2),
      
      radioButtons(ns("q_self_compare"), "Б2. Я хорошо разбираюсь в финансовых вопросах по сравнению с другими инвесторами:",
                   choices = c("1 — Абсолютно не согласен" = 1, 
                               "2 — Не согласен" = 2,
                               "3 — Нейтрально" = 3, 
                               "4 — Согласен" = 4, 
                               "5 — Полностью согласен" = 5), 
                   selected = 3),
      
      radioButtons(ns("q_herd_effect"), "Б3. Когда я вижу, что многие покупают какую-то акцию, я тоже склонен её купить:",
                   choices = c("1 — Абсолютно не согласен" = 1, 
                               "2 — Не согласен" = 2,
                               "3 — Нейтрально" = 3, 
                               "4 — Согласен" = 4, 
                               "5 — Полностью согласен" = 5), 
                   selected = 3)
    ),
    
    wellPanel(
      h4("В. Инвестиционные предпочтения"),
      
      radioButtons(ns("q_risk_aversion"), "В1. Какой тип портфеля больше всего соответствует вашему подходу к инвестированию?",
                   choices = c("Консервативный – главное сохранить капитал, готов мириться с невысокой доходностью" = 1, 
                               "Умеренный – готов к умеренному риску ради дохода выше депозита" = 2,
                               "Агрессивный – готов к высокому риску ради потенциально высокой доходности" = 3, 
                               "Спекулятивный – готов к очень высокому риску, возможны значительные потери" = 4), 
                   selected = 2)
    ),
    
    br(),
    actionButton(ns("calculate"), "Определить инвестиционный профиль", 
                 class = "btn-primary", style = "width: 100%"),
    br(), br(),
    uiOutput(ns("result"))
  )
}

# -----------------------------------------------------------------------------
# 4. СЕРВЕРНАЯ ЛОГИКА МОДУЛЯ
# -----------------------------------------------------------------------------
risk_profile_server <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    # Реактивное значение для хранения последнего рассчитанного профиля
    last_profile <- reactiveVal(NULL)
    
    observeEvent(input$calculate, {
      # Сбор и валидация ответов
      user_data <- validate_and_prepare(as.list(input))
      
      # Расчёт доли рисковых активов
      risk_share <- calculate_risk_share(user_data)
      
      # Определение профиля
      profile_info <- determine_profile(risk_share)
      
      # Сохранение в реактивные значения (для доступа из других модулей)
      values$gamma <- profile_info$gamma
      values$risk_profile <- profile_info$profile
      values$risk_share <- risk_share
      values$recommended_strategy <- profile_info$recommended
      
      last_profile(list(
        risk_share = risk_share,
        profile = profile_info$profile,
        gamma = profile_info$gamma,
        recommended = profile_info$recommended
      ))
      
      # Рендер результатов
      output$result <- renderUI({
        prof <- last_profile()
        
        tagList(
          hr(),
          fluidRow(
            column(4,
                   valueBox(
                     value = paste0(prof$risk_share, "%"), 
                     subtitle = "Рекомендуемая доля рисковых активов",
                     icon = icon("percent"), 
                     color = "blue", 
                     width = 12
                   )
            ),
            column(4,
                   valueBox(
                     value = prof$profile, 
                     subtitle = "Ваш инвестиционный профиль",
                     icon = icon("user-tag"), 
                     color = "green", 
                     width = 12
                   )
            ),
            column(4,
                   valueBox(
                     value = paste0("γ = ", prof$gamma), 
                     subtitle = "Коэффициент неприятия риска",
                     icon = icon("balance-scale"), 
                     color = "purple", 
                     width = 12
                   )
            )
          ),
          box(
            title = paste0("📊 Ваш профиль: ", prof$profile),
            status = "success",
            solidHeader = TRUE,
            width = 12,
            p(strong("Рекомендуемая стратегия управления портфелем:"), 
              code(prof$recommended)),
            p(strong("Коэффициент γ =", prof$gamma), 
              "— этот параметр будет использован при оптимизации портфеля"),
            hr(),
            h5("💡 Рекомендации:"),
            tags$ul(
              tags$li("Регулярно пересматривайте состав портфеля (1-2 раза в год)"),
              tags$li("По мере приближения к пенсионному возрасту снижайте долю рисковых активов"),
              tags$li("Диверсификация между разными классами активов снижает риски"),
              tags$li("Повышайте финансовую грамотность — это помогает принимать взвешенные решения")
            ),
            hr(),
            p(em("Расчёт выполнен на основе эконометрической модели (Тобит-регрессия), 
                  оценённой по данным опроса частных инвесторов (N=133)."),
              style = "font-size: 12px; color: gray;")
          )
        )
      })
    })
  })
}

# -----------------------------------------------------------------------------
# 5. ПРИМЕР ИСПОЛЬЗОВАНИЯ В ОСНОВНОМ APP
# -----------------------------------------------------------------------------
# В файле app.R нужно:
#
# source("risk_profile.R")
#
# ui <- fluidPage(
#   risk_profile_ui("riskProfile")
# )
#
# server <- function(input, output, session) {
#   values <- reactiveValues(gamma = 8, risk_profile = "Сбалансированный")
#   risk_profile_server("riskProfile", values)
#   
#   # В других модулях можно использовать values$gamma
#   # Например, в модуле оптимизации портфеля
# }
#
# shinyApp(ui, server)