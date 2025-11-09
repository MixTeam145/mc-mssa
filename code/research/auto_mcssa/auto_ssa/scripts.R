source("auto_ssa.R", TRUE)
source("mcmssa_utils.R", TRUE)
# source("ic.R")
#source("ic_AR1noise.R")

estimate_rank <- function(time_series,
                          L = length(time_series) %/% 2,
                          method = "svd",
                          rank_range=0:30,
                          check_stat=FALSE,
                          ...){
  # Метод для автоматической оценки ранга временного ряда
  # 
  # :param time_series: временной ряд
  # :param L: длина окна
  # :param method: название метода, 
  #                'svd', 'bic', 'trmat', 'trmat_fast', 'mgn', 
  #                   'samos', 'ester', 'mdl' - для белого шума; 
  #                'svd_ar1', 'bic_ar1', 'trmat_ar1', 'mdl_ar1' - для красного шума; 
  #                'svd_ar2', 'bic_ar2', 'trmat_ar2', 'mdl_ar2' - для AR(2)
  # :param check_stat: параметр для методов AR(1), AR(2)
  h <- calculate_ic_factory(method = method,
                            r_range = rank_range,
                            L = L,
                            check_stat=check_stat,
                            ...)
  r_range <- rank_range
  result <- h(time_series[!is.na(time_series)])
  return(r_range[which.max(result$bic)])
}

auto_trend_model <- function(ssa_obj,
                             method = c("iossa.post_grouping", "iossa.manual", "iossa.auto_grouping.ssa",
                                        "iossa.auto_grouping.fossa", "eossa.auto", "basic.auto", "fossa.auto"),
                             signal_rank,
                             base='series',
                             ...){
  # Input: methods --- vector of names of ssa methods:
  #           "iossa.post_grouping"       --- IOSSA with separate initial grouping, using 
  #                                           automated trend identification after IOSSA. 
  #                                           Tolerance parameter is  provided in 
  #                                           parameters["iossa_tolerance"]. Maxiter  parameter
  #                                           is in parameters["iossa_maxiter"].
  #           "iossa.manual"              --- IOSSA with manual initial grouping, grouping is 
  #                                           provided in parameters["manual_grouping"]. 
  #                                           Tolerance parameter is provided in 
  #                                           parameters["iossa_tolerance"]. Maxiter parameter 
  #                                           is in parameters["iossa_maxiter"].
  #           "iossa.auto_grouping.ssa"   --- IOSSA with auto initial grouping, Tolerance 
  #                                           parameter is provided in 
  #                                           parameters["iossa_tolerance"]. Maxiter
  #                                           parameter is in parameters["iossa_maxiter"].
  #           "iossa.auto_grouping.fossa" --- same as iossa.auto_grouping.ssa but using
  #                                           FOSSA for initial grouping.
  #           "eossa.auto"                --- EOSSA with auto grouping.
  #           "basic.auto"                --- SSA with auto grouping.
  #           "fossa.auto"                --- FOSSA with auto grouping.
  #       parameters --- list with method parameters:
  #           "series_length"             --- time series length.
  #           "signal_rank"               --- rank of signal for nested decomposition.
  #           "window_length"             --- window length in SSA.
  #           "auto_trend_freq"           --- trend max frequency for auto grouping.
  #           "auto_threshold"            --- threshold for auto grouping.
  #           "iossa_maxiter"             --- max num iterations in IOSSA.
  #           "iossa_tolerance"           --- tolerance in IOSSA.
  #           "manual_grouping"           --- nested grouping, list of vectors, in IOSSA.
  # Output: function creates and returns an object of class auto_trend,
  #         object contains result of decomposition and ssa information.
  
  if (!prod(method %in% c("iossa.post_grouping", "iossa.manual", "iossa.auto_grouping.ssa",
                          "iossa.auto_grouping.fossa", "eossa.auto", "basic.auto", "fossa.auto"))){
    stop("Wrong method parameter.")
  }
  method <- c(method)
  
  kwargs <- list(...)
  
  this <- list(method = method)
  attr(this, ".env") <- new.env();
  
  if (length(method) > 1){
    for (m in method){
      this[[m]] <- list()
    }
  }
  
  results <- list()
  
  class(this) <- c("auto_trend")
  if ("iossa.post_grouping" %in% method){
    res <- iossa_post_grouping_trend(ssa_obj = ssa_obj, 
                                     signal_rank = signal_rank, 
                                     auto_trend_freq = kwargs$auto_trend_freq, 
                                     auto_threshold = kwargs$auto_threshold,
                                     iossa_maxiter = kwargs$iossa_maxiter, 
                                     iossa_tolerance = kwargs$iossa_tolerance,
                                     kappa = kwargs$kappa)
    results[["iossa.post_grouping"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else {
      this[["iossa.post_grouping"]]$iterations = res$iterations
    }
  }
  if ("iossa.manual" %in% method){
    res <- iossa_manual_trend(ssa_obj = ssa_obj, 
                              signal_rank = signal_rank, 
                              auto_trend_freq = auto_trend_freq, 
                              auto_threshold = auto_threshold,
                              iossa_maxiter = iossa_maxiter, 
                              iossa_tolerance = iossa_tolerance, 
                              manual_grouping = manual_grouping,
                              kappa = kappa)
    results[["iossa.manual"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else{
      this[["iossa.manual"]]$iterations = res$iterations
    }
  }
  if ("iossa.auto_grouping.ssa" %in% method){
    res <- iossa_auto_grouping_ssa_trend(ssa_obj = ssa_obj, 
                                         signal_rank = signal_rank, 
                                         auto_trend_freq = auto_trend_freq, 
                                         auto_threshold = auto_threshold,
                                         auto_initial_threshold = auto_initial_threshold, 
                                         iossa_tolerance = iossa_tolerance,
                                         kappa = kappa)
    results[["iossa.auto_grouping.ssa"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else{
      this[["iossa.auto_grouping.ssa"]]$iterations = res$iterations
    }
  }
  if ("iossa.auto_grouping.fossa" %in% method){
    res <- iossa_auto_grouping_fossa_trend(ssa_obj = ssa_obj, 
                                           signal_rank = signal_rank, 
                                           auto_trend_freq = kwargs$auto_trend_freq, 
                                           auto_threshold = kwargs$auto_threshold,
                                           iossa_tolerance = kwargs$iossa_tolerance,
                                           iossa_maxiter = kwargs$iossa_maxiter,
                                           kappa = kwargs$iossa_kappa)
    results[["iossa.auto_grouping.fossa"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else{
      this[["iossa.auto_grouping.fossa"]]$iterations = res$iterations
    }
  }
  if ("eossa.auto" %in% method){
    res <- eossa_auto_trend(ssa_obj, 
                            signal_rank = signal_rank, 
                            clust_type = kwargs$clust_type,
                            auto_trend_freq = kwargs$auto_trend_freq, 
                            auto_threshold = kwargs$auto_threshold,
                            delta = kwargs$delta,
                            base=base)
    results[["eossa.auto"]] <- res
  }
  if ("basic.auto" %in% method){
    res <- basic_auto_trend(ssa_obj = ssa_obj, 
                            signal_rank = signal_rank, 
                            auto_trend_freq = kwargs$auto_trend_freq, 
                            auto_threshold = kwargs$auto_threshold)
    results[["basic.auto"]] <- res
  }
  if ("fossa.auto" %in% method){
    res <- fossa_auto_trend(ssa_obj, 
                            signal_rank = signal_rank, 
                            auto_trend_freq = auto_trend_freq, 
                            auto_threshold = auto_threshold)
    results[["fossa.auto"]] <- res
  }
  
  this$original <- ssa_obj$F
  this$signal_rank <- signal_rank
  
  if (length(method) == 1){
    this$model_object <- res$model_object
    this$grouping <- res$grouping
    this$trend <- res$trend
    this$periodics <- res$periodics
    this$residuals <- res$residuals
  }
  else{
    for (m in method){
      this[[m]]$model_object <- results[[m]]$model_object
      this[[m]]$grouping <- results[[m]]$grouping
      this[[m]]$trend <- results[[m]]$trend
      this[[m]]$periodics <- results[[m]]$periodics
      this[[m]]$residuals <- results[[m]]$residuals
    }
  }
  
  this
}

auto_periodic_model <- function(ssa_obj,
                                groups,
                                method = "frequency",
                                ...) {
  # Input:
  #       ssa_obj --- ssa object.
  #       groups  --- list with groups for periodic identification.
  #       method --- string with name of periodic identification method. Must be one of:
  #           "frequency"                 --- frequency method (default).
  #           "angle_reg"                 --- angle regularity method.
  #     *extra parameters*:
  #       s0                              --- parameter for frequency method 
  #                                           (frequencies proximity), default = 1
  #       rho0                            --- threshold for rho measure, default = 0.9
  
  allowed_methods <- c("frequency", "angle_reg")
  
  if (!(method %in% allowed_methods))
    stop(paste0("Wrong method parameter, must be one of: ", allowed_methods))
  
  this <- list(method = method)
  attr(this, ".env") <- new.env();
  class(this) <- c("auto_periodic_model")
  
  if (method == "frequency"){
    res <- periodic_grouping_freq(ssa_obj = ssa_obj,
                                  groups = groups,
                                  ...)
  }
  if (method == "angle_reg"){
    res <- periodic_grouping_angle_reg(ssa_obj = ssa_obj,
                                       groups = groups,
                                       ...)
    this$two_el_tau <- res$two_el_tau
    this$one_el_gamma <- res$one_el_gamma
  }
  
  this$series <- list()
  I_1_harmonics_num <- length(res$I_1) %/% 2
  
  index <- 1
  while (index <= I_1_harmonics_num){
    series_name <- paste0("P", index)
    this$series[[series_name]] <- list(series = reconstruct(ssa_obj, 
                                                            groups = list(c(res$I_1[index], 
                                                                            res$I_1[index + 1])))$F1,
                                       freq = res$I_1_freqs)
    index <- index + 1
  }
  
  index <- 1
  while (index <= length(res$I_2)){
    series_name <- paste0("P", index + I_1_harmonics_num)
    this$series[[series_name]] <- list(series = reconstruct(ssa_obj, 
                                                            groups = list(c(res$I_2[index])))$F1,
                                       freq = 0.5)
    index <- index + 1
  }
  
  
  this$two_el_indices <- res$I_1
  this$two_el_freqs <- res$I_1_freqs
  this$one_el_indices <- res$I_2
  this$undetermined_indices <- setdiff(unique(groups), 
                                       c(res$I_1, res$I_2))
  if (length(this$undetermined_indices) > 0){
    this[["undetermined_series"]] <- reconstruct(ssa_obj, 
                                                 groups = list(this$undetermined_indices))$F1
  } else {
    this[["undetermined_series"]] <- numeric(0)
  }
  
  this
}

noise_hypothesis <- function(time_series,
                             method='mc',
                             sign_level=0.05,
                             noise_type='red',
                             L=NULL,
                             fixed_freqs=FALSE,
                             freq.range=NULL){
  # Функция для проверки гипотезы о шуме
  # 
  # :param time_series: ряд, для которого проверяется гипотеза
  # :param method: критерий, для noise_type=='white' поддерживаются 'mc', 'box',
  #                'wavelet'; для noise_type=='red' поддерживается 'mc'
  # :param sign_level: уровень значимости
  # :param noise_type: тип шума, поддерживаются 'white', 'red'
  # :param L: длина окна, только для method=='mc'
  # :param fixed_freqs: флаг для оценки параметров шума по части спектра (только
  #                     для method=='mc')
  # :param freq.range: диапазон частот шума, по которым вычисляется оценка 
  #                    параметров, только для fixed_freqs==TRUE
  # 
  # :return: list с ключами
  #   :key result: TRUE, если гипотеза не отвергается, FALSE - если отвергается
  #   :key pval: pvalue
  
  if (method == "mc"){
    if (fixed_freqs)
      pval <- MonteCarloSSA(
        f = time_series,
        freq.range = freq.range,
        L = L,
        model = NULL,
        basis = 'sin0',
        kind = "ev",
        D = 1,
        G = 1000,
        level.conf = NULL,
        composite = TRUE,
        noise_type=noise_type
      )$p.value
    else
      pval <- MonteCarloSSA(
        f = time_series,
        L = L,
        model = NULL,
        basis = 'sin0',
        kind = "ev",
        D = 1,
        G = 1000,
        level.conf = NULL,
        composite = TRUE,
        noise_type=noise_type
      )$p.value
  } else if (method == "box"){
    pval <- Box.test(time_series)$p.value
  } else if (method == "wavelet"){
    pval <- genwwn.test(time_series, filter.number = 10)$p.value
  }
  
  return(list('result'=!(pval<sign_level),
              'pval'=pval))
}

get_periodics_trend <- function(detrend_series,
                                model_obj,
                                periodic_indices,
                                L,
                                method,
                                basis=NULL,
                                sign_level=0.05,
                                tau_threshold=0.05,
                                p_0=0.03,
                                stop_flag=TRUE,
                                trace=FALSE,
                                noise_type="white",
                                auto_trend_freq=NULL,
                                fixed_freqs=TRUE,
                                check_hypothesis=TRUE){
  # Метод автоматической идентификации периодик после выделения тренда
  # 
  # :param detrend_series: временной ряд, из которого извлекли тренд
  # :param model_obj: объект разложения
  # :param periodic_indices: индексы нетрендовых компонент в сигнале
  # :param L: длина окна
  # :param method: критерий для проверки гипотезы о шуме
  # :param trend_method: метод автоматической идентификации тренда, 
  #                      поддерживаются 'eossa.auto', 'iossa.auto_grouping.fossa'
  #                      'iossa.auto_grouping.ssa', 'fossa.auto', 'basic.auto',
  #                      'iossa.post_grouping'
  # :param basis: базис для критерия MC-SSA, поддерживается 'sin0', 'sin1'
  # :param sign_level: уровень значимости для проверки гипотез о шуме
  #                    в алгоритме автоматической идентификации гармоник
  # :param tau_threshold: порог для метода регулярности углов, гармоники
  #                       с мерой tau меньше этого порога отбираются в сигнал
  # :param auto_trend_freq: порог низких частот, частоты в разложении Фурье,
  #                         которые меньше этой частоты, считаем трендовыми
  # :param auto_threshold: порог для доли низких частот, компоненты с долей 
  #                        низких частот больше этого порога считаем трендовыми
  # :param trace: флаг для подробного вывода всех этапов вычислений
  # :param p_0: порог для метода идентификации периодик с периодом 2 (порог для 
  #             доли пар без смены знака), если доля пар без смены знака меньше 
  #             этого порога, то компоненту относим к сигналу
  # :param stop_flag: флаг остановки проверки гипотез при не отвержении
  # :param noise_type: тип шума в гипотезе (внутри алгоритма идентификации 
  #                    периодик), 'white' или 'red'
  # :param fixed_freqs: флаг для оценки параметров шума по части спектра
  # :param check_hypothesis: флаг для проверки гипотез
  # 
  # :return: list с ключами
  #   :key raw_two_el_indices: Индексы парных компонент, отобранных методом
  #                            регулярности углов
  #   :key true_two_el_indices: Индексы периодик с частотой < 0.5, отобранных в
  #                             сигнал после проверки гипотез
  #   :key raw_one_el_indices: Индексы парных компонент, отобранных методом
  #                            чередования знака
  #   :key true_one_el_indices: Индексы периодик с частотой = 0.5, отобранных в
  #                             сигнал после проверки гипотез
  #   :key result_df: Таблица, в которой хранятся меры и pval
  #   :key periodics: периодическая компонента
  
  result_df <- data.frame(matrix(nrow=0, ncol=0))
  result_df[1, "periodic_index"] <- 0
  result_df[1, "measure"] <- -1
  result_df[1, "measure_type"] <- -1
  result_df[1, "mss"] <- -1
  result_df[1, "period"] <- 0
  result_df[1, "index1"] <- 0
  result_df[1, "index2"] <- 0
  
  
  noise_estimate <- detrend_series
  true_two_el_indices <- numeric(0)
  true_one_el_indices <- numeric(0)
  
  if (check_hypothesis){
    res <- noise_hypothesis(time_series=noise_estimate,
                            method=method,
                            sign_level=sign_level,
                            noise_type=noise_type,
                            L=L,
                            fixed_freqs=fixed_freqs,
                            freq.range=c(auto_trend_freq, 0.5))
    
    result_df[1, "pval"] <- res$pval
    
    end_flag <- res$result 
  }
  
  j <- 1
  
  model <- auto_periodic_model(ssa_obj=model_obj,
                               groups=periodic_indices,
                               method="angle_reg",
                               threshold=tau_threshold,
                               p_0=p_0)
  
  indices <- c(model$two_el_indices, model$one_el_indices)
  raw_two_el_indices <- model$two_el_indices
  raw_one_el_indices <- model$one_el_indices
  measure <- c(rep(model$two_el_tau, each=2), model$one_el_gamma)
  measure_type <- c(rep(1, length(raw_two_el_indices)), rep(2, length(raw_one_el_indices)))
  
  if (!check_hypothesis){
    
    all_indices <- c(raw_two_el_indices, raw_one_el_indices)
    periodic_estimate <- reconstruct(model_obj, groups=list(all_indices))$F1

    return (list("raw_two_el_indices"=raw_two_el_indices,
                 "true_two_el_indices" = NULL,
                 "raw_one_el_indices" = raw_one_el_indices,
                 "true_one_el_indices" = NULL,
                 "result_df"=NULL,
                 "periodics"=periodic_estimate))
  }
  
  if (trace){
    print("All periodic indices: ")
    print(indices)
    
    print("All identified two el periodics: ")
    print(raw_two_el_indices)
    
    print("All identified one el periodics: ")
    print(raw_one_el_indices)
    
    print("All measure values: ")
    print(measure)
    
    print("Measure type: ")
    print(measure_type) 
  }
  
  i <- 1
  mss <- numeric(0)
  while (i < length(raw_two_el_indices)){
    mss <- c(mss, mean(reconstruct(model_obj, groups=list(c(raw_two_el_indices[i], raw_two_el_indices[i + 1])))$F1 ^ 2))
    i <- i + 2
  }
  mss <- rep(mss, each=2)
  
  i <- 1
  while(i <= length(raw_one_el_indices)){
    mss <- c(mss, mean(reconstruct(model_obj, groups=list(c(raw_one_el_indices[i])))$F1 ^ 2))
    i <- i + 1
  }
  
  if (trace){
    print("MSS: ")
    print(mss)
    
    print("Start filtering...")
  }
  
  combined_matrix <- cbind(indices, mss, measure, measure_type)
  
  sorted_matrix <- combined_matrix[order(combined_matrix[, 2], decreasing = TRUE), ]
  
  if (is.null(nrow(sorted_matrix))){
    indices <- sorted_matrix[1]
    mss <- sorted_matrix[2]
    measure <- sorted_matrix[3]
    measure_type <- sorted_matrix[4]
  } else {
    indices <- sorted_matrix[, 1]
    mss <- sorted_matrix[, 2]
    measure <- sorted_matrix[, 3]
    measure_type <- sorted_matrix[, 4]
  }

  i <- 1
  j <- 1
  
  periodic_estimate <- rep(0, length(detrend_series))
  
  while (i <= length(indices)){
    if (trace & (!end_flag)){
      if (measure_type[i] == 1)
        print(paste0("Indices: ", i, ", ", i + 1))
      else
        print(paste0("Index: ", i))
      if (measure_type[i] == 1)
        print(paste0("Decomposition indices: ", indices[i], ", ", indices[i + 1]))
      else
        print(paste0("Decomposition index: ", indices[i]))
      print(paste0("Measure: ", measure_type[i], ", value: ", measure[i], ", MSS: ", mss[i], ", periodic index: ", j))
    }
    if (!end_flag){
      if (measure_type[i] == 1) {
        true_two_el_indices <- c(true_two_el_indices, indices[i], indices[i + 1])
      } else {
        true_one_el_indices <- c(true_one_el_indices, indices[i])
      }
    }
    else{
      if (stop_flag)
        break
    }
    
    result_df[j + 1, "periodic_index"] <- j
    result_df[j + 1, "measure"] <- measure[i]
    result_df[j + 1, "measure_type"] <- measure_type[i]
    result_df[j + 1, "mss"] <- mss[i]
    
    if (measure_type[i] == 1) {
      periodic <- reconstruct(model_obj, groups=list(c(indices[i],
                                                       indices[i + 1])))$F1
      result_df[j + 1, "period"] <- parestimate(model_obj, groups=list(c(indices[i],
                                                                         indices[i + 1])))$period[1]
      result_df[j + 1, "index1"] <- indices[i]
      result_df[j + 1, "index2"] <- indices[i + 1]
      i <- i + 2
    } else {
      periodic <- reconstruct(model_obj, groups=list(c(indices[i])))$F1
      result_df[j + 1, "period"] <- 2
      result_df[j + 1, "index1"] <- indices[i]
      i <- i + 1
    }
    
    periodic_estimate <- periodic_estimate + periodic
    noise_estimate <- noise_estimate - periodic
    
    res <- noise_hypothesis(time_series=noise_estimate,
                            method=method,
                            sign_level=sign_level,
                            noise_type=noise_type,
                            L=L,
                            fixed_freqs=fixed_freqs,
                            freq.range=c(auto_trend_freq, 0.5))
    
    result_df[j + 1, "pval"] <- res$pval
    
    end_flag <- res$result
    j <- j + 1
  }
  
  if ((length(true_two_el_indices) > 0) & trace)
    print(parestimate(model_obj, groups=list(true_two_el_indices)))
  
  return (list("raw_two_el_indices"=raw_two_el_indices,
               "true_two_el_indices" = true_two_el_indices,
               "raw_one_el_indices" = raw_one_el_indices,
               "true_one_el_indices" = true_one_el_indices,
               "result_df"=result_df,
               "periodics"=periodic_estimate))
}

procedure <- function(
    time_series,
    L,
    signal_rank,
    method='mc',
    trend_method="eossa.auto",
    basis='sin0',
    sign_level=0.05,
    tau_threshold=0.05,
    auto_trend_freq=1/24,
    auto_threshold=0.5,
    p_0=0.03,
    base='series',
    trace=FALSE,
    noise_type="white",
    fixed_freqs=TRUE,
    check_hypothesis=TRUE,
    rank_method='trmat',
    max_rank=30,
    eossa_delta=1e-3,
    ...
  ){
  # Метод автоматической идентификации компонент
  # :param time_series: временной ряд, по которому будет построено разложение
  # :param L: длина окна
  # :param signal_rank: ранг сигнала, целое число, либо 'auto'
  # :param method: критерий для проверки гипотезы о шуме
  # :param trend_method: метод автоматической идентификации тренда, 
  #                      поддерживаются 'eossa.auto', 'iossa.auto_grouping.fossa'
  #                      'iossa.auto_grouping.ssa', 'fossa.auto', 'basic.auto',
  #                      'iossa.post_grouping'
  # :param basis: базис для критерия MC-SSA, поддерживается 'sin0', 'sin1'
  # :param sign_level: уровень значимости для проверки гипотез о шуме
  #                    в алгоритме автоматической идентификации гармоник
  # :param tau_threshold: порог для метода регулярности углов, гармоники
  #                       с мерой tau меньше этого порога отбираются в сигнал
  # :param auto_trend_freq: порог низких частот, частоты в разложении Фурье,
  #                         которые меньше этой частоты, считаем трендовыми
  # :param auto_threshold: порог для доли низких частот, компоненты с долей 
  #                        низких частот больше этого порога считаем трендовыми
  # :param base: тип ряда для автоматической группировки, 'series' или 'eigen'
  # :param trace: флаг для подробного вывода всех этапов вычислений
  # :param p_0: порог для метода идентификации периодик с периодом 2 (порог для 
  #             доли пар без смены знака), если доля пар без смены знака меньше 
  #             этого порога, то компоненту относим к сигналу
  # :param noise_type: тип шума в гипотезе (внутри алгоритма идентификации 
  #                    периодик), 'white' или 'red'
  # :param fixed_freqs: флаг для оценки параметров шума по части спектра
  # :param check_hypothesis: флаг для проверки гипотез
  # :param rank_method: метод автоматической оценки ранга сигнала, только для 
  #                     signal_rank=='auto'; значения:
  #                     'svd', 'bic', 'trmat', 'mdl' - для белого и красного шума;
  #                     'trmat_fast', 'mgn', 'samos', 'ester' - для белого шума
  # :param max_rank: максимальная оценка ранга сигнала, только для signal_rank=='auto'
  # 
  # :return: list с ключами
  #   :key residuals: оценка шума
  #   :key trend: оценка тренда
  #   :key periodics: оценка суммы периодических компонент
  #   :key signal_rank: оценка ранга сигнала
  #   :key true_two_el_indices: индексы истинных двуэлементных периодических 
  #        компонент
  #   :key true_one_el_indices: индексы истинных периодик с периодом 2
  #   :key raw_two_el_indices: индексы всех двуэлементных периодических
  #        компонент, отобранных методом регулярности углов
  #   :key raw_one_el_indices: индексы всех периодик с периодом 2,
  #        отобранных методом регулярности углов
  #   :key trend_indices: индексы трендовых компонент
  #   :key model_obj: объект SSA (вложенное разложение)
  
  if (signal_rank == 'auto'){
    if (noise_type == 'red'){
      if (rank_method %in% c('svd', 'bic', 'trmat', 'mdl'))
        rank_method <- paste0(rank_method, '_ar1')
      else
        stop(paste0('Unknown method for ar1 noise: ', rank_method))
    }

    signal_rank <- estimate_rank(time_series=time_series,
                                 L=L,
                                 rank_range=0:max_rank,
                                 method=rank_method)
  }
  
  s <- ssa(time_series, L=L, svd.method="svd")
  trend_model <- auto_trend_model(ssa_obj=s,
                                  method=trend_method,
                                  signal_rank=signal_rank,
                                  clust_type="distance",
                                  auto_trend_freq=auto_trend_freq,
                                  auto_threshold=auto_threshold,
                                  delta=eossa_delta,
                                  base=base,
                                  ...)
  
  if (trace){
    plot(trend_model$model_obj)
    plot(trend_model$model_obj, type="paired")
    print("Trend indices")
    print(trend_model$grouping$Trend)
  }
  
  detrend_series <- trend_model$periodics + trend_model$residuals
  
  res <- get_periodics_trend(detrend_series,
                             model_obj=trend_model$model_obj,
                             periodic_indices=trend_model$grouping$Periodics,
                             L=L,
                             noise_type=noise_type,
                             method=method,
                             basis=basis,
                             sign_level=sign_level,
                             tau_threshold=tau_threshold,
                             p_0=0.03,
                             trace=trace,
                             fixed_freqs=fixed_freqs,
                             auto_trend_freq=auto_trend_freq,
                             check_hypothesis=check_hypothesis)
  
  res[["residuals"]] <- time_series - trend_model$trend - res$periodics
  res[["trend"]] <- trend_model$trend
  res[["trend_indices"]] <- trend_model$grouping$Trend
  res[["model_obj"]] <- trend_model$model_obj
  res[['signal_rank']] <- length(trend_model$grouping$Trend) + 
    length(res$true_one_el_indices) + length(res$true_two_el_indices)
  
  return(res)
}

auto_ssa <- function(
  time_series,
  L,
  signal_rank,
  method='mc',
  trend_method="eossa.auto",
  basis='sin0',
  sign_level=0.001,
  tau_threshold=0.05,
  auto_trend_freq=1/24,
  auto_threshold=0.5,
  p_0=0.03,
  base='series',
  trace=FALSE,
  noise_type="white",
  fixed_freqs=TRUE,
  check_hypothesis=TRUE,
  rank_method='trmat',
  max_rank=30,
  eossa_delta=1e-4,
  nstages='auto',
  ...){
  # Метод автоматической идентификации компонент
  # :param time_series: временной ряд, по которому будет построено разложение
  # :param L: длина окна
  # :param signal_rank: ранг сигнала, целое число, либо 'auto'
  # :param method: критерий для проверки гипотезы о шуме
  # :param trend_method: метод автоматической идентификации тренда, 
  #                      поддерживаются 'eossa.auto', 'iossa.auto_grouping.fossa'
  #                      'iossa.auto_grouping.ssa', 'fossa.auto', 'basic.auto',
  #                      'iossa.post_grouping'
  # :param basis: базис для критерия MC-SSA, поддерживается 'sin0', 'sin1'
  # :param sign_level: уровень значимости для проверки гипотез о шуме
  #                    в алгоритме автоматической идентификации гармоник
  # :param tau_threshold: порог для метода регулярности углов, гармоники
  #                       с мерой tau меньше этого порога отбираются в сигнал
  # :param auto_trend_freq: порог низких частот, частоты в разложении Фурье,
  #                         которые меньше этой частоты, считаем трендовыми
  # :param auto_threshold: порог для доли низких частот, компоненты с долей 
  #                        низких частот больше этого порога считаем трендовыми
  # :param base: тип ряда для автоматической группировки, 'series' или 'eigen'
  # :param trace: флаг для подробного вывода всех этапов вычислений
  # :param p_0: порог для метода идентификации периодик с периодом 2 (порог для 
  #             доли пар без смены знака), если доля пар без смены знака меньше 
  #             этого порога, то компоненту относим к сигналу
  # :param noise_type: тип шума в гипотезе (внутри алгоритма идентификации 
  #                    периодик), 'white' или 'red'
  # :param fixed_freqs: флаг для оценки параметров шума по части спектра
  # :param check_hypothesis: флаг для проверки гипотез
  # :param rank_method: метод автоматической оценки ранга сигнала, только для 
  #                     signal_rank=='auto'; значения:
  #                     'svd', 'bic', 'trmat', 'mdl' - для белого и красного шума;
  #                     'trmat_fast', 'mgn', 'samos', 'ester' - для белого шума
  # :param max_rank: максимальная оценка ранга сигнала, только для signal_rank=='auto'
  # :param nstages: количество этапов алгоритма. 'auto' или int. nstages > 1 можно использовать для
  #                 более точной оценки ранга сигнала
  # 
  # :return: list с ключами
  #   :key residuals: оценка шума
  #   :key trend: оценка тренда
  #   :key periodics: оценка суммы периодических компонент
  #   :key signal_rank: оценка ранга сигнала
  #   :key true_two_el_indices: индексы истинных двуэлементных периодических 
  #        компонент
  #   :key true_one_el_indices: индексы истинных периодик с периодом 2
  #   :key raw_two_el_indices: индексы всех двуэлементных периодических
  #        компонент, отобранных методом регулярности углов
  #   :key raw_one_el_indices: индексы всех периодик с периодом 2,
  #        отобранных методом регулярности углов
  #   :key trend_indices: индексы трендовых компонент
  #   :key model_obj: объект SSA (вложенное разложение)
  
  if (nstages == 'auto')
    nstages_max <- 30
  else
    nstages_max <- nstages
  current_signal_rank <- signal_rank
  for (stage in 1:nstages_max){
    dec <- procedure(
      time_series=time_series,
      L=L,
      signal_rank=current_signal_rank,
      method='mc',
      trend_method="eossa.auto",
      basis='sin0',
      sign_level=sign_level,
      tau_threshold=tau_threshold,
      auto_trend_freq=auto_trend_freq,
      auto_threshold=auto_threshold,
      p_0=p_0,
      base=base,
      trace=trace,
      noise_type=noise_type,
      fixed_freqs=fixed_freqs,
      check_hypothesis = check_hypothesis,
      rank_method=rank_method,
      max_rank=max_rank,
      eossa_delta=eossa_delta
    )
    if (((nstages == 'auto') && (current_signal_rank == dec$signal_rank)))
      break
    current_signal_rank <- dec$signal_rank
  }
  return(dec)
}