library("lattice")
library("svd")
library("forecast")
library(ggplot2)
source("eossa_new.R", local = TRUE)

iossa_post_grouping <- function(ssa_obj,
                                signal_rank,
                                auto_trend_freq = 1/24,
                                auto_threshold = 0.5,
                                iossa_maxiter = 10,
                                iossa_tolerance = 1e-3,
                                kappa = 2){
  ioss <- iossa(ssa_obj, nested.groups = as.list(c(1:signal_rank)), tol = iossa_tolerance, 
                maxiter = iossa_maxiter, kappa = kappa)
  g <- grouping.auto(ioss, base = "series", groups = list(c(1:signal_rank)),
                     freq.bins = list(Trend = auto_trend_freq), 
                     threshold = auto_threshold)
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  rec <- reconstruct(ioss, groups = list(Trend = trend_comp, Periodics = res_comp))
  
  return(list(model_object = ioss, 
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              iterations = ioss$iossa.result$iter,
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

iossa_manual <- function(ssa_obj,
                         signal_rank,
                         auto_trend_freq = 1 / 24,
                         auto_threshold = 0.5,
                         iossa_maxiter = 10,
                         iossa_tolerance = 1e-3,
                         manual_grouping,
                         kappa = 2){
  ioss <- iossa(ssa_obj, nested.groups = manual_grouping, tol = iossa_tolerance,
                maxiter = iossa_maxiter, kappa = kappa)
  g <- grouping.auto(ioss, base = "series", groups = manual_grouping,
                     freq.bins = list(Trend = auto_trend_freq), 
                     threshold = auto_threshold)
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  rec <- reconstruct(ioss, groups = list(Trend = trend_comp, Periodics = res_comp))
  return(list(model_object = ioss,
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              iterations = ioss$iossa.result$iter,
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

iossa_auto_grouping_ssa <- function(ssa_obj,
                                    signal_rank,
                                    auto_trend_freq = 1 / 24,
                                    auto_threshold = 0.5,
                                    auto_initial_threshold = 0.2,
                                    iossa_tolerance = 1e-3,
                                    iossa_maxiter = 10,
                                    kappa = 2){
  auto_grouping <- grouping.auto(ssa_obj, base = "series", groups = list(c(1:signal_rank)),
                                 freq.bins = list(Trend = auto_trend_freq), 
                                 threshold = auto_initial_threshold)
  
  trend_comp <- auto_grouping$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  result <- list()
  
  if ((length(trend_comp) > 0) & (length(res_comp) > 0)) {
    ioss <- iossa(ssa_obj, nested.groups = list(trend_comp, res_comp), tol = iossa_tolerance,
                  maxiter = iossa_maxiter, kappa = kappa)
    g <- grouping.auto(ioss, base = "series", groups = list(c(1:signal_rank)),
                       freq.bins = list(Trend = auto_trend_freq), 
                       threshold = auto_threshold)
    
    trend_comp <- g$Trend
    signal_indices <- 1:signal_rank
    res_comp <- signal_indices[!signal_indices %in% trend_comp]
    
    result[["model_object"]] <- ioss
    result[["iterations"]] <- ioss$iossa.result$iter
  } else {
    result[["model_object"]] <- ssa_obj
    result[["iterations"]] <- 1
  }
  
  rec <- reconstruct(result[["model_object"]], groups = list(Trend = trend_comp, Periodics = res_comp))
  
  result[["grouping"]] <- list(Trend = trend_comp, Periodics = res_comp)
  result[["trend"]] <- rec$Trend
  result[["periodics"]] <- rec$Periodics
  result[["residuals"]] <- attr(rec, "residuals")
  
  return(result)
}

iossa_auto_grouping_fossa <- function(ssa_obj,
                                      signal_rank,
                                      auto_trend_freq = 1 / 24,
                                      auto_threshold = 0.5,
                                      auto_initial_threshold = 0.2,
                                      iossa_tolerance = 1e-3,
                                      iossa_maxiter = 1,
                                      kappa = kappa){
  foss <- fossa(ssa_obj, nested.groups = 1:signal_rank)
  auto_grouping <- grouping.auto(foss, base = "series", groups = list(c(1:signal_rank)),
                                 freq.bins = list(Trend = auto_trend_freq), 
                                 threshold = auto_initial_threshold)
  
  trend_comp <- auto_grouping$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  result <- list()
  
  if ((length(trend_comp) > 0) & (length(res_comp) > 0)) {
    ioss <- iossa(foss, nested.groups = list(trend_comp, res_comp), tol = iossa_tolerance,
                  maxiter = iossa_maxiter, kappa = kappa)
    g <- grouping.auto(ioss, base = "series", groups = list(c(1:signal_rank)),
                       freq.bins = list(Trend = auto_trend_freq), 
                       threshold = auto_threshold)
    
    trend_comp <- g$Trend
    signal_indices <- 1:signal_rank
    res_comp <- signal_indices[!signal_indices %in% trend_comp]
    
    result[["model_object"]] <- ioss
    result[["iterations"]] <- ioss$iossa.result$iter
  } else {
    result[["model_object"]] <- foss
    result[["iterations"]] <- 1
  }
  
  rec <- reconstruct(result[["model_object"]], groups = list(Trend = trend_comp, Periodics = res_comp))
  
  result[["grouping"]] <- list(Trend = trend_comp, Periodics = res_comp)
  result[["trend"]] <- rec$Trend
  result[["periodics"]] <- rec$Periodics
  result[["residuals"]] <- attr(rec, "residuals")
  
  return(result)
}

eossa_auto <- function(ssa_obj,
                        signal_rank,
                        clust_type = c("distance", "frequency", "hierarchial"),
                        k = 2,
                        delta = 1e-3,
                        auto_trend_freq = 1 / 24,
                        auto_threshold = 0.5){
  eoss <- eossa_new(ssa_obj, 
                    nested.groups = list(1:signal_rank), 
                    clust_type = clust_type, 
                    k = k, 
                    delta = delta,
                    auto_trend_freq = auto_trend_freq)
  
  g <- grouping.auto(eoss, base = "series", groups = eoss$iossa.groups,
                     freq.bins = list(Trend = auto_trend_freq), 
                     threshold = auto_threshold)
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  rec <- reconstruct(eoss, groups = list(Trend = trend_comp, Periodics = res_comp))
  return(list(model_object = eoss,
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

basic_auto <- function(ssa_obj,
                       signal_rank,
                       auto_trend_freq = 1 / 24,
                       auto_threshold = 0.5){
  g <- grouping.auto(ssa_obj, base = "series", groups = list(1:signal_rank),
                     freq.bins = list(Trend = auto_trend_freq), 
                     threshold = auto_threshold)
  
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  rec <- reconstruct(ssa_obj, groups = list(Trend = trend_comp, Periodics = res_comp))
  
  return(list(model_object = ssa_obj,
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

fossa_auto <- function(ssa_obj,
                       signal_rank,
                       auto_trend_freq = 1 / 24,
                       auto_threshold = 0.5){
  foss <- fossa(ssa_obj, nested.groups = 1:signal_rank)
  g <- grouping.auto(foss, base = "series", groups = list(1:signal_rank),
                     freq.bins = list(Trend = auto_trend_freq), 
                     threshold = auto_threshold)
  
  
  trend_comp <- g$Trend
  signal_indices <- 1:signal_rank
  res_comp <- signal_indices[!signal_indices %in% trend_comp]
  
  rec <- reconstruct(foss, groups = list(Trend = trend_comp, Periodics = res_comp))
  return(list(model_object = foss,
              grouping = list(Trend = trend_comp, Periodics = res_comp),
              trend = rec$Trend,
              periodics = rec$Periodics,
              residuals = attr(rec, "residuals")))
}

auto_trend_model <- function(ssa_obj,
                             method = c("iossa.post_grouping", "iossa.manual", "iossa.auto_grouping.ssa",
                                        "iossa.auto_grouping.fossa", "eossa.auto", "basic.auto", "fossa.auto"),
                             signal_rank,
                             auto_trend_freq = 1 / 24,
                             auto_threshold = 0.5,
                             auto_initial_threshold = 0.2,
                             iossa_maxiter = 10,
                             iossa_tolerance = 1e-3,
                             kappa = 2,
                             manual_grouping = list(),
                             eossa_delta = 1e-3,
                             eossa_clustering = "distance"){
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
    res <- iossa_post_grouping(ssa_obj = ssa_obj, 
                               signal_rank = signal_rank, 
                               auto_trend_freq = auto_trend_freq, 
                               auto_threshold = auto_threshold,
                               iossa_maxiter = iossa_maxiter, 
                               iossa_tolerance = iossa_tolerance,
                               kappa = kappa)
    results[["iossa.post_grouping"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else {
      this[["iossa.post_grouping"]]$iterations = res$iterations
    }
  }
  if ("iossa.manual" %in% method){
    res <- iossa_manual(ssa_obj = ssa_obj, 
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
    res <- iossa_auto_grouping_ssa(ssa_obj = ssa_obj, 
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
    res <- iossa_auto_grouping_fossa(ssa_obj = ssa_obj, 
                                     signal_rank = signal_rank, 
                                     auto_trend_freq = auto_trend_freq, 
                                     auto_threshold = auto_threshold,
                                     auto_initial_threshold = auto_initial_threshold, 
                                     iossa_tolerance = iossa_tolerance,
                                     kappa = kappa)
    results[["iossa.auto_grouping.fossa"]] <- res
    if (length(method) == 1)
      this$iterations <- res$iterations
    else{
      this[["iossa.auto_grouping.fossa"]]$iterations = res$iterations
    }
  }
  if ("eossa.auto" %in% method){
    res <- eossa_auto(ssa_obj, 
                        signal_rank = signal_rank, 
                        clust_type = eossa_clustering,
                        auto_trend_freq = auto_trend_freq, 
                        auto_threshold = auto_threshold,
                        delta = eossa_delta)
    results[["eossa.auto"]] <- res
  }
  if ("basic.auto" %in% method){
    res <- basic_auto(ssa_obj = ssa_obj, 
                      signal_rank = signal_rank, 
                      auto_trend_freq = auto_trend_freq, 
                      auto_threshold = auto_threshold)
    results[["basic.auto"]] <- res
  }
  if ("fossa.auto" %in% method){
    res <- fossa_auto(ssa_obj, 
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

get_weight <- function(index, N, window_length){
  # Input: index          --- index of series element
  #        N              --- series length
  #        window_length  --- length of window in SSA
  # Result: returns weight in weighted correlation
  K <- N - window_length + 1
  L_star <- min(window_length, K)
  K_star <- max(window_length, K)
  return(ifelse((index <= L_star), index, ifelse(index <= K_star, L_star, N - index + 1)))
}

wcor_series <- function(comp1, comp2, window_length) {
  # Input: comp1         --- series
  #        comp2         --- series
  #        window_length --- length of window in SSA
  if (length(comp1) != length(comp2))
    stop("Lengths of series are not equal.")
  N <- length(comp1)
  wcor_sum <- 0
  fr_norm1 <- 0
  fr_norm2 <- 0
  for (i in 1:N){
    wcor_sum <- wcor_sum + get_weight(i, N, window_length) * comp1[i] * comp2[i]
    fr_norm1 <- fr_norm1 + get_weight(i, N, window_length) * comp1[i] * comp1[i]
    fr_norm2 <- fr_norm2 + get_weight(i, N, window_length) * comp2[i] * comp2[i]
  }
  fr_norm1 <- sqrt(fr_norm1)
  fr_norm2 <- sqrt(fr_norm2)
  return(wcor_sum / fr_norm1 / fr_norm2)
}

plot.auto_trend <- function(x, type = c("series", "vectors", "values", "elementary"),
                            filename = NA,
                            file_width = 1000,
                            file_height = 600,
                            size = 0.7,
                            compare = FALSE,
                            actual_trend = numeric(0),
                            compare_only = FALSE,
                            plot.contrib = FALSE,
                            original_series = numeric(0),
                            color_list,
                            methods_map,
                            indices = NA
){
  if (any(is.na(indices)))
    indices <- 1:x$signal_rank
  
  type = match.arg(type)
  if (type == "series"){
    if (compare){
      if (length(actual_trend) == 0)
        stop("Actual trend must be passed to plot.")
      if (length(x$method) == 1){
        if (!is.na(filename))
          pdf(file = filename, width = file_width, height = file_height)
        
        plot_df <- data.frame("Reconstructed trend" = x$trend, 
                              "Index" = 0:(length(x$original)-1), 
                              "Actual trend" = actual_trend)
        
        if (length(original_series) > 0)
          plot_df[,"Original series"] <- original_series
        
        p <- ggplot(plot_df, aes(Index)) + 
          geom_line(aes(y = Actual.trend, colour = "Actual trend"), size = size) + 
          geom_line(aes(y = Reconstructed.trend, colour = "Reconstructed trend"), size = size)
        
        if (length(original_series) > 0)
          p <- p  + geom_line(aes(y = Actual.trend, colour = "Original series"), size = size)
        
        p <- p + 
          labs(x = "Index", y = "Value", color = "Trends") + 
          scale_color_manual(values=c("#07973D", "#EE3D3D"))
        
        print(p + theme_classic())
        
        if (!is.na(filename))
          dev.off()
      }
      else{
        index <- 1
        if (!compare_only){
          for (m in x$method){
            if (!is.na(filename)){
              if (length(c(filename)) == 1)
                png(filename = paste0(filename, '.', m, '.png'), width = file_width, height = file_height)
              else
                png(filename = filename[index], width = file_width, height = file_height)
            }
            
            plot_df <- data.frame("Reconstructed trend" = x[[m]]$trend, 
                                  "Index" = 0:(length(x$original)-1), 
                                  "Actual trend" = actual_trend)
            
            p <- ggplot(plot_df, aes(Index)) + 
              geom_line(aes(y = Actual.trend, colour = "Actual trend"), size = size) + 
              geom_line(aes(y = Reconstructed.trend, colour = "Reconstructed trend"), size = size) + 
              labs(x = "Index", y = "Value", color = "Series") + 
              scale_color_manual(values=c("#EE3D3D", "#0805B5"))
            
            print(p)
            
            if (!is.na(filename))
              dev.off()
            
            index <- index + 1
          }
        }
        
        if (!is.na(filename)){
          if (length(c(filename)) == 1)
            pdf(file = paste0(filename, '.', 'compare', '.pdf'), width = file_width, height = file_height)
          else
            pdf(file = filename[index], width = file_width, height = file_height)
        }
        
        plot_df <- data.frame("Index" = 0:(length(x$original)-1), 
                              "Actual trend" = actual_trend)
        
        if (length(original_series) > 0)
          plot_df$Original.series <- original_series
        for (m in x$method){
          plot_df[m] <- x[[m]]$trend
        }
        
        index <- 1
        p <- ggplot(plot_df, aes(Index)) + 
          geom_line(aes_string(y = "Actual.trend", colour = as.factor(index)), size = size)
        index <- 2
        colors <- c(color_list$actual.trend)
        
        if (length(original_series) > 0) {
          p <- p + geom_line(aes_string(y = "Original.series", colour = as.factor(index)), size = size)
          index <- 3
        }
        labels <- c("Actual trend")
        for (m in x$method){
          p <- p + geom_line(aes_string(y = m, colour = as.factor(index)), size = size)
          colors <- c(colors, color_list[[m]])
          labels <- c(labels, methods_map[[m]])
          index <- index + 1
        }
        
        p <- p + scale_colour_manual(name = "Trends", labels = labels, values = colors) +
                     labs(title = "", x = "Index", y = "Value")
        
        print(p + theme_classic())
        
        if (!is.na(filename))
          dev.off()
      }
    }
    else{
      if (length(x$method) == 1){
        if (!is.na(filename))
          pdf(file = filename, width = file_width, height = file_height)
        
        plot_df <- data.frame("Original series" = x$original, "Reconstructed trend" = x$trend, 
                              "Reconstructed periodics" = x$periodics, "Index" = 0:(length(x$original)-1), 
                              "Reconstructed residuals" = x$residuals)
        
        p <- ggplot(plot_df, aes(Index)) + 
          geom_line(aes(y = Original.series, colour = "Original series"), size = size) + 
          geom_line(aes(y = Reconstructed.trend, colour = "Reconstructed trend"), size = size) + 
          geom_line(aes(y = Reconstructed.periodics, colour = "Reconstructed periodics"), size = size) + 
          geom_line(aes(y = Reconstructed.residuals, colour = "Reconstructed residuals"), size = size) + 
          labs(title = "", x = "Index", y = "Value", color = "Series") + 
          scale_color_manual(values=c(color_list[["original.series"]], 
                                      color_list[["reconstructed.periodics"]], 
                                      color_list[["reconstructed.residuals"]], 
                                      color_list[[x$method]]))
        
        print(p + theme_classic())
        
        if (!is.na(filename))
          dev.off()
      }
      else{
        index <- 1
        if (!compare_only){
          for (m in x$method){
            if (!is.na(filename)){
              if (length(filename) == 1)
                png(filename = paste0(filename, '.', m, '.png'), width = file_width, height = file_height)
              else
                png(filename = filename[index], width = file_width, height = file_height)
            }
            
            plot_df <- data.frame("Original series" = x$original, "Trend" = x[[m]]$trend, 
                                  "Periodics" = x[[m]]$periodics, "Index" = 0:(length(x$original)-1), 
                                  "Residuals" = x[[m]]$residuals)
            
            p <- ggplot(plot_df, aes(Index)) + 
              geom_line(aes(y = Original.series, colour = "Original series"), size = size) + 
              geom_line(aes(y = Trend, colour = "Trend"), size = size) + 
              geom_line(aes(y = Periodics, colour = "Periodics"), size = size) + 
              geom_line(aes(y = Residuals, colour = "Residuals"), size = size) + 
              labs(title = paste0(m, ' components'), x = "Index", y = "Value", color = "Series") + 
              scale_color_manual(values=c("#52FA5E", "#EE3D3D", "#645502", "#0805B5"))
            
            print(p)
            
            if (!is.na(filename))
              dev.off()
            
            index <- index + 1
          } 
        }
        
        if (!is.na(filename)){
          if (length(c(filename)) == 1)
            pdf(file = paste0(filename, '.', 'compare', '.pdf'), width = file_width, height = file_height)
          else
            pdf(file = filename[index], width = file_width, height = file_height)
        }
        
        plot_df <- data.frame("Index" = 0:(length(x$original)-1))
        
        if (length(original_series) > 0)
          plot_df$Original.series <- original_series
        
        for (m in x$method){
          plot_df[m] <- x[[m]]$trend
        }
        
        p <- ggplot(plot_df, aes(Index))
        
        index <- 0
        colors <- numeric(0)
        labels <- numeric(0)
        if (length(original_series) > 0) {
          p <- p + geom_line(aes_string(y = "Original.series", colour = as.factor(index)), size = size)
          labels <- c("Original series")
          colors <- c(color_list$original.series)
        }
        
        index <- 1
        for (m in x$method){
          p <- p + geom_line(aes_string(y = m, colour = as.factor(index)), size = size)
          colors <- c(colors, color_list[[m]])
          labels <- c(labels, methods_map[[m]])
          index <- index + 1
        }
        
        
        p <- p + scale_colour_manual(labels = labels, values = colors) +
                     labs(x = "Index", y = "Value", color = "")
        
        print(p + theme_classic())
        
        if (!is.na(filename))
          dev.off()
      }
    }
  }
  if (type == "vectors"){
    if (length(x$method) == 1){
      if (!is.na(filename))
        png(filename = filename, width = file_width, height = file_height)
      
      print(plot(x$model_object, type = "vectors", idx = 1:x$signal_rank, main = paste0("Eigenvectors, ", x$method)))
      
      if (!is.na(filename))
        dev.off() 
    }
    else {
      index <- 1
      for (m in x$method){
        if (!is.na(filename)){
          if (length(filename) == 1)
            png(filename = paste0(filename,'.', m, '.png'), width = file_width, height = file_height)
          else
            png(filename = filename[index], width = file_width, height = file_height)
        }
        
        print(plot(x[[m]]$model_object, type = "vectors", idx = 1:x$signal_rank, main = paste0("Eigenvectors, ", m)))
        
        if (!is.na(filename))
          dev.off() 
        index <- index + 1
      }
    }
  }
  if (type == "values"){
    if (length(x$method) == 1){
      if (!is.na(filename))
        png(filename = filename, width = file_width, height = file_height)
      
      print(plot(x$model_object, type = "values", main = paste0("Component norms, ", x$method)))
      
      if (!is.na(filename))
        dev.off() 
    }
    else {
      index <- 1
      for (m in x$method){
        if (!is.na(filename)){
          if (length(filename) == 1)
            png(filename = paste0(filename,'.', m, '.png'), width = file_width, height = file_height)
          else
            png(filename = filename[index], width = file_width, height = file_height)
        }
        
        print(plot(x[[m]]$model_object, type = "values", main = paste0("Component norms, ", m)))
        
        if (!is.na(filename))
          dev.off() 
        index <- index + 1
      }
    }
  }
  if (type == "elementary"){
    if (length(x$method) == 1){
      if (!is.na(filename))
        pdf(file = filename, width = file_width, height = file_height)
      
      print(plot(x$model_object, type = "series", groups = indices, 
                 main = "",
                 plot.contrib = plot.contrib))
      
      if (!is.na(filename))
        dev.off() 
    }
    else {
      index <- 1
      for (m in x$method){
        if (!is.na(filename)){
          if (length(filename) == 1)
            png(filename = paste0(filename,'.', m, '.png'), width = file_width, height = file_height)
          else
            png(filename = filename[index], width = file_width, height = file_height)
        }
        
        print(plot(x[[m]]$model_object, type = "series", groups = 1:x$signal_rank, 
                   main = paste0("Reconstructed series, ", m),
                   plot.contrib = plot.contrib))
        
        if (!is.na(filename))
          dev.off() 
        index <- index + 1
      }
    }
    
  }
}