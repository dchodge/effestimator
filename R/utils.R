check_events <- function(events, KM_est) {

  N_events <- nrow(events)
  events$at_risk_before <- sum(events$n) - cumsum(append(0, events$n[1:(N_events - 1)]))

  t_data <- c()
  n_data <- c()
  j <- 0 # Index for KM_data, where "I" events take up two rows

  for (i in seq(1, N_events, 1)) {

    j <- j + 1

    if (events$Type[i] == "C") {
      t_data[i] <- KM_est[j, 1]
      n_data[i] <- NA
    } else if (events$Type[i] == "I") {
      t_data[i] <- mean(KM_est[j:(j + 1), 1])
      n_data[i] <- -diff(KM_est[j:(j + 1), 2])*events$at_risk_before[i]
      j <- j + 1
    }
  }

  check <- data.frame(t_data = t_data, t = events$t, dt = events$t - t_data, 
                      n_data = n_data, n = events$n, dn = events$n - n_data)

  return(check)
}

# Work out KM estimator from events data
KM_est <- function(events) {
  
  N_events <- nrow(events)
  events$at_risk_before <- sum(events$n) - cumsum(append(0, events$n[1:(N_events - 1)]))
   
  KM_est <- data.frame(t = 0, Proportion_Free = 1)
  j <- 0
  
  for (i in seq(1, N_events, 1)) {
    
    time <- events$t[i]
    j <- j + 1
    
    if (events$Type[i] == "C") {
      KM_est <- rbind(KM_est, c(time, KM_est[j, 2]))
    } else if (events$Type[i] == "I") {
      KM_est <- rbind(KM_est, c(time, KM_est[j, 2]))
      KM_est <- rbind(KM_est, c(time, KM_est[j, 2]*(1 - events$n[i]/events$at_risk_before[i])))
      j <- j + 1
    }
  }
  
  events$pt_elapsed <- cumsum(events$at_risk_before*diff(append(0, events$t)))
  
  data <- list(Events = events, KM_estimator = KM_est)
  
  return(data)
}


events_agg <- function(data, ts) {
  
  n_ts <- length(ts)
  
  events <- data$Events
  KM_est <- data$KM_estimator

  time_between_events <- diff(append(0, events$t))
  time_between_events[1] <- time_between_events[1] + 1
  events$pt_elapsed <- cumsum(events$at_risk_before*time_between_events) # up to and including day of event itself
  
  at_risk_f <- events$at_risk_before[1]
  
  output <- data.frame()
  
  for (i in seq(1, n_ts - 1, 1)) {
    
    t_min <- ts[i]
    t_max <- ts[i + 1]
    
    data <- subset(events, events$t >= t_min & events$t < t_max)

    if (nrow(data) == 0) {

      pt <- at_risk_f*(t_max - t_min)
      n <- 0
      
    } else {
      
      pt_elapsed_i_before_event <- head(data, 1)$pt_elapsed - head(data, 1)$at_risk_before # Second term = person-time on day of event
      at_risk_i <- head(data, 1)$at_risk_before
      t_before <- head(data, 1)$t - t_min
      pt_start <- pt_elapsed_i_before_event - at_risk_i*t_before
      
      pt_elapsed_f_after_event <- tail(data, 1)$pt_elapsed
      at_risk_f <- tail(data, 1)$at_risk_before - tail(data, 1)$n
      t_after <- t_max - 1 - tail(data, 1)$t
      pt_end <- pt_elapsed_f_after_event + at_risk_f*t_after
      
      pt <- pt_end - pt_start
    
      n <-  sum(data[which(data$Type == "I"), ]$n)
    
    }
    
    output <- rbind(output, c(paste("[", t_min, ", ", t_max, ")", sep = ""),
                              (t_min + t_max)/2, pt, n))
  }
  
  names(output) <- c("t_interval", "t", "person_time", "n")
  
  delta <- rep(1e-9, n_ts)
  delta[1] <- 0
  proportion_free <- approx(KM_est$t, KM_est$Proportion_Free, ts - delta, ties = "ordered")$y
  output$proportion_inf <- -diff(proportion_free)/proportion_free[1:(n_ts - 1)]
  
  # Why is this necessary?
  output$t <- as.numeric(output$t)
  output$person_time <- as.numeric(output$person_time)
  output$n <- as.numeric(output$n)
  output$proportion_inf <- as.numeric(output$proportion_inf)
  
  return(output)
}

get_mean_ci <- function(x) {
    data.frame(
        mean = mean(x),
        lb_95 = quantile(x, 0.025),
        ub_95 = quantile(x, 0.975)
    )
}
