
source(here::here("R", "utils.R") )
get_nmab_data <- function() {
    data_pla <- read.csv(here::here("data", "Placebo_pooled.csv"), header = F)
    var <- c("t", "Proportion_Free")
    names(data_pla) <- var
    data_vac <- read.csv(here::here("data", "Nirsevimab_pooled.csv"), header = F)
    names(data_vac) <- var

    # Assumed event data (to be checked below)
    events_pla <- data.frame(Type = c("C", "I", "I", "I", "I", "I", "I", "I", "C",
                                    "I", "I", "C", "I", "I", "I", "I", "I", "I", "I",
                                    "I", "I", "I", "I", "I", "I", "I", "I", "C", "I", "I", "I",
                                    "I", "I", "I", "C", "I", "I",
                                    "I", "I", "I", "C", "C", "C", "C"),
                            t = c(1, 7, 8, 15, 17, 19, 22, 23, 26,
                                30, 33, 33, 37, 48, 49, 50, 51, 53, 57,
                                63, 64, 66, 69, 73, 75, 77, 78, 79, 80, 82, 84,
                                91, 92, 111, 113, 115, 117,
                                121, 122, 132, 144, 149, 150, 151),
                            n = c(5, 1, 1, 1, 1, 2, 1, 1, 1,
                                1, 2, 1, 3, 4, 1, 1, 1, 1, 1,
                                2, 2, 1, 1, 1, 3, 1, 2, 1, 1, 3, 1,
                                1, 2, 1, 1, 2, 1,
                                1, 1, 1, 1, 1, 0, 724))
    events_vac <- data.frame(Type = c("C", "I", "I", "I",
                                    "C", "C", "C", "I", "I", "I",
                                    "C", "C", "C", "I", "C", "C", "I", "C",
                                    "I", "I", "C", "C", "I", "I", "I", "I", "I",
                                    "C", "I", "C", "I", "I", "C", "C", "C", "C"),
                            t = c(1, 7, 21, 28,
                                35, 38, 45, 51, 54, 57,
                                68, 71, 73, 74, 74, 75, 82, 85,
                                92, 93, 97, 99, 102, 103, 108, 109, 112, 123, 127, 127,
                                136, 138, 140, 143, 146, 151),
                            n = c(8, 1, 1, 1,
                                1, 1, 1, 2, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 2, 2, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1519))

    check_pla <- check_events(events_pla, data_pla)
    check_vac <- check_events(events_vac, data_vac)

    data_pla_clean <- KM_est(events_pla)
    data_vac_clean <- KM_est(events_vac)

    dt <- 1
    ts <- seq(0, 150, dt)
    agg_pla <- events_agg(data_pla_clean, ts) 
    agg_vac <- events_agg(data_vac_clean, ts) 
    list(placebo = agg_pla, vaccine = agg_vac)
}

# Maternal vaccination
get_matvac_data <- function() {
    dt <- 30
    ts <- seq(15, 165, dt)
    at_risk_pla <- c(3480, 3288, 2964, 2879, 2804, 2738, 2700)
    cumul_case_pla <- c(0, 15, 38, 56, 81, 99, 117)
    at_risk_vac <- c(3495, 3348, 3035, 2968, 2898, 2845, 2792)
    cumul_case_vac <- c(0, 2, 14, 24, 35, 47, 57)
    l <- length(at_risk_pla)
    agg_pla <- data.frame(
                        tp = seq_len(length(ts)),
                        t = ts,
                        person_time = (at_risk_pla[1:(l - 1)] + at_risk_pla[2:l])*dt/2,
                        n = diff(cumul_case_pla))
    agg_vac <- data.frame(
                        tp = seq_len(length(ts)),
                        t = ts,
                        person_time = (at_risk_vac[1:(l - 1)] + at_risk_vac[2:l])*dt/2,
                        n = diff(cumul_case_vac))
    list(placebo = agg_pla, vaccine = agg_vac)

}


get_oa_walsh_data <- function() {

    ts <- c(0, 35, 72, 104, 147, 176)
    at_risk_pla <- rep(17069, 6)
    cumul_case_pla <- c(0, 7, 12, 17, 27, 29)
    at_risk_vac <- rep(17215, 6)
    cumul_case_vac <- c(0, 3, 5, 8, 8, 10)
    l <- length(at_risk_pla)

    agg_pla <- data.frame(
                        tp = seq_len(length(ts) - 1),
                        t = ts[2:length(ts)],
                        person_time = (at_risk_pla[1:(l - 1)] + at_risk_pla[2:(l)]) * diff(ts) / 2,
                        n = diff(cumul_case_pla))
    agg_vac <- data.frame(t = ts[2:length(ts)],
                        tp = seq_len(length(ts) - 1),
                        person_time = (at_risk_vac[1:(l - 1)] + at_risk_vac[2:(l)]) * diff(ts) / 2,
                        n = diff(cumul_case_vac))
    list(placebo = agg_pla, vaccine = agg_vac)
}


# OA GSK
get_oa_papi_data <- function() {

    dt <- 30
    ts <- seq(15, 315, dt)
    at_risk_pla <- c(12494, 12390, 12268, 11853, 11597, 10973, 8255, 5441, 2697, 554, 2, 0)
    cumul_case_pla <- c(0, 22, 43, 62, 76, 86, 90, 95, 95, 95, 95, 95)
    at_risk_vac <- c(12466, 12390, 12282, 11881, 11641, 11029, 8305, 5481, 2717, 570, 2, 0)
    cumul_case_vac <- c(0, 3, 7, 15, 19, 23, 24, 26, 27, 27, 27, 27)
    l <- length(at_risk_pla)
    agg_pla <- data.frame(
                        tp = seq_len(length(ts)),
                        t = ts,
                        person_time = (at_risk_pla[1:(l - 1)] + at_risk_pla[2:l])*dt/2,
                        n = diff(cumul_case_pla))
    agg_vac <- data.frame(
                        tp = seq_len(length(ts)),
                        t = ts,
                        person_time = (at_risk_vac[1:(l - 1)] + at_risk_vac[2:l])*dt/2,
                        n = diff(cumul_case_vac))
    list(placebo = agg_pla, vaccine = agg_vac)
}
