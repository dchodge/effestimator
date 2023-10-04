# This function calls stan and fits eveything using HMC 

fitted_functions <- function(data_vac, name, lower_bounds) {

    # define the data
    data_list <- list(
        N_pla = nrow(data_vac$placebo),
        tp_pla = data_vac$placebo$tp,
        t_pla = data_vac$placebo$t,
        persons_pla =  data_vac$placebo$person_time,
        n_pla = data_vac$placebo$n,

        N_vac = nrow(data_vac$vaccine),
        tp_vac = data_vac$vaccine$tp,
        t_vac = data_vac$vaccine$t,
        persons_vac = data_vac$vaccine$person_time,
        n_vac = data_vac$vaccine$n
    )

    T <- max(data_vac$placebo$t)
    mod <- cmdstan_model(here::here("src", "eff_est.stan"), compile = TRUE)
    mod_old <- cmdstan_model(here::here("src", "eff_est_old.stan"), compile = TRUE)

    fit_stan_exp <- mod$sample(
        data = c(model_id = 1, data_list, lower_bound = lower_bounds[1]), 
        seed = 123, 
        chains = 4, 
        parallel_chains = 4,
        refresh = 500 # print update every 500 iters
    )
    fit_stan_er2 <- mod$sample(
      data = c(model_id = 2, data_list, lower_bound = lower_bounds[2]), 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500 # print update every 500 iters
    )
    fit_stan_er3 <- mod$sample(
      data = c(model_id = 3, data_list, lower_bound = lower_bounds[3]), 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500 # print update every 500 iters
    )
    fit_stan_old <- mod_old$sample(
      data = c(data_list, lower_bounds = list(lower_bounds)), 
      seed = 123, 
      chains = 4, 
      parallel_chains = 4,
      refresh = 500 # print update every 500 iters
    )

    if (!dir.exists(here::here("outputs", name))) {
        # create the folder if it doesn't exist
        dir.create(here::here("outputs", name))
    }
    
    fit_stan_exp$save_object(file = here::here("outputs", name, "fits_exp.rda"))
    fit_stan_er2$save_object(file = here::here("outputs", name, "fits_er2.rda"))
    fit_stan_er3$save_object(file = here::here("outputs", name, "fits_er3.rda"))
    fit_stan_old$save_object(file = here::here("outputs", name, "fits_old.rda"))

    fit_stan_exp %>% as_draws_df %>% spread_draws(wane_a, wane_b) %>% 
      rename(wane_a_exp = wane_a, wane_b_exp = wane_b) -> posteriors_wane_exp
    fit_stan_er2 %>% as_draws_df %>% spread_draws(wane_a, wane_b) %>% 
      rename(wane_a_er2 = wane_a, wane_b_er2 = wane_b) -> posteriors_wane_er2
    fit_stan_er3 %>% as_draws_df %>% spread_draws(wane_a, wane_b) %>% 
      rename(wane_a_er3 = wane_a, wane_b_er3 = wane_b) -> posteriors_wane_er3
    
    posteriors_wane_exp %>% 
      left_join(posteriors_wane_er2, by = c(".chain", ".iteration", ".draw")) %>% 
      left_join(posteriors_wane_er3, by = c(".chain", ".iteration", ".draw")) -> posteriors_wane

    save(posteriors_wane, file = here::here("outputs", name, paste0(name, "_post_wane.RData") ))

    df_data_pla <- data.frame(
        t = data_list$t_pla, 
        inci = (data_list$n_pla) / (data_list$persons_pla)
    )

    df_data_vac <- data.frame(
        t = data_list$t_vac, 
        inci = (data_list$n_vac) / (data_list$persons_vac)
    )

    # Technically this fit (p1) is different for each waning fn. Just show exp fit for now
    p1 <- fit_stan_exp %>% as_draws_df %>% spread_draws(incidence_pla[t], time_p[t]) %>%
        ggplot() + 
            stat_lineribbon(aes(time_p, incidence_pla), .width = 0.95, fill = "blue", alpha = 0.5) + 
            geom_point(data = df_data_pla, aes(t, inci), shape = 21, fill = "red") +
            theme_bw() + labs(x = "Time post-vaccination (days)", y = "Incidence (placebo group)")

    p2 <- fit_stan_exp %>% as_draws_df %>% spread_draws(waning[t]) %>% 
        ggplot() + 
            stat_lineribbon(aes(t, waning), .width = 0.95, fill = "red", alpha = 0.5) + theme_bw() +
            labs(x = "Time post-vaccination (days)", y = "Estimated efficacy")
    p3 <- fit_stan_er2 %>% as_draws_df %>% spread_draws(waning[t]) %>% 
        ggplot() + 
            stat_lineribbon(aes(t, waning), .width = 0.95, fill = "red", alpha = 0.5) + theme_bw() +
            labs(x = "Time post-vaccination (days)", y = "Estimated efficacy")
    p4 <- fit_stan_er3 %>% as_draws_df %>% spread_draws(waning[t]) %>% 
        ggplot() + 
            stat_lineribbon(aes(t, waning), .width = 0.95, fill = "red", alpha = 0.5) + theme_bw() +
            labs(x = "Time post-vaccination (days)", y = "Estimated efficacy")
    (p1) / (p2 + p3 + p4)
    ggsave(here::here("outputs", name, "summary_fig.png"), height = 10, width = 10)

}