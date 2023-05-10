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
        n_vac = data_vac$vaccine$n,
        lower_bounds = lower_bounds
    )

    T <- max(data_vac$placebo$t)
    mod <- cmdstan_model(here::here("src", "eff_est.stan"))

    fit_stan <- mod$sample(
        data = data_list, 
        seed = 123, 
        chains = 4, 
        parallel_chains = 4,
        refresh = 500 # print update every 500 iters
    )

    if (!dir.exists(here::here("outputs", name))) {
        # create the folder if it doesn't exist
        dir.create(here::here("outputs", name))
    }


    fit_stan$save_object(file = here::here("outputs", name, "fits.rda"))

    posteriors_wane <- fit_stan %>% as_draws_df %>%
        spread_draws(wane_a_exp, wane_b_exp, wane_a_er2, wane_b_er2, wane_a_er3, wane_b_er3)
    save(posteriors_wane, file = here::here("outputs", name, paste0(name, "_post_wane.RData") ))

    df_data_pla <- data.frame(
        t = data_list$t_pla, 
        inci = (data_list$n_pla) / (data_list$persons_pla)
    )

    df_data_vac <- data.frame(
        t = data_list$t_vac, 
        inci = (data_list$n_vac) / (data_list$persons_vac)
    )

    p1 <- fit_stan %>% as_draws_df %>% spread_draws(incidence_pla[t], time_p[t]) %>%
        ggplot() + 
            stat_lineribbon(aes(time_p, incidence_pla), .width = 0.95, fill = "blue", alpha = 0.5) + 
            geom_point(data = df_data_pla, aes(t, inci), shape = 21, fill = "red") +
            theme_bw() + labs(x = "Time post-vaccination (days)", y = "Incidence (placebo group)")

    p2 <- fit_stan %>% as_draws_df %>% spread_draws(waning_exp[t]) %>% 
        ggplot() + 
            stat_lineribbon(aes(t, waning_exp), .width = 0.95, fill = "red", alpha = 0.5) + theme_bw() +
            labs(x = "Time post-vaccination (days)", y = "Estimated efficacy")
    p3 <- fit_stan %>% as_draws_df %>% spread_draws(waning_er2[t]) %>% 
        ggplot() + 
            stat_lineribbon(aes(t, waning_er2), .width = 0.95, fill = "red", alpha = 0.5) + theme_bw() +
            labs(x = "Time post-vaccination (days)", y = "Estimated efficacy")
    p4 <- fit_stan %>% as_draws_df %>% spread_draws(waning_er2[t]) %>% 
        ggplot() + 
            stat_lineribbon(aes(t, waning_er3), .width = 0.95, fill = "red", alpha = 0.5) + theme_bw() +
            labs(x = "Time post-vaccination (days)", y = "Estimated efficacy")
    (p1) / (p2 + p3 + p4)
    ggsave(here::here("outputs", name, "summary_fig.png"), height = 10, width = 10)

}