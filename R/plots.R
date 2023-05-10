# COMPARE fits for the MATERNAL VACCINE

plot_outputs




plot_outputs <- function(data_eff, modelname) { 
    post_mat_get <- load(file = here::here("outputs", modelname, paste0(modelname, "_post_wane.Rdata")))
    efficacy_mat <- get(post_mat_get) 
    fits_mat <- read_rds(here::here("outputs", modelname, "fits.rda") )



    p1 <- fits_mat %>% as_draws_df %>% spread_draws(waning_exp[t], waning_er2[t], waning_er3[t]) %>%
            pivot_longer(c(waning_exp, waning_er2, waning_er3), names_to = "waning_dist", values_to = "eff") %>%
            ggplot() + 
                xlim(0, 365) + 
                    stat_lineribbon(aes(t, eff, fill = waning_dist, color = waning_dist), .width = 0.95, alpha = 0.5, size = 3 ) + 
                    theme_bw() + 
                    labs(x = "Days after birth", y = "Estimated efficacy per day", color = "Type") +
                    theme(text = element_text(size = 18))


    model_efficacy <- data_eff$days %>%
        map_df(
            function(x) {
                data.frame(
                    days = x,
                    efficacy_exp = 1:4000 %>% map_dbl(~efficacy_mat$wane_a_exp[.x] * exp( 1:x * - 1/efficacy_mat$wane_b_exp[.x]) %>% mean ),
                    efficacy_er2 = 1:4000 %>% map_dbl(~efficacy_mat$wane_a_er2[.x] * (1 - pgamma(1:x, 2, 1/efficacy_mat$wane_b_er2[.x])) %>% mean ),
                    efficacy_er3 = 1:4000 %>% map_dbl(~efficacy_mat$wane_a_er3[.x] * (1 - pgamma(1:x, 3, 1/efficacy_mat$wane_b_er3[.x])) %>% mean )
                )
            }
        )

    efficacy_mat_comparison <- bind_rows(
        model_efficacy %>% group_by(days) %>% summarise(get_mean_ci(efficacy_exp)) %>% mutate(type = "Fitted efficacy (exp)"),
        model_efficacy %>% group_by(days) %>% summarise(get_mean_ci(efficacy_er2)) %>% mutate(type = "Fitted efficacy (erlang2)"),
        model_efficacy %>% group_by(days) %>% summarise(get_mean_ci(efficacy_er3)) %>% mutate(type = "Fitted efficacy (erlang3)"),
        data_eff
    )

    p2 <- efficacy_mat_comparison %>% mutate(days = factor(days, levels = data_eff$days )) %>%
        ggplot() + 
            geom_linerange(aes(days, ymin = lb_95, ymax = ub_95 , color = type),
                position = position_dodge(width = 0.5), alpha = 0.6, size = 3) +
            geom_point(aes(days, y = mean, color = type), size = 6,
                position = position_dodge(width = 0.5)) + 
            geom_line(aes(days, y = mean, color = type, group = type), size = 3,
                position = position_dodge(width = 0.5)) + 
                theme_bw() + labs(x = "Days after birth", y = "Estimated efficacy (cumulative)", color = "Type") +
                theme(text = element_text(size = 18))
    p1 / p2
    ggsave(here::here("outputs", modelname, "fitscompare.png"),height = 10, width = 10)
    efficacy_mat_comparison
}








