/// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6120589/pdf/nihms1502750.pdf
functions{
    // functions assocated with waning
    vector exp_wane(real wane_a, real wane_b, vector t_vac) {
        return wane_a * exp(- (wane_b) * t_vac);
    }
    vector erlangk_wane(real wane_a, real wane_b, vector t_vac, int k) {
        int N = size(t_vac);
        vector[N] gamma_vec;
        for (i in 1:N) {
            gamma_vec[i] = gamma_cdf(t_vac[i]| k, wane_b);
        }
        return wane_a * (1 - gamma_vec);
    }
}
data {
    int N_pla;
    array[N_pla] int tp_pla;
    vector[N_pla] t_pla; // number of time points in a vector (for gaussian process)
    vector[N_pla] persons_pla;
    array[N_pla] int n_pla;

    int N_vac;
    array[N_pla] int tp_vac;
    vector[N_vac] t_vac;
    vector[N_vac] persons_vac; // number of time points in a vector (for gaussian process)
    array[N_vac] int  n_vac;
    vector[3] lower_bounds;
} 
transformed data {
    real delta = 1e-10; // just a thing to stabilise gp
}
parameters{
    real<lower = 0, upper = 1> wane_a_exp; // wane a rate for exponential waning
    real<lower = lower_bounds[1], upper = 0.1> wane_b_exp; // wane b rate for exponential waning
 
    real<lower = 0, upper = 1> wane_a_er2; // wane a rate for erlang-2 waning
    real<lower = lower_bounds[2], upper = 0.1> wane_b_er2; // wane b rate for erlang-2 waning

    real<lower = 0, upper = 1> wane_a_er3; // wane a rate for erlang-3 waning
    real<lower = lower_bounds[3], upper = 0.1> wane_b_er3; // wane b rate for erlang-3 waning

    // add in gaussian processes to model the incidence over time
    real<lower = 0> alpha; // step size
    real<lower = 0> rho; // length
    vector[N_pla] eta; // non-centralisation parameters

}
transformed parameters {
    vector[N_pla] inci_gp;
    matrix[N_pla, N_pla] K;
    matrix[N_pla, N_pla] L_K;

    // This is the Gaussian process stuff, K, L_K are covariance matrices
    K = gp_exp_quad_cov(tp_pla, alpha, rho) + diag_matrix(rep_vector(delta, N_pla));
    L_K = cholesky_decompose(K);
    inci_gp = L_K * eta;

}
model {
    // exponential of inci_gp to ensure positivity 
    n_pla ~ poisson(persons_pla .* exp(inci_gp));
    n_vac ~ poisson(persons_vac .* exp(inci_gp) .* (1 - exp_wane(wane_a_exp, wane_b_exp, t_vac) ) );
    n_vac ~ poisson(persons_vac .* exp(inci_gp) .* (1 - erlangk_wane(wane_a_er2, wane_b_er2, t_vac, 2) ) );
    n_vac ~ poisson(persons_vac .* exp(inci_gp) .* (1 - erlangk_wane(wane_a_er3, wane_b_er3, t_vac, 3) ) );

    // GP priors 
    rho ~ inv_gamma(5, 5);
    alpha ~ normal(0, 1);
    eta ~ std_normal(); 

    // waning priors
    wane_a_exp ~ uniform(0, 1);
    wane_b_exp ~ uniform(lower_bounds[1], 0.1);
    wane_a_er2 ~ uniform(0, 1);
    wane_b_er2 ~ uniform(lower_bounds[2], 0.1);
    wane_a_er3 ~ uniform(0, 1);
    wane_b_er3 ~ uniform(lower_bounds[3], 0.1);
}
generated quantities {
    // this section just gives useful outputs
    vector[N_pla] incidence_pla;
    vector[N_pla] incidence_vac_exp;
    vector[N_pla] incidence_vac_er2;
    vector[N_pla] incidence_vac_er3;

    vector[730] waning_exp;
    vector[730] waning_er2;
    vector[730] waning_er3;

    vector[730] time_points = linspaced_vector(730, 0, 730 - 1);
    // Waning efficacy over 730 days
    waning_exp = exp_wane(wane_a_exp, wane_b_exp, time_points);
    waning_er2 = erlangk_wane(wane_a_er2, wane_b_er2, time_points, 2);
    waning_er3 = erlangk_wane(wane_a_er3, wane_b_er3, time_points, 3);

    // Waning efficacy at time points in data
    vector[N_pla] waning_exp_inc;
    vector[N_pla] waning_er2_inc;
    vector[N_pla] waning_er3_inc;

    vector[N_pla] time_p;
    waning_exp_inc = exp_wane(wane_a_exp, wane_b_exp, t_vac);
    waning_er2_inc = erlangk_wane(wane_a_er2, wane_b_er2, t_vac, 2);
    waning_er3_inc = erlangk_wane(wane_a_er3, wane_b_er3, t_vac, 3);

    // Estimated incidence at time points in data, for comparison with data
    for(i in 1:N_pla) {
        time_p[i] = t_pla[i];
        incidence_pla[i] = exp(inci_gp[i]);
        incidence_vac_exp[i] = exp(inci_gp[i]) * (1 - waning_exp_inc[i]);
        incidence_vac_er2[i] = exp(inci_gp[i]) * (1 - waning_er2_inc[i]);
        incidence_vac_er3[i] = exp(inci_gp[i]) * (1 - waning_er3_inc[i]);
    }
}
