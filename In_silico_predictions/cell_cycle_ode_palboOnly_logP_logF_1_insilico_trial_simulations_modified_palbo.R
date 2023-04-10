library(ggplot2)
library(ggsci)
library(rstan)
#library(bayesplot)
library(reshape2)
library(dplyr)
library(tidyr)
library(loo)
library(deSolve)
library(foreach)
library(doMC)
registerDoMC(cores = 5)

source('cell_cycle_posterior_predictions_modified.R')
#source('/home/srs54/ERBreastCancer/code/cell_cycle_ode_combine_datasets.R')

load("test_singleparameter_24ODES_even_1000iters_3rates.RData")
load("Rescale_ND_singleparameter_24ODES_even_1000iters_3rates.RData")

#palbo_cell_cycle_dat_ND_wide_noDay0 = palbo_cell_cycle_dat_ND_wide %>% filter(Day > 0.5) %>% mutate(Day = Day - 1)
#abema_cell_cycle_dat_ND_wide_noDay0 = abema_cell_cycle_dat_ND_wide %>% filter(Day > 0.5) %>% mutate(Day = Day - 1)

#palbo_cell_cycle_dat_Dox_wide_noDay0 = palbo_cell_cycle_dat_Dox_wide %>% filter(Day > 0.5) %>% mutate(Day = Day - 1)
#abema_cell_cycle_dat_Dox_wide_noDay0 = abema_cell_cycle_dat_Dox_wide %>% filter(Day > 0.5) %>% mutate(Day = Day - 1)

#palbo_cell_total_data_ND_noDay0 = palbo_cell_total_data_ND %>% filter(Day > 0.5) %>% mutate(Day = Day - 1)
#abema_cell_total_data_ND_noDay0 = abema_cell_total_data_ND %>% filter(Day > 0.5) %>% mutate(Day = Day - 1)

#palbo_cell_total_data_Dox_noDay0 = palbo_cell_total_data_Dox %>% filter(Day > 0.5) %>% mutate(Day = Day - 1)
#abema_cell_total_data_Dox_noDay0 = abema_cell_total_data_Dox %>% filter(Day > 0.5) %>% mutate(Day = Day - 1)

mu1 = 25.9
sd1 = 25.9*29/100
mu2 = 0.385
sd2 = 0.385*48/100
mu3 = 14000
sd3 = 14000*38/100
mu4 = 31
sd4 = 31*4.7/100

pk1_mean = log(mu1^2/(sqrt(mu1^2 + sd1^2)))
pk2_mean = log(mu2^2/(sqrt(mu2^2 + sd2^2)))
pk3_mean = log(mu3^2/(sqrt(mu3^2 + sd3^2)))
pk4_mean = log(mu4^2/(sqrt(mu4^2 + sd4^2)))

pk1_sd = sqrt(log(1 + (sd1^2)/(mu1^2)))
pk2_sd = sqrt(log(1 + (sd2^2)/(mu2^2)))
pk3_sd = sqrt(log(1 + (sd3^2)/(mu3^2)))
pk4_sd = sqrt(log(1 + (sd4^2)/(mu4^2)))

pk1 = rlnorm(3000, meanlog = pk1_mean, sdlog = pk1_sd)
pk2 = rlnorm(3000, meanlog = pk2_mean, sdlog = pk2_sd)
pk3 = rlnorm(3000, meanlog = pk3_mean, sdlog = pk3_sd)
pk4 = rlnorm(3000, meanlog = pk4_mean, sdlog = pk4_sd)

ND_posterior_params = as.data.frame(Rescale_ND_singleparameter_24ODES_even_1000iters_3rates,
                                    pars = c("alpha", "gamma", "delta", "b_P","b_F","c_P","c_F", "alpha_max", "a_FP")) %>% mutate(a_FA = 0, c_A = 1, b_A = 1,
                                                                                                                                  t_half_P=pk1, ka_F=pk2, v1_F=pk3, cl1_F=pk4 )
ND_init_palbo = as.data.frame(Rescale_ND_singleparameter_24ODES_even_1000iters_3rates, 
                              pars = c("y0_ct_p"))

#sigmoid_effectiveDose_ND_noDay0_2_posterior_init_abema = as.data.frame(fit_sigmoid_effectiveDose_ND_noDay0_2, pars = c('y0_cc_a[1]','y0_cc_a[2]','y0_cc_a[3]'))

colnames(ND_init_palbo) = c("G0G1","S","G2M")

ND_init_palbo$G0G1_1 = ND_init_palbo$G0G1/8
ND_init_palbo$G0G1_2 = ND_init_palbo$G0G1/8
ND_init_palbo$G0G1_3 = ND_init_palbo$G0G1/8
ND_init_palbo$G0G1_4 = ND_init_palbo$G0G1/8
ND_init_palbo$G0G1_5 = ND_init_palbo$G0G1/8
ND_init_palbo$G0G1_6 = ND_init_palbo$G0G1/8
ND_init_palbo$G0G1_7 = ND_init_palbo$G0G1/8
ND_init_palbo$G0G1_8 = ND_init_palbo$G0G1/8

ND_init_palbo$S_1 = ND_init_palbo$S/8
ND_init_palbo$S_2 = ND_init_palbo$S/8
ND_init_palbo$S_3 = ND_init_palbo$S/8
ND_init_palbo$S_4 = ND_init_palbo$S/8
ND_init_palbo$S_5 = ND_init_palbo$S/8
ND_init_palbo$S_6 = ND_init_palbo$S/8
ND_init_palbo$S_7 = ND_init_palbo$S/8
ND_init_palbo$S_8 = ND_init_palbo$S/8

ND_init_palbo$G2M_1 = ND_init_palbo$G2M/8
ND_init_palbo$G2M_2 = ND_init_palbo$G2M/8
ND_init_palbo$G2M_3 = ND_init_palbo$G2M/8
ND_init_palbo$G2M_4 = ND_init_palbo$G2M/8
ND_init_palbo$G2M_5 = ND_init_palbo$G2M/8
ND_init_palbo$G2M_6 = ND_init_palbo$G2M/8
ND_init_palbo$G2M_7 = ND_init_palbo$G2M/8
ND_init_palbo$G2M_8 = ND_init_palbo$G2M/8


ND_init_palbo = subset(ND_init_palbo, select = -c(1,2,3) )

ND_init_palbo$PP = 0
ND_init_palbo$FF = 0
ND_init_palbo$AA = 0

ND_init_palbo$PP_conc = 0
ND_init_palbo$FF_dose = 0
ND_init_palbo$FF_conc = 0

#odefit_test = ode(parms=ND_posterior_params[1,],
#                  func = ode_sigmoid_effectiveDose_g1s_transition_2,
#                  y = c(G0G1 = ND_init_palbo[1,1], 
#			S = ND_init_palbo[1,2], 
#			G2M = ND_init_palbo[1,3]),
#                  times = seq(0, 10, by = 1), drugA = 0, drugP = 0.005, drugF = 0)

Dox_posterior_params = as.data.frame(test_singleparameter_24ODES_even_1000iters_3rates,
                                     pars = c("alpha", "gamma", "delta", "b_P","b_F","c_P","c_F", "alpha_max", "a_FP")) %>% mutate(a_FA = 0, c_A = 1, b_A = 1,
                                                                                                                                   t_half_P=pk1, ka_F=pk2, v1_F=pk3, cl1_F=pk4)
Dox_init_palbo = as.data.frame(test_singleparameter_24ODES_even_1000iters_3rates,
                               pars = c("y0_ct_p"))
colnames(Dox_init_palbo) = c("G0G1","S","G2M")

Dox_init_palbo$G0G1_1 = Dox_init_palbo$G0G1/8
Dox_init_palbo$G0G1_2 = Dox_init_palbo$G0G1/8
Dox_init_palbo$G0G1_3 = Dox_init_palbo$G0G1/8
Dox_init_palbo$G0G1_4 = Dox_init_palbo$G0G1/8
Dox_init_palbo$G0G1_5 = Dox_init_palbo$G0G1/8
Dox_init_palbo$G0G1_6 = Dox_init_palbo$G0G1/8
Dox_init_palbo$G0G1_7 = Dox_init_palbo$G0G1/8
Dox_init_palbo$G0G1_8 = Dox_init_palbo$G0G1/8

Dox_init_palbo$S_1 = Dox_init_palbo$S/8
Dox_init_palbo$S_2 = Dox_init_palbo$S/8
Dox_init_palbo$S_3 = Dox_init_palbo$S/8
Dox_init_palbo$S_4 = Dox_init_palbo$S/8
Dox_init_palbo$S_5 = Dox_init_palbo$S/8
Dox_init_palbo$S_6 = Dox_init_palbo$S/8
Dox_init_palbo$S_7 = Dox_init_palbo$S/8
Dox_init_palbo$S_8 = Dox_init_palbo$S/8

Dox_init_palbo$G2M_1 = Dox_init_palbo$G2M/8
Dox_init_palbo$G2M_2 = Dox_init_palbo$G2M/8
Dox_init_palbo$G2M_3 = Dox_init_palbo$G2M/8
Dox_init_palbo$G2M_4 = Dox_init_palbo$G2M/8
Dox_init_palbo$G2M_5 = Dox_init_palbo$G2M/8
Dox_init_palbo$G2M_6 = Dox_init_palbo$G2M/8
Dox_init_palbo$G2M_7 = Dox_init_palbo$G2M/8
Dox_init_palbo$G2M_8 = Dox_init_palbo$G2M/8

Dox_init_palbo = subset(Dox_init_palbo, select = -c(1,2,3) )

Dox_init_palbo$PP = 0
Dox_init_palbo$FF = 0
Dox_init_palbo$AA = 0

Dox_init_palbo$PP_conc = 0
Dox_init_palbo$FF_dose = 0
Dox_init_palbo$FF_conc = 0

#Dox_posterior_init_abema = as.data.frame(fit_Dox, pars = c('y0_cc_a[1]','y0_cc_a[2]','y0_cc_a[3]'))

posterior_params = bind_rows(gather(ND_posterior_params, "parameter","value") %>% group_by(parameter) %>% mutate(treatment = "ND", id = row_number()),
			     gather(Dox_posterior_params, "parameter","value") %>% group_by(parameter) %>% mutate(treatment = "Dox", id = row_number()))

drug_schedules = list()
drug_schedules[["palbo_3w_1w_125mg"]] = data.frame(dose_times = c(seq(1, 21, by = 1), seq(29, 49, by = 1), seq(57, 77, by = 1)),
						   drug_amount = rep(52, length(c(seq(1, 21, by = 1), seq(29, 49, by = 1), seq(57, 77, by = 1)))))
drug_schedules[["palbo_daily_100mg"]] = data.frame(dose_times = seq(1, 84, by = 1), drug_amount = rep(47, 84))
drug_schedules[["palbo_daily_75mg"]] = data.frame(dose_times = seq(1, 84, by = 1), drug_amount = rep(29, 84))
drug_schedules[["palbo_BID_50mg_50mg"]] = data.frame(dose_times = seq(1, 84, by = 0.5), drug_amount = rep(21, 167))
drug_schedules[["palbo_BID_50mg_25mg"]] = data.frame(dose_times = seq(1, 84.5, by = 0.5), drug_amount = rep(c(21, 10), 84))
drug_schedules[["fulv_500mg_with_loading"]] = data.frame(dose_times = c(1.1, 15.1, 29.1, 57.1), drug_amount = rep(500000, 4))
drug_schedules[["fulv_500mg_only"]] = data.frame(dose_times = seq(1, 84, by = 1), drug_amount = rep(0, 84))


combine_drug_schedule = function(palbo_schedule, fulv_schedule){
	palbo_schedule = palbo_schedule %>% mutate(drug_type = "P")
	fulv_schedule = fulv_schedule %>% mutate(drug_type = "F")
	combo_schedule = rbind(palbo_schedule, fulv_schedule)
	combo_schedule = combo_schedule %>% arrange(dose_times)
	return(combo_schedule)
}
 
simulation_results_24 = list()
simulation_results_24[["fulv_500mg_only"]] = list()
simulation_results_24$fulv_500mg_only$combo_schedule = combine_drug_schedule(drug_schedules$fulv_500mg_only, drug_schedules$fulv_500mg_with_loading)
simulation_results_24$fulv_500mg_only$simulation_ND = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_ND(params = ND_posterior_params[i,],
                                                                         init = ND_init_palbo[i,],
                                                                         tmax = 100,
                                                                         dose_times = simulation_results_24$fulv_500mg_only$combo_schedule$dose_times,
                                                                         drug_type = simulation_results_24$fulv_500mg_only$combo_schedule$drug_type,
                                                                         drug_amount = simulation_results_24$fulv_500mg_only$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "ND")
}


simulation_results_24$fulv_500mg_only$simulation_Dox = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_Dox(params = Dox_posterior_params[i,],
                                                                          init = Dox_init_palbo[i,],
                                                                          tmax = 100,
                                                                          dose_times = simulation_results_24$fulv_500mg_only$combo_schedule$dose_times,
                                                                          drug_type = simulation_results_24$fulv_500mg_only$combo_schedule$drug_type,
                                                                          drug_amount = simulation_results_24$fulv_500mg_only$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "Dox")
}


simulation_results_24[["palbo_3w_1w_125mg"]] = list()
simulation_results_24$palbo_3w_1w_125mg$combo_schedule = combine_drug_schedule(drug_schedules$palbo_3w_1w_125mg, drug_schedules$fulv_500mg_with_loading)
simulation_results_24$palbo_3w_1w_125mg$simulation_ND = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_ND(params = ND_posterior_params[i,],
                                               init = ND_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_3w_1w_125mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_3w_1w_125mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_3w_1w_125mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "ND")
}


simulation_results_24$palbo_3w_1w_125mg$simulation_Dox = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_Dox(params = Dox_posterior_params[i,],
                                               init = Dox_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_3w_1w_125mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_3w_1w_125mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_3w_1w_125mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "Dox")
}


simulation_results_24[["palbo_daily_100mg"]] = list()
simulation_results_24$palbo_daily_100mg$combo_schedule = combine_drug_schedule(drug_schedules$palbo_daily_100mg, drug_schedules$fulv_500mg_with_loading)
simulation_results_24$palbo_daily_100mg$simulation_ND = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_ND(params = ND_posterior_params[i,],
                                               init = ND_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_daily_100mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_daily_100mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_daily_100mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "ND")
}
simulation_results_24$palbo_daily_100mg$simulation_Dox = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_Dox(params = Dox_posterior_params[i,],
                                                                    init = Dox_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_daily_100mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_daily_100mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_daily_100mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "Dox")
}

simulation_results_24[["palbo_daily_75mg"]] = list()
simulation_results_24$palbo_daily_75mg$combo_schedule = combine_drug_schedule(drug_schedules$palbo_daily_75mg, drug_schedules$fulv_500mg_with_loading)
simulation_results_24$palbo_daily_75mg$simulation_ND = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_ND(params = ND_posterior_params[i,],
                                                                    init = ND_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_daily_75mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_daily_75mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_daily_75mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "ND")
}

simulation_results_24$palbo_daily_75mg$simulation_Dox = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_Dox(params = Dox_posterior_params[i,],
                                                                    init = Dox_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_daily_75mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_daily_75mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_daily_75mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "Dox")
}

simulation_results_24[["palbo_BID_50mg_50mg"]] = list()
simulation_results_24$palbo_BID_50mg_50mg$combo_schedule = combine_drug_schedule(drug_schedules$palbo_BID_50mg_50mg, drug_schedules$fulv_500mg_with_loading)
simulation_results_24$palbo_BID_50mg_50mg$simulation_ND = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_ND(params = ND_posterior_params[i,],
                                                                    init = ND_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_BID_50mg_50mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_BID_50mg_50mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_BID_50mg_50mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "ND")
}
simulation_results_24$palbo_BID_50mg_50mg$simulation_Dox = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_Dox(params = Dox_posterior_params[i,],
                                                                    init = Dox_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_BID_50mg_50mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_BID_50mg_50mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_BID_50mg_50mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "Dox")
}

simulation_results_24[["palbo_BID_50mg_25mg"]] = list()
simulation_results_24$palbo_BID_50mg_25mg$combo_schedule = combine_drug_schedule(drug_schedules$palbo_BID_50mg_25mg, drug_schedules$fulv_500mg_with_loading)
simulation_results_24$palbo_BID_50mg_25mg$simulation_ND = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_ND(params = ND_posterior_params[i,],
                                                                    init = ND_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_BID_50mg_25mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_BID_50mg_25mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_BID_50mg_25mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "ND")
}
simulation_results_24$palbo_BID_50mg_25mg$simulation_Dox = foreach(i=1:1500, .combine = rbind) %dopar% {
  df = ode_sigmoid_effectiveDose_g1s_transition_1_invivo_pk_wrapper_Dox(params = Dox_posterior_params[i,],
                                                                    init = Dox_init_palbo[i,],
                                               tmax = 100,
                                               dose_times = simulation_results_24$palbo_BID_50mg_25mg$combo_schedule$dose_times,
                                               drug_type = simulation_results_24$palbo_BID_50mg_25mg$combo_schedule$drug_type,
                                               drug_amount = simulation_results_24$palbo_BID_50mg_25mg$combo_schedule$drug_amount)
  df = as.data.frame(df)
  df = df %>% mutate(id = i, total = G0G1_1 + G0G1_2 +G0G1_3 +G0G1_4 +G0G1_5 +G0G1_6 +G0G1_7 + 
                       G0G1_8 + G2M_1 +G2M_2 +G2M_3 +G2M_4 +G2M_5 +G2M_6 +G2M_7 +G2M_8 + 
                       S_1+ S_2+ S_3+ S_4+ S_5+ S_6+ S_7+ S_8, treatment = "Dox")
}


save.image("cell_cycle_ode_noLongtermData_palboOnly_logF_logP_insilico_trial_simulations_ct_24_FulOnly_3rates_100days_pk.Rdata")

#save.image("cell_cycle_ode_noLongtermData_palboOnly_logF_logP_insilico_trial_simulations_ct_24_FulOnly_3rates_100days.Rdata")

#save.image("cell_cycle_ode_noLongtermData_palboOnly_logF_logP_insilico_trial_simulations_ct.Rdata")

#save.image("cell_cycle_ode_noLongtermData_palboOnly_logF_logP_insilico_trial_simulations.Rdata")


