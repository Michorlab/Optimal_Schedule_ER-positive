setwd("~/Desktop/PostDoc_Research/2021/Modified Code/Cluster")

source("posterior_prediction_functions_validation2.R")
load("stan_input_include_longterm_data.RData")
load("stan_input_include_longterm_data_2.RData")
load("stan_input_include_longterm_data_new.RData")
load("stan_input_include_longterm_data_modified.RData")
load("stan_input_include_longterm_data_modified_reduced.RData")

source("post_processing_script.R")

library(rstan)
library(abind)

load("cell_cycle_ode_effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_neutral.RData")

posterior_prediction_plots_no_abema(model_data = effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_neutral", dox = "Dox", skip_pass = TRUE)

model_qc(data = effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_1, name = "effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_neutral",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))















load("test_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_15ODES_even_1000iters, name = "test_singleparameter_15ODES_even_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("test_singleparameter_15ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)



load("test_singleparameter_15ODES_even_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_newdata", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_15ODES_even_1000iters, name = "test_singleparameter_15ODES_even_newdata",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("test_singleparameter_15ODES_even_test_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_test_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_test_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_15ODES_even_test1_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_test1_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_test1_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_15ODES_even_test2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_test2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_test2_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_15ODES_even_test3_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_test3_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_test3_newdata", dox = "Dox", skip_pass = TRUE)


looic_15_even_test3_newdata_logP_logF = calculate_looic(test_singleparameter_15ODES_even_test3_newdata)
print(looic_15_even_test3_newdata_logP_logF)

load("test_singleparameter_15ODES_even_ful_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_ful_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_ful_G2_newdata", dox = "Dox", skip_pass = TRUE)

looic_15_even_ful_G2_newdata_logP_logF = calculate_looic(test_singleparameter_15ODES_even_ful_G2_newdata)
print(looic_15_even_ful_G2_newdata_logP_logF)

load("test_singleparameter_15ODES_even_ful_S_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_ful_S_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_ful_S_newdata", dox = "Dox", skip_pass = TRUE)

looic_15_even_ful_S_newdata_logP_logF = calculate_looic(test_singleparameter_15ODES_even_ful_S_newdata)
print(looic_15_even_ful_S_newdata_logP_logF)

load("test_singleparameter_15ODES_even_Loewe_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_Loewe_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_Loewe_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_15ODES_even_multiple_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_multiple_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_multiple_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_15ODES_even_syn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_syn_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_15ODES_even_syn_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_16ODES_G2_modified.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_modified, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_modified", dox = "Dox", skip_pass = TRUE)


model_qc(data = test_singleparameter_16ODES_G2_modified, name = "test_singleparameter_16ODES_G2_modified",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F", "alpha_max"),
         pars_acf1 = c("alpha", "b_P"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("test_singleparameter_16ODES_G2_syn_modified.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_syn_modified, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_syn_modified", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_16ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_16ODES_G2_syn_modified, name = "test_singleparameter_16ODES_G2_syn_modified",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("test_singleparameter_16ODES_G2_both_G2_nonsyn_modified.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_both_G2_nonsyn_modified, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_both_G2_nonsyn_modified", dox = "Dox", skip_pass = TRUE)


model_qc(data = test_singleparameter_16ODES_G2_both_G2_nonsyn_modified, name = "test_singleparameter_16ODES_G2_both_G2_nonsyn_modified",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F", "b_P_2", "b_F_2", "c_P_2","c_F_2", "alpha_max"),
         pars_acf1 = c("alpha","b_P","b_F","c_P","c_F"),
         pars_acf2 = c("alpha","b_P_2","b_F_2","c_P_2","c_F_2"))



load("Rescale_ND_singleparameter_16ODES_G2_modified.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_modified, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_modified", dox = "ND", skip_pass = TRUE)
model_qc(data = Rescale_ND_singleparameter_16ODES_G2_modified, name = "Rescale_ND_singleparameter_16ODES_G2_modified",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","alpha_max"),
         pars_acf1 = c("alpha","b_P"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Rescale_ND_singleparameter_16ODES_G2_syn_modified.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_syn_modified, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_syn_modified", dox = "ND", skip_pass = TRUE)
model_qc(data = Rescale_ND_singleparameter_16ODES_G2_syn_modified, name = "Rescale_ND_singleparameter_16ODES_G2_syn_modified",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_modified.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_modified, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_modified", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_modified, name = "Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_modified",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F", "b_P_2", "b_F_2", "c_P_2","c_F_2", "alpha_max"),
         pars_acf1 = c("alpha","b_P","b_F","c_P","c_F"),
         pars_acf2 = c("alpha","b_P_2","b_F_2","c_P_2","c_F_2"))





load("test_singleparameter_16ODES_G2_apop_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_apop_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_apop_newdata", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_21ODES_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_21ODES_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_21ODES_G2_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_21ODES_even_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_21ODES_even_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_21ODES_even_newdata", dox = "Dox", skip_pass = TRUE)




load("test_singleparameter_10ODES_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_10ODES_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_10ODES_G2_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_11ODES_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_11ODES_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_11ODES_G2_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_13ODES_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_13ODES_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_13ODES_G2_newdata", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_16ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_1000iters_3rates", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_17ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_17ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_17ODES_G2_1000iters_3rates", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_18ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_18ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_18ODES_G2_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_19ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_19ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_19ODES_G2_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_20ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_20ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_2ODES_G2_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_21ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_21ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_21DES_G2_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_22ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_22ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_22DES_G2_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_23ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_23ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_23DES_G2_1000iters_3rates", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_14ODES_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_14ODES_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_14ODES_G2_newdata", dox = "Dox", skip_pass = TRUE)


looic_14_G2_newdata_logP_logF = calculate_looic(test_singleparameter_14ODES_G2_newdata)
print(looic_14_G2_newdata_logP_logF)

load("test_singleparameter_14ODES_G2_ful_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_14ODES_G2_ful_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_14ODES_G2_ful_G2_newdata", dox = "Dox", skip_pass = TRUE)



looic_14_G2_ful_G2_newdata_logP_logF = calculate_looic(test_singleparameter_14ODES_G2_ful_G2_newdata)
print(looic_14_G2_ful_G2_newdata_logP_logF)

load("test_singleparameter_14ODES_G2_ful_S_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_14ODES_G2_ful_S_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_14ODES_G2_ful_S_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_14ODES_G2_palbo_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_14ODES_G2_palbo_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_14ODES_G2_palbo_G2_newdata", dox = "Dox", skip_pass = TRUE)


looic_14_G2_palbo_G2_newdata_logP_logF = calculate_looic(test_singleparameter_14ODES_G2_palbo_G2_newdata)
print(looic_14_G2_palbo_G2_newdata_logP_logF)

load("test_singleparameter_14ODES_G2_palbo_S_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_14ODES_G2_palbo_S_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_14ODES_G2_palbo_S_newdata", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_14ODES_G2_both_G2_syn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_14ODES_G2_both_G2_syn_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_14ODES_G2_both_G2_syn_newdata", dox = "Dox", skip_pass = TRUE)

looic_14_G2_both_G2_syn_newdata_logP_logF = calculate_looic(test_singleparameter_14ODES_G2_both_G2_syn_newdata)
print(looic_14_G2_both_G2_syn_newdata_logP_logF)


load("test_singleparameter_14ODES_G2_both_G2_syn3_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_14ODES_G2_both_G2_syn3_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_14ODES_G2_both_G2_syn3_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_14ODES_G2_both_G2_nonsyn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_14ODES_G2_both_G2_nonsyn_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_14ODES_G2_both_G2_nonsyn_newdata", dox = "Dox", skip_pass = TRUE)

looic_14_G2_both_G2_nonsyn_newdata_logP_logF = calculate_looic(test_singleparameter_14ODES_G2_both_G2_nonsyn_newdata)
print(looic_14_G2_both_G2_nonsyn_newdata_logP_logF)

load("test_singleparameter_16ODES_G2_ful_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_ful_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_ful_G2_newdata", dox = "Dox", skip_pass = TRUE)

looic_16_G2_ful_G2_newdata_logP_logF = calculate_looic(test_singleparameter_16ODES_G2_ful_G2_newdata)
print(looic_16_G2_ful_G2_newdata_logP_logF)

load("test_singleparameter_16ODES_G2_both_G2_nonsyn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_both_G2_nonsyn_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_both_G2_nonsyn_newdata", dox = "Dox", skip_pass = TRUE)

looic_16_G2_both_G2_nonsyn_newdata_logP_logF = calculate_looic(test_singleparameter_16ODES_G2_both_G2_nonsyn_newdata)
print(looic_16_G2_both_G2_nonsyn_newdata_logP_logF)

load("test_singleparameter_16ODES_G2_both_G2_syn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_both_G2_syn_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_both_G2_syn_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_16ODES_G2_both_G2_syn1_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_both_G2_syn1_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_both_G2_syn1_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_16ODES_G2_both_G2_syn2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_both_G2_syn2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_both_G2_syn2_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_16ODES_G2_both_G2_syn3_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_both_G2_syn3_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_16ODES_G2_both_G2_syn3_newdata", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_20ODES_G2_both_G2_nonsyn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_20ODES_G2_both_G2_nonsyn_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_20ODES_G2_both_G2_nonsyn_newdata", dox = "Dox", skip_pass = TRUE)




load("test_singleparameter_20ODES_G2_ful_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_20ODES_G2_ful_G2_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_20ODES_G2_ful_G2_newdata", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_20ODES_G2_both_G2_syn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_20ODES_G2_both_G2_syn_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_20ODES_G2_both_G2_syn_newdata", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_20ODES_G2_both_G2_nonsyn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_20ODES_G2_both_G2_nonsyn_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_20ODES_G2_both_G2_nonsyn_newdata", dox = "Dox", skip_pass = TRUE)

load("cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_logF_Dox_1_modified.RData")

posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1_modified, stan_data = stan_input_include_longterm_data_2,
                                    name = "effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1_modified", dox = "Dox", skip_pass = TRUE)

looic_logic_newdata_logP_logF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1_modified)
print(looic_logic_newdata_logP_logF)



load("cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_logF_ND_1_modified.RData")

posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1_modified, stan_data = stan_input_include_longterm_data_2,
                                    name = "effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1_modified", dox = "ND", skip_pass = TRUE)

looic_logic_ND_newdata_logP_logF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1_modified)
print(looic_logic_ND_newdata_logP_logF)

looic_logic_Dox_newdata_logP_logF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1_modified)
print(looic_logic_Dox_newdata_logP_logF)


looic_logic_Dox_ND_newdata_logP_logF = calculate_looic_2(effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1_modified, effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1_modified)
print(looic_logic_Dox_ND_newdata_logP_logF)


load("Rescale_ND_singleparameter_15ODES_even_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_15ODES_even_newdata", dox = "ND", skip_pass = TRUE)

looic_ND_15_even_newdata_logP_logF = calculate_looic(Rescale_ND_singleparameter_15ODES_even_newdata)
print(looic_ND_15_even_newdata_logP_logF)


load("Rescale_ND_singleparameter_16ODES_G2_1000iters_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_1000iters_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_1000iters_newdata", dox = "ND", skip_pass = TRUE)

looic_ND_16_G2_newdata_logP_logF = calculate_looic(Rescale_ND_singleparameter_16ODES_G2_1000iters_newdata)
print(looic_ND_16_G2_newdata_logP_logF)

load("Rescale_ND_singleparameter_16ODES_G2_ful_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_ful_G2_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_ful_G2_newdata", dox = "ND", skip_pass = TRUE)




load("Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_newdata", dox = "ND", skip_pass = TRUE)

looic_ND_16_G2_both_G2_nonsyn_newdata_logP_logF = calculate_looic(Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_newdata)
print(looic_ND_16_G2_both_G2_nonsyn_newdata_logP_logF)

load("Rescale_ND_singleparameter_16ODES_G2_both_G2_syn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_both_G2_syn_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_both_G2_syn_newdata", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_16ODES_G2_both_G2_syn1_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_both_G2_syn1_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_both_G2_syn1_newdata", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_16ODES_G2_both_G2_syn2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_both_G2_syn2_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_both_G2_syn2_newdata", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_16ODES_G2_both_G2_syn3_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_both_G2_syn3_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_both_G2_syn3_newdata", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_20ODES_G2_1000iters_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_1000iters_newdata, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_1000iters_newdata", dox = "ND", skip_pass = TRUE)



load("Rescale_ND_singleparameter_20ODES_G2_ful_G2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_ful_G2_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_ful_G2_newdata", dox = "ND", skip_pass = TRUE)




load("Rescale_ND_singleparameter_20ODES_G2_both_G2_nonsyn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_both_G2_nonsyn_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_both_G2_nonsyn_newdata", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_20ODES_G2_both_G2_syn_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_both_G2_syn_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_both_G2_syn_newdata", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_20ODES_G2_both_G2_syn1_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_both_G2_syn1_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_both_G2_syn1_newdata", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_20ODES_G2_both_G2_syn2_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_both_G2_syn2_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_both_G2_syn2_newdata", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_20ODES_G2_both_G2_syn3_newdata.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_both_G2_syn3_newdata, stan_data = stan_input_include_longterm_data_new,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_both_G2_syn3_newdata", dox = "ND", skip_pass = TRUE)

load("test_singleparameter_15ODES_even_multiple.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_multiple, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_15ODES_even_multiple", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_15ODES_even_multiple, name = "test_singleparameter_15ODES_even_multiple",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max", "Beta", "Gamma"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("test_singleparameter_15ODES_even_IC.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_IC, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_15ODES_even_IC", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_15ODES_even_multiple, name = "test_singleparameter_15ODES_even_IC",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max", "Beta", "Gamma"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("test_singleparameter_15ODES_even_syn.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_syn, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_15ODES_even_syn", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_15ODES_even_Bliss.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_Bliss, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_15ODES_even_Bliss", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_15ODES_even_Loewe.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_Loewe, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_15ODES_even_Loewe", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_15ODES_even_Loewe_2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_Loewe_2, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_15ODES_even_Loewe_2", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_15ODES_even_Loewe_3.RData")


load("test_singleparameter_15ODES_even_ful_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_15ODES_even_ful_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_15ODES_even_ful_G2_1000iters", dox = "Dox", skip_pass = TRUE)


model_qc(data = test_singleparameter_15ODES_even_ful_G2_1000iters, name = "test_singleparameter_15ODES_even_ful_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))




load("test_singleparameter_16ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_16ODES_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_16ODES_G2_1000iters, name = "test_singleparameter_16ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("test_singleparameter_16ODES_G2_ful_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_16ODES_G2_ful_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_16ODES_G2_ful_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_16ODES_G2_ful_G2_1000iters, name = "test_singleparameter_16ODES_G2_ful_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("test_singleparameter_17ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_17ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_17ODES_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_17ODES_G2_1000iters, name = "test_singleparameter_17ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("test_singleparameter_17ODES_G2_ful_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_17ODES_G2_ful_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_17ODES_G2_ful_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_17ODES_G2_ful_G2_1000iters, name = "test_singleparameter_17ODES_G2_ful_G2_1000iters",
        params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
        pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
        pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))



load("test_singleparameter_18ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_18ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_18ODES_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_18ODES_G2_1000iters, name = "test_singleparameter_18ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("test_singleparameter_18ODES_G2_ful_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_18ODES_G2_ful_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_18ODES_G2_ful_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_18ODES_G2_ful_G2_1000iters, name = "test_singleparameter_18ODES_G2_ful_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("test_singleparameter_19ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_19ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_19ODES_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_19ODES_G2_1000iters, name = "test_singleparameter_19ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("test_singleparameter_19ODES_G2_ful_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_19ODES_G2_ful_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_19ODES_G2_ful_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_19ODES_G2_ful_G2_1000iters, name = "test_singleparameter_19ODES_G2_ful_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("test_singleparameter_20ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_20ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_20ODES_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_20ODES_G2_1000iters, name = "test_singleparameter_20ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))



load("test_singleparameter_20ODES_G2_ful_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_20ODES_G2_ful_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_20ODES_G2_ful_G2_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_20ODES_G2_ful_G2_1000iters, name = "test_singleparameter_20ODES_G2_ful_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))




load("test_singleparameter_18ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_18ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_18ODES_G2", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_18ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_18ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_18ODES_G2_ful_G2", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_19ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_19ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_19ODES_G2", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_19ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_19ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_19ODES_G2_ful_G2", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_20ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_20ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_20ODES_G2", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_20ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_20ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_20ODES_G2_ful_G2", dox = "Dox", skip_pass = TRUE)









load("ND_singleparameter_15ODES_even.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_15ODES_even, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_15ODES_even", dox = "ND", skip_pass = TRUE)

load("ND_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_15ODES_even_1000iters", dox = "ND", skip_pass = TRUE)



load("ND_singleparameter_15ODES_even_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_15ODES_even_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_15ODES_even_ful_G2", dox = "ND", skip_pass = TRUE)

model_qc(data = ND_singleparameter_15ODES_even_ful_G2, name = "ND_singleparameter_15ODES_even_ful_G2",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("ND_singleparameter_15ODES_even_temp.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_15ODES_even_temp, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_15ODES_even_temp", dox = "ND", skip_pass = TRUE)




load("ND_singleparameter_16ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_16ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_16ODES_G2", dox = "ND", skip_pass = TRUE)

model_qc(data = ND_singleparameter_16ODES_G2, name = "ND_singleparameter_16ODES_G2",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("ND_singleparameter_16ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_16ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_16ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)

model_qc(data = ND_singleparameter_16ODES_G2_1000iters, name = "ND_singleparameter_16ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("ND_singleparameter_16ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_16ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_16ODES_G2_ful_G2", dox = "ND", skip_pass = TRUE)



load("ND_singleparameter_17ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_17ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_17ODES_G2", dox = "ND", skip_pass = TRUE)



load("ND_singleparameter_17ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_17ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_17ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)
model_qc(data = ND_singleparameter_17ODES_G2_1000iters, name = "ND_singleparameter_17ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_17ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_17ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_17ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)



load("ND_singleparameter_17ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_17ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_17ODES_G2_ful_G2", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_17ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_17ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_17ODES_G2_ful_G2", dox = "ND", skip_pass = TRUE)



load("Rescale_ND_singleparameter_17ODES_G2_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_17ODES_G2_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_17ODES_G2_ful_G2_temp", dox = "ND", skip_pass = TRUE)



load("ND_singleparameter_18ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_18ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_18ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)

model_qc(data = ND_singleparameter_18ODES_G2_1000iters, name = "ND_singleparameter_18ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_18ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_18ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_18ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)



load("Rescale_ND_singleparameter_18ODES_G2_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_18ODES_G2_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_18ODES_G2_ful_G2_temp", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_19ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_19ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_19ODES_G2_ful_G2", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_20ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_ful_G2", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_20ODES_G2_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_ful_G2_temp", dox = "ND", skip_pass = TRUE)




load("ND_singleparameter_19ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_19ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_19ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)

model_qc(data = ND_singleparameter_19ODES_G2_1000iters, name = "ND_singleparameter_19ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_19ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_19ODES_G2_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_19ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_19ODES_G2_1000iters, name = "Rescale_ND_singleparameter_19ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp", dox = "ND", skip_pass = TRUE)

load("ND_singleparameter_20ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = ND_singleparameter_20ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "ND_singleparameter_20ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_20ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_20ODES_G2_1000iters, name = "Rescale_ND_singleparameter_20ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_20ODES_G2_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_ful_G2_temp", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_15ODES_even.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even_1000iters", dox = "ND", skip_pass = TRUE)


model_qc(data = Rescale_ND_singleparameter_15ODES_even_1000iters, name = "Rescale_ND_singleparameter_15ODES_even_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))



load("Rescale_ND_singleparameter_15ODES_even_multiple.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_multiple, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even_multiple", dox = "ND", skip_pass = TRUE)


model_qc(data = Rescale_ND_singleparameter_15ODES_even_multiple, name = "Rescale_ND_singleparameter_15ODES_even_multiple",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_15ODES_even_IC.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_IC, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even_IC", dox = "ND", skip_pass = TRUE)


model_qc(data = Rescale_ND_singleparameter_15ODES_even_IC, name = "Rescale_ND_singleparameter_15ODES_even_IC",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))



load("Rescale_ND_singleparameter_15ODES_even_syn.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_syn, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even_syn", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_15ODES_even_Bliss.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_Bliss, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even_Bliss", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_15ODES_even_Loewe.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_Loewe, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even_Loewe", dox = "ND", skip_pass = TRUE)

load("Asymp_Rescale_ND_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_Rescale_ND_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_Rescale_ND_singleparameter_15ODES_even_1000iters", dox = "ND", skip_pass = TRUE)


model_qc(data = Asymp_Rescale_ND_singleparameter_15ODES_even_1000iters, name = "Asymp_Rescale_ND_singleparameter_15ODES_even_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))




load("Rescale_ND_singleparameter_15ODES_even_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even_ful_G2", dox = "ND", skip_pass = TRUE)



load("Rescale_ND_singleparameter_15ODES_even_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even_ful_G2_temp", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_15ODES_even_ful_G2_temp, name = "Rescale_ND_singleparameter_15ODES_even_ful_G2_temp",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Rescale_ND_singleparameter_15ODES_even_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_15ODES_even_ful_G2", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_15ODES_even_ful_G2, name = "Rescale_ND_singleparameter_15ODES_even_ful_G2",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_16ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)


model_qc(data = Rescale_ND_singleparameter_16ODES_G2_1000iters, name = "Rescale_ND_singleparameter_16ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))




load("Asymp_Rescale_ND_singleparameter_16ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_Rescale_ND_singleparameter_16ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_Rescale_ND_singleparameter_16ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)


model_qc(data = Asymp_Rescale_ND_singleparameter_16ODES_G2_1000iters, name = "Asymp_Rescale_ND_singleparameter_16ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Asymp_Rescale_ND_singleparameter_17ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_Rescale_ND_singleparameter_17ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_Rescale_ND_singleparameter_17ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)



load("Asymp_Rescale_ND_singleparameter_18ODES_G2_1000iters.RData")

load("Asymp_Rescale_ND_singleparameter_19ODES_G2_1000iters.RData")


load("Asymp_Rescale_ND_singleparameter_20ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_Rescale_ND_singleparameter_20ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_Rescale_ND_singleparameter_20ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)

model_qc(data = Asymp_Rescale_ND_singleparameter_20ODES_G2_1000iters, name = "Asymp_Rescale_ND_singleparameter_20ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))



load("Asymp_Rescale_ND_singleparameter_15ODES_even_ful_G2.RData")

load("Asymp_Rescale_ND_singleparameter_16ODES_G2_ful_G2.RData")

load("Asymp_Rescale_ND_singleparameter_17ODES_G2_ful_G2.RData")

load("Asymp_Rescale_ND_singleparameter_18ODES_G2_ful_G2.RData")

load("Asymp_Rescale_ND_singleparameter_19ODES_G2_ful_G2.RData")

load("Asymp_Rescale_ND_singleparameter_20ODES_G2_ful_G2.RData")


load("Rescale_ND_singleparameter_16ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_1000iters_3rates", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_17ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_17ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_17ODES_G2_1000iters_3rates", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_18ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_18ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_18ODES_G2_1000iters_3rates", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_19ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_19ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_19ODES_G2_1000iters_3rates", dox = "ND", skip_pass = TRUE)



load("Rescale_ND_singleparameter_20ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_20ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_20ODES_G2_1000iters_3rates", dox = "ND", skip_pass = TRUE)



load("Rescale_ND_singleparameter_21ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_21ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_21ODES_G2_1000iters_3rates", dox = "ND", skip_pass = TRUE)



load("Rescale_ND_singleparameter_22ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_22ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_22ODES_G2_1000iters_3rates", dox = "ND", skip_pass = TRUE)



load("Rescale_ND_singleparameter_23ODES_G2_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_23ODES_G2_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_23ODES_G2_1000iters_3rates", dox = "ND", skip_pass = TRUE)




load("Rescale_ND_singleparameter_16ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_ful_G2", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_16ODES_G2_ful_G2, name = "Rescale_ND_singleparameter_16ODES_G2_ful_G2",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Rescale_ND_singleparameter_16ODES_G2_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_16ODES_G2_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_16ODES_G2_ful_G2_temp", dox = "ND", skip_pass = TRUE)



load("Rescale_ND_singleparameter_17ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_17ODES_G2_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_17ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)


model_qc(data = Rescale_ND_singleparameter_17ODES_G2_1000iters, name = "Rescale_ND_singleparameter_17ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Asymp_Rescale_ND_singleparameter_16ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_Rescale_ND_singleparameter_16ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_Rescale_ND_singleparameter_16ODES_G2_ful_G2", dox = "ND", skip_pass = TRUE)



load("Asymp_Rescale_ND_singleparameter_17ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_Rescale_ND_singleparameter_17ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_Rescale_ND_singleparameter_17ODES_G2_ful_G2", dox = "ND", skip_pass = TRUE)




load("Rescale_ND_singleparameter_18ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_18ODES_G2_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_18ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)


model_qc(data = Rescale_ND_singleparameter_18ODES_G2_1000iters, name = "Rescale_ND_singleparameter_18ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))



load("Rescale_ND_singleparameter_18ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_18ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_18ODES_G2_ful_G2", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_18ODES_G2_ful_G2, name = "Rescale_ND_singleparameter_18ODES_G2_ful_G2",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max","b_F_2","c_F_2","alpha_max_2"),
         pars_acf1 = c("b_F_2","c_F_2","alpha_max_2","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Rescale_ND_singleparameter_18ODES_G2_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_18ODES_G2_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_18ODES_G2_ful_G2_temp", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp, stan_data = stan_input_include_longterm_data,
                                    name = "Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp", dox = "ND", skip_pass = TRUE)

load("cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_logF_Dox_1.RData")
load("cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_logF_ND_1.RData")


load("cell_cycle_ode_effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_1.RData.RData")


posterior_prediction_plots_no_abema(model_data = effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_1", dox = "Dox", skip_pass = TRUE)


load("cell_cycle_ode_effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_ful.RData")


posterior_prediction_plots_no_abema(model_data = effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_ful, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_ful", dox = "Dox", skip_pass = TRUE)


load("cell_cycle_ode_effectiveDose_multistage_singleRate_noLongtermData_palboOnly_ND_logP_logF_1.RData")


posterior_prediction_plots_no_abema(model_data = effectiveDose_multistage_singleRate_noLongtermData_palboOnly_ND_logP_logF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose_multistage_singleRate_noLongtermData_palboOnly_ND_logP_logF_1", dox = "ND", skip_pass = TRUE)


load("cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_logF_Dox_1.RData")

load("cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_logF_ND_1.RData")
model_qc(data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1, name = "effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","b_max"),
         pars_acf1 = c("alpha","b_max"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))



load("Weibull_test_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Weibull_test_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Weibull_test_singleparameter_15ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = Weibull_test_singleparameter_15ODES_even_1000iters, name = "Weibull_test_singleparameter_15ODES_even_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Weibull_Rescale_ND_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Weibull_Rescale_ND_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Weibull_Rescale_ND_singleparameter_15ODES_even_1000iters", dox = "ND", skip_pass = TRUE)


model_qc(data = Rescale_ND_singleparameter_15ODES_even_1000iters, name = "Rescale_ND_singleparameter_15ODES_even_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Weibull_test_singleparameter_16ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = Weibull_test_singleparameter_16ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Weibull_test_singleparameter_16ODES_G2", dox = "Dox", skip_pass = TRUE)

model_qc(data = Weibull_test_singleparameter_16ODES_G2, name = "Weibull_test_singleparameter_16ODES_G2",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Weibull_Rescale_ND_singleparameter_16ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Weibull_Rescale_ND_singleparameter_16ODES_G2_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Weibull_Rescale_ND_singleparameter_16ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)


model_qc(data = Weibull_Rescale_ND_singleparameter_16ODES_G2_1000iters, name = "Weibull_Rescale_ND_singleparameter_16ODES_G2_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))



load("Asymp_Rescale_ND_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_Rescale_ND_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_Rescale_ND_singleparameter_15ODES_even_1000iters", dox = "ND", skip_pass = TRUE)



load("Asymp_test_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_test_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_test_singleparameter_15ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)


load("Asymp_test_singleparameter_15ODES_even_ful_G2_1000iters.RData")

load("Asymp_test_singleparameter_16ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_test_singleparameter_16ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_test_singleparameter_16ODES_G2", dox = "Dox", skip_pass = TRUE)

load("Asymp_test_singleparameter_16ODES_G2_ful_G2.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_test_singleparameter_16ODES_G2_ful_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_test_singleparameter_16ODES_G2_ful_G2", dox = "Dox", skip_pass = TRUE)

load("Asymp_test_singleparameter_17ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_test_singleparameter_17ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_test_singleparameter_17ODES_G2", dox = "Dox", skip_pass = TRUE)

load("Asymp_test_singleparameter_17ODES_G2_ful_G2.RData")


load("Asymp_test_singleparameter_18ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = Asymp_test_singleparameter_18ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "Asymp_test_singleparameter_18ODES_G2", dox = "Dox", skip_pass = TRUE)


load("Asymp_test_singleparameter_18ODES_G2_ful_G2.RData")


load("Asymp_test_singleparameter_19ODES_G2_ful_G2.RData")

load("Asymp_test_singleparameter_19ODES_G2.RData")

load("Asymp_test_singleparameter_20ODES_G2_ful_G2.RData")


load("Asymp_test_singleparameter_20ODES_G2.RData")


load("test_singleparameter_3ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_3ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_3ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_3ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_3ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_3ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_6ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_6ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_6ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_6ODES_even_1000iters, name = "test_singleparameter_6ODES_even_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("test_singleparameter_6ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_6ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_6ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_6ODES_even_1000iters_3rates, name = "test_singleparameter_6ODES_even_1000iters_3rates",
         params_to_plot = c("sigma", "alpha", "gamma", "delta", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha", "gamma", "delta", "a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("test_singleparameter_9ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_9ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_9ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_9ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_9ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_9ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_12ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_12ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_12ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_12ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_12ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_12ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_18ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_18ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_18ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_18ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_18ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_18ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_21ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_21ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_21ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_21ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_21ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_21ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_24ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_24ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_24ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_24ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_24ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_24ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_27ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_27ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_27ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_27ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_27ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_27ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)


load("Rescale_ND_singleparameter_3ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_3ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_3ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_3ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_3ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_3ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_6ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_6ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_6ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_6ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_6ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_6ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_9ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_9ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_9ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_9ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_9ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_9ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_12ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_12ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_12ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_12ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_12ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_12ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_15ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_15ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_15ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_18ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_18ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_18ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_18ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_18ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_18ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_21ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_21ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_21ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_21ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_21ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_21ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_24ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_24ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_24ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_24ODES_even_1000iters_3rates.RData")



posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_24ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_24ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_24ODES_even_1000iters_3rates, name = "Rescale_ND_singleparameter_24ODES_even_1000iters_3rates",
         params_to_plot = c("sigma", "alpha", "gamma", "delta", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha", "gamma", "delta", "a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_27ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_27ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_27ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_27ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_27ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_27ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_30ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_30ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_30ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_30ODES_even_1000iters, name = "Rescale_ND_singleparameter_30ODES_even_1000iters",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("Rescale_ND_singleparameter_30ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_30ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_30ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_60ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_60ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_60ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_90ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_90ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_90ODES_even_1000iters_3rates", dox = "ND", skip_pass = TRUE)


load("test_singleparameter_3ODES_even_1000iters_noMax.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_3ODES_even_1000iters_noMax, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_3ODES_even_1000iters_noMax", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_3ODES_even_1000iters_noMax, name = "test_singleparameter_3ODES_even_1000iters_noMax",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("test_singleparameter_6ODES_even_1000iters_noMax.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_6ODES_even_1000iters_noMax, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_6ODES_even_1000iters_noMax", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_6ODES_even_1000iters_noMax, name = "test_singleparameter_6ODES_even_1000iters_noMax",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))


load("test_singleparameter_9ODES_even_1000iters_noMax.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_9ODES_even_1000iters_noMax, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_9ODES_even_1000iters_noMax", dox = "Dox", skip_pass = TRUE)

model_qc(data = test_singleparameter_9ODES_even_1000iters_noMax, name = "test_singleparameter_9ODES_even_1000iters_noMax",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP"),
         pars_acf1 = c("alpha","a_FP"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("test_singleparameter_12ODES_even_1000iters_noMax.RData")

load("test_singleparameter_12ODES_even_1000iters_noMax_test.RData")


load("test_singleparameter_12ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_12ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_12ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_18ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_18ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_18ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_21ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_21ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_21ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_24ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_24ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_24ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_27ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_27ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_27ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_30ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_30ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_30ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_30ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_30ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_30ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_60ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_60ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_60ODES_even_1000iters", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_60ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_60ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_60ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)


load("test_singleparameter_90ODES_even_1000iters_3rates.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_90ODES_even_1000iters_3rates, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_90ODES_even_1000iters_3rates", dox = "Dox", skip_pass = TRUE)


load("Rescale_ND_singleparameter_3ODES_even_1000iters_noMax.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_3ODES_even_1000iters_noMax, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_3ODES_even_1000iters_noMax", dox = "ND", skip_pass = TRUE)

model_qc(data = Rescale_ND_singleparameter_3ODES_even_1000iters_noMax, name = "Rescale_ND_singleparameter_3ODES_even_1000iters_noMax",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F", "a_FP"),
         pars_acf1 = c("alpha","b_P"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_6ODES_even_1000iters_noMax.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_6ODES_even_1000iters_noMax, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_6ODES_even_1000iters_noMax", dox = "ND", skip_pass = TRUE)


model_qc(data = Rescale_ND_singleparameter_6ODES_even_1000iters_noMax, name = "Rescale_ND_singleparameter_6ODES_even_1000iters_noMax",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F", "a_FP"),
         pars_acf1 = c("alpha","b_P"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_9ODES_even_1000iters_noMax.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_9ODES_even_1000iters_noMax, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_9ODES_even_1000iters_noMax", dox = "ND", skip_pass = TRUE)


model_qc(data = Rescale_ND_singleparameter_9ODES_even_1000iters_noMax, name = "Rescale_ND_singleparameter_9ODES_even_1000iters_noMax",
         params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F", "a_FP"),
         pars_acf1 = c("alpha","b_P"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))

load("Rescale_ND_singleparameter_12ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_12ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_12ODES_even_1000iters", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_15ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_15ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_15ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_18ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_18ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_18ODES_even_1000iters", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_21ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_21ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_21ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_24ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_24ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_24ODES_even_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_27ODES_even_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_27ODES_even_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_27ODES_even_1000iters", dox = "ND", skip_pass = TRUE)



load("test_singleparameter_21ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_21ODES_G2, stan_data = stan_input_include_longterm_data_2,
                                    name = "test_singleparameter_21ODES_G2", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_22ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_22ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_22ODES_G2", dox = "Dox", skip_pass = TRUE)

load("test_singleparameter_23ODES_G2.RData")

posterior_prediction_plots_no_abema(model_data = test_singleparameter_23ODES_G2, stan_data = stan_input_include_longterm_data,
                                    name = "test_singleparameter_23ODES_G2", dox = "Dox", skip_pass = TRUE)


load("Rescale_ND_singleparameter_21ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_21ODES_G2_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_21ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)

load("Rescale_ND_singleparameter_22ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_22ODES_G2_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_22ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)


load("Rescale_ND_singleparameter_23ODES_G2_1000iters.RData")

posterior_prediction_plots_no_abema(model_data = Rescale_ND_singleparameter_23ODES_G2_1000iters, stan_data = stan_input_include_longterm_data_2,
                                    name = "Rescale_ND_singleparameter_23ODES_G2_1000iters", dox = "ND", skip_pass = TRUE)


calculate_looic = function(stan_model){
  log_lik_cc <- extract_log_lik(stan_model,parameter_name = "log_lik_cc_p", merge_chains = FALSE)
  log_lik_ct <- extract_log_lik(stan_model,parameter_name = "log_lik_ct_p", merge_chains = FALSE)
  log_lik = abind(log_lik_cc,  log_lik_ct, along = 3)
  r_eff <- relative_eff(exp(log_lik), cores = 2)
  loo <- loo(log_lik, r_eff = r_eff, cores = 2)
  return(loo)
}


calculate_looic_2 = function(stan_model_1, stan_model_2){
  log_lik_cc_1 <- extract_log_lik(stan_model_1,parameter_name = "log_lik_cc_p", merge_chains = FALSE)
  log_lik_ct_1 <- extract_log_lik(stan_model_1,parameter_name = "log_lik_ct_p", merge_chains = FALSE)
  log_lik_cc_2 <- extract_log_lik(stan_model_2,parameter_name = "log_lik_cc_p", merge_chains = FALSE)
  log_lik_ct_2 <- extract_log_lik(stan_model_2,parameter_name = "log_lik_ct_p", merge_chains = FALSE)
  log_lik = abind(log_lik_cc_1,  log_lik_ct_1, log_lik_cc_2, log_lik_ct_2,  along = 3)
  r_eff <- relative_eff(exp(log_lik), cores = 2)
  loo <- loo(log_lik, r_eff = r_eff, cores = 2)
  return(loo)
}

looic_16_Dox_modified_logP_logF = calculate_looic(test_singleparameter_16ODES_G2_modified)
print(looic_16_Dox_modified_logP_logF)

looic_16_Dox_nonsyn_modified_logP_logF = calculate_looic(test_singleparameter_16ODES_G2_both_G2_nonsyn_modified)
print(looic_16_Dox_nonsyn_modified_logP_logF)

looic_16_ND_modified_logP_logF = calculate_looic(Rescale_ND_singleparameter_16ODES_G2_modified)
print(looic_16_ND_modified_logP_logF)

looic_16_ND_nonsyn_modified_logP_logF = calculate_looic(Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_modified)
print(looic_16_ND_nonsyn_modified_logP_logF)



looic_dox_ND_logP_logF = calculate_looic_2(effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1, effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1)
print(looic_dox_ND_logP_logF)

looic_dox_ND_15_logP_logF = calculate_looic_2(test_singleparameter_15ODES_even_1000iters, Rescale_ND_singleparameter_15ODES_even_1000iters)
print(looic_dox_ND_15_logP_logF)



looic_dox_ND_16_logP_logF = calculate_looic_2(test_singleparameter_16ODES_G2_syn_modified, Rescale_ND_singleparameter_16ODES_G2_syn_modified)
print(looic_dox_ND_16_logP_logF)

looic_dox_ND_17_logP_logF = calculate_looic_2(test_singleparameter_17ODES_G2_1000iters, Rescale_ND_singleparameter_17ODES_G2_1000iters)
print(looic_dox_ND_17_logP_logF)

looic_dox_ND_18_logP_logF = calculate_looic_2(test_singleparameter_18ODES_G2_1000iters, Rescale_ND_singleparameter_18ODES_G2_1000iters)
print(looic_dox_ND_18_logP_logF)

looic_dox_ND_19_logP_logF = calculate_looic_2(test_singleparameter_19ODES_G2_1000iters, Rescale_ND_singleparameter_19ODES_G2_1000iters)
print(looic_dox_ND_19_logP_logF)

looic_dox_ND_20_logP_logF = calculate_looic_2(test_singleparameter_20ODES_G2_1000iters, Rescale_ND_singleparameter_20ODES_G2_1000iters)
print(looic_dox_ND_20_logP_logF)


looic_dox_ND_15_ful_logP_logF = calculate_looic_2(test_singleparameter_15ODES_even_ful_G2_1000iters, Rescale_ND_singleparameter_15ODES_even_ful_G2)
print(looic_dox_ND_15_ful_logP_logF)

looic_dox_ND_16_ful_logP_logF = calculate_looic_2(test_singleparameter_16ODES_G2_ful_G2_1000iters, Rescale_ND_singleparameter_16ODES_G2_ful_G2)
print(looic_dox_ND_16_ful_logP_logF)

looic_dox_ND_17_ful_logP_logF = calculate_looic_2(test_singleparameter_17ODES_G2_ful_G2_1000iters, Rescale_ND_singleparameter_17ODES_G2_ful_G2)
print(looic_dox_ND_17_ful_logP_logF)

looic_dox_ND_18_ful_logP_logF = calculate_looic_2(test_singleparameter_18ODES_G2_ful_G2_1000iters, Rescale_ND_singleparameter_18ODES_G2_ful_G2)
print(looic_dox_ND_18_ful_logP_logF)

looic_dox_ND_19_ful_logP_logF = calculate_looic_2(test_singleparameter_19ODES_G2_ful_G2_1000iters, Rescale_ND_singleparameter_19ODES_G2_ful_G2)
print(looic_dox_ND_19_ful_logP_logF)

looic_dox_ND_20_ful_logP_logF = calculate_looic_2(test_singleparameter_20ODES_G2_ful_G2_1000iters, Rescale_ND_singleparameter_20ODES_G2_ful_G2)
print(looic_dox_ND_20_ful_logP_logF)

looic_dox_logP_logF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1)
print(looic_dox_logP_logF)

looic_ND_logP_logF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1)
print(looic_ND_logP_logF)

looic_test_logP_logF = calculate_looic(effectiveDose_multistage_singleRate_noLongtermData_palboOnly_Dox_logP_logF_1)
print(looic_test_logP_logF)

looic_15_syn_logP_logF = calculate_looic(test_singleparameter_15ODES_even_syn)
print(looic_15_syn_logP_logF)

looic_15_Bliss_logP_logF = calculate_looic(test_singleparameter_15ODES_even_Bliss)
print(looic_15_Bliss_logP_logF)

looic_15_Loewe_logP_logF = calculate_looic(test_singleparameter_15ODES_even_Loewe)
print(looic_15_Loewe_logP_logF)


looic_15_logP_logF = calculate_looic(test_singleparameter_15ODES_even_1000iters)
print(looic_15_logP_logF)

looic_15_asymp_logP_logF = calculate_looic(Asymp_test_singleparameter_15ODES_even_1000iters)
print(looic_15_asymp_logP_logF)

looic_15_ful_logP_logF = calculate_looic(test_singleparameter_15ODES_even_ful_G2_1000iters)
print(looic_15_ful_logP_logF)

looic_15_ND_1000_logP_logF = calculate_looic(ND_singleparameter_15ODES_even_1000iters)
print(looic_15_ND_1000_logP_logF)

looic_15_rescale_ND_logP_logF = calculate_looic(Rescale_ND_singleparameter_15ODES_even_1000iters)
print(looic_15_rescale_ND_logP_logF)

looic_15_rescale_ND_IC_logP_logF = calculate_looic(Rescale_ND_singleparameter_15ODES_even_IC)
print(looic_15_rescale_ND_IC_logP_logF)

looic_15_rescale_ND_ful_logP_logF = calculate_looic(Rescale_ND_singleparameter_15ODES_even_ful_G2)
print(looic_15_rescale_ND_ful_logP_logF)

looic_15_rescale_ND_ful_temp_logP_logF = calculate_looic(Rescale_ND_singleparameter_15ODES_even_ful_G2_temp)
print(looic_15_rescale_ND_ful_temp_logP_logF)

looic_15_ND_ful_logP_logF = calculate_looic(ND_singleparameter_15ODES_even_ful_G2)
print(looic_15_ND_ful_logP_logF)


looic_15_multi_logP_logF = calculate_looic(test_singleparameter_15ODES_even_multiple)
print(looic_15_multi_logP_logF)

looic_15_multi_ND_logP_logF = calculate_looic(Rescale_ND_singleparameter_15ODES_even_multiple)
print(looic_15_multi_ND_logP_logF)

looic_15_IC_logP_logF = calculate_looic(test_singleparameter_15ODES_even_IC)
print(looic_15_IC_logP_logF)

looic_15_IC_rescale_ND_logP_logF = calculate_looic(Rescale_ND_singleparameter_15ODES_even_IC)
print(looic_15_IC_rescale_ND_logP_logF)


looic_16_logP_logF = calculate_looic(test_singleparameter_16ODES_G2_1000iters)
print(looic_16_logP_logF)


looic_16_ND_logP_logF = calculate_looic(ND_singleparameter_16ODES_G2)
print(looic_16_ND_logP_logF)

looic_16_ND_1000_logP_logF = calculate_looic(ND_singleparameter_16ODES_G2_1000iters)
print(looic_16_ND_1000_logP_logF)


looic_16_ful_logP_logF = calculate_looic(test_singleparameter_16ODES_G2_ful_G2_1000iters)
print(looic_16_ful_logP_logF)

looic_16_rescale_ND_logP_logF = calculate_looic(Rescale_ND_singleparameter_16ODES_G2_1000iters)
print(looic_16_rescale_ND_logP_logF)

looic_16_rescale_ND_ful_logP_logF = calculate_looic(Rescale_ND_singleparameter_16ODES_G2_ful_G2)
print(looic_16_rescale_ND_ful_logP_logF)


looic_16_asymp_ND_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_16ODES_G2_1000iters)
print(looic_16_asymp_ND_logP_logF)




looic_17_logP_logF = calculate_looic(test_singleparameter_17ODES_G2_1000iters)
print(looic_17_logP_logF)

looic_17_ND_1000_logP_logF = calculate_looic(ND_singleparameter_17ODES_G2_1000iters)
print(looic_17_ND_1000_logP_logF)

looic_17_rescale_ND_logP_logF = calculate_looic(Rescale_ND_singleparameter_17ODES_G2_1000iters)
print(looic_17_rescale_ND_logP_logF)

looic_17_ful_logP_logF = calculate_looic(test_singleparameter_17ODES_G2_ful_G2_1000iters)
print(looic_17_ful_logP_logF)

looic_17_rescale_ND_ful_logP_logF = calculate_looic(Rescale_ND_singleparameter_17ODES_G2_ful_G2)
print(looic_17_rescale_ND_ful_logP_logF)

looic_17_asymp_ND_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_17ODES_G2_1000iters)
print(looic_17_asymp_ND_logP_logF)

looic_18_logP_logF = calculate_looic(test_singleparameter_18ODES_G2_1000iters)
print(looic_18_logP_logF)

looic_18_ND_1000_logP_logF = calculate_looic(ND_singleparameter_18ODES_G2_1000iters)
print(looic_18_ND_1000_logP_logF)

looic_18_ful_logP_logF = calculate_looic(test_singleparameter_18ODES_G2_ful_G2_1000iters)
print(looic_18_ful_logP_logF)

looic_18_rescale_ND_logP_logF = calculate_looic(Rescale_ND_singleparameter_18ODES_G2_1000iters)
print(looic_18_rescale_ND_logP_logF)

looic_18_rescale_ND_ful_logP_logF = calculate_looic(Rescale_ND_singleparameter_18ODES_G2_ful_G2)
print(looic_18_rescale_ND_ful_logP_logF)

looic_18_rescale_ND_ful_temp_logP_logF = calculate_looic(Rescale_ND_singleparameter_18ODES_G2_ful_G2_temp)
print(looic_18_rescale_ND_ful_temp_logP_logF)

looic_18_asymp_ND_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_18ODES_G2_1000iters)
print(looic_18_asymp_ND_logP_logF)


looic_19_logP_logF = calculate_looic(test_singleparameter_19ODES_G2_1000iters)
print(looic_19_logP_logF)

looic_19_ful_logP_logF = calculate_looic(test_singleparameter_19ODES_G2_ful_G2_1000iters)
print(looic_19_ful_logP_logF)

looic_19_ND_1000_logP_logF = calculate_looic(ND_singleparameter_19ODES_G2_1000iters)
print(looic_19_ND_1000_logP_logF)

looic_19_rescale_ND_logP_logF = calculate_looic(Rescale_ND_singleparameter_19ODES_G2_1000iters)
print(looic_19_rescale_ND_logP_logF)

looic_19_rescale_ND_ful_logP_logF = calculate_looic(Rescale_ND_singleparameter_19ODES_G2_ful_G2)
print(looic_19_rescale_ND_ful_logP_logF)

looic_19_rescale_ND_ful_temp_logP_logF = calculate_looic(Rescale_ND_singleparameter_19ODES_G2_ful_G2_temp)
print(looic_19_rescale_ND_ful_temp_logP_logF)

looic_19_asymp_ND_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_19ODES_G2_1000iters)
print(looic_19_asymp_ND_logP_logF)

looic_20_logP_logF = calculate_looic(test_singleparameter_20ODES_G2_1000iters)
print(looic_20_logP_logF)

looic_20_ful_logP_logF = calculate_looic(test_singleparameter_20ODES_G2_ful_G2_1000iters)
print(looic_20_ful_logP_logF)

looic_20_ND_1000_logP_logF = calculate_looic(ND_singleparameter_20ODES_G2_1000iters)
print(looic_20_ND_1000_logP_logF)

looic_20_rescale_ND_logP_logF = calculate_looic(Rescale_ND_singleparameter_20ODES_G2_1000iters)
print(looic_20_rescale_ND_logP_logF)

looic_20_rescale_ND_ful_logP_logF = calculate_looic(Rescale_ND_singleparameter_20ODES_G2_ful_G2_temp)
print(looic_20_rescale_ND_ful_logP_logF)


looic_15_asymp_ND_ful_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_15ODES_even_ful_G2)
print(looic_15_asymp_ND_ful_logP_logF)

looic_16_asymp_ND_ful_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_16ODES_G2_ful_G2)
print(looic_16_asymp_ND_ful_logP_logF)

looic_17_asymp_ND_ful_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_17ODES_G2_ful_G2)
print(looic_17_asymp_ND_ful_logP_logF)

looic_18_asymp_ND_ful_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_18ODES_G2_ful_G2)
print(looic_18_asymp_ND_ful_logP_logF)


looic_19_asymp_ND_ful_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_19ODES_G2_ful_G2)
print(looic_19_asymp_ND_ful_logP_logF)


looic_20_asymp_ND_ful_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_20ODES_G2_ful_G2)
print(looic_20_asymp_ND_ful_logP_logF)




looic_15_rescale_ND_test_logP_logF = calculate_looic(Asymp_Rescale_ND_singleparameter_15ODES_even_1000iters)
print(looic_15_rescale_ND_test_logP_logF)

looic_15_a_test_logP_logF = calculate_looic(Asymp_test_singleparameter_15ODES_even_1000iters)
print(looic_15_a_test_logP_logF)

looic_15_a_test_ful_logP_logF = calculate_looic(Asymp_test_singleparameter_15ODES_even_ful_G2_1000iters)
print(looic_15_a_test_ful_logP_logF)

looic_16_a_test_logP_logF = calculate_looic(Asymp_test_singleparameter_16ODES_G2)
print(looic_16_a_test_logP_logF)

looic_16_a_test_ful_logP_logF = calculate_looic(Asymp_test_singleparameter_16ODES_G2_ful_G2)
print(looic_16_a_test_ful_logP_logF)

looic_17_a_test_logP_logF = calculate_looic(Asymp_test_singleparameter_17ODES_G2)
print(looic_17_a_test_logP_logF)

looic_17_a_test_ful_logP_logF = calculate_looic(Asymp_test_singleparameter_17ODES_G2_ful_G2)
print(looic_17_a_test_ful_logP_logF)

looic_18_a_test_logP_logF = calculate_looic(Asymp_test_singleparameter_18ODES_G2)
print(looic_18_a_test_logP_logF)


looic_18_a_test_ful_logP_logF = calculate_looic(Asymp_test_singleparameter_18ODES_G2_ful_G2)
print(looic_18_a_test_ful_logP_logF)

looic_19_a_test_logP_logF = calculate_looic(Asymp_test_singleparameter_19ODES_G2)
print(looic_19_a_test_logP_logF)

looic_19_a_test_ful_logP_logF = calculate_looic(Asymp_test_singleparameter_19ODES_G2_ful_G2)
print(looic_19_a_test_ful_logP_logF)

looic_20_a_test_logP_logF = calculate_looic(Asymp_test_singleparameter_20ODES_G2)
print(looic_20_a_test_logP_logF)

looic_20_a_test_ful_logP_logF = calculate_looic(Asymp_test_singleparameter_20ODES_G2_ful_G2)
print(looic_20_a_test_ful_logP_logF)

looic_dox_ND_3_e_logP_logF_noMax = calculate_looic_2(test_singleparameter_3ODES_even_1000iters_noMax, Rescale_ND_singleparameter_3ODES_even_1000iters_noMax)
print(looic_dox_ND_3_e_logP_logF_noMax)

looic_dox_ND_6_e_logP_logF_noMax = calculate_looic_2(test_singleparameter_6ODES_even_1000iters_noMax, Rescale_ND_singleparameter_6ODES_even_1000iters_noMax)
print(looic_dox_ND_6_e_logP_logF_noMax)


looic_dox_ND_9_e_logP_logF_noMax = calculate_looic_2(test_singleparameter_9ODES_even_1000iters_noMax, Rescale_ND_singleparameter_9ODES_even_1000iters_noMax)
print(looic_dox_ND_9_e_logP_logF_noMax)



looic_dox_ND_6_e_logP_logF = calculate_looic_2(test_singleparameter_6ODES_even_1000iters, Rescale_ND_singleparameter_6ODES_even_1000iters)
print(looic_dox_ND_6_e_logP_logF)



looic_dox_ND_9_e_logP_logF = calculate_looic_2(test_singleparameter_9ODES_even_1000iters, Rescale_ND_singleparameter_9ODES_even_1000iters)
print(looic_dox_ND_9_e_logP_logF)



looic_dox_ND_12_e_logP_logF = calculate_looic_2(test_singleparameter_12ODES_even_1000iters, Rescale_ND_singleparameter_12ODES_even_1000iters)
print(looic_dox_ND_12_e_logP_logF)

looic_dox_ND_15_e_logP_logF = calculate_looic_2(test_singleparameter_15ODES_even_1000iters, Rescale_ND_singleparameter_15ODES_even_1000iters)
print(looic_dox_ND_15_e_logP_logF)


looic_dox_ND_18_e_logP_logF = calculate_looic_2(test_singleparameter_18ODES_even_1000iters, Rescale_ND_singleparameter_18ODES_even_1000iters)
print(looic_dox_ND_18_e_logP_logF)


looic_dox_ND_21_e_logP_logF = calculate_looic_2(test_singleparameter_21ODES_even_1000iters, Rescale_ND_singleparameter_21ODES_even_1000iters)
print(looic_dox_ND_21_e_logP_logF)


looic_dox_ND_24_e_logP_logF = calculate_looic_2(test_singleparameter_24ODES_even_1000iters, Rescale_ND_singleparameter_24ODES_even_1000iters)
print(looic_dox_ND_24_e_logP_logF)


looic_dox_ND_27_e_logP_logF = calculate_looic_2(test_singleparameter_27ODES_even_1000iters, Rescale_ND_singleparameter_27ODES_even_1000iters)
print(looic_dox_ND_27_e_logP_logF)

looic_dox_ND_30_e_logP_logF = calculate_looic_2(test_singleparameter_30ODES_even_1000iters, Rescale_ND_singleparameter_30ODES_even_1000iters)
print(looic_dox_ND_30_e_logP_logF)


looic_dox_ND_3_e_logP_logF_noMax = calculate_looic_2(test_singleparameter_3ODES_even_1000iters_noMax, Rescale_ND_singleparameter_3ODES_even_1000iters_noMax)
print(looic_dox_ND_3_e_logP_logF_noMax)

looic_dox_ND_6_e_logP_logF_noMax = calculate_looic_2(test_singleparameter_6ODES_even_1000iters_noMax, Rescale_ND_singleparameter_6ODES_even_1000iters_noMax)
print(looic_dox_ND_6_e_logP_logF_noMax)


looic_dox_ND_9_e_logP_logF_noMax = calculate_looic_2(test_singleparameter_9ODES_even_1000iters_noMax, Rescale_ND_singleparameter_9ODES_even_1000iters_noMax)
print(looic_dox_ND_9_e_logP_logF_noMax)


looic_ND_3_e_logP_logF = calculate_looic(Rescale_ND_singleparameter_3ODES_even_1000iters)
print(looic_ND_3_e_logP_logF)

looic_ND_3_e_logP_logF_3rates = calculate_looic(Rescale_ND_singleparameter_3ODES_even_1000iters_3rates)
print(looic_ND_3_e_logP_logF_3rates)

looic_ND_6_e_logP_logF = calculate_looic(Rescale_ND_singleparameter_6ODES_even_1000iters)
print(looic_ND_6_e_logP_logF)

looic_ND_6_e_logP_logF_3rates = calculate_looic(Rescale_ND_singleparameter_6ODES_even_1000iters_3rates)
print(looic_ND_6_e_logP_logF_3rates)


looic_ND_9_e_logP_logF = calculate_looic(Rescale_ND_singleparameter_9ODES_even_1000iters)
print(looic_ND_9_e_logP_logF)

looic_ND_9_e_logP_logF_3rates = calculate_looic(Rescale_ND_singleparameter_9ODES_even_1000iters_3rates)
print(looic_ND_9_e_logP_logF_3rates)


looic_ND_12_e_logP_logF = calculate_looic(Rescale_ND_singleparameter_12ODES_even_1000iters)
print(looic_ND_12_e_logP_logF)

looic_ND_12_e_logP_logF_3rates = calculate_looic(Rescale_ND_singleparameter_12ODES_even_1000iters_3rates)
print(looic_ND_12_e_logP_logF_3rates)

looic_ND_15_e_logP_logF = calculate_looic(Rescale_ND_singleparameter_15ODES_even_1000iters)
print(looic_ND_15_e_logP_logF)

looic_ND_15_e_logP_logF_3rates = calculate_looic(Rescale_ND_singleparameter_15ODES_even_1000iters_3rates)
print(looic_ND_15_e_logP_logF_3rates)

looic_ND_18_e_logP_logF = calculate_looic(Rescale_ND_singleparameter_18ODES_even_1000iters)
print(looic_ND_18_e_logP_logF)

looic_ND_18_e_logP_logF_3rates = calculate_looic(Rescale_ND_singleparameter_18ODES_even_1000iters_3rates)
print(looic_ND_18_e_logP_logF_3rates)


looic_ND_21_e_logP_logF = calculate_looic(Rescale_ND_singleparameter_21ODES_even_1000iters)
print(looic_ND_21_e_logP_logF)

looic_ND_21_e_logP_logF_3rates = calculate_looic(Rescale_ND_singleparameter_21ODES_even_1000iters_3rates)
print(looic_ND_21_e_logP_logF_3rates)

looic_ND_24_e_logP_logF = calculate_looic(Rescale_ND_singleparameter_24ODES_even_1000iters)
print(looic_ND_24_e_logP_logF)

looic_ND_24_e_logP_logF_3rates = calculate_looic(Rescale_ND_singleparameter_24ODES_even_1000iters_3rates)
print(looic_ND_24_e_logP_logF_3rates)

looic_ND_27_e_logP_logF = calculate_looic( Rescale_ND_singleparameter_27ODES_even_1000iters)
print(looic_ND_27_e_logP_logF)

looic_ND_27_e_logP_logF_3rates = calculate_looic(Rescale_ND_singleparameter_27ODES_even_1000iters_3rates)
print(looic_ND_27_e_logP_logF_3rates)

looic_ND_30_e_logP_logF = calculate_looic( Rescale_ND_singleparameter_30ODES_even_1000iters)
print(looic_ND_30_e_logP_logF)

looic_ND_30_e_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_30ODES_even_1000iters_3rates)
print(looic_ND_30_e_logP_logF_3rates)


looic_ND_60_e_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_60ODES_even_1000iters_3rates)
print(looic_ND_60_e_logP_logF_3rates)

looic_ND_90_e_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_90ODES_even_1000iters_3rates)
print(looic_ND_90_e_logP_logF_3rates)

looic_dox_3_e_logP_logF = calculate_looic(test_singleparameter_3ODES_even_1000iters)
print(looic_dox_3_e_logP_logF)

looic_dox_3_e_logP_logF_3rates = calculate_looic(test_singleparameter_3ODES_even_1000iters_3rates)
print(looic_dox_3_e_logP_logF_3rates)

looic_dox_6_e_logP_logF = calculate_looic(test_singleparameter_6ODES_even_1000iters)
print(looic_dox_6_e_logP_logF)

looic_dox_6_e_logP_logF_3rates = calculate_looic(test_singleparameter_6ODES_even_1000iters_3rates)
print(looic_dox_6_e_logP_logF_3rates)


looic_dox_9_e_logP_logF = calculate_looic(test_singleparameter_9ODES_even_1000iters)
print(looic_dox_9_e_logP_logF)

looic_dox_9_e_logP_logF_3rates = calculate_looic(test_singleparameter_9ODES_even_1000iters_3rates)
print(looic_dox_9_e_logP_logF_3rates)

looic_dox_12_e_logP_logF = calculate_looic(test_singleparameter_12ODES_even_1000iters)
print(looic_dox_12_e_logP_logF)

looic_dox_12_e_logP_logF_3rates = calculate_looic(test_singleparameter_12ODES_even_1000iters_3rates)
print(looic_dox_12_e_logP_logF_3rates)

looic_dox_15_e_logP_logF = calculate_looic(test_singleparameter_15ODES_even_1000iters)
print(looic_dox_15_e_logP_logF)

looic_dox_15_e_logP_logF_3rates = calculate_looic(test_singleparameter_15ODES_even_1000iters_3rates)
print(looic_dox_15_e_logP_logF_3rates)

looic_dox_18_e_logP_logF = calculate_looic(test_singleparameter_18ODES_even_1000iters)
print(looic_dox_18_e_logP_logF)

looic_dox_18_e_logP_logF_3rates = calculate_looic(test_singleparameter_18ODES_even_1000iters_3rates)
print(looic_dox_18_e_logP_logF_3rates)



looic_dox_21_e_logP_logF = calculate_looic(test_singleparameter_21ODES_even_1000iters)
print(looic_dox_21_e_logP_logF)

looic_dox_21_e_logP_logF_3rates = calculate_looic(test_singleparameter_21ODES_even_1000iters_3rates)
print(looic_dox_21_e_logP_logF_3rates)



looic_dox_24_e_logP_logF = calculate_looic(test_singleparameter_24ODES_even_1000iters)
print(looic_dox_24_e_logP_logF)

looic_dox_24_e_logP_logF_3rates = calculate_looic(test_singleparameter_24ODES_even_1000iters_3rates)
print(looic_dox_24_e_logP_logF_3rates)

looic_dox_27_e_logP_logF = calculate_looic(test_singleparameter_27ODES_even_1000iters)
print(looic_dox_27_e_logP_logF)

looic_dox_27_e_logP_logF_3rates = calculate_looic(test_singleparameter_27ODES_even_1000iters_3rates)
print(looic_dox_27_e_logP_logF_3rates)

looic_dox_30_e_logP_logF = calculate_looic(test_singleparameter_30ODES_even_1000iters)
print(looic_dox_30_e_logP_logF)

looic_dox_30_e_logP_logF_3rates = calculate_looic(test_singleparameter_30ODES_even_1000iters_3rates)
print(looic_dox_30_e_logP_logF_3rates)

looic_dox_60_e_logP_logF = calculate_looic(test_singleparameter_60ODES_even_1000iters)
print(looic_dox_60_e_logP_logF)

looic_dox_60_e_logP_logF_3rates = calculate_looic(test_singleparameter_60ODES_even_1000iters_3rates)
print(looic_dox_60_e_logP_logF_3rates)

looic_dox_90_e_logP_logF_3rates = calculate_looic(test_singleparameter_90ODES_even_1000iters_3rates)
print(looic_dox_90_e_logP_logF_3rates)


looic_dox_ND_16_logP_logF = calculate_looic_2(test_singleparameter_16ODES_G2_syn_modified, Rescale_ND_singleparameter_16ODES_G2_syn_modified)
print(looic_dox_ND_16_logP_logF)

looic_dox_ND_17_g2_logP_logF = calculate_looic_2(test_singleparameter_17ODES_G2_newdata, Rescale_ND_singleparameter_17ODES_G2_1000iters)
print(looic_dox_ND_17_g2_logP_logF)

looic_dox_ND_18_g2_logP_logF = calculate_looic_2(test_singleparameter_18ODES_G2_newdata, Rescale_ND_singleparameter_18ODES_G2_1000iters)
print(looic_dox_ND_18_g2_logP_logF)

looic_dox_ND_19_g2_logP_logF = calculate_looic_2(test_singleparameter_19ODES_G2_newdata, Rescale_ND_singleparameter_19ODES_G2_1000iters)
print(looic_dox_ND_19_g2_logP_logF)

looic_dox_ND_20_g2_logP_logF = calculate_looic_2(test_singleparameter_20ODES_G2_newdata, Rescale_ND_singleparameter_20ODES_G2_1000iters_newdata)
print(looic_dox_ND_20_g2_logP_logF)



looic_dox_ND_21_g2_logP_logF = calculate_looic_2(test_singleparameter_21ODES_G2, Rescale_ND_singleparameter_21ODES_G2_1000iters)
print(looic_dox_ND_21_g2_logP_logF)

looic_dox_ND_22_g2_logP_logF = calculate_looic_2(test_singleparameter_22ODES_G2_newdata, Rescale_ND_singleparameter_22ODES_G2_1000iters)
print(looic_dox_ND_22_g2_logP_logF)

looic_dox_ND_23_g2_logP_logF = calculate_looic_2(test_singleparameter_23ODES_G2_newdata, Rescale_ND_singleparameter_23ODES_G2_1000iters)
print(looic_dox_ND_23_g2_logP_logF)

looic_ND_16_logP_logF = calculate_looic( Rescale_ND_singleparameter_16ODES_G2_syn_modified)
print(looic_ND_16_logP_logF)

looic_ND_17_g2_logP_logF = calculate_looic( Rescale_ND_singleparameter_17ODES_G2_1000iters)
print(looic_ND_17_g2_logP_logF)

looic_ND_18_g2_logP_logF = calculate_looic( Rescale_ND_singleparameter_18ODES_G2_1000iters)
print(looic_ND_18_g2_logP_logF)

looic_ND_19_g2_logP_logF = calculate_looic( Rescale_ND_singleparameter_19ODES_G2_1000iters)
print(looic_ND_19_g2_logP_logF)

looic_ND_20_g2_logP_logF = calculate_looic( Rescale_ND_singleparameter_20ODES_G2_1000iters_newdata)
print(looic_ND_20_g2_logP_logF)



looic_ND_21_g2_logP_logF = calculate_looic( Rescale_ND_singleparameter_21ODES_G2_1000iters)
print(looic_ND_21_g2_logP_logF)

looic_ND_22_g2_logP_logF = calculate_looic( Rescale_ND_singleparameter_22ODES_G2_1000iters)
print(looic_ND_22_g2_logP_logF)

looic_ND_23_g2_logP_logF = calculate_looic( Rescale_ND_singleparameter_23ODES_G2_1000iters)
print(looic_ND_23_g2_logP_logF)


looic_ND_16_g2_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_16ODES_G2_1000iters_3rates)
print(looic_ND_16_g2_logP_logF_3rates)

looic_ND_17_g2_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_17ODES_G2_1000iters_3rates)
print(looic_ND_17_g2_logP_logF_3rates)

looic_ND_18_g2_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_18ODES_G2_1000iters_3rates)
print(looic_ND_18_g2_logP_logF_3rates)

looic_ND_19_g2_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_19ODES_G2_1000iters_3rates)
print(looic_ND_19_g2_logP_logF_3rates)

looic_ND_20_g2_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_20ODES_G2_1000iters_3rates)
print(looic_ND_20_g2_logP_logF_3rates)

looic_ND_21_g2_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_21ODES_G2_1000iters_3rates)
print(looic_ND_21_g2_logP_logF_3rates)

looic_ND_22_g2_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_22ODES_G2_1000iters_3rates)
print(looic_ND_22_g2_logP_logF_3rates)

looic_ND_23_g2_logP_logF_3rates = calculate_looic( Rescale_ND_singleparameter_23ODES_G2_1000iters_3rates)
print(looic_ND_23_g2_logP_logF_3rates)


looic_dox_16_logP_logF = calculate_looic(test_singleparameter_16ODES_G2_syn_modified)
print(looic_dox_16_logP_logF)

looic_dox_17_g2_logP_logF = calculate_looic(test_singleparameter_17ODES_G2_newdata)
print(looic_dox_17_g2_logP_logF)

looic_dox_18_g2_logP_logF = calculate_looic(test_singleparameter_18ODES_G2_newdata)
print(looic_dox_18_g2_logP_logF)

looic_dox_19_g2_logP_logF = calculate_looic(test_singleparameter_19ODES_G2_newdata)
print(looic_dox_19_g2_logP_logF)

looic_dox_20_g2_logP_logF = calculate_looic(test_singleparameter_20ODES_G2_newdata)
print(looic_dox_20_g2_logP_logF)


looic_dox_21_g2_logP_logF = calculate_looic(test_singleparameter_21ODES_G2)
print(looic_dox_21_g2_logP_logF)

looic_dox_22_g2_logP_logF = calculate_looic(test_singleparameter_22ODES_G2_newdata)
print(looic_dox_22_g2_logP_logF)

looic_dox_23_g2_logP_logF = calculate_looic(test_singleparameter_23ODES_G2_newdata)
print(looic_dox_23_g2_logP_logF)




looic_dox_16_g2_logP_logF_3rates = calculate_looic(test_singleparameter_16ODES_G2_1000iters_3rates)
print(looic_dox_16_g2_logP_logF_3rates)

looic_dox_17_g2_logP_logF_3rates = calculate_looic(test_singleparameter_17ODES_G2_1000iters_3rates)
print(looic_dox_17_g2_logP_logF_3rates)

looic_dox_18_g2_logP_logF_3rates = calculate_looic(test_singleparameter_18ODES_G2_1000iters_3rates)
print(looic_dox_18_g2_logP_logF_3rates)

looic_dox_19_g2_logP_logF_3rates = calculate_looic(test_singleparameter_19ODES_G2_1000iters_3rates)
print(looic_dox_19_g2_logP_logF_3rates)

looic_dox_20_g2_logP_logF_3rates = calculate_looic(test_singleparameter_20ODES_G2_1000iters_3rates)
print(looic_dox_20_g2_logP_logF_3rates)


looic_dox_21_g2_logP_logF_3rates = calculate_looic(test_singleparameter_21ODES_G2_1000iters_3rates)
print(looic_dox_21_g2_logP_logF_3rates)

looic_dox_22_g2_logP_logF_3rates = calculate_looic(test_singleparameter_22ODES_G2_1000iters_3rates)
print(looic_dox_22_g2_logP_logF_3rates)

looic_dox_23_g2_logP_logF_3rates = calculate_looic(test_singleparameter_23ODES_G2_1000iters_3rates)
print(looic_dox_23_g2_logP_logF_3rates)

diff <- loo_compare(looic_dox_ND_15_e_logP_logF, looic_dox_ND_16_logP_logF, looic_dox_ND_17_g2_logP_logF, looic_dox_ND_18_g2_logP_logF,looic_dox_ND_19_g2_logP_logF,looic_dox_ND_20_g2_logP_logF, looic_dox_ND_21_g2_logP_logF, looic_dox_ND_22_g2_logP_logF, looic_dox_ND_23_g2_logP_logF, looic_dox_ND_24_e_logP_logF, looic_dox_ND_27_e_logP_logF)
print(diff)


diff_g2_ND <- loo_compare(looic_ND_15_e_logP_logF_3rates, looic_ND_16_g2_logP_logF_3rates, looic_ND_17_g2_logP_logF_3rates, looic_ND_18_g2_logP_logF_3rates,looic_ND_19_g2_logP_logF_3rates,looic_ND_20_g2_logP_logF_3rates, looic_ND_21_g2_logP_logF_3rates, looic_ND_22_g2_logP_logF_3rates, looic_ND_23_g2_logP_logF_3rates, looic_ND_24_e_logP_logF_3rates)
print(diff_g2_ND)

diff_g2_Dox <- loo_compare(looic_dox_15_e_logP_logF_3rates, looic_dox_16_g2_logP_logF_3rates, looic_dox_17_g2_logP_logF_3rates, looic_dox_18_g2_logP_logF_3rates,looic_dox_19_g2_logP_logF_3rates,looic_dox_20_g2_logP_logF_3rates, looic_dox_21_g2_logP_logF_3rates, looic_dox_22_g2_logP_logF_3rates, looic_dox_23_g2_logP_logF_3rates, looic_dox_24_e_logP_logF_3rates)
print(diff_g2_Dox)


diff_even_logistic <- loo_compare(looic_logic_Dox_ND_newdata_logP_logF, looic_dox_ND_3_e_logP_logF, looic_dox_ND_6_e_logP_logF, looic_dox_ND_9_e_logP_logF, looic_dox_ND_12_e_logP_logF, looic_dox_ND_15_e_logP_logF, looic_dox_ND_18_e_logP_logF, looic_dox_ND_21_e_logP_logF, looic_dox_ND_24_e_logP_logF, looic_dox_ND_27_e_logP_logF)
print(diff_even_logistic)

diff_even_logistic_Dox <- loo_compare(looic_logic_Dox_newdata_logP_logF, looic_dox_3_e_logP_logF_3rates, looic_dox_6_e_logP_logF_3rates, looic_dox_9_e_logP_logF_3rates, looic_dox_12_e_logP_logF_3rates, looic_dox_15_e_logP_logF_3rates, looic_dox_18_e_logP_logF_3rates, looic_dox_21_e_logP_logF_3rates, looic_dox_24_e_logP_logF_3rates, looic_dox_27_e_logP_logF_3rates, looic_dox_30_e_logP_logF_3rates, looic_dox_60_e_logP_logF_3rates)
print(diff_even_logistic_Dox)

diff_even_logistic_ND <- loo_compare(looic_logic_ND_newdata_logP_logF, looic_ND_3_e_logP_logF_3rates, looic_ND_6_e_logP_logF_3rates, looic_ND_9_e_logP_logF_3rates, looic_ND_12_e_logP_logF_3rates, looic_ND_15_e_logP_logF_3rates, looic_ND_18_e_logP_logF_3rates, looic_ND_21_e_logP_logF_3rates, looic_ND_24_e_logP_logF_3rates, looic_ND_27_e_logP_logF_3rates, looic_ND_30_e_logP_logF_3rates, looic_ND_60_e_logP_logF_3rates)
print(diff_even_logistic_ND)

diff_even_Dox <- loo_compare(looic_dox_3_e_logP_logF_3rates, looic_dox_6_e_logP_logF_3rates, looic_dox_9_e_logP_logF_3rates, looic_dox_12_e_logP_logF_3rates, looic_dox_15_e_logP_logF_3rates, looic_dox_18_e_logP_logF_3rates, looic_dox_21_e_logP_logF_3rates, looic_dox_24_e_logP_logF_3rates, looic_dox_27_e_logP_logF_3rates, looic_dox_30_e_logP_logF_3rates, looic_dox_60_e_logP_logF_3rates)
print(diff_even_Dox)

diff_even_ND <- loo_compare(looic_ND_3_e_logP_logF_3rates, looic_ND_6_e_logP_logF_3rates, looic_ND_9_e_logP_logF_3rates, looic_ND_12_e_logP_logF_3rates, looic_ND_15_e_logP_logF_3rates, looic_ND_18_e_logP_logF_3rates, looic_ND_21_e_logP_logF_3rates, looic_ND_24_e_logP_logF_3rates, looic_ND_27_e_logP_logF_3rates, looic_ND_30_e_logP_logF_3rates, looic_ND_60_e_logP_logF_3rates)
print(diff_even_ND)


diff_even <- loo_compare(looic_dox_ND_3_e_logP_logF, looic_dox_ND_6_e_logP_logF, looic_dox_ND_9_e_logP_logF, looic_dox_ND_12_e_logP_logF, looic_dox_ND_15_e_logP_logF, looic_dox_ND_18_e_logP_logF, looic_dox_ND_21_e_logP_logF, looic_dox_ND_24_e_logP_logF, looic_dox_ND_27_e_logP_logF)
print(diff_even)

diff_total <- loo_compare(looic_dox_ND_3_e_logP_logF, looic_dox_ND_6_e_logP_logF, looic_dox_ND_9_e_logP_logF, looic_dox_ND_12_e_logP_logF, looic_dox_ND_16_logP_logF, looic_dox_ND_18_e_logP_logF, looic_dox_ND_21_g2_logP_logF, looic_dox_ND_21_e_logP_logF, looic_dox_ND_22_g2_logP_logF, looic_dox_ND_23_g2_logP_logF, looic_dox_ND_24_e_logP_logF)
print(diff_total)


df <- data.frame(elpd_diff = diff[,1],
                 se_diff = diff[,2],
                 Multistages=c("8:8:7 (23)", "8:8:8 (24)", "8:8:6 (22)", "8:8:5 (21)", "8:8:4 (20)", "8:8:3 (19)", "5:5:5 (15)", "8:8:2 (18)", "7:7:2 (16)", "7:7:3 (17)"))

df_even <- data.frame(elpd_diff = diff_even[,1],
                 se_diff = diff_even[,2],
                 Multistages=c( "8:8:8 (24)", "9:9:9 (27)",  "7:7:7 (21)", "6:6:6 (18)", "5:5:5 (15)", "4:4:4 (12)", "3:3:3 (9)", "2:2:2 (6)", "1:1:1 (3)"    ))

df_even_logistic <- data.frame(elpd_diff = diff_even_logistic[,1],
                      se_diff = diff_even_logistic[,2],
                      Multistages=c( "8:8:8 (24)", "9:9:9 (27)",  "7:7:7 (21)", "6:6:6 (18)", "5:5:5 (15)", "4:4:4 (12)", "3:3:3 (9)", "1:1:1 (logistic)", "2:2:2 (6)", "1:1:1 (3)"    ))

df_even_logistic_Dox <- data.frame(elpd_diff = diff_even_logistic_Dox[,1],
                               se_diff = diff_even_logistic_Dox[,2],
                               Multistages=c( "2:2:2 (6)", "3:3:3 (9)", "1:1:1 (3)", "4:4:4 (12)", "5:5:5 (15)", "6:6:6 (18)", "7:7:7 (21)", "8:8:8 (24)", "9:9:9 (27)", "10:10:10 (30)", "20:20:20 (60)", "1:1:1 (logistic)"))

df_even_ND <- data.frame(elpd_diff = diff_even_ND[,1],
                                   se_diff = diff_even_ND[,2],
                                  Multistages=c("8:8:8 (24)", "9:9:9 (27)", "7:7:7 (21)", "10:10:10 (30)", "6:6:6 (18)", "5:5:5 (15)", "4:4:4 (12)", "20:20:20 (60)", "3:3:3 (9)",  "2:2:2 (6)", "1:1:1 (3)"    ))

df_even_Dox <- data.frame(elpd_diff = diff_even_Dox[,1],
                                   se_diff = diff_even_Dox[,2],
                                   Multistages=c( "2:2:2 (6)", "3:3:3 (9)", "1:1:1 (3)", "4:4:4 (12)", "5:5:5 (15)", "6:6:6 (18)", "7:7:7 (21)", "8:8:8 (24)", "9:9:9 (27)", "10:10:10 (30)", "20:20:20 (60)"))

df_even_logistic_ND <- data.frame(elpd_diff = diff_even_logistic_ND[,1],
                                  se_diff = diff_even_logistic_ND[,2],
                                  Multistages=c("8:8:8 (24)", "9:9:9 (27)", "7:7:7 (21)", "10:10:10 (30)", "6:6:6 (18)", "5:5:5 (15)", "1:1:1 (logistic)", "4:4:4 (12)", "20:20:20 (60)", "3:3:3 (9)",  "2:2:2 (6)", "1:1:1 (3)"      ))

df_total <- data.frame(elpd_diff = diff_total[,1],
                      se_diff = diff_total[,2],
                      Multistages=c( "8:8:7 (23)", "8:8:8 (24)", "8:8:6 (22)", "7:7:7 (21)", "8:8:5 (21)",  "6:6:6 (18)", "7:7:2 (16)",  "4:4:4 (12)", "3:3:3 (9)", "2:2:2 (6)", "1:1:1 (3)"    ))

ggplot(df, aes(x=Multistages, y=elpd_diff)) + 
  geom_point()+
  geom_errorbar(aes(ymin=elpd_diff-se_diff, ymax=elpd_diff+se_diff), width=.2,
                position=position_dodge(0.05)) 

ggplot(df_even, aes(x=Multistages, y=elpd_diff)) + 
  geom_point()+
  geom_errorbar(aes(ymin=elpd_diff-se_diff, ymax=elpd_diff+se_diff), width=.2,
                position=position_dodge(0.05)) 

ggplot(df_total, aes(x=Multistages, y=elpd_diff)) + 
  geom_point()+
  geom_errorbar(aes(ymin=elpd_diff-se_diff, ymax=elpd_diff+se_diff), width=.2,
                position=position_dodge(0.05)) 

ggplot(df_even_logistic, aes(x=Multistages, y=elpd_diff)) + 
  geom_point()+
  geom_errorbar(aes(ymin=elpd_diff-se_diff, ymax=elpd_diff+se_diff), width=.2,
                position=position_dodge(0.05)) 

ggplot(df_even_logistic_Dox, aes(x=Multistages, y=elpd_diff)) + 
  geom_point()+
  geom_errorbar(aes(ymin=elpd_diff-se_diff, ymax=elpd_diff+se_diff), width=.2,
                position=position_dodge(0.05))

ggplot(df_even_logistic_ND, aes(x=Multistages, y=elpd_diff)) + 
  geom_point()+
  geom_errorbar(aes(ymin=elpd_diff-se_diff, ymax=elpd_diff+se_diff), width=.2,
                position=position_dodge(0.05))

Number_of_subphases <- factor(df_even_logistic_Dox$Multistages, level = c("1:1:1 (logistic)", "1:1:1 (3)", "2:2:2 (6)", "3:3:3 (9)", "4:4:4 (12)", "5:5:5 (15)", "6:6:6 (18)", "7:7:7 (21)", "8:8:8 (24)", "9:9:9 (27)", "10:10:10 (30)", "20:20:20 (60)") )
  

ggplot(df_even_logistic_Dox, aes(x=Number_of_subphases, y=elpd_diff)) + 
  geom_point()+
  geom_errorbar(aes(ymin=elpd_diff-se_diff, ymax=elpd_diff+se_diff), width=.2,
                position=position_dodge(0.05))

ggplot(df_even_logistic_ND, aes(x=Number_of_subphases, y=elpd_diff)) + 
  geom_point()+
  geom_errorbar(aes(ymin=elpd_diff-se_diff, ymax=elpd_diff+se_diff), width=.2,
                position=position_dodge(0.05))

calculate_waic_2 = function(stan_model_1, stan_model_2){
  log_lik_cc_1 <- extract_log_lik(stan_model_1,parameter_name = "log_lik_cc_p", merge_chains = FALSE)
  log_lik_ct_1 <- extract_log_lik(stan_model_1,parameter_name = "log_lik_ct_p", merge_chains = FALSE)
  log_lik_cc_2 <- extract_log_lik(stan_model_2,parameter_name = "log_lik_cc_p", merge_chains = FALSE)
  log_lik_ct_2 <- extract_log_lik(stan_model_2,parameter_name = "log_lik_ct_p", merge_chains = FALSE)
  log_lik = abind(log_lik_cc_1,  log_lik_ct_1, log_lik_cc_2, log_lik_ct_2,  along = 3)
  #r_eff <- relative_eff(exp(log_lik), cores = 2)
  w <- WAIC(log_lik, cores=2)
  return(w)
}

calculate_waic_2(test_singleparameter_23ODES_G2_newdata, Rescale_ND_singleparameter_23ODES_G2_1000iters)




source("posterior_prediction_functions_validation2.R")
source("post_processing_script.R")


validation_schedules_10day_96well_palbo_only_15_even(validation_fn = "validation_growth_study_12-2020.csv", 
                                             model_data = test_singleparameter_15ODES_even_test3_newdata,
                                             name = "test_singleparameter_15ODES_even_test3_newdata",
                                             params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max"),
                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)

validation_schedules_10day_96well_palbo_only_15_even_ND(validation_fn = "validation_growth_study_12-2020.csv", 
                                                     model_data = Rescale_ND_singleparameter_15ODES_even_newdata,
                                                     name = "Rescale_ND_singleparameter_15ODES_even_newdata",
                                                     params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max"),
                                                     num_iters = 1500, dox = "ND", passage_day = 5.1, passage_number = 3000)


validation_schedules_10day_96well_palbo_only_15_even_ful_G2(validation_fn = "validation_growth_study_12-2020.csv", 
                                                 model_data = test_singleparameter_15ODES_even_ful_G2_newdata,
                                                 name = "test_singleparameter_15ODES_even_ful_G2_newdata",
                                                 params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                 num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)

validation_schedules_10day_96well_palbo_only_15_even_ful_S(validation_fn = "validation_growth_study_12-2020.csv", 
                                                            model_data = test_singleparameter_15ODES_even_ful_S_newdata,
                                                            name = "test_singleparameter_15ODES_even_ful_S_newdata",
                                                            params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                            num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)

validation_schedules_10day_96well_palbo_only_14_G2(validation_fn = "validation_growth_study_12-2020.csv", 
                                                   model_data = test_singleparameter_14ODES_G2_newdata,
                                                   name = "test_singleparameter_14ODES_G2_newdata",
                                                   params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max"),
                                                   num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)



validation_schedules_10day_96well_palbo_only_14_G2_palbo(validation_fn = "validation_growth_study_12-2020.csv", 
                                                       model_data = test_singleparameter_14ODES_G2_palbo_G2_newdata,
                                                       name = "test_singleparameter_14ODES_G2_palbo_G2_newdata",
                                                       params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max", "b_P_2", "c_P_2", "alpha_max_2"),
                                                       num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)



validation_schedules_10day_96well_palbo_only_14_G2_ful(validation_fn = "validation_growth_study_12-2020.csv", 
                                                       model_data = test_singleparameter_14ODES_G2_ful_G2_newdata,
                                                       name = "test_singleparameter_14ODES_G2_ful_G2_newdata",
                                                       params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                       num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)


validation_schedules_10day_96well_palbo_only_14_G2_both_G2_nonsyn(validation_fn = "validation_growth_study_12-2020.csv", 
                                                                  model_data = test_singleparameter_14ODES_G2_both_G2_nonsyn_newdata,
                                                                  name = "test_singleparameter_14ODES_G2_both_G2_nonsyn_newdata",
                                                                  params = c("alpha","b_P","b_F","c_P","c_F", "alpha_max", "b_P_2", "b_F_2", "c_P_2", "c_F_2", "alpha_max_2"),
                                                                  num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)

validation_schedules_10day_96well_palbo_only_14_G2_both_G2_syn(validation_fn = "validation_growth_study_12-2020.csv", 
                                                                  model_data = test_singleparameter_14ODES_G2_both_G2_syn_newdata,
                                                                  name = "test_singleparameter_14ODES_G2_both_G2_syn_newdata",
                                                                  params = c("alpha","b_P","b_F","c_P","c_F", "a_FP", "a_FP_2", "alpha_max", "b_P_2", "b_F_2", "c_P_2", "c_F_2", "alpha_max_2"),
                                                                  num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)


validation_schedules_10day_96well_palbo_only_16_G2_ND(validation_fn = "validation_growth_study_12-2020.csv", 
                                                       model_data = Rescale_ND_singleparameter_16ODES_G2_1000iters_newdata,
                                                       name = "Rescale_ND_singleparameter_16ODES_G2_1000iters_newdata",
                                                       params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max"),
                                                       num_iters = 1500, dox = "ND", passage_day = 5.1, passage_number = 3000)


validation_schedules_10day_96well_palbo_only_16_G2_ful(validation_fn = "validation_growth_study_12-2020.csv", 
                                                       model_data = test_singleparameter_16ODES_G2_ful_G2_newdata,
                                                       name = "test_singleparameter_16ODES_G2_ful_G2_newdata",
                                                       params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                       num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)

validation_schedules_10day_96well_palbo_only_16_G2_both_G2_nonsyn(validation_fn = "validation_growth_study_12-2020.csv", 
                                                                  model_data = test_singleparameter_16ODES_G2_both_G2_nonsyn_newdata,
                                                                  name = "test_singleparameter_16ODES_G2_both_G2_nonsyn_newdata",
                                                                  params = c("alpha","b_P","b_F","c_P","c_F", "alpha_max", "b_P_2", "b_F_2", "c_P_2", "c_F_2", "alpha_max_2"),
                                                                  num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)


validation_schedules_10day_96well_palbo_only_16_G2_both_G2_nonsyn_ND(validation_fn = "validation_growth_study_12-2020.csv", 
                                                                  model_data = Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_newdata,
                                                                  name = "Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_newdata",
                                                                  params = c("alpha","b_P","b_F","c_P","c_F", "alpha_max", "b_P_2", "b_F_2", "c_P_2", "c_F_2", "alpha_max_2"),
                                                                  num_iters = 1500, dox = "ND", passage_day = 5.1, passage_number = 3000)

validation_schedules_10day_96well_palbo_only_16_G2_both_G2_nonsyn_ND(validation_fn = "validation_growth_study_12-2020.csv", 
                                                                     model_data = Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_newdata,
                                                                     name = "Rescale_ND_singleparameter_16ODES_G2_both_G2_nonsyn_newdata_2",
                                                                     params = c("alpha","b_P","b_F","c_P","c_F", "alpha_max", "b_P_2", "b_F_2", "c_P_2", "c_F_2", "alpha_max_2"),
                                                                     num_iters = 1500, dox = "ND", passage_day = 5.1, passage_number = 3000)




validation_schedules_10day_96well_palbo_only_mod(validation_fn = "validation_growth_study_12-2020.csv", 
                                             model_data = test_singleparameter_15ODES_even_1000iters,
                                             name = "test_singleparameter_15ODES_even_1000iters_mod",
                                             params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max"),
                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 2000)

validation_schedules_10day_96well_palbo_only_multiple(validation_fn = "validation_growth_study_12-2020.csv", 
                                             model_data = test_singleparameter_15ODES_even_multiple,
                                             name = "test_singleparameter_15ODES_even_multiple",
                                             params = c("alpha", "Beta", "Gamma", "b_P","b_F","c_P","c_F","a_FP", "alpha_max"),
                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000)


validation_schedules_10day_96well_palbo_only_IC(validation_fn = "validation_growth_study_12-2020.csv", 
                                             model_data = test_singleparameter_15ODES_even_IC,
                                             name = "test_singleparameter_15ODES_even_IC",
                                             params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max"),
                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000)

validation_schedules_10day_96well_palbo_only(validation_fn = "validation_growth_study_12-2020.csv", 
                                             model_data = test_singleparameter_15ODES_even_1000iters,
                                             name = "test_singleparameter_15ODES_even_1000iters",
                                             params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max"),
                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 3000)

validation_schedules_10day_96well_palbo_only_ful(validation_fn = "validation_growth_study_12-2020.csv", 
                                             model_data = test_singleparameter_15ODES_even_ful_G2_1000iters,
                                             name = "test_singleparameter_15ODES_even_ful_G2_1000iters",
                                             params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000)





validation_schedules_10day_96well_palbo_only_20_G2(validation_fn = "validation_growth_study_12-2020.csv", 
                                             model_data = test_singleparameter_20ODES_G2_1000iters,
                                             name = "test_singleparameter_20ODES_G2_1000iters",
                                             params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max"),
                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000)


validation_schedules_10day_96well_palbo_only_20_G2_ful(validation_fn = "validation_growth_study_12-2020.csv", 
                                                 model_data = test_singleparameter_20ODES_G2_ful_G2_1000iters,
                                                 name = "test_singleparameter_20ODES_G2_ful_G2_1000iters",
                                                 params = c("alpha","b_P","b_F","c_P","c_F","a_FP", "alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                 num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000)
