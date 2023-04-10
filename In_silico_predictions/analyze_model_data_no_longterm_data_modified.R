source("from_cluster_include_longterm_data/posterior_prediction_functions_validation2.R")
#load("stan_input_include_longterm_data.RData")
source("from_cluster_include_longterm_data/post_processing_script.R")

install.packages("ggpubr")

library("ggpubr")
library(abind)

library(usethis) 
usethis::edit_r_environ()

calculate_looic = function(stan_model){
  log_lik_cc <- extract_log_lik(stan_model,parameter_name = "log_lik_cc_p", merge_chains = FALSE)
  log_lik_ct <- extract_log_lik(stan_model,parameter_name = "log_lik_ct_p", merge_chains = FALSE)
  log_lik = abind(log_lik_cc,  log_lik_ct, along = 3)
  r_eff <- relative_eff(exp(log_lik), cores = 2)
  loo <- loo(log_lik, r_eff = r_eff, cores = 2)
  return(loo)
}

## log_log_ND
load("/Users/shaynastein/Google Drive/Harvard/MichorLab/BreastCancerCombination/code/BreastCancerModeling/from_cluster_include_longterm_data/model_fits/cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_logF_ND_1.RData")
model_qc(data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1, name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-ND-1",
         params_to_plot = c("sigma","alpha","b_P","b_F","c_P","c_F","a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf1 = c("a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))
posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-ND-1", dox = "ND", skip_pass = TRUE)

looic_ND_logP_logF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1)
print(looic_ND_logP_logF)

validation_schedules_24well_palbo_only(validation_fn = "../../combination_data/validation_1_24wells_include_tx_passaging_info.csv", 
                                                   model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1, 
                                                   name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-ND-1",
                                                   params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_max","b_min","Gamma","delta","K"),
                                                   num_iters = 1500, dox = "ND")

validation_schedules_24well_palbo_only(validation_fn = "../../combination_data/validation_1_24wells_include_tx_passaging_info.csv", 
                                                   model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1, 
                                                   name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-Dox-1",
                                                   params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_max","b_min","Gamma","delta","K"),
                                                   num_iters = 1500, dox = "Dox")

validation_schedules_10day_96well_palbo_only(validation_fn = "../../combination_data/validation_growth_study_12-2020.csv", 
                                             model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1,
                                             name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-Dox-1",
                                             params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_min","Gamma","delta","K"),
                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000)

validation_schedules_10day_96well_palbo_only(validation_fn = "../../combination_data/validation_growth_study_12-2020.csv", 
                                             model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1,
                                             name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-ND-1",
                                             params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_min","Gamma","delta","K"),
                                             num_iters = 1500, dox = "ND", passage_day = 5.1, passage_number = 5000)

validation_predictions_nopassage_ND = validation_schedules_10day_96well_palbo_only_nopassaging(validation_fn = "../../combination_data/validation_growth_study_12-2020.csv", 
                                             model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1,
                                             name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-ND-1",
                                             params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_min","Gamma","delta","K"),
                                             num_iters = 1500, dox = "ND")

validation_predictions_nopassage_Dox = validation_schedules_10day_96well_palbo_only_nopassaging(validation_fn = "../../combination_data/validation_growth_study_12-2020.csv", 
                                                                                               model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1,
                                                                                               name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-Dox-1",
                                                                                               params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_min","Gamma","delta","K"),
                                                                                               num_iters = 1500, dox = "Dox")

validation_predictions_nopassage = bind_rows(validation_predictions_nopassage_ND %>% mutate(Dox = "ND"), validation_predictions_nopassage_Dox %>% mutate(Dox = "Dox"))

pdf("from_cluster_include_longterm_data/plots/posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-no_passaging-10days-96wells.pdf", width = 10, height = 4)
validation_predictions_nopassage %>% ggplot(aes(x = time, y = total_median, color = Experiment)) + geom_line() + 
  geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL, fill = Experiment, color = NULL), alpha = 0.25) +
  facet_wrap(Dox, ncol = 2) + scale_y_log10() +  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") + theme(legend.position="bottom")
dev.off()

invitro_predictions_ND = fourweek_schedule_prediction_no_passaging_24well_palbo_only(validation_fn = "../../combination_data/validation_1_24wells_include_tx_passaging_info.csv", 
                                       model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1, 
                                       name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-ND-1",
                                       params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_max","b_min","Gamma","delta","K"),
                                       num_iters = 1500, dox = "ND")

invitro_predictions_Dox = fourweek_schedule_prediction_no_passaging_24well_palbo_only(validation_fn = "../../combination_data/validation_1_24wells_include_tx_passaging_info.csv", 
                                                                                     model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1, 
                                                                                     name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-Dox-1",
                                                                                     params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_max","b_min","Gamma","delta","K"),
                                                                                     num_iters = 1500, dox = "Dox")

p1<-invitro_predictions_ND %>% ggplot(aes(x = time, y = total_median, color = Experiment)) + geom_line() + 
  geom_ribbon(aes(ymin = total_low, ymax = total_high, fill = Experiment, color = NULL), alpha = 0.25) +
  scale_y_log10() + scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") + ylab("Predicted total cell count") +
  xlab("Time (days)") + ggtitle("Predicted 4 week in vitro schedules - No Dox")

p2<-invitro_predictions_Dox %>% ggplot(aes(x = time, y = total_median, color = Experiment)) + geom_line() + 
  geom_ribbon(aes(ymin = total_low, ymax = total_high, fill = Experiment, color = NULL), alpha = 0.25) +
  scale_y_log10() + scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") + ylab("Predicted total cell count") +
  xlab("Time (days)") + ggtitle("Predicted 4 week in vitro schedules - Dox")

pdf("from_cluster_include_longterm_data/plots/posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-no_passaging-24wells.pdf")
ggarrange(p1, p2, common.legend = T, legend = "bottom")
dev.off()

as.data.frame(effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1, 
              pars = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_min","b_max","Gamma","delta","K")) %>% 
  gather(factor_key = TRUE) %>% group_by(key) %>% summarize(median = median(value))

simulate_validation_schedules_24well_oneSample(validation_fn = "../../combination_data/validation_1_24wells_include_tx_passaging_info.csv",
                                               model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_logF_1,
                                               params = data.frame(alpha = 4.29, b_P = 0.992, b_F = .165, b_A = .7,c_P = 0.0122, c_F = .0001, c_A = .004,
                                                                   a_FP = -1, a_FA = 4, b = 5,Gamma = 1.28,delta = 0.06,K = 6),
                                               dox = "ND")



load("cell_cycle_ode_noLongtermData_palboOnly_logF_logP_insilico_trial_simulations_ct_6_FulOnly_3rates_100days_pk.Rdata")
load("cell_cycle_ode_noLongtermData_palboOnly_logF_logP_insilico_trial_simulations_ct_24_FulOnly_3rates_100days_pk.Rdata")
palbo_pk_df = rbind(simulation_results_24$palbo_3w_1w_125mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "3 weeks on, 1 week off, 125mg"), 
                    simulation_results_24$palbo_daily_100mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "daily, 100mg"),
                    simulation_results_24$palbo_daily_75mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "daily, 75mg"),
                    simulation_results_24$palbo_BID_50mg_50mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "BID, 50mg"),
                    simulation_results_24$palbo_BID_50mg_25mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "BID, 50mg in morning, 25mg at night"))

#palbo_pk_df = rbind(simulation_results$palbo_3w_1w_125mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "3 weeks on, 1 week off, 125mg"), 
#                    simulation_results$palbo_daily_100mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "daily, 100mg"),
#                    simulation_results$palbo_daily_75mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "daily, 75mg"),
#                    simulation_results$palbo_BID_50mg_50mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "BID, 50mg"),
#                    simulation_results$palbo_BID_50mg_25mg$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "BID, 50mg in morning, 25mg at night"),
#                    simulation_results$fulv_500mg_only$simulation_ND %>% filter(id == 1) %>% mutate(schedule = "fulv_500mg_only"))

pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_palbociclib-pk-insilico-trial_ct_combined_3rates_logistic_pk_test.pdf")
ggplot(palbo_pk_df, aes(time, PP_conc)) + 
  geom_line() + theme_bw() + facet_wrap(~ schedule, nrow = 5) + 
  ylab("Palbociclib concentration (ug/L)") + xlab("Time (days)")
dev.off()

test_palbo_100d = filter(palbo_pk_df, schedule =="daily, 100mg")
palbo_100d = mean(test_palbo_100d$PP_conc)

test_palbo_100b = filter(palbo_pk_df, schedule =="BID, 50mg")
palbo_100b = mean(test_palbo_100b$PP_conc)

test_palbo_75d = filter(palbo_pk_df, schedule =="daily, 75mg")
palbo_75d = mean(test_palbo_75d$PP_conc)

test_palbo_75b = filter(palbo_pk_df, schedule =="BID, 50mg in morning, 25mg at night")
palbo_75b = mean(test_palbo_75b$PP_conc)

ggplot(simulation_results_24$palbo_3w_1w_125mg$simulation_ND %>% filter(id == 1), 
       aes(time, FF_conc)) + geom_line() + theme_bw() + 
  ylab("Fulvestrant concentration (ug/L)") + xlab("Time (days)")

simulation_results_24$palbo_3w_1w_125mg$summary_ND = simulation_results_24$palbo_3w_1w_125mg$simulation_ND %>% 
  group_by(time) %>%
  summarize(total_low = quantile(total, 0.025, na.rm = T), total_median = median(total, na.rm = T), 
            total_high = quantile(total, 0.975, na.rm = T))

simulation_results_6$palbo_3w_1w_125mg$summary_Dox = simulation_results_6$palbo_3w_1w_125mg$simulation_Dox %>% 
  group_by(time) %>%
  summarize(total_low = quantile(total, 0.025, na.rm = T), total_median = median(total, na.rm = T), 
            total_high = quantile(total, 0.975, na.rm = T))



ggplot(simulation_results_24$palbo_3w_1w_125mg$summary_ND %>% filter(time <= 86), aes(time, total_median)) + 
  theme_bw() + scale_y_log10() + 
  geom_rect(xmin = 1, xmax = 21, ymin = 0, 
            ymax = max(simulation_results_24$palbo_3w_1w_125mg$summary_ND$total_median), alpha = 0.005, fill = "yellow") +
  geom_rect(xmin = 29, xmax = 49, ymin = 0, 
            ymax = max(simulation_results_24$palbo_3w_1w_125mg$summary_ND$total_median), alpha = 0.005, fill = "yellow") +
  geom_rect(xmin = 57, xmax = 77, ymin = 0, 
            ymax = max(simulation_results_24$palbo_3w_1w_125mg$summary_ND$total_median), alpha = 0.005, fill = "yellow") +
  geom_ribbon(aes(ymin = total_low, ymax = total_high), fill = "grey70") + geom_line() +
  ylab("Total cell count") + xlab("Time (days)") + ggtitle("3 weeks on, 1 week off, 125mg (DOX-)")

ggplot(simulation_results_6$palbo_3w_1w_125mg$summary_Dox %>% filter(time <= 86), aes(time, total_median)) + 
  theme_bw() + scale_y_log10() + 
  geom_rect(xmin = 1, xmax = 21, ymin = 0, 
            ymax = max(simulation_results_6$palbo_3w_1w_125mg$summary_Dox$total_median), alpha = 0.005, fill = "yellow") +
  geom_rect(xmin = 29, xmax = 49, ymin = 0, 
            ymax = max(simulation_results_6$palbo_3w_1w_125mg$summary_Dox$total_median), alpha = 0.005, fill = "yellow") +
  geom_rect(xmin = 57, xmax = 77, ymin = 0, 
            ymax = max(simulation_results_6$palbo_3w_1w_125mg$summary_Dox$total_median), alpha = 0.005, fill = "yellow") +
  geom_ribbon(aes(ymin = total_low, ymax = total_high), fill = "grey70") + geom_line() +
  ylab("Total cell count") + xlab("Time (days)") + ggtitle("3 weeks on, 1 week off, 125mg (DOX+)")

insilico_ND = rbind(simulation_results_24$palbo_3w_1w_125mg$simulation_ND  %>% mutate(schedule = "3 weeks on, 1 week off, 125mg"), 
                    simulation_results_24$palbo_daily_100mg$simulation_ND  %>% mutate(schedule = "daily, 100mg"),
                    simulation_results_24$palbo_daily_75mg$simulation_ND  %>% mutate(schedule = "daily, 75mg"),
                    simulation_results_24$palbo_BID_50mg_50mg$simulation_ND  %>% mutate(schedule = "BID, 50mg"),
                    simulation_results_24$palbo_BID_50mg_25mg$simulation_ND  %>% mutate(schedule = "BID, 50mg in morning, 25mg at night"))

insilico_Dox = rbind(simulation_results_6$palbo_3w_1w_125mg$simulation_Dox  %>% mutate(schedule = "3 weeks on, 1 week off, 125mg"), 
                     simulation_results_6$palbo_daily_100mg$simulation_Dox  %>% mutate(schedule = "daily, 100mg"),
                     simulation_results_6$palbo_daily_75mg$simulation_Dox  %>% mutate(schedule = "daily, 75mg"),
                     simulation_results_6$palbo_BID_50mg_50mg$simulation_Dox  %>% mutate(schedule = "BID, 50mg"),
                     simulation_results_6$palbo_BID_50mg_25mg$simulation_Dox  %>% mutate(schedule = "BID, 50mg in morning, 25mg at night"))

insilico_growth_ND = insilico_ND %>% group_by(time, schedule) %>% summarize(total_low = quantile(total, 0.025, na.rm = T),
                                                       total_median = median(total, na.rm = T),
                                                       total_high = quantile(total, 0.975, na.rm = T)) %>%
  ggplot(aes(x = time, y = total_median, color = schedule)) + geom_line() + scale_y_log10() +
  geom_ribbon(aes(ymin = total_low, ymax = total_high, fill = schedule, color = NULL), alpha = 0.15) + theme(aspect.ratio=4/3)+
  xlab("Time (days)") + ylab("Predicted total cell count")

insilico_growth_Dox = insilico_Dox %>% group_by(time, schedule) %>% summarize(total_low = quantile(total, 0.025, na.rm = T),
                                                       total_median = median(total, na.rm = T),
                                                       total_high = quantile(total, 0.975, na.rm = T)) %>%
  ggplot(aes(x = time, y = total_median, color = schedule)) + geom_line() + scale_y_log10() +
  geom_ribbon(aes(ymin = total_low, ymax = total_high, fill = schedule, color = NULL), alpha = 0.15) + theme(aspect.ratio=4/3)+
  xlab("Time (days)") + ylab("Predicted total cell count")



pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-cell-count-insilico-trial_ct_combined_3rates_logistic_pk_test.pdf")
ggarrange(insilico_growth_ND + ggtitle("-DOX"), insilico_growth_Dox + ggtitle("+DOX"),legend = "right", common.legend = T)+
  xlab("Dose (mg)") + ylab("Teeth length")
dev.off()

insilico_boxplot_ND = insilico_ND %>% group_by(schedule) %>% filter(time == max(time)) %>%
  ggplot(aes(x = schedule, y = total, color = schedule)) + geom_boxplot() + scale_y_log10() + theme(axis.text.x = element_blank()) +
  stat_compare_means(comparisons = list( c("daily, 100mg","BID, 50mg"),  c("daily, 75mg", "BID, 50mg in morning, 25mg at night"), c("daily, 100mg", "BID, 50mg in morning, 25mg at night")), method = "wilcox.test", size=2) + theme(aspect.ratio=4/3)+
 ylab("Total cell count at 100 days") + theme(plot.title = element_text(size =11), axis.text=element_text(size=11),
                                                        axis.title=element_text(size=11))
insilico_boxplot_Dox = insilico_Dox %>% group_by(schedule) %>% filter(time == max(time)) %>%
  ggplot(aes(x = schedule, y = total, color = schedule)) + geom_boxplot() + scale_y_log10() + theme(axis.text.x = element_blank()) +
  stat_compare_means(comparisons = list(c("daily, 100mg","BID, 50mg"),  c("daily, 75mg", "BID, 50mg in morning, 25mg at night"), c("daily, 100mg", "BID, 50mg in morning, 25mg at night")), method = "wilcox.test", size=2) + theme(aspect.ratio=4/3)+
  ylab("Total cell count at 100 days")+ theme(plot.title = element_text(size =11), axis.text=element_text(size=11),
                                                        axis.title=element_text(size=11))

pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-boxplot-insilico-trial_ct_combined_3rates_logistic_pk_test.pdf")
ggarrange(insilico_boxplot_ND + ggtitle("-DOX cells"), insilico_boxplot_Dox + ggtitle("+DOX cells"),legend = "right", common.legend = T)
dev.off()

# orginal code: group by schedule then sample 30 from 1500 simulations (different cohort of 30 patients)
waterfall_ND = insilico_ND  %>%  group_by(schedule) %>% filter(time == max(time), id %in% sample(1:1500, 30)) %>% unite(patient,id, schedule, sep = ": " , remove = F) %>%
  ggplot(aes(x = reorder(patient, -total), y = total, fill = schedule)) + geom_bar(stat="identity", width = 0.8) + scale_y_log10()+ theme(axis.text.x = element_blank())+
  xlab("Simulated patients") + ylab("Total cell count at 100 days")+ theme(aspect.ratio=2/5)  + theme(plot.title = element_text(size =14), axis.text=element_text(size=14),
                                                                                                                                        axis.title=element_text(size=14))

waterfall_Dox = insilico_Dox %>% group_by(schedule) %>% filter(time == max(time), id %in% sample(1:1500, 30))   %>% unite(patient,id, schedule, sep = ": " , remove = F) %>%
  ggplot(aes(x = reorder(patient, -total), y = total, fill = schedule)) + geom_bar(stat="identity", wdith = 0.8) + scale_y_log10() + theme(axis.text.x = element_blank())+
  xlab("Simulated patients") + ylab("Total cell count at 100 days") + theme(aspect.ratio=2/5) + theme(plot.title = element_text(size =14), axis.text=element_text(size=14),
                                                                                                      axis.title=element_text(size=14))




pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-waterfall-ND-insilico-trial_combined_3rates_logistic_pk_test.pdf")
ggarrange(waterfall_ND + ggtitle("-DOX cells"),legend = "right", common.legend = T, ncol = 1)
dev.off()

pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-waterfall-DOX-insilico-trial_combined_3rates_logistic_pk_test.pdf")
ggarrange(waterfall_Dox + ggtitle("+DOX cells"),legend = "right", common.legend = T, ncol = 1)
dev.off()

insilico_ND_2 = full_join(insilico_ND %>% filter(time %% 0.5 == 0), 
                          simulation_results_24$palbo_3w_1w_125mg$simulation_ND %>% filter(time %% 0.5 == 0) %>% dplyr::select(time, id, total),
                          by = c("time", "id"))

insilico_Dox_2 = full_join(insilico_Dox %>% filter(time %% 0.5 == 0), 
                           simulation_results_6$palbo_3w_1w_125mg$simulation_Dox %>% filter(time %% 0.5 == 0) %>% dplyr::select(time, id, total),
                           by = c("time", "id"))

insilico_ND_2 = insilico_ND_2 %>% mutate(percent_of_standard = total.x/total.y)
insilico_Dox_2 = insilico_Dox_2 %>% mutate(percent_of_standard = total.x/total.y)

## daily 100mg versus BID 50mg/25mg
#a = filter(insilico_Dox_2, schedule=="daily, 100mg", time == max(time))
#median(a$percent_of_standard)
#quantile(a$percent_of_standard, 0.95, na.rm = T)

b = filter(insilico_Dox_2, schedule=="BID, 50mg in morning, 25mg at night", time == max(time))
median(b$percent_of_standard)
quantile(b$percent_of_standard, 0.442, na.rm = T)


insilico_ND_2_sum = insilico_ND_2 %>% group_by(time, schedule) %>% 
  summarize(percent_low = quantile(percent_of_standard,0.05, na.rm = T),
            percent_median = median(percent_of_standard, na.rm = T),
            percent_high = quantile(percent_of_standard,0.95, na.rm = T))
insilico_Dox_2_sum = insilico_Dox_2 %>% group_by(time, schedule) %>% 
  summarize(percent_low = quantile(percent_of_standard,0.05, na.rm = T),
            percent_median = median(percent_of_standard, na.rm = T),
            percent_high = quantile(percent_of_standard,0.95, na.rm = T))

p1 <- ggplot(insilico_ND_2_sum, aes(time, percent_median, color = schedule)) + 
  geom_ribbon(aes(ymin = percent_low, ymax = percent_high, fill = schedule, color = NULL), alpha = 0.3) +
  geom_line() + theme_bw() + scale_y_log10() + xlab("Time (days)") + ylab("Percent of the standard schedule") +
  ggtitle("-DOX cells") + theme(plot.title = element_text(size =14), axis.text=element_text(size=14),
                                axis.title=element_text(size=14)) + theme(legend.text = element_text(size = 7))+ theme(aspect.ratio=5/4)
p2 <- ggplot(insilico_Dox_2_sum, aes(time, percent_median, color = schedule)) + 
  geom_ribbon(aes(ymin = percent_low, ymax = percent_high, fill = schedule, color = NULL), alpha = 0.3) +
  geom_line() + theme_bw() + scale_y_log10() + xlab("Time (days)") + ylab("Percent of the standard schedule") +
  ggtitle("+DOX cells") + theme(plot.title = element_text(size =14), axis.text=element_text(size=14),
                                axis.title=element_text(size=14))+ theme(legend.text = element_text(size = 7))+ theme(aspect.ratio=5/4)

pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-insilico-trial_ct_combined_3rates_logistic_pk_test.pdf")
ggarrange(p1,p2, ncol = 2, nrow = 1,common.legend = T, legend = "bottom")
dev.off()



## Fulvestrant_only

insilico_ND_FulOnly = rbind(simulation_results_24$palbo_3w_1w_125mg$simulation_ND  %>% mutate(schedule = "3 weeks on, 1 week off, 125mg"), 
                            simulation_results_24$palbo_daily_100mg$simulation_ND  %>% mutate(schedule = "daily, 100mg"),
                            simulation_results_24$palbo_daily_75mg$simulation_ND  %>% mutate(schedule = "daily, 75mg"),
                            simulation_results_24$palbo_BID_50mg_50mg$simulation_ND  %>% mutate(schedule = "BID, 50mg"),
                            simulation_results_24$palbo_BID_50mg_25mg$simulation_ND  %>% mutate(schedule = "BID, 50mg in morning, 25mg at night"),
                            simulation_results_24$fulv_500mg_only$simulation_ND  %>% mutate(schedule = "fulv_500mg_only") )

insilico_Dox_FulOnly = rbind(simulation_results_6$palbo_3w_1w_125mg$simulation_Dox  %>% mutate(schedule = "3 weeks on, 1 week off, 125mg"), 
                             simulation_results_6$palbo_daily_100mg$simulation_Dox  %>% mutate(schedule = "daily, 100mg"),
                             simulation_results_6$palbo_daily_75mg$simulation_Dox  %>% mutate(schedule = "daily, 75mg"),
                             simulation_results_6$palbo_BID_50mg_50mg$simulation_Dox  %>% mutate(schedule = "BID, 50mg"),
                             simulation_results_6$palbo_BID_50mg_25mg$simulation_Dox  %>% mutate(schedule = "BID, 50mg in morning, 25mg at night"),
                             simulation_results_6$fulv_500mg_only$simulation_Dox  %>% mutate(schedule = "fulv_500mg_only"))

insilico_growth_ND_FulOnly = insilico_ND_FulOnly %>% group_by(time, schedule) %>% summarize(total_low = quantile(total, 0.025, na.rm = T),
                                                                            total_median = median(total, na.rm = T),
                                                                            total_high = quantile(total, 0.975, na.rm = T)) %>%
  ggplot(aes(x = time, y = total_median, color = schedule)) + geom_line() + scale_y_log10() +
  geom_ribbon(aes(ymin = total_low, ymax = total_high, fill = schedule, color = NULL), alpha = 0.15)

insilico_growth_Dox_FulOnly = insilico_Dox_FulOnly %>% group_by(time, schedule) %>% summarize(total_low = quantile(total, 0.025, na.rm = T),
                                                                              total_median = median(total, na.rm = T),
                                                                              total_high = quantile(total, 0.975, na.rm = T)) %>%
  ggplot(aes(x = time, y = total_median, color = schedule)) + geom_line() + scale_y_log10() +
  geom_ribbon(aes(ymin = total_low, ymax = total_high, fill = schedule, color = NULL), alpha = 0.15)


pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-cell-count-insilico-trial_ct_combined_FulOnly_3rates_logistic_pk.pdf")
ggarrange(insilico_growth_ND_FulOnly + ggtitle("-DOX"), insilico_growth_Dox_FulOnly + ggtitle("+DOX"),legend = "right", common.legend = T)
dev.off()


insilico_boxplot_ND_FulOnly = insilico_ND_FulOnly %>% group_by(schedule) %>% filter(time == max(time)) %>%
  ggplot(aes(x = schedule, y = total, color = schedule)) + geom_boxplot() + scale_y_log10() + theme(axis.text.x = element_blank()) +
  stat_compare_means(comparisons = list(c("fulv_500mg_only","3 weeks on, 1 week off, 125mg"), c("daily, 100mg","3 weeks on, 1 week off, 125mg"), c("daily, 100mg","BID, 50mg"), c("daily, 75mg", "BID, 50mg in morning, 25mg at night"),
                                        c("daily, 100mg", "daily, 75mg"), c("BID, 50mg", "BID, 50mg in morning, 25mg at night")), method = "wilcox.test", size=2)
insilico_boxplot_Dox_FulOnly = insilico_Dox_FulOnly %>% group_by(schedule) %>% filter(time == max(time)) %>%
  ggplot(aes(x = schedule, y = total, color = schedule)) + geom_boxplot() + scale_y_log10() + theme(axis.text.x = element_blank()) +
  stat_compare_means(comparisons = list(c("fulv_500mg_only","3 weeks on, 1 week off, 125mg"), c("daily, 100mg","3 weeks on, 1 week off, 125mg"), c("daily, 100mg","BID, 50mg"), c("daily, 75mg", "BID, 50mg in morning, 25mg at night"),
                                        c("daily, 100mg", "daily, 75mg"), c("BID, 50mg", "BID, 50mg in morning, 25mg at night")), method = "wilcox.test", size=2)

pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-boxplot-insilico-trial_ct_combined_FulOnly_3rates_logistic_pk.pdf")
ggarrange(insilico_boxplot_ND_FulOnly + ggtitle("-DOX"), insilico_boxplot_Dox_FulOnly + ggtitle("+DOX"),legend = "right", common.legend = T)
dev.off()

waterfall_ND_FulOnly = insilico_ND_FulOnly %>% group_by(schedule) %>% filter(time == max(time), id %in% sample(1:1500, 30)) %>% unite(patient,id, schedule, sep = ": " , remove = F) %>% 
  ggplot(aes(x = reorder(patient, -total), y = total, fill = schedule)) + geom_bar(stat="identity", width = 0.8) + scale_y_log10() + theme(axis.text.x = element_blank()) 

waterfall_Dox_FulOnly = insilico_Dox_FulOnly %>% group_by(schedule) %>% filter(time == max(time), id %in% sample(1:1500, 30)) %>% unite(patient,id, schedule, sep = ": " , remove = F) %>% 
  ggplot(aes(x = reorder(patient, -total), y = total, fill = schedule)) + geom_bar(stat="identity", wdith = 0.8) + scale_y_log10() + theme(axis.text.x = element_blank()) 

pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-waterfall-insilico-trial_ct_combined_FulOnly_3rates_logistic_pk.pdf")
ggarrange(waterfall_ND_FulOnly + ggtitle("-DOX"), waterfall_Dox_FulOnly + ggtitle("+DOX"),legend = "right", common.legend = T, ncol = 1)
dev.off()

insilico_ND_2_FulOnly = full_join(insilico_ND_FulOnly %>% filter(time %% 0.5 == 0), 
                                  simulation_results_24$palbo_3w_1w_125mg$simulation_ND %>% filter(time %% 0.5 == 0) %>% dplyr::select(time, id, total),
                          by = c("time", "id"))

insilico_Dox_2_FulOnly = full_join(insilico_Dox_FulOnly %>% filter(time %% 0.5 == 0), 
                           simulation_results_6$palbo_3w_1w_125mg$simulation_Dox %>% filter(time %% 0.5 == 0) %>% dplyr::select(time, id, total),
                           by = c("time", "id"))

insilico_ND_2_FulOnly = insilico_ND_2_FulOnly %>% mutate(percent_of_standard = total.x/total.y)
insilico_Dox_2_FulOnly = insilico_Dox_2_FulOnly %>% mutate(percent_of_standard = total.x/total.y)

insilico_ND_2_sum_FulOnly = insilico_ND_2_FulOnly %>% group_by(time, schedule) %>% 
  summarize(percent_low = quantile(percent_of_standard,0.05, na.rm = T),
            percent_median = median(percent_of_standard, na.rm = T),
            percent_high = quantile(percent_of_standard,0.95, na.rm = T))
insilico_Dox_2_sum_FulOnly = insilico_Dox_2_FulOnly %>% group_by(time, schedule) %>% 
  summarize(percent_low = quantile(percent_of_standard,0.05, na.rm = T),
            percent_median = median(percent_of_standard, na.rm = T),
            percent_high = quantile(percent_of_standard,0.95, na.rm = T))


p1 <- ggplot(insilico_ND_2_sum_FulOnly, aes(time, percent_median, color = schedule)) + 
  geom_ribbon(aes(ymin = percent_low, ymax = percent_high, fill = schedule, color = NULL), alpha = 0.3) +
  geom_line() + theme_bw() + scale_y_log10() + xlab("Time (days)") + ylab("Percent of current standard (3 weeks on, 1 week off, 125mg)") +
  ggtitle("In silico trial schedule comparison (DOX-)") + theme(plot.title = element_text(size = 10)) + theme(legend.text = element_text(size = 6))


p2 <- ggplot(insilico_Dox_2_sum_FulOnly, aes(time, percent_median, color = schedule)) + 
  geom_ribbon(aes(ymin = percent_low, ymax = percent_high, fill = schedule, color = NULL), alpha = 0.3) +
  geom_line() + theme_bw() + scale_y_log10() + xlab("Time (days)") + ylab("Percent of current standard (3 weeks on, 1 week off, 125mg)") +
  ggtitle("In silico trial schedule comparison (DOX+)") + theme(plot.title = element_text(size = 10)) + theme(legend.text = element_text(size = 6))


pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-insilico-trial_ct_combined_FulOnly_2_3rates_logistic.pdf")
ggarrange(p1,p2, ncol = 2, nrow = 1,common.legend = T, legend = "bottom")
dev.off()


insilico_ND_3_FulOnly = full_join(insilico_ND_FulOnly %>% filter(time %% 0.5 == 0), 
                                  simulation_results_24$fulv_500mg_only$simulation_ND %>% filter(time %% 0.5 == 0) %>% dplyr::select(time, id, total),
                                  by = c("time", "id"))

insilico_Dox_3_FulOnly = full_join(insilico_Dox_FulOnly %>% filter(time %% 0.5 == 0), 
                                   simulation_results_6$fulv_500mg_only$simulation_Dox %>% filter(time %% 0.5 == 0) %>% dplyr::select(time, id, total),
                                   by = c("time", "id"))

insilico_ND_3_FulOnly = insilico_ND_3_FulOnly %>% mutate(percent_of_standard = total.x/total.y)
insilico_Dox_3_FulOnly = insilico_Dox_3_FulOnly %>% mutate(percent_of_standard = total.x/total.y)

insilico_ND_3_sum_FulOnly = insilico_ND_3_FulOnly %>% group_by(time, schedule) %>% 
  summarize(percent_low = quantile(percent_of_standard,0.05, na.rm = T),
            percent_median = median(percent_of_standard, na.rm = T),
            percent_high = quantile(percent_of_standard,0.95, na.rm = T))
insilico_Dox_3_sum_FulOnly = insilico_Dox_3_FulOnly %>% group_by(time, schedule) %>% 
  summarize(percent_low = quantile(percent_of_standard,0.05, na.rm = T),
            percent_median = median(percent_of_standard, na.rm = T),
            percent_high = quantile(percent_of_standard,0.95, na.rm = T))


p1_3 <- ggplot(insilico_ND_3_sum_FulOnly, aes(time, percent_median, color = schedule)) + 
  geom_ribbon(aes(ymin = percent_low, ymax = percent_high, fill = schedule, color = NULL), alpha = 0.3) +
  geom_line() + theme_bw() + scale_y_log10() + xlab("Time (days)") + ylab("Percent of Fulvestrant Only") +
  ggtitle("In silico trial schedule comparison (DOX-)") + theme(plot.title = element_text(size = 10)) + theme(legend.text = element_text(size = 6))

p2_3 <- ggplot(insilico_Dox_3_sum_FulOnly, aes(time, percent_median, color = schedule)) + 
  geom_ribbon(aes(ymin = percent_low, ymax = percent_high, fill = schedule, color = NULL), alpha = 0.3) +
  geom_line() + theme_bw() + scale_y_log10() + xlab("Time (days)") + ylab("Percent of Fulvestrant Only") +
  ggtitle("In silico trial schedule comparison (DOX+)") + theme(plot.title = element_text(size = 10)) + theme(legend.text = element_text(size = 6))

pdf("posterior_predictions_effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF_posterior-predictions-insilico-trial_ct_combined_FulOnly_3rates_logistic.pdf")
ggarrange(p1_3,p2_3, ncol = 2, nrow = 1,common.legend = T, legend = "bottom")
dev.off()


## log_log Dox
load("/Users/shaynastein/Google Drive/Harvard/MichorLab/BreastCancerCombination/code/BreastCancerModeling/from_cluster_include_longterm_data/model_fits/cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_logF_Dox_1.RData")
model_qc(data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1, name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-Dox-1",
         params_to_plot = c("sigma","alpha","b_P","b_F","c_P","c_F","a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf1 = c("a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))
posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-logF-Dox-1", dox = "Dox", skip_pass = TRUE)
looic_Dox_logP_logF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_logF_1)
print(looic_Dox_logP_logF)

## log_weibull_ND
load("/Users/shaynastein/Google Drive/Harvard/MichorLab/BreastCancerCombination/code/BreastCancerModeling/from_cluster_include_longterm_data/model_fits/cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_weibullF_ND_1.RData")
model_qc(data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_weibullF_1, name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-weibullF-ND-1",
         params_to_plot = c("sigma","alpha","b_P","b_F","c_P","c_F","a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf1 = c("a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))
posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_weibullF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-weibullF-ND-1", dox = "ND", skip_pass = TRUE)

looic_ND_logP_weibullF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_weibullF_1)
print(looic_ND_logP_weibullF)

## log_weibull_Dox
load("/Users/shaynastein/Google Drive/Harvard/MichorLab/BreastCancerCombination/code/BreastCancerModeling/from_cluster_include_longterm_data/model_fits/cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_logP_weibullF_Dox_1.RData")
model_qc(data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_logP_weibullF_1, name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-weibullF-Dox-1",
         params_to_plot = c("sigma","alpha","b_P","b_F","c_P","c_F","a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf1 = c("a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))
posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_weibullF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-logP-weibullF-Dox-1", dox = "Dox", skip_pass = TRUE)

looic_Dox_logP_weibullF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_logP_weibullF_1)
print(looic_Dox_logP_weibullF)

## weibull_log_ND
load("/Users/shaynastein/Google Drive/Harvard/MichorLab/BreastCancerCombination/code/BreastCancerModeling/from_cluster_include_longterm_data/model_fits/cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_weibullP_logF_ND_1.RData")
model_qc(data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_weibullP_logF_1, name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-weibullP-logF-ND-1",
         params_to_plot = c("sigma","alpha","b_P","b_F","c_P","c_F","a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf1 = c("a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))
posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_weibullP_logF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-weibullP-logF-ND-1", dox = "ND", skip_pass = TRUE)

looic_ND_weibullP_logF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_weibullP_logF_1)
print(looic_ND_weibullP_logF)

## weibull_log_Dox
load("/Users/shaynastein/Google Drive/Harvard/MichorLab/BreastCancerCombination/code/BreastCancerModeling/from_cluster_include_longterm_data/model_fits/cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_weibullP_logF_Dox_1.RData")
model_qc(data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_weibullP_logF_1, name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-weibullP-logF-Dox-1",
         params_to_plot = c("sigma","alpha","b_P","b_F","c_P","c_F","a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf1 = c("a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))
posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_weibullP_logF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-weibullP-logF-Dox-1", dox = "Dox", skip_pass = TRUE)

looic_Dox_weibullP_logF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_weibullP_logF_1)
print(looic_Dox_weibullP_logF)

## weibull_weibull_ND
load("/Users/shaynastein/Google Drive/Harvard/MichorLab/BreastCancerCombination/code/BreastCancerModeling/from_cluster_include_longterm_data/model_fits/cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_weibullP_weibullF_ND_1.RData")
model_qc(data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_weibullP_weibullF_1, name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-weibullP-weibullF-ND-1",
         params_to_plot = c("sigma","alpha","b_P","b_F","c_P","c_F","a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf1 = c("a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))
posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_weibullP_weibullF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-weibullP-lweibullF-ND-1", dox = "ND", skip_pass = TRUE)

looic_ND_weibullP_weibullF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_ND_weibullP_weibullF_1)
print(looic_ND_weibullP_weibullF)

## weibull_weibull_Dox
load("/Users/shaynastein/Google Drive/Harvard/MichorLab/BreastCancerCombination/code/BreastCancerModeling/from_cluster_include_longterm_data/model_fits/cell_cycle_ode_effectiveDose_logisticGrowth_noLongtermData_palboOnly_weibullP_weibullF_Dox_1.RData")
model_qc(data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_weibullP_weibullF_1, name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-weibullP-weibullF-Dox-1",
         params_to_plot = c("sigma","alpha","b_P","b_F","c_P","c_F","a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf1 = c("a_FP","b_min","b_max","Gamma","delta","K"),
         pars_acf2 = c("alpha","b_P","b_F","c_P","c_F"))
posterior_prediction_plots_no_abema(model_data = effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_weibullP_weibullF_1, stan_data = stan_input_include_longterm_data,
                                    name = "effectiveDose-logisticGrowth-noLongtermData-palboOnly-weibullP-lweibullF-Dox-1", dox = "Dox", skip_pass = TRUE)

looic_Dox_weibullP_weibullF = calculate_looic(effectiveDose_logisticGrowth_noLongtermData_palboOnly_Dox_weibullP_weibullF_1)
print(looic_Dox_weibullP_weibullF)
