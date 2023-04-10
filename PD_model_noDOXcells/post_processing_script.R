library(ggplot2)
library(ggsci)
library(rstan)
library(bayesplot)
library(reshape2)
library(dplyr)
library(tidyr)
library(loo)
library(deSolve)

theme_set(theme_classic())

compare_parameters = function(data1, data2, name1, name2, params_in_both = c("sigma","alpha","b_A","b_P","b_F","c_A","c_P","c_F","b","a_FP","a_FA","b","Gamma","delta")){
  params1 = as.data.frame(data1, pars = c("alpha","b_P","b_F","c_P","c_F","b","a_FP","b","gamma","delta")) %>%
    gather(factor_key = TRUE) %>%
    mutate(key= str_replace(key, "gamma","Gamma"))
  
  params2 = as.data.frame(data2, pars = c("alpha","b_P","b_F","c_P","c_F","b","a_FP","b","Gamma","delta"))  %>%
    gather(factor_key = TRUE) %>%
    mutate(key= str_replace(key, "gamma","Gamma"))
  
  compare_params = bind_rows(list(model1 = params1, model2 = params2), .id = "Model")
  
  pdf(paste("from_cluster_include_longterm_data/plots/compare_", name2,"_vs_",name1,".pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(compare_params, aes(x = key, y = value, color = Model)) + geom_boxplot() + theme_bw())
  dev.off()  
}

model_qc = function(data, name, params_to_plot = c("sigma", "alpha", "b_P", "b_F", "c_P","c_F","a_FP","alpha_max"),
                    pars_acf1 = c("sigma[1]","sigma[2]","sigma[3]","sigma[4]","alpha","a_FP"), 
                    pars_acf2 = c("b_P","b_F","c_P","c_F")){
  print("Plot Model QC")
  pdf(paste("model-qc_",name,"_plots.pdf", sep = "", collapse = ""), useDingbats = F)
  print(rstan::traceplot(data, pars = params_to_plot))
  print(mcmc_rhat(rhat(data, pars = params_to_plot)) + yaxis_text(hjust = 1))
  print(rstan::traceplot(data, pars = c("lp__"))  + theme(axis.text.x = element_text(size=6, angle=45)))
  
  ratios <- neff_ratio(data,pars = params_to_plot)
  color_scheme_set("brightblue")
  print(mcmc_neff(ratios) + yaxis_ticks(on = TRUE) + yaxis_text(on = TRUE))
  print(mcmc_acf(data, pars_acf1))
  print(mcmc_acf(data, pars_acf2))
  print(mcmc_parcoord(as.array(data, pars = pars_acf1), np = nuts_params(data)))
  print(mcmc_parcoord(as.array(data, pars = pars_acf2), np = nuts_params(data)))
  print(mcmc_parcoord(as.array(data, pars = pars_acf1), np = nuts_params(data)) + theme_bw())
  print(mcmc_parcoord(as.array(data, pars = pars_acf2), np = nuts_params(data)) + theme_bw())
  dev.off()
}

posterior_prediction_plots_differentRates_palbociclib = function(model_data, stan_data, name, dox = "Dox"){
  predictions_cc =  as.data.frame(model_data, pars = c('y_pred_cc_p') )%>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lowerBound = quantile(value, probs = 0.025),
              median = quantile(value, probs = 0.5),
              upperBound = quantile(value, probs = 0.975)) %>% 
    mutate(data = c(stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$`G0/G1`,
                    stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$S,
                    stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$`G2/M`),
           Stage = rep(c("G0/G1","S","G2/M"), each = stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$T_obs),
           Day = rep(stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$t_obs,3),
           drug_dose = rep(stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$experiment_obs, 3),
           Experiment = rep("Palbociclib", 3*stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$T_obs)
    )  %>% separate(drug_dose, c("Fulvestrant","Palbociclib"), "_") %>%
    mutate( Fulvestrant = as.numeric(Fulvestrant), Palbociclib = as.numeric(Palbociclib)) %>% 
    mutate(Palbociclib = case_when(Palbociclib == 0 ~ 0,
                                   Palbociclib == 1 ~ 1.25e-08,
                                   Palbociclib == 2 ~ 2.5e-08,
                                   Palbociclib == 3 ~ 5e-08,
                                   Palbociclib == 4 ~ 1e-07,
                                   Palbociclib == 5 ~ 2e-07,
                                   Palbociclib == 6 ~ 4e-07),
           Fulvestrant = case_when(Fulvestrant == 0 ~ 0,
                                   Fulvestrant == 1 ~ 6.5e-10,
                                   Fulvestrant == 2 ~ 1.3e-09,
                                   Fulvestrant == 3 ~ 2.6e-09,
                                   Fulvestrant == 4 ~ 5.2e-09,
                                   Fulvestrant == 5 ~ 1.04e-08,
                                   Fulvestrant == 6 ~ 2.08e-08))
  
  predictions_ct <- as.data.frame(model_data, pars = c('y_pred_ct_p')) %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lowerBound = quantile(value, probs = 0.025),
              median = quantile(value, probs = 0.5),
              upperBound = quantile(value, probs = 0.975)) %>%
    mutate(data = stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$y_obs,
           Day = stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$t_obs,
           drug_dose = stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$experiment_obs,
           Experiment = rep("Palbociclib", stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$T_obs)
    ) %>% separate(drug_dose, c("Fulvestrant","Palbociclib"), "_") %>%
    mutate(Palbociclib = as.numeric(Palbociclib),
           Fulvestrant = as.numeric(Fulvestrant)) %>%
    mutate(Palbociclib = case_when(Palbociclib == 0 ~ 0,
                                   Palbociclib == 1 ~ 1.25e-08,
                                   Palbociclib == 2 ~ 2.5e-08,
                                   Palbociclib == 3 ~ 5e-08,
                                   Palbociclib == 4 ~ 1e-07,
                                   Palbociclib == 5 ~ 2e-07,
                                   Palbociclib == 6 ~ 4e-07),
           Fulvestrant = case_when(Fulvestrant == 0 ~ 0,
                                   Fulvestrant == 1 ~ 6.5e-10,
                                   Fulvestrant == 2 ~ 1.3e-09,
                                   Fulvestrant == 3 ~ 2.6e-09,
                                   Fulvestrant == 4 ~ 5.2e-09,
                                   Fulvestrant == 5 ~ 1.04e-08,
                                   Fulvestrant == 6 ~ 2.08e-08))
  
  pdf(paste("posterior_predictions_",name,"_plots.pdf", sep = "", collapse = ""), useDingbats = F)
  ## Cell total (Palbociclib)
  print(ggplot(predictions_ct %>% filter( Experiment == "Palbociclib"),
               aes(x = Day, y = data)) + geom_point() +
          geom_line(aes(x= Day, y = median)) +
          geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, color = NULL), alpha = 0.25) +
          ylab("Cell Count") +
          theme(legend.position = "bottom") +
          facet_grid(Palbociclib ~ Fulvestrant, labeller = label_both, scales = "free_y") +
          theme(axis.text.x = element_text(size=6), strip.text.x = element_text(size = 6)) + scale_y_log10()+
          ggtitle(" Palbociclib experiment") + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
  ## Cell cycle (Fulvestrant alone - Palbo exp)
  print(ggplot(subset(predictions_cc, Palbociclib == 0 & Experiment == "Palbociclib"),
               aes(x = Day, y = data, color = Stage)) + geom_point() +
          geom_line(aes(x= Day, y = median)) +
          geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
          ylab("Cell Count") +
          theme(legend.position = "bottom") +
          facet_grid(Stage~ Fulvestrant, labeller = label_both, scales = "free_y") + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10()+
          ggtitle("Fulvestrant alone -- Palbociclib experiment"))
  print(ggplot(subset(predictions_cc, Fulvestrant == 0  & Experiment == "Palbociclib"),
               aes(x = Day, y = data, color = Stage)) + geom_point() +
          geom_line(aes(x= Day, y = median)) +
          geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
          ylab("Cell Count") +
          theme(legend.position = "bottom") +
          facet_grid(Stage~ Palbociclib, labeller = label_both, scales = "free_y") + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10()+
          ggtitle("Palbociclib alone"))
  print(ggplot(subset(predictions_cc, Fulvestrant >0 & Palbociclib > 0 & Experiment == "Palbociclib"),
               aes(x = Day, y = data, color = Stage)) + geom_point() +
          geom_line(aes(x= Day, y = median)) +
          geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
          ylab("Cell Count") +
          theme(legend.position = "bottom") +
          facet_grid(Stage ~ Palbociclib + Fulvestrant, labeller = label_both, scales = "free_y") + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10() +
          ggtitle("Fulvestrant + Palbociclib"))
  dev.off()
}

posterior_prediction_plots = function(model_data, stan_data, name, dox = "Dox", skip_pass = FALSE){
  predictions_cc =  as.data.frame(model_data, pars = c('y_pred_cc_p', 'y_pred_cc_a') )%>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lowerBound = quantile(value, probs = 0.025),
              median = quantile(value, probs = 0.5),
              upperBound = quantile(value, probs = 0.975)) %>%
    mutate(data = c(stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$`G0/G1`,
                    stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$S,
                    stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$`G2/M`,
                    stan_data[[paste("abema_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$`G0/G1`,
                    stan_data[[paste("abema_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$S,
                    stan_data[[paste("abema_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$`G2/M`),
           Stage = c(rep(c("G0/G1","S","G2/M"), each = stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$T_obs),
                     rep(c("G0/G1","S","G2/M"), each = stan_data[[paste("abema_cc_",dox,"_input", sep = "", collapse = "")]]$T_obs)),
           Day = c(rep(stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$t_obs,3),
                   rep(stan_data[[paste("abema_cc_",dox,"_input", sep = "", collapse = "")]]$t_obs,3)),
           drug_dose = c(rep(stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$experiment_obs, 3),
                         rep(stan_data[[paste("abema_cc_",dox,"_input", sep = "", collapse = "")]]$experiment_obs, 3)),
           Experiment = c(rep("Palbociclib", 3*stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$T_obs),
                          rep("Abemaciclib", 3*stan_data[[paste("abema_cc_",dox,"_input", sep = "", collapse = "")]]$T_obs))
    )  %>% separate(drug_dose, c("Fulvestrant","Drug2"), "_") %>%
    mutate(Palbociclib = as.numeric(case_when(Experiment == "Palbociclib" ~ Drug2,
                                              TRUE ~ "0")),
           Abemaciclib = as.numeric(case_when(Experiment == "Abemaciclib" ~ Drug2,
                                              TRUE ~ "0")),
           Fulvestrant = as.numeric(Fulvestrant)) %>% select(-Drug2) %>%
    mutate(Palbociclib = case_when(Palbociclib == 0 ~ 0,
                                   Palbociclib == 1 ~ 1.25e-08,
                                   Palbociclib == 2 ~ 2.5e-08,
                                   Palbociclib == 3 ~ 5e-08,
                                   Palbociclib == 4 ~ 1e-07,
                                   Palbociclib == 5 ~ 2e-07,
                                   Palbociclib == 6 ~ 4e-07),
           Abemaciclib = case_when(Abemaciclib == 0 ~ 0,
                                   Abemaciclib == 1 ~ 8.02e-10,
                                   Abemaciclib == 2 ~ 1.6e-09,
                                   Abemaciclib == 3 ~ 3.21e-09,
                                   Abemaciclib == 4 ~ 6.41e-09,
                                   Abemaciclib == 5 ~ 1.28e-08,
                                   Abemaciclib == 6 ~ 2.57e-08),
           Fulvestrant = case_when(Fulvestrant == 0 ~ 0,
                                   Fulvestrant == 1 ~ 6.5e-10,
                                   Fulvestrant == 2 ~ 1.3e-09,
                                   Fulvestrant == 3 ~ 2.6e-09,
                                   Fulvestrant == 4 ~ 5.2e-09,
                                   Fulvestrant == 5 ~ 1.04e-08,
                                   Fulvestrant == 6 ~ 2.08e-08))
  
  predictions_ct <- as.data.frame(model_data, pars = c('y_pred_ct_p','y_pred_ct_a')) %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lowerBound = quantile(value, probs = 0.025),
              median = quantile(value, probs = 0.5),
              upperBound = quantile(value, probs = 0.975)) %>%
    mutate(data = c(stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$y_obs,
                    stan_data[[paste("abema_ct_",dox,"_input", sep = "", collapse = "")]]$y_obs),
           Day = c(stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$t_obs,
                   stan_data[[paste("abema_ct_",dox,"_input", sep = "", collapse = "")]]$t_obs),
           drug_dose = c(stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$experiment_obs,
                         stan_data[[paste("abema_ct_",dox,"_input", sep = "", collapse = "")]]$experiment_obs),
           Experiment = c(rep("Palbociclib", stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$T_obs),
                          rep("Abemaciclib", stan_data[[paste("abema_ct_",dox,"_input", sep = "", collapse = "")]]$T_obs))
    ) %>% separate(drug_dose, c("Fulvestrant","Drug2"), "_") %>%
    mutate(Palbociclib = as.numeric(case_when(Experiment == "Palbociclib" ~ Drug2,
                                              TRUE ~ "0")),
           Abemaciclib = as.numeric(case_when(Experiment == "Abemaciclib" ~ Drug2,
                                              TRUE ~ "0")),
           Fulvestrant = as.numeric(Fulvestrant)) %>% select(-Drug2) %>%
    mutate(Palbociclib = case_when(Palbociclib == 0 ~ 0,
                                   Palbociclib == 1 ~ 1.25e-08,
                                   Palbociclib == 2 ~ 2.5e-08,
                                   Palbociclib == 3 ~ 5e-08,
                                   Palbociclib == 4 ~ 1e-07,
                                   Palbociclib == 5 ~ 2e-07,
                                   Palbociclib == 6 ~ 4e-07),
           Abemaciclib = case_when(Abemaciclib == 0 ~ 0,
                                   Abemaciclib == 1 ~ 8.02e-10,
                                   Abemaciclib == 2 ~ 1.6e-09,
                                   Abemaciclib == 3 ~ 3.21e-09,
                                   Abemaciclib == 4 ~ 6.41e-09,
                                   Abemaciclib == 5 ~ 1.28e-08,
                                   Abemaciclib == 6 ~ 2.57e-08),
           Fulvestrant = case_when(Fulvestrant == 0 ~ 0,
                                   Fulvestrant == 1 ~ 6.5e-10,
                                   Fulvestrant == 2 ~ 1.3e-09,
                                   Fulvestrant == 3 ~ 2.6e-09,
                                   Fulvestrant == 4 ~ 5.2e-09,
                                   Fulvestrant == 5 ~ 1.04e-08,
                                   Fulvestrant == 6 ~ 2.08e-08)
    )
  if(!skip_pass){
    predictions_pass <- as.data.frame(model_data, pars = c('y_pred_pass_p')) %>%
      gather(factor_key = TRUE) %>%
      group_by(key) %>%
      summarize(lowerBound = quantile(value, probs = 0.025),
                median = quantile(value, probs = 0.5),
                upperBound = quantile(value, probs = 0.975)) %>%
      mutate(data = stan_data[[paste("validation_data_1_",dox,"_input", sep = "", collapse = "")]]$y_obs,
             Day = stan_data[[paste("validation_data_1_",dox,"_input", sep = "", collapse = "")]]$t_obs,
             Experiment = stan_data[[paste("validation_data_1_",dox,"_input", sep = "", collapse = "")]]$experiment_obs)
  }
  pdf(paste("posterior_predictions_",name,"_plots.pdf", sep = "", collapse = ""), useDingbats = F)
  ## Cell total (Palbociclib)
  print(ggplot(predictions_ct %>% filter( Experiment == "Palbociclib"),
         aes(x = Day, y = data)) + geom_point() +
    geom_line(aes(x= Day, y = median)) +
    geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, color = NULL), alpha = 0.25) +
    ylab("Cell Count") +
    theme(legend.position = "bottom") +
    facet_grid(Palbociclib ~ Fulvestrant, labeller = label_both, scales = "free_y") +
    theme(axis.text.x = element_text(size=6), strip.text.x = element_text(size = 6)) + scale_y_log10()+
    ggtitle(" Palbociclib experiment") + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
  ## Cell total (Abemaciclib)
  print(ggplot(predictions_ct %>% filter( Experiment == "Abemaciclib"),
         aes(x = Day, y = data)) + geom_point() +
    geom_line(aes(x= Day, y = median)) +
    geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, color = NULL), alpha = 0.25) +
    ylab("Cell Count") +
    theme(legend.position = "bottom") +
    facet_grid(Abemaciclib ~ Fulvestrant, labeller = label_both, scales = "free_y") +
    theme(axis.text.x = element_text(size=6), strip.text.x = element_text(size = 6)) + scale_y_log10()+
    ggtitle(" Abemaciclib experiment") + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
  ## Cell cycle (Fulvestrant alone - Palbo exp)
  print(ggplot(subset(predictions_cc, Palbociclib == 0 & Abemaciclib == 0 & Experiment == "Palbociclib"),
         aes(x = Day, y = data, color = Stage)) + geom_point() +
    geom_line(aes(x= Day, y = median)) +
    geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
    ylab("Cell Count") +
    theme(legend.position = "bottom") +
    facet_grid(Stage~ Fulvestrant, labeller = label_both, scales = "free_y") + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10()+
    ggtitle("Fulvestrant alone -- Palbociclib experiment"))
  ## cell cycle (Fulvestrant alone - abema exp)
  print(ggplot(subset(predictions_cc, Palbociclib == 0 & Abemaciclib == 0 & Experiment == "Abemaciclib"),
         aes(x = Day, y = data, color = Stage)) + geom_point() +
    geom_line(aes(x= Day, y = median)) +
    geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
    ylab("Cell Count") +
    theme(legend.position = "bottom") +
    facet_grid(Stage~ Fulvestrant, labeller = label_both, scales = "free_y") + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10() +
    ggtitle("Fulvestrant alone -- Abemaciclib experiment"))
  print(ggplot(subset(predictions_cc, Fulvestrant == 0 & Abemaciclib == 0 & Experiment == "Palbociclib"),
         aes(x = Day, y = data, color = Stage)) + geom_point() +
    geom_line(aes(x= Day, y = median)) +
    geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
    ylab("Cell Count") +
    theme(legend.position = "bottom") +
    facet_grid(Stage~ Palbociclib, labeller = label_both, scales = "free_y") + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10()+
    ggtitle("Palbociclib alone"))
  print(ggplot(subset(predictions_cc, Fulvestrant == 0 & Palbociclib == 0 & Experiment == "Abemaciclib"),
         aes(x = Day, y = data, color = Stage)) + geom_point() +
    geom_line(aes(x= Day, y = median)) +
    geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
    ylab("Cell Count") +
    theme(legend.position = "bottom") +
    facet_grid(Stage~ Abemaciclib, labeller = label_both, scales = "free_y") + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10()+
    ggtitle("Abemaciclib alone"))
  print(ggplot(subset(predictions_cc, Fulvestrant >0 & Palbociclib > 0 & Experiment == "Palbociclib"),
         aes(x = Day, y = data, color = Stage)) + geom_point() +
    geom_line(aes(x= Day, y = median)) +
    geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
    ylab("Cell Count") +
    theme(legend.position = "bottom") +
    facet_grid(Stage ~ Palbociclib + Fulvestrant, labeller = label_both, scales = "free_y") + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10() +
    ggtitle("Fulvestrant + Palbociclib"))
  print(ggplot(subset(predictions_cc, Fulvestrant >0 & Abemaciclib > 0 & Experiment == "Abemaciclib"),
         aes(x = Day, y = data, color = Stage)) + geom_point() +
    geom_line(aes(x= Day, y = median)) +
    geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
    ylab("Cell Count") +
    theme(legend.position = "bottom") +
    facet_grid(Stage ~ Abemaciclib + Fulvestrant, labeller = label_both, scales = "free_y") + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    scale_y_log10() +
    ggtitle("Fulvestrant + Abemaciclib"))
  if(!skip_pass){
    print(ggplot(predictions_pass, aes(x = Day, y = data)) + geom_point() +
      geom_line(aes(x = Day, y = median)) + 
      geom_ribbon(aes(ymin=lowerBound, ymax = upperBound, color = NULL, alpha = 0.25)) + 
      ylab("Cell Count") +
      theme(legend.position = "bottom") +
      facet_wrap(~Experiment, labeller = label_both, scales = "free_y") + theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10() +
      ggtitle("Longterm data") + theme(position = NULL))
    print(ggplot(predictions_pass %>% mutate(Day = round(Day)), aes(x = Day, y = data)) + geom_point(size=2) +
      geom_errorbar(data = predictions_pass %>% mutate(Day = round(Day)) %>% group_by(Day, Experiment) %>% mutate(lowerBound = min(lowerBound), upperBound = max(upperBound), median = median(median)), 
                    aes(x = Day, ymin = lowerBound, ymax = upperBound)) +
      facet_wrap(~Experiment, labeller = label_both, scales = "free_y") + theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      ggtitle("Longterm data") )
  }
  dev.off()
}


predictions_cc =  as.data.frame(test_singleparameter_6ODES_even_1000iters_3rates, pars = c('y_pred_cc_p') )%>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lowerBound = quantile(value, probs = 0.025),
            median = quantile(value, probs = 0.5),
            upperBound = quantile(value, probs = 0.975))
p_cc_mean <- mean(predictions_cc$median)
p_cc_lower <- mean(predictions_cc$lowerBound)
p_cc_upper <- mean(predictions_cc$upperBound)
p_cc_upper - p_cc_lower


predictions_ct =  as.data.frame(test_singleparameter_6ODES_even_1000iters_3rates, pars = c('y_pred_ct_p') )%>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lowerBound = quantile(value, probs = 0.025),
            median = quantile(value, probs = 0.5),
            upperBound = quantile(value, probs = 0.975))
p_ct_mean <- mean(predictions_ct$median)
p_ct_lower <- mean(predictions_ct$lowerBound)
p_ct_upper <- mean(predictions_ct$upperBound)
p_ct_upper - p_ct_lower

posterior_prediction_plots_no_abema = function(model_data, stan_data, name, dox = "Dox", skip_pass = FALSE){
  predictions_cc =  as.data.frame(model_data, pars = c('y_pred_cc_p') )%>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lowerBound = quantile(value, probs = 0.025),
              median = quantile(value, probs = 0.5),
              upperBound = quantile(value, probs = 0.975)) %>%
    mutate(data = c(stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$`G0/G1`,
                    stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$S,
                    stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$y_obs$`G2/M`),
           Stage = c(rep(c("G0/G1","S","G2/M"), each = stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$T_obs)),
           Day = c(rep(stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$t_obs,3)),
           drug_dose = c(rep(stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$experiment_obs, 3)),
           Experiment = c(rep("Palbociclib", 3*stan_data[[paste("palbo_cc_",dox,"_input", sep = "", collapse = "")]]$T_obs))
    )  %>% separate(drug_dose, c("Fulvestrant","Palbociclib"), "_") %>%
    mutate(Fulvestrant = as.numeric(Fulvestrant),
           Palbociclib = as.numeric(Palbociclib)) %>%
    mutate(Palbociclib = case_when(Palbociclib == 0 ~ 0,
                                   Palbociclib == 1 ~ 12.50,
                                   Palbociclib == 2 ~ 25.00,
                                   Palbociclib == 3 ~ 50.00,
                                   Palbociclib == 4 ~ 100.00,
                                   Palbociclib == 5 ~ 200.00,
                                   Palbociclib == 6 ~ 400.00),
           Fulvestrant = case_when(Fulvestrant == 0 ~ 0,
                                   Fulvestrant == 1 ~ 1.25,
                                   Fulvestrant == 2 ~ 2.50,
                                   Fulvestrant == 3 ~ 5.00,
                                   Fulvestrant == 4 ~ 10.00,
                                   Fulvestrant == 5 ~ 20.00,
                                   Fulvestrant == 6 ~ 40.00))
    
  
  predictions_ct <- as.data.frame(model_data, pars = c('y_pred_ct_p')) %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lowerBound = quantile(value, probs = 0.025),
              median = quantile(value, probs = 0.5),
              upperBound = quantile(value, probs = 0.975)) %>%
    mutate(data = c(stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$y_obs),
           Day = c(stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$t_obs),
           drug_dose = c(stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$experiment_obs),
           Experiment = c(rep("Palbociclib", stan_data[[paste("palbo_ct_",dox,"_input", sep = "", collapse = "")]]$T_obs))
    ) %>% separate(drug_dose, c("Fulvestrant","Palbociclib"), "_") %>%
    mutate(Fulvestrant = as.numeric(Fulvestrant),
           Palbociclib = as.numeric(Palbociclib)) %>%
    mutate(Palbociclib = case_when(Palbociclib == 0 ~ 0,
                                   Palbociclib == 1 ~ 12.50,
                                   Palbociclib == 2 ~ 25.50,
                                   Palbociclib == 3 ~ 50.10,
                                   Palbociclib == 4 ~ 100.00,
                                   Palbociclib == 5 ~ 200.00,
                                   Palbociclib == 6 ~ 401.00),
           Fulvestrant = case_when(Fulvestrant == 0 ~ 0,
                                   Fulvestrant == 1 ~ 0.65,
                                   Fulvestrant == 2 ~ 1.30,
                                   Fulvestrant == 3 ~ 2.60,
                                   Fulvestrant == 4 ~ 5.20,
                                   Fulvestrant == 5 ~ 10.40,
                                   Fulvestrant == 6 ~ 20.80)
    )
  if(!skip_pass){
    predictions_pass <- as.data.frame(model_data, pars = c('y_pred_pass_p')) %>%
      gather(factor_key = TRUE) %>%
      group_by(key) %>%
      summarize(lowerBound = quantile(value, probs = 0.025),
                median = quantile(value, probs = 0.5),
                upperBound = quantile(value, probs = 0.975)) %>%
      mutate(data = stan_data[[paste("validation_data_1_",dox,"_input", sep = "", collapse = "")]]$y_obs,
             Day = stan_data[[paste("validation_data_1_",dox,"_input", sep = "", collapse = "")]]$t_obs,
             Experiment = stan_data[[paste("validation_data_1_",dox,"_input", sep = "", collapse = "")]]$experiment_obs)
  }
  pdf(paste("posterior_predictions_",name,"_plots.pdf", sep = "", collapse = ""), useDingbats = F)
  ## Cell total (Palbociclib)
  print(ggplot(predictions_ct %>% filter( Experiment == "Palbociclib"),
               aes(x = Day, y = data)) + geom_point() +
          geom_line(aes(x= Day, y = median)) +
          geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, color = NULL), alpha = 0.25) +
          ylab("Cell Count") +
          theme(legend.position = "bottom") +
          facet_grid(Palbociclib ~ Fulvestrant, labeller = label_both, scales = "free_y") +
          theme(axis.text.x = element_text(size=6), strip.text.x = element_text(size = 6)) + scale_y_log10()+
          ggtitle(" Palbociclib experiment") + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
  ## Cell cycle (Fulvestrant alone - Palbo exp)
  print(ggplot(subset(predictions_cc, Palbociclib == 0  & Experiment == "Palbociclib"),
               aes(x = Day, y = data, color = Stage)) + geom_point() +
          geom_line(aes(x= Day, y = median)) +
          geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
          ylab("Cell Count") +
          theme(legend.position = "bottom") +
          facet_grid(Stage~ Fulvestrant, labeller = label_both, scales = "free_y") + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10()+
          ggtitle("Fulvestrant alone -- Palbociclib experiment"))
  print(ggplot(subset(predictions_cc, Fulvestrant == 0  & Experiment == "Palbociclib"),
               aes(x = Day, y = data, color = Stage)) + geom_point() +
          geom_line(aes(x= Day, y = median)) +
          geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
          ylab("Cell Count") +
          theme(legend.position = "bottom") +
          facet_grid(Stage~ Palbociclib, labeller = label_both, scales = "free_y") + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10()+
          ggtitle("Palbociclib alone"))
  print(ggplot(subset(predictions_cc, Fulvestrant >0 & Palbociclib > 0 & Experiment == "Palbociclib"),
               aes(x = Day, y = data, color = Stage)) + geom_point() +
          geom_line(aes(x= Day, y = median)) +
          geom_ribbon(aes(ymin = lowerBound, ymax = upperBound, fill = Stage, color = NULL), alpha = 0.25) +
          ylab("Cell Count") +
          theme(legend.position = "bottom") +
          facet_grid(Stage ~ Palbociclib + Fulvestrant, labeller = label_both, scales = "free_y") + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10() +
          ggtitle("Fulvestrant + Palbociclib"))
  if(!skip_pass){
    print(ggplot(predictions_pass, aes(x = Day, y = data)) + geom_point() +
            geom_line(aes(x = Day, y = median)) + 
            geom_ribbon(aes(ymin=lowerBound, ymax = upperBound, color = NULL, alpha = 0.25)) + 
            ylab("Cell Count") +
            theme(legend.position = "bottom") +
            facet_wrap(~Experiment, labeller = label_both, scales = "free_y") + theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_log10() +
            ggtitle("Longterm data") + theme(position = NULL))
    print(ggplot(predictions_pass %>% mutate(Day = round(Day)), aes(x = Day, y = data)) + geom_point(size=2) +
            geom_errorbar(data = predictions_pass %>% mutate(Day = round(Day)) %>% group_by(Day, Experiment) %>% mutate(lowerBound = min(lowerBound), upperBound = max(upperBound), median = median(median)), 
                          aes(x = Day, ymin = lowerBound, ymax = upperBound)) +
            facet_wrap(~Experiment, labeller = label_both, scales = "free_y") + theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
            ggtitle("Longterm data") )
  }
  dev.off()
}

validation_schedules_96well = function(validation_fn, model_data, name,
                                       params= c("alpha","b_P","b_F","b_A","c_P","c_F","c_A","a_FP","a_FA","b","Gamma","delta","K"),
                                       num_iters = 2250, dox = "Dox"){
  valid2 = read.csv(validation_fn)
  valid2_tx = valid2 %>% filter(Measurement == "Treatment", Dox == dox) %>% arrange(Experiment, Day) %>% 
    mutate(Fulvestrant = Fulvestrant * 1e6, Palbociclib = Palbociclib * 1e6, Abemaciclib = Abemaciclib*1e6)
  valid2_obs = valid2 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)

  params_for_prediction = as.data.frame(model_data, pars = params) 
  
  init_conditions = as.data.frame(model_data, pars = c("y0_pass_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  ## with passaging
  prediction_results = data.frame()
  for (exp in unique(valid2_tx$Experiment)){
    pred = lapply(1:num_iters, 
                  function(i){
                    ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper(params = params_for_prediction[i,], 
                                                                                                    init = init_conditions[i,], 
                                                                                                    tmax = 28, 
                                                                                                    dose_times = (valid2_tx %>% filter(Experiment == exp))$Day,
                                                                                                    drug_amount = valid2_tx %>% filter(Experiment == exp) %>% 
                                                                                                      select (Fulvestrant,Palbociclib,Abemaciclib,Passage)) })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% 
      summarize(total_mean = mean(total), total_low = quantile(total, 0.025, na.rm = T), 
                total_high = quantile(total, 0.975, na.rm = T), total_median = median(total, na.rm = T),
                AA = mean(AA), PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
   
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name,"_posterior-predictions-passaging.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
    geom_point(data = valid2_obs, aes(x = Day, y = total_count))  +
    geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + facet_wrap(Experiment ~. , ncol = 5))
  dev.off()

  prediction_results_no_passaging = data.frame()
  for (exp in unique(valid2_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper_no_passaging(params = params_for_prediction[i,], 
                                                                                                   init = init_conditions[i,], 
                                                                                                   tmax = 28, 
                                                                                                   dose_times = (valid2_tx %>% filter(Experiment == exp))$Day,
                                                                                                   drug_amount = valid2_tx %>% filter(Experiment == exp) %>% 
                                                                                                     select (Fulvestrant,Palbociclib,Abemaciclib,Passage))})
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>%
      mutate(total = G0G1 + S + G2M) %>% 
      summarize(total_mean = mean(total, na.rm = T), total_low = quantile(total, 0.025, na.rm = T), 
                total_high = quantile(total, 0.975, na.rm = T), total_median = median(total, na.rm = T),
                AA = mean(AA), PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results_no_passaging = bind_rows(prediction_results_no_passaging, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_", name, "_posterior-predictions-noPassaging.pdf",sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results_no_passaging, aes(x = time, y = total_median)) + geom_line()  + geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + facet_wrap(Experiment ~. , ncol = 5) +scale_y_log10())
  print(ggplot(prediction_results_no_passaging, aes(x = time, y = total_median, color = Experiment)) + geom_line() + scale_y_log10())
  dev.off()
  
  write.csv(prediction_results_no_passaging %>% filter(time ==28) %>% arrange(-total_median),
            paste("from_cluster_include_longterm_data/plots/posterior_predictions_", name, "_posterior-predictions_noPassaging_ranks.csv", sep = "", collapse = ""),
            quote = FALSE, row.names = FALSE)
  
}


validation_schedules_24well = function(validation_fn, model_data, name,
                                       params = c("alpha","b_P","b_F","b_A","c_P","c_F","c_A","a_FP","a_FA","b","Gamma","delta","K"),
                                       num_iters = 2250, dox = "Dox"){
  valid1_24wells = read.csv(validation_fn)
  valid1_24wells_tx = valid1_24wells %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e6, Palbociclib = Palbociclib * 1e6, Abemaciclib = Abemaciclib*1e6)
  valid1_24wells_obs = valid1_24wells %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions_24wells = init_conditions %>% 
    mutate(total = G0G1 + S + G2M, G0G1 = G0G1/total, S = S/total, G2M = G2M/total) %>% 
    mutate(total = mean((valid1_24wells_obs %>% filter(Day == 0, Experiment == "Control"))$total_count)) %>% 
    mutate(G0G1 = G0G1 * total,
           S = S * total,
           G2M = G2M * total) %>%
    select(-total)
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  params_for_prediction_24wells = params_for_prediction %>% mutate(K = 6)
  
  prediction_results_24wells = data.frame()
  for (exp in unique(valid1_24wells_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
                    ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper(params = params_for_prediction_24wells[i,], 
                                                                                                    init = init_conditions_24wells[i,], 
                                                                                                    tmax = 28, 
                                                                                                    dose_times = (valid1_24wells_tx %>% filter(Experiment == exp))$Day,
                                                                                                    drug_amount = valid1_24wells_tx %>% filter(Experiment == exp) %>% 
                                                                                                      select (Fulvestrant,Palbociclib,Abemaciclib,Passage))
                  })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                   total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                   total_median = median(total, na.rm = T),
                                                   AA = mean(AA), PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results_24wells = bind_rows(prediction_results_24wells, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name, "_posterior-predictions-passaging-24wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results_24wells, aes(x = time, y = total_median)) + geom_line() + 
    geom_point(data = valid1_24wells_obs, aes(x = Day, y = total_count))  + 
      geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
      facet_wrap(Experiment ~. , ncol = 3)) 
  
  print(ggplot(prediction_results_24wells %>% filter(Experiment == "Control"), aes(x = time, y = total_median)) + 
          geom_line() + 
    geom_point(data = valid1_24wells_obs %>% filter(Experiment == "Control"), aes(x = Day, y = total_count))  + 
      geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
  
}

validation_schedules_24well_palbo_only = function(validation_fn, model_data, name,
                                       params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_min","Gamma","delta","K"),
                                       num_iters = 1500, dox = "Dox"){
  valid1_24wells = read.csv(validation_fn)
  valid1_24wells_tx = valid1_24wells %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e10, Palbociclib = Palbociclib * 1e8, Abemaciclib = Abemaciclib*1e6)
  valid1_24wells_obs = valid1_24wells %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions_24wells = init_conditions %>% 
    mutate(total = G0G1 + S + G2M, G0G1 = G0G1/total, S = S/total, G2M = G2M/total) %>% 
    mutate(total = mean((valid1_24wells_obs %>% filter(Day == 0, Experiment == "Control"))$total_count)) %>% 
    mutate(G0G1 = G0G1 * total,
           S = S * total,
           G2M = G2M * total) %>%
    select(-total)
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  params_for_prediction_24wells = params_for_prediction %>% mutate(K = 6)
  params_for_prediction_24wells = params_for_prediction_24wells %>% mutate(b = b_min, c_A = 1, b_A = 1, a_FA = 0)
  
  prediction_results_24wells = data.frame()
  for (exp in unique(valid1_24wells_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper(params = params_for_prediction_24wells[i,], 
                                                                                      init = init_conditions_24wells[i,], 
                                                                                      tmax = 28, 
                                                                                      dose_times = (valid1_24wells_tx %>% filter(Experiment == exp))$Day,
                                                                                      drug_amount = valid1_24wells_tx %>% filter(Experiment == exp) %>% 
                                                                                        select (Fulvestrant,Palbociclib,Abemaciclib,Passage))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                   total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                   total_median = median(total, na.rm = T),
                                                   PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results_24wells = bind_rows(prediction_results_24wells, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name, "_posterior-predictions-passaging-24wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results_24wells, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_24wells_obs, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 3)) 
  
  print(ggplot(prediction_results_24wells %>% filter(Experiment == "Control"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_24wells_obs %>% filter(Experiment == "Control"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
  
}


validation_schedules_10day_96well_palbo_only_mod = function(validation_fn, model_data, name,
                                                        params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max"),
                                                        num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  init_conditions = init_conditions/5
  colnames(init_conditions) = c("G0G1_1","S_1","G2M_1")
  init_conditions$G0G1_2 = init_conditions$G0G1_1
  init_conditions$G0G1_3 = init_conditions$G0G1_1
  init_conditions$G0G1_4 = init_conditions$G0G1_1
  init_conditions$G0G1_5 = init_conditions$G0G1_1
  init_conditions$S_2 = init_conditions$S_1
  init_conditions$S_3 = init_conditions$S_1
  init_conditions$S_4 = init_conditions$S_1
  init_conditions$S_5 = init_conditions$S_1
  init_conditions$G2M_2 = init_conditions$G2M_1
  init_conditions$G2M_3 = init_conditions$G2M_1
  init_conditions$G2M_4 = init_conditions$G2M_1
  init_conditions$G2M_5 = init_conditions$G2M_1
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_mod(params = params_for_prediction[i,], 
                                                                                            init = init_conditions[i,], 
                                                                                            tmax = 10, 
                                                                                            passage_day = passage_day,
                                                                                            passage_number = passage_number,
                                                                                            dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                            drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                              select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + S_1 + S_2 + S_3 + S_4 + S_5  + G2M_1 + G2M_2 + G2M_3 + G2M_4 + G2M_5) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_multiple = function(validation_fn, model_data, name,
                                                        params = c("alpha", "Beta", "Gamma", "b_P","b_F","c_P","c_F","a_FP","alpha_max"),
                                                        num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  init_conditions = init_conditions/5
  colnames(init_conditions) = c("G0G1_1","S_1","G2M_1")
  init_conditions$G0G1_2 = init_conditions$G0G1_1
  init_conditions$G0G1_3 = init_conditions$G0G1_1
  init_conditions$G0G1_4 = init_conditions$G0G1_1
  init_conditions$G0G1_5 = init_conditions$G0G1_1
  init_conditions$S_2 = init_conditions$S_1
  init_conditions$S_3 = init_conditions$S_1
  init_conditions$S_4 = init_conditions$S_1
  init_conditions$S_5 = init_conditions$S_1
  init_conditions$G2M_2 = init_conditions$G2M_1
  init_conditions$G2M_3 = init_conditions$G2M_1
  init_conditions$G2M_4 = init_conditions$G2M_1
  init_conditions$G2M_5 = init_conditions$G2M_1
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_multiple(params = params_for_prediction[i,], 
                                                                                            init = init_conditions[i,], 
                                                                                            tmax = 10, 
                                                                                            passage_day = passage_day,
                                                                                            passage_number = passage_number,
                                                                                            dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                            drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                              select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + S_1 + S_2 + S_3 + S_4 + S_5  + G2M_1 + G2M_2 + G2M_3 + G2M_4 + G2M_5) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_IC = function(validation_fn, model_data, name,
                                                        params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max"),
                                                        num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_pp")) 
  colnames(init_conditions) = c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_1", "S_2", "S_3", "S_4", "S_5", "G2M_1", "G2M_2", "G2M_3", "G2M_4", "G2M_5")
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_IC(params = params_for_prediction[i,], 
                                                                                            init = init_conditions[i,], 
                                                                                            tmax = 10, 
                                                                                            passage_day = passage_day,
                                                                                            passage_number = passage_number,
                                                                                            dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                            drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                              select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + S_1 + S_2 + S_3 + S_4 + S_5  + G2M_1 + G2M_2 + G2M_3 + G2M_4 + G2M_5) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}

validation_schedules_10day_96well_palbo_only_15_even = function(validation_fn, model_data, name,
                                                        params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max"),
                                                        num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  init_conditions = init_conditions/5
  colnames(init_conditions) = c("G0G1_1","S_1","G2M_1")
  init_conditions$G0G1_2 = init_conditions$G0G1_1
  init_conditions$G0G1_3 = init_conditions$G0G1_1
  init_conditions$G0G1_4 = init_conditions$G0G1_1
  init_conditions$G0G1_5 = init_conditions$G0G1_1
  init_conditions$S_2 = init_conditions$S_1
  init_conditions$S_3 = init_conditions$S_1
  init_conditions$S_4 = init_conditions$S_1
  init_conditions$S_5 = init_conditions$S_1
  init_conditions$G2M_2 = init_conditions$G2M_1
  init_conditions$G2M_3 = init_conditions$G2M_1
  init_conditions$G2M_4 = init_conditions$G2M_1
  init_conditions$G2M_5 = init_conditions$G2M_1
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_15_even(params = params_for_prediction[i,], 
                                                                                            init = init_conditions[i,], 
                                                                                            tmax = 10, 
                                                                                            passage_day = passage_day,
                                                                                            passage_number = passage_number,
                                                                                            dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                            drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                              select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + S_1 + S_2 + S_3 + S_4 + S_5  + G2M_1 + G2M_2 + G2M_3 + G2M_4 + G2M_5) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells_15_even.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_15_even_ND = function(validation_fn, model_data, name,
                                                                params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max"),
                                                                num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e14, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  init_conditions = init_conditions/5
  colnames(init_conditions) = c("G0G1_1","S_1","G2M_1")
  init_conditions$G0G1_2 = init_conditions$G0G1_1
  init_conditions$G0G1_3 = init_conditions$G0G1_1
  init_conditions$G0G1_4 = init_conditions$G0G1_1
  init_conditions$G0G1_5 = init_conditions$G0G1_1
  init_conditions$S_2 = init_conditions$S_1
  init_conditions$S_3 = init_conditions$S_1
  init_conditions$S_4 = init_conditions$S_1
  init_conditions$S_5 = init_conditions$S_1
  init_conditions$G2M_2 = init_conditions$G2M_1
  init_conditions$G2M_3 = init_conditions$G2M_1
  init_conditions$G2M_4 = init_conditions$G2M_1
  init_conditions$G2M_5 = init_conditions$G2M_1
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_15_even_ND(params = params_for_prediction[i,], 
                                                                                                    init = init_conditions[i,], 
                                                                                                    tmax = 10, 
                                                                                                    passage_day = passage_day,
                                                                                                    passage_number = passage_number,
                                                                                                    dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                    drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                      select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + S_1 + S_2 + S_3 + S_4 + S_5  + G2M_1 + G2M_2 + G2M_3 + G2M_4 + G2M_5) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells_15_even.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}

validation_schedules_10day_96well_palbo_only_15_even_ful_G2 = function(validation_fn, model_data, name,
                                                            params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                            num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  init_conditions = init_conditions/5
  colnames(init_conditions) = c("G0G1_1","S_1","G2M_1")
  init_conditions$G0G1_2 = init_conditions$G0G1_1
  init_conditions$G0G1_3 = init_conditions$G0G1_1
  init_conditions$G0G1_4 = init_conditions$G0G1_1
  init_conditions$G0G1_5 = init_conditions$G0G1_1
  init_conditions$S_2 = init_conditions$S_1
  init_conditions$S_3 = init_conditions$S_1
  init_conditions$S_4 = init_conditions$S_1
  init_conditions$S_5 = init_conditions$S_1
  init_conditions$G2M_2 = init_conditions$G2M_1
  init_conditions$G2M_3 = init_conditions$G2M_1
  init_conditions$G2M_4 = init_conditions$G2M_1
  init_conditions$G2M_5 = init_conditions$G2M_1
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_15_even_ful_G2(params = params_for_prediction[i,], 
                                                                                                init = init_conditions[i,], 
                                                                                                tmax = 10, 
                                                                                                passage_day = passage_day,
                                                                                                passage_number = passage_number,
                                                                                                dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                  select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + S_1 + S_2 + S_3 + S_4 + S_5  + G2M_1 + G2M_2 + G2M_3 + G2M_4 + G2M_5) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}

validation_schedules_10day_96well_palbo_only_15_even_ful_S = function(validation_fn, model_data, name,
                                                                       params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                                       num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  init_conditions = init_conditions/5
  colnames(init_conditions) = c("G0G1_1","S_1","G2M_1")
  init_conditions$G0G1_2 = init_conditions$G0G1_1
  init_conditions$G0G1_3 = init_conditions$G0G1_1
  init_conditions$G0G1_4 = init_conditions$G0G1_1
  init_conditions$G0G1_5 = init_conditions$G0G1_1
  init_conditions$S_2 = init_conditions$S_1
  init_conditions$S_3 = init_conditions$S_1
  init_conditions$S_4 = init_conditions$S_1
  init_conditions$S_5 = init_conditions$S_1
  init_conditions$G2M_2 = init_conditions$G2M_1
  init_conditions$G2M_3 = init_conditions$G2M_1
  init_conditions$G2M_4 = init_conditions$G2M_1
  init_conditions$G2M_5 = init_conditions$G2M_1
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_15_even_ful_S(params = params_for_prediction[i,], 
                                                                                                           init = init_conditions[i,], 
                                                                                                           tmax = 10, 
                                                                                                           passage_day = passage_day,
                                                                                                           passage_number = passage_number,
                                                                                                           dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                           drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                             select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + S_1 + S_2 + S_3 + S_4 + S_5  + G2M_1 + G2M_2 + G2M_3 + G2M_4 + G2M_5) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_14_G2 = function(validation_fn, model_data, name,
                                                              params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max"),
                                                              num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/6
  init_conditions$G0G1_2 = init_conditions$G0G1/6
  init_conditions$G0G1_3 = init_conditions$G0G1/6
  init_conditions$G0G1_4 = init_conditions$G0G1/6
  init_conditions$G0G1_5 = init_conditions$G0G1/6
  init_conditions$G0G1_6 = init_conditions$G0G1/6

  
  init_conditions$S_1 = init_conditions$S/6
  init_conditions$S_2 = init_conditions$S/6
  init_conditions$S_3 = init_conditions$S/6
  init_conditions$S_4 = init_conditions$S/6
  init_conditions$S_5 = init_conditions$S/6
  init_conditions$S_6 = init_conditions$S/6

  
  init_conditions$G2M_1 = init_conditions$G2M/2
  init_conditions$G2M_2 = init_conditions$G2M/2

  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2(params = params_for_prediction[i,], 
                                                                                                  init = init_conditions[i,], 
                                                                                                  tmax = 10, 
                                                                                                  passage_day = passage_day,
                                                                                                  passage_number = passage_number,
                                                                                                  dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                  drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                    select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6 + S_1 + S_2 + S_3 + S_4 + S_5  + S_6 + G2M_1 + G2M_2) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                                                           total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                                                           total_median = median(total, na.rm = T),
                                                                                                                                                                                           PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_14_G2_palbo = function(validation_fn, model_data, name,
                                                                  params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max", "b_P_2", "c_P_2", "alpha_max_2"),
                                                                  num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/6
  init_conditions$G0G1_2 = init_conditions$G0G1/6
  init_conditions$G0G1_3 = init_conditions$G0G1/6
  init_conditions$G0G1_4 = init_conditions$G0G1/6
  init_conditions$G0G1_5 = init_conditions$G0G1/6
  init_conditions$G0G1_6 = init_conditions$G0G1/6

  
  init_conditions$S_1 = init_conditions$S/6
  init_conditions$S_2 = init_conditions$S/6
  init_conditions$S_3 = init_conditions$S/6
  init_conditions$S_4 = init_conditions$S/6
  init_conditions$S_5 = init_conditions$S/6
  init_conditions$S_6 = init_conditions$S/6

  
  init_conditions$G2M_1 = init_conditions$G2M/2
  init_conditions$G2M_2 = init_conditions$G2M/2

  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2_palbo(params = params_for_prediction[i,], 
                                                                                                      init = init_conditions[i,], 
                                                                                                      tmax = 10, 
                                                                                                      passage_day = passage_day,
                                                                                                      passage_number = passage_number,
                                                                                                      dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                      drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                        select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6  + S_1 + S_2 + S_3 + S_4 + S_5  + S_6  + G2M_1 + G2M_2 ) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                                                           total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                                                           total_median = median(total, na.rm = T),
                                                                                                                                                                                           PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_14_G2_ful = function(validation_fn, model_data, name,
                                                                    params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                                    num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/6
  init_conditions$G0G1_2 = init_conditions$G0G1/6
  init_conditions$G0G1_3 = init_conditions$G0G1/6
  init_conditions$G0G1_4 = init_conditions$G0G1/6
  init_conditions$G0G1_5 = init_conditions$G0G1/6
  init_conditions$G0G1_6 = init_conditions$G0G1/6
  
  
  init_conditions$S_1 = init_conditions$S/6
  init_conditions$S_2 = init_conditions$S/6
  init_conditions$S_3 = init_conditions$S/6
  init_conditions$S_4 = init_conditions$S/6
  init_conditions$S_5 = init_conditions$S/6
  init_conditions$S_6 = init_conditions$S/6
  
  
  init_conditions$G2M_1 = init_conditions$G2M/2
  init_conditions$G2M_2 = init_conditions$G2M/2
  
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2_ful(params = params_for_prediction[i,], 
                                                                                                        init = init_conditions[i,], 
                                                                                                        tmax = 10, 
                                                                                                        passage_day = passage_day,
                                                                                                        passage_number = passage_number,
                                                                                                        dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                        drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                          select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6  + S_1 + S_2 + S_3 + S_4 + S_5  + S_6  + G2M_1 + G2M_2 ) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                total_median = median(total, na.rm = T),
                                                                                                                                                PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_16_G2_ND = function(validation_fn, model_data, name,
                                                                  params = c("alpha","b_P","b_F","c_P","c_F","a_FP"),
                                                                  num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e14, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/7
  init_conditions$G0G1_2 = init_conditions$G0G1/7
  init_conditions$G0G1_3 = init_conditions$G0G1/7
  init_conditions$G0G1_4 = init_conditions$G0G1/7
  init_conditions$G0G1_5 = init_conditions$G0G1/7
  init_conditions$G0G1_6 = init_conditions$G0G1/7
  init_conditions$G0G1_7 = init_conditions$G0G1/7
  
  
  init_conditions$S_1 = init_conditions$S/7
  init_conditions$S_2 = init_conditions$S/7
  init_conditions$S_3 = init_conditions$S/7
  init_conditions$S_4 = init_conditions$S/7
  init_conditions$S_5 = init_conditions$S/7
  init_conditions$S_6 = init_conditions$S/7
  init_conditions$S_7 = init_conditions$S/7
  
  init_conditions$G2M_1 = init_conditions$G2M/2
  init_conditions$G2M_2 = init_conditions$G2M/2
  
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_16_G2_ND(params = params_for_prediction[i,], 
                                                                                                      init = init_conditions[i,], 
                                                                                                      tmax = 10, 
                                                                                                      passage_day = passage_day,
                                                                                                      passage_number = passage_number,
                                                                                                      dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                      drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                        select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6 + G0G1_7  + S_1 + S_2 + S_3 + S_4 + S_5  + S_6 + S_7  + G2M_1 + G2M_2 ) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                               total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                               total_median = median(total, na.rm = T),
                                                                                                                                                               PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}




validation_schedules_10day_96well_palbo_only_16_G2_ful = function(validation_fn, model_data, name,
                                                                  params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                                  num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/7
  init_conditions$G0G1_2 = init_conditions$G0G1/7
  init_conditions$G0G1_3 = init_conditions$G0G1/7
  init_conditions$G0G1_4 = init_conditions$G0G1/7
  init_conditions$G0G1_5 = init_conditions$G0G1/7
  init_conditions$G0G1_6 = init_conditions$G0G1/7
  init_conditions$G0G1_7 = init_conditions$G0G1/7
  
  
  init_conditions$S_1 = init_conditions$S/7
  init_conditions$S_2 = init_conditions$S/7
  init_conditions$S_3 = init_conditions$S/7
  init_conditions$S_4 = init_conditions$S/7
  init_conditions$S_5 = init_conditions$S/7
  init_conditions$S_6 = init_conditions$S/7
  init_conditions$S_7 = init_conditions$S/7
  
  init_conditions$G2M_1 = init_conditions$G2M/2
  init_conditions$G2M_2 = init_conditions$G2M/2
  
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_16_G2_ful(params = params_for_prediction[i,], 
                                                                                                      init = init_conditions[i,], 
                                                                                                      tmax = 10, 
                                                                                                      passage_day = passage_day,
                                                                                                      passage_number = passage_number,
                                                                                                      dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                      drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                        select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6 + G0G1_7  + S_1 + S_2 + S_3 + S_4 + S_5  + S_6 + S_7  + G2M_1 + G2M_2 ) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                total_median = median(total, na.rm = T),
                                                                                                                                                PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}



validation_schedules_10day_96well_palbo_only_14_G2_both_G2_nonsyn = function(validation_fn, model_data, name,
                                                                    params = c("alpha","b_P","b_F","c_P","c_F", "alpha_max", "b_P_2", "b_F_2", "c_P_2", "c_F_2", "alpha_max_2"),
                                                                    num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/6
  init_conditions$G0G1_2 = init_conditions$G0G1/6
  init_conditions$G0G1_3 = init_conditions$G0G1/6
  init_conditions$G0G1_4 = init_conditions$G0G1/6
  init_conditions$G0G1_5 = init_conditions$G0G1/6
  init_conditions$G0G1_6 = init_conditions$G0G1/6
  
  
  init_conditions$S_1 = init_conditions$S/6
  init_conditions$S_2 = init_conditions$S/6
  init_conditions$S_3 = init_conditions$S/6
  init_conditions$S_4 = init_conditions$S/6
  init_conditions$S_5 = init_conditions$S/6
  init_conditions$S_6 = init_conditions$S/6
  
  
  init_conditions$G2M_1 = init_conditions$G2M/2
  init_conditions$G2M_2 = init_conditions$G2M/2
  
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2_both_G2_nonsyn(params = params_for_prediction[i,], 
                                                                                                        init = init_conditions[i,], 
                                                                                                        tmax = 10, 
                                                                                                        passage_day = passage_day,
                                                                                                        passage_number = passage_number,
                                                                                                        dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                        drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                          select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6  + S_1 + S_2 + S_3 + S_4 + S_5  + S_6  + G2M_1 + G2M_2 ) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                total_median = median(total, na.rm = T),
                                                                                                                                                PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_16_G2_both_G2_nonsyn = function(validation_fn, model_data, name,
                                                                             params = c("alpha","b_P","b_F","c_P","c_F", "alpha_max", "b_P_2", "b_F_2", "c_P_2", "c_F_2", "alpha_max_2"),
                                                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/7
  init_conditions$G0G1_2 = init_conditions$G0G1/7
  init_conditions$G0G1_3 = init_conditions$G0G1/7
  init_conditions$G0G1_4 = init_conditions$G0G1/7
  init_conditions$G0G1_5 = init_conditions$G0G1/7
  init_conditions$G0G1_6 = init_conditions$G0G1/7
  init_conditions$G0G1_7 = init_conditions$G0G1/7
  
  init_conditions$S_1 = init_conditions$S/7
  init_conditions$S_2 = init_conditions$S/7
  init_conditions$S_3 = init_conditions$S/7
  init_conditions$S_4 = init_conditions$S/7
  init_conditions$S_5 = init_conditions$S/7
  init_conditions$S_6 = init_conditions$S/7
  init_conditions$S_7 = init_conditions$S/7
  
  init_conditions$G2M_1 = init_conditions$G2M/2
  init_conditions$G2M_2 = init_conditions$G2M/2
  
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_16_G2_both_G2_nonsyn(params = params_for_prediction[i,], 
                                                                                                                 init = init_conditions[i,], 
                                                                                                                 tmax = 10, 
                                                                                                                 passage_day = passage_day,
                                                                                                                 passage_number = passage_number,
                                                                                                                 dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                                 drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                                   select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6 + G0G1_7  + S_1 + S_2 + S_3 + S_4 + S_5  + S_6 + S_7  + G2M_1 + G2M_2 ) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                total_median = median(total, na.rm = T),
                                                                                                                                                PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_16_G2_both_G2_nonsyn_ND = function(validation_fn, model_data, name,
                                                                             params = c("alpha","b_P","b_F","c_P","c_F", "alpha_max", "b_P_2", "b_F_2", "c_P_2", "c_F_2", "alpha_max_2"),
                                                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e14, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/7
  init_conditions$G0G1_2 = init_conditions$G0G1/7
  init_conditions$G0G1_3 = init_conditions$G0G1/7
  init_conditions$G0G1_4 = init_conditions$G0G1/7
  init_conditions$G0G1_5 = init_conditions$G0G1/7
  init_conditions$G0G1_6 = init_conditions$G0G1/7
  init_conditions$G0G1_7 = init_conditions$G0G1/7
  
  init_conditions$S_1 = init_conditions$S/7
  init_conditions$S_2 = init_conditions$S/7
  init_conditions$S_3 = init_conditions$S/7
  init_conditions$S_4 = init_conditions$S/7
  init_conditions$S_5 = init_conditions$S/7
  init_conditions$S_6 = init_conditions$S/7
  init_conditions$S_7 = init_conditions$S/7
  
  init_conditions$G2M_1 = init_conditions$G2M/2
  init_conditions$G2M_2 = init_conditions$G2M/2
  
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_16_G2_both_G2_nonsyn_ND(params = params_for_prediction[i,], 
                                                                                                                 init = init_conditions[i,], 
                                                                                                                 tmax = 10, 
                                                                                                                 passage_day = passage_day,
                                                                                                                 passage_number = passage_number,
                                                                                                                 dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                                 drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                                   select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6 + G0G1_7  + S_1 + S_2 + S_3 + S_4 + S_5  + S_6 + S_7  + G2M_1 + G2M_2 ) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                               total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                               total_median = median(total, na.rm = T),
                                                                                                                                                               PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_14_G2_both_G2_syn = function(validation_fn, model_data, name,
                                                                             params = c("alpha","b_P","b_F","c_P","c_F", "a_FP", "a_FP_2", "alpha_max", "b_P_2", "b_F_2", "c_P_2", "c_F_2", "alpha_max_2"),
                                                                             num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/6
  init_conditions$G0G1_2 = init_conditions$G0G1/6
  init_conditions$G0G1_3 = init_conditions$G0G1/6
  init_conditions$G0G1_4 = init_conditions$G0G1/6
  init_conditions$G0G1_5 = init_conditions$G0G1/6
  init_conditions$G0G1_6 = init_conditions$G0G1/6
  
  
  init_conditions$S_1 = init_conditions$S/6
  init_conditions$S_2 = init_conditions$S/6
  init_conditions$S_3 = init_conditions$S/6
  init_conditions$S_4 = init_conditions$S/6
  init_conditions$S_5 = init_conditions$S/6
  init_conditions$S_6 = init_conditions$S/6
  
  
  init_conditions$G2M_1 = init_conditions$G2M/2
  init_conditions$G2M_2 = init_conditions$G2M/2
  
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2_both_G2_syn(params = params_for_prediction[i,], 
                                                                                                                 init = init_conditions[i,], 
                                                                                                                 tmax = 10, 
                                                                                                                 passage_day = passage_day,
                                                                                                                 passage_number = passage_number,
                                                                                                                 dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                                 drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                                   select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6  + S_1 + S_2 + S_3 + S_4 + S_5  + S_6  + G2M_1 + G2M_2 ) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                total_median = median(total, na.rm = T),
                                                                                                                                                PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}



validation_schedules_10day_96well_palbo_only = function(validation_fn, model_data, name,
                                                        params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max"),
                                                        num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  init_conditions = init_conditions/5
  colnames(init_conditions) = c("G0G1_1","S_1","G2M_1")
  init_conditions$G0G1_2 = init_conditions$G0G1_1
  init_conditions$G0G1_3 = init_conditions$G0G1_1
  init_conditions$G0G1_4 = init_conditions$G0G1_1
  init_conditions$G0G1_5 = init_conditions$G0G1_1
  init_conditions$S_2 = init_conditions$S_1
  init_conditions$S_3 = init_conditions$S_1
  init_conditions$S_4 = init_conditions$S_1
  init_conditions$S_5 = init_conditions$S_1
  init_conditions$G2M_2 = init_conditions$G2M_1
  init_conditions$G2M_3 = init_conditions$G2M_1
  init_conditions$G2M_4 = init_conditions$G2M_1
  init_conditions$G2M_5 = init_conditions$G2M_1
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper(params = params_for_prediction[i,], 
                                                                                            init = init_conditions[i,], 
                                                                                            tmax = 10, 
                                                                                            passage_day = passage_day,
                                                                                            passage_number = passage_number,
                                                                                            dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                            drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                              select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + S_1 + S_2 + S_3 + S_4 + S_5  + G2M_1 + G2M_2 + G2M_3 + G2M_4 + G2M_5) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_ful = function(validation_fn, model_data, name,
                                                  params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                  num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  init_conditions = init_conditions/5
  colnames(init_conditions) = c("G0G1_1","S_1","G2M_1")
  init_conditions$G0G1_2 = init_conditions$G0G1_1
  init_conditions$G0G1_3 = init_conditions$G0G1_1
  init_conditions$G0G1_4 = init_conditions$G0G1_1
  init_conditions$G0G1_5 = init_conditions$G0G1_1
  init_conditions$S_2 = init_conditions$S_1
  init_conditions$S_3 = init_conditions$S_1
  init_conditions$S_4 = init_conditions$S_1
  init_conditions$S_5 = init_conditions$S_1
  init_conditions$G2M_2 = init_conditions$G2M_1
  init_conditions$G2M_3 = init_conditions$G2M_1
  init_conditions$G2M_4 = init_conditions$G2M_1
  init_conditions$G2M_5 = init_conditions$G2M_1
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_ful(params = params_for_prediction[i,], 
                                                                                      init = init_conditions[i,], 
                                                                                      tmax = 10, 
                                                                                      passage_day = passage_day,
                                                                                      passage_number = passage_number,
                                                                                      dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                      drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                        select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + S_1 + S_2 + S_3 + S_4 + S_5  + G2M_1 + G2M_2 + G2M_3 + G2M_4 + G2M_5) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                   total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                   total_median = median(total, na.rm = T),
                                                   PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}



validation_schedules_10day_96well_palbo_only_20_G2 = function(validation_fn, model_data, name,
                                                        params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max"),
                                                        num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/8
  init_conditions$G0G1_2 = init_conditions$G0G1/8
  init_conditions$G0G1_3 = init_conditions$G0G1/8
  init_conditions$G0G1_4 = init_conditions$G0G1/8
  init_conditions$G0G1_5 = init_conditions$G0G1/8
  init_conditions$G0G1_6 = init_conditions$G0G1/8
  init_conditions$G0G1_7 = init_conditions$G0G1/8
  init_conditions$G0G1_8 = init_conditions$G0G1/8

  init_conditions$S_1 = init_conditions$S/8
  init_conditions$S_2 = init_conditions$S/8
  init_conditions$S_3 = init_conditions$S/8
  init_conditions$S_4 = init_conditions$S/8
  init_conditions$S_5 = init_conditions$S/8
  init_conditions$S_6 = init_conditions$S/8
  init_conditions$S_7 = init_conditions$S/8
  init_conditions$S_8 = init_conditions$S/8

  init_conditions$G2M_1 = init_conditions$G2M/4
  init_conditions$G2M_2 = init_conditions$G2M/4
  init_conditions$G2M_3 = init_conditions$G2M/4
  init_conditions$G2M_4 = init_conditions$G2M/4
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_20_G2(params = params_for_prediction[i,], 
                                                                                            init = init_conditions[i,], 
                                                                                            tmax = 10, 
                                                                                            passage_day = passage_day,
                                                                                            passage_number = passage_number,
                                                                                            dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                            drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                              select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6 + G0G1_7 + G0G1_8 + S_1 + S_2 + S_3 + S_4 + S_5  + S_6 + S_7 + S_8 + G2M_1 + G2M_2 + G2M_3 + G2M_4) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}


validation_schedules_10day_96well_palbo_only_20_G2_ful = function(validation_fn, model_data, name,
                                                            params = c("alpha","b_P","b_F","c_P","c_F","a_FP","alpha_max", "b_F_2", "c_F_2", "alpha_max_2"),
                                                            num_iters = 1500, dox = "Dox", passage_day = 5.1, passage_number = 5000){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e9, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$G0G1_1 = init_conditions$G0G1/8
  init_conditions$G0G1_2 = init_conditions$G0G1/8
  init_conditions$G0G1_3 = init_conditions$G0G1/8
  init_conditions$G0G1_4 = init_conditions$G0G1/8
  init_conditions$G0G1_5 = init_conditions$G0G1/8
  init_conditions$G0G1_6 = init_conditions$G0G1/8
  init_conditions$G0G1_7 = init_conditions$G0G1/8
  init_conditions$G0G1_8 = init_conditions$G0G1/8
  
  init_conditions$S_1 = init_conditions$S/8
  init_conditions$S_2 = init_conditions$S/8
  init_conditions$S_3 = init_conditions$S/8
  init_conditions$S_4 = init_conditions$S/8
  init_conditions$S_5 = init_conditions$S/8
  init_conditions$S_6 = init_conditions$S/8
  init_conditions$S_7 = init_conditions$S/8
  init_conditions$S_8 = init_conditions$S/8
  
  init_conditions$G2M_1 = init_conditions$G2M/4
  init_conditions$G2M_2 = init_conditions$G2M/4
  init_conditions$G2M_3 = init_conditions$G2M/4
  init_conditions$G2M_4 = init_conditions$G2M/4
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions = subset(init_conditions, select = -c(1,2,3) )
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction = params_for_prediction %>% mutate(b = alpha_max)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_20_G2_ful(params = params_for_prediction[i,], 
                                                                                                init = init_conditions[i,], 
                                                                                                tmax = 10, 
                                                                                                passage_day = passage_day,
                                                                                                passage_number = passage_number,
                                                                                                dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                                drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                                  select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1_1 + G0G1_2 + G0G1_3 + G0G1_4 + G0G1_5 + G0G1_6 + G0G1_7 + G0G1_8 + S_1 + S_2 + S_3 + S_4 + S_5  + S_6 + S_7 + S_8 + G2M_1 + G2M_2 + G2M_3 + G2M_4) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                                                                                                                      total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                                                                                                                      total_median = median(total, na.rm = T),
                                                                                                                                                      PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("posterior_predictions_",name, "_posterior-predictions-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_obs_long, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 4)) 
  
  print(ggplot(prediction_results %>% filter(Experiment == "DMSO"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_obs_long %>% filter(Experiment == "DMSO"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
}



validation_schedules_10day_96well_palbo_only_nopassaging = function(validation_fn, model_data, name,
                                                        params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_min","Gamma","delta","K"),
                                                        num_iters = 1500, dox = "Dox"){
  valid1 = read.csv(validation_fn)
  valid1_tx = valid1 %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e10, Palbociclib = Palbociclib * 1e8)
  valid1_obs = valid1 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  valid1_obs_long = valid1_obs %>% pivot_longer(cols = contains("count"), names_to = "replicate", values_to = "total_count")
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  params_for_prediction = params_for_prediction %>% mutate(b = b_min, K = 1e10)
  
  prediction_results = data.frame()
  for (exp in unique(valid1_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_no_passaging_validation_wrapper(params = params_for_prediction[i,], 
                                                                                            init = init_conditions[i,], 
                                                                                            tmax = 10, 
                                                                                            dose_times = (valid1_tx %>% filter(Experiment == exp))$Day,
                                                                                            drug_amount = valid1_tx %>% filter(Experiment == exp) %>% 
                                                                                              select(Fulvestrant,Palbociclib))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                   total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                   total_median = median(total, na.rm = T),
                                                   PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name, "_posterior-predictions-no-passaging-10day-96wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median, color = Experiment)) + geom_line() + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL, fill = Experiment), alpha = 0.25)) 
  dev.off()
  return(prediction_results)
}


fourweek_schedule_prediction_no_passaging_24well_palbo_only = function(validation_fn, model_data, name,
                                                  params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_min","Gamma","delta","K"),
                                                  num_iters = 1500, dox = "Dox"){
  valid1_24wells = read.csv(validation_fn)
  valid1_24wells_tx = valid1_24wells %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e10, Palbociclib = Palbociclib * 1e8, Abemaciclib = Abemaciclib*1e6)
  valid1_24wells_obs = valid1_24wells %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions_24wells = init_conditions %>% 
    mutate(total = G0G1 + S + G2M, G0G1 = G0G1/total, S = S/total, G2M = G2M/total) %>% 
    mutate(total = mean((valid1_24wells_obs %>% filter(Day == 0, Experiment == "Control"))$total_count)) %>% 
    mutate(G0G1 = G0G1 * total,
           S = S * total,
           G2M = G2M * total) %>%
    select(-total)
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  params_for_prediction_24wells = params_for_prediction %>% mutate(K = 6)
  params_for_prediction_24wells = params_for_prediction_24wells %>% mutate(b = b_min, c_A = 1, b_A = 1, a_FA = 0)
  
  prediction_results_24wells = data.frame()
  for (exp in unique(valid1_24wells_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper_no_passaging(params = params_for_prediction_24wells[i,], 
                                                                                      init = init_conditions_24wells[i,], 
                                                                                      tmax = 28, 
                                                                                      dose_times = (valid1_24wells_tx %>% filter(Experiment == exp))$Day,
                                                                                      drug_amount = valid1_24wells_tx %>% filter(Experiment == exp) %>% 
                                                                                        select (Fulvestrant,Palbociclib,Abemaciclib,Passage))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                   total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                   total_median = median(total, na.rm = T),
                                                   PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results_24wells = bind_rows(prediction_results_24wells, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name, "_posterior-predictions-no_passaging-24wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results_24wells, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_24wells_obs, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 3)) 
  
  print(ggplot(prediction_results_24wells %>% filter(Experiment == "Control"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_24wells_obs %>% filter(Experiment == "Control"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
  return(prediction_results_24wells)
}


validation_schedules_loglog_bmin_24well_palbo_only = function(validation_fn, model_data, name,
                                       params = c("alpha","b_P","b_F","c_P","c_F","a_FP","b_max","b_min","Gamma","delta","K"),
                                       num_iters = 1500, dox = "Dox"){
  valid1_24wells = read.csv(validation_fn)
  valid1_24wells_tx = valid1_24wells %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e10, Palbociclib = Palbociclib * 1e8)
  valid1_24wells_obs = valid1_24wells %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  #init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions_24wells = init_conditions %>% 
    mutate(total = G0G1 + S + G2M, G0G1 = G0G1/total, S = S/total, G2M = G2M/total) %>% 
    mutate(total = mean((valid1_24wells_obs %>% filter(Day == 0, Experiment == "Control"))$total_count)) %>% 
    mutate(G0G1 = G0G1 * total,
           S = S * total,
           G2M = G2M * total) %>%
    select(-total)
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  params_for_prediction_24wells = params_for_prediction %>% mutate(K = 6)
  
  prediction_results_24wells = data.frame()
  for (exp in unique(valid1_24wells_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_bmin_palboOnly_validation_wrapper(params = params_for_prediction_24wells[i,], 
                                                                                      init = init_conditions_24wells[i,], 
                                                                                      tmax = 28, 
                                                                                      dose_times = (valid1_24wells_tx %>% filter(Experiment == exp))$Day,
                                                                                      drug_amount = valid1_24wells_tx %>% filter(Experiment == exp) %>% 
                                                                                        select (Fulvestrant,Palbociclib,Passage))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                   total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                   total_median = median(total, na.rm = T),
                                                   PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results_24wells = bind_rows(prediction_results_24wells, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name, "_posterior-predictions-passaging-24wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results_24wells, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_24wells_obs, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 3)) 
  
  print(ggplot(prediction_results_24wells %>% filter(Experiment == "Control"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_24wells_obs %>% filter(Experiment == "Control"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
  
}


simulate_validation_schedules_24well = function(validation_fn, model_data, name,
                                       params,
                                       num_iters = 2250, dox = "Dox"){
  valid1_24wells = read.csv(validation_fn)
  valid1_24wells_tx = valid1_24wells %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e6, Palbociclib = Palbociclib * 1e6, Abemaciclib = Abemaciclib*1e6)
  valid1_24wells_obs = valid1_24wells %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions_24wells = init_conditions %>% 
    mutate(total = G0G1 + S + G2M, G0G1 = G0G1/total, S = S/total, G2M = G2M/total) %>% 
    mutate(total = mean((valid1_24wells_obs %>% filter(Day == 0, Experiment == "Control"))$total_count)) %>% 
    mutate(G0G1 = G0G1 * total,
           S = S * total,
           G2M = G2M * total) %>%
    select(-total)
  
  #params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction_24wells = params_for_prediction %>% mutate(K = 6)
  
  prediction_results_24wells = data.frame()
  for (exp in unique(valid1_24wells_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper(params = params[i,], 
                                                                                      init = init_conditions_24wells[i,], 
                                                                                      tmax = 28, 
                                                                                      dose_times = (valid1_24wells_tx %>% filter(Experiment == exp))$Day,
                                                                                      drug_amount = valid1_24wells_tx %>% filter(Experiment == exp) %>% 
                                                                                        select (Fulvestrant,Palbociclib,Abemaciclib,Passage))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                   total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                   total_median = median(total, na.rm = T),
                                                   AA = mean(AA), PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results_24wells = bind_rows(prediction_results_24wells, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name, "_posterior-predictions-passaging-24wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results_24wells, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_24wells_obs, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 3)) 
  
  print(ggplot(prediction_results_24wells %>% filter(Experiment == "Control"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_24wells_obs %>% filter(Experiment == "Control"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
  
  return(prediction_results_24wells)
}

simulate_validation_schedules_24well_oneSample = function(validation_fn, model_data,
                                                params,
                                                dox = "Dox"){
  valid1_24wells = read.csv(validation_fn)
  valid1_24wells_tx = valid1_24wells %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e6, Palbociclib = Palbociclib * 1e6, Abemaciclib = Abemaciclib*1e6)
  valid1_24wells_obs = valid1_24wells %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions_24wells = init_conditions %>% 
    mutate(total = G0G1 + S + G2M, G0G1 = G0G1/total, S = S/total, G2M = G2M/total) %>% 
    mutate(total = mean((valid1_24wells_obs %>% filter(Day == 0, Experiment == "Control"))$total_count)) %>% 
    mutate(G0G1 = G0G1 * total,
           S = S * total,
           G2M = G2M * total) %>%
    select(-total)
  
  #params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction_24wells = params_for_prediction %>% mutate(K = 6)
  
  prediction_results_24wells = data.frame()
  for (exp in unique(valid1_24wells_tx$Experiment)){
    pred = ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper(params = params, 
                                                                                      init = init_conditions_24wells[1,], 
                                                                                      tmax = 28, 
                                                                                      dose_times = (valid1_24wells_tx %>% filter(Experiment == exp))$Day,
                                                                                      drug_amount = valid1_24wells_tx %>% filter(Experiment == exp) %>% 
                                                                                        select (Fulvestrant,Palbociclib,Abemaciclib,Passage))
    pred_summary = bind_rows(pred, .id = "rep") %>% 
      mutate(total = G0G1 + S + G2M)  %>% mutate(Experiment = exp)
    prediction_results_24wells = bind_rows(prediction_results_24wells, pred_summary)
  }
  print(ggplot(prediction_results_24wells, aes(x = time, y = total)) + geom_line() + 
          geom_point(data = valid1_24wells_obs, aes(x = Day, y = total_count))  + 
          facet_wrap(Experiment ~. , ncol = 3)) 
  
  #return(prediction_results_24wells)
}

simulate_validation_schedules_96well_test_prior_means = function(validation_fn, model_data,
                                                          params,
                                                          dox = "Dox",
                                                          drug2 = "Palbociclib"){
  valid1_24wells = read.csv(validation_fn)
  if(!("Palbociclib" %in% names(valid1_24wells))){
    valid1_24wells = valid1_24wells %>% mutate(Palbociclib = 0)
  } else{
    valid1_24wells = valid1_24wells %>% mutate(Abemaciclib = 0)
  }
  valid1_24wells_tx = valid1_24wells %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e6, Palbociclib = Palbociclib * 1e6, Abemaciclib = Abemaciclib*1e6)
  valid1_24wells_obs = valid1_24wells %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  prediction_results_24wells = data.frame()
  for (exp in unique(valid1_24wells_tx$Experiment)){
    pred = ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper(params = params, 
                                                                                           init = init_conditions[1,], 
                                                                                           tmax = 6, 
                                                                                           dose_times = (valid1_24wells_tx %>% filter(Experiment == exp))$Day,
                                                                                           drug_amount = valid1_24wells_tx %>% filter(Experiment == exp) %>% 
                                                                                             select (Fulvestrant,Palbociclib,Abemaciclib,Passage))
    pred_summary = bind_rows(pred, .id = "rep") %>% 
      mutate(total = G0G1 + S + G2M)  %>% mutate(Experiment = exp)
    prediction_results_24wells = bind_rows(prediction_results_24wells, pred_summary)
  }
  prediction_results_24wells = prediction_results_24wells %>% separate(Experiment, into = c("Fulvestrant", drug2))
  valid1_24wells_obs = valid1_24wells_obs %>% select(-Fulvestrant, -Abemaciclib, -Palbociclib) %>% separate(Experiment, into = c("Fulvestrant", drug2))
  print(ggplot(prediction_results_24wells, aes(x = time, y = total)) + geom_line() + 
          geom_point(data = valid1_24wells_obs, aes(x = Day, y = total_count))  + 
          facet_grid(Palbociclib ~ Fulvestrant , labeller = label_both, scales = "free_y" )) 
  
  #return(prediction_results_24wells)
}

validation_schedules_96well_differentIC50 = function(validation_fn, model_data, name,
                                       params= c("alpha","b_P","b_F","b_A","c_P","c_F","c_A","c_P_eff","c_A_eff","a_FP","a_FA","b","Gamma","delta","K"),
                                       num_iters = 2250, dox = "Dox"){
  valid2 = read.csv(validation_fn)
  valid2_tx = valid2 %>% filter(Measurement == "Treatment", Dox == dox) %>% arrange(Experiment, Day) %>% 
    mutate(Fulvestrant = Fulvestrant * 1e6, Palbociclib = Palbociclib * 1e6, Abemaciclib = Abemaciclib*1e6)
  valid2_obs = valid2 %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  
  init_conditions = as.data.frame(model_data, pars = c("y0_pass_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  ## with passaging
  prediction_results = data.frame()
  for (exp in unique(valid2_tx$Experiment)){
    pred = lapply(1:num_iters, 
                  function(i){
                    ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_differentIC50_validation_wrapper(params = params_for_prediction[i,], 
                                                                                                    init = init_conditions[i,], 
                                                                                                    tmax = 28, 
                                                                                                    dose_times = (valid2_tx %>% filter(Experiment == exp))$Day,
                                                                                                    drug_amount = valid2_tx %>% filter(Experiment == exp) %>% 
                                                                                                      select (Fulvestrant,Palbociclib,Abemaciclib,Passage)) })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% 
      summarize(total_mean = mean(total), total_low = quantile(total, 0.025, na.rm = T), 
                total_high = quantile(total, 0.975, na.rm = T), total_median = median(total, na.rm = T),
                AA = mean(AA), PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results = bind_rows(prediction_results, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name,"_posterior-predictions-passaging.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid2_obs, aes(x = Day, y = total_count))  +
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + facet_wrap(Experiment ~. , ncol = 5))
  dev.off()
  
  prediction_results_no_passaging = data.frame()
  for (exp in unique(valid2_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_differentIC50_validation_wrapper(params = params_for_prediction[i,], 
                                                                                                   init = init_conditions[i,], 
                                                                                                   tmax = 28, 
                                                                                                   dose_times = (valid2_tx %>% filter(Experiment == exp))$Day,
                                                                                                   drug_amount = valid2_tx %>% filter(Experiment == exp) %>% 
                                                                                                     select (Fulvestrant,Palbociclib,Abemaciclib,Passage))})
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>%
      mutate(total = G0G1 + S + G2M) %>% 
      summarize(total_mean = mean(total, na.rm = T), total_low = quantile(total, 0.025, na.rm = T), 
                total_high = quantile(total, 0.975, na.rm = T), total_median = median(total, na.rm = T),
                AA = mean(AA), PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results_no_passaging = bind_rows(prediction_results_no_passaging, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_", name, "_posterior-predictions-noPassaging.pdf",sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results_no_passaging, aes(x = time, y = total_median)) + geom_line()  + geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + facet_wrap(Experiment ~. , ncol = 5) +scale_y_log10())
  print(ggplot(prediction_results_no_passaging, aes(x = time, y = total_median, color = Experiment)) + geom_line() + scale_y_log10())
  dev.off()
  
  write.csv(prediction_results_no_passaging %>% filter(time ==28) %>% arrange(-total_median),
            paste("from_cluster_include_longterm_data/plots/posterior_predictions_", name, "_posterior-predictions_noPassaging_ranks.csv", sep = "", collapse = ""),
            quote = FALSE, row.names = FALSE)
  
}


validation_schedules_24well_differentIC50 = function(validation_fn, model_data, name,
                                       params = c("alpha","b_P","b_F","b_A","c_P","c_F","c_A","c_P_eff","c_A_eff","a_FP","a_FA","b","Gamma","delta","K"),
                                       num_iters = 2250, dox = "Dox"){
  valid1_24wells = read.csv(validation_fn)
  valid1_24wells_tx = valid1_24wells %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e6, Palbociclib = Palbociclib * 1e6, Abemaciclib = Abemaciclib*1e6)
  valid1_24wells_obs = valid1_24wells %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions_24wells = init_conditions %>% 
    mutate(total = G0G1 + S + G2M, G0G1 = G0G1/total, S = S/total, G2M = G2M/total) %>% 
    mutate(total = mean((valid1_24wells_obs %>% filter(Day == 0, Experiment == "Control"))$total_count)) %>% 
    mutate(G0G1 = G0G1 * total,
           S = S * total,
           G2M = G2M * total) %>%
    select(-total)
  
  params_for_prediction = as.data.frame(model_data, pars = params) 
  names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  params_for_prediction_24wells = params_for_prediction %>% mutate(K = 6)
  
  prediction_results_24wells = data.frame()
  for (exp in unique(valid1_24wells_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_differentIC50_validation_wrapper(params = params_for_prediction_24wells[i,], 
                                                                                      init = init_conditions_24wells[i,], 
                                                                                      tmax = 28, 
                                                                                      dose_times = (valid1_24wells_tx %>% filter(Experiment == exp))$Day,
                                                                                      drug_amount = valid1_24wells_tx %>% filter(Experiment == exp) %>% 
                                                                                        select (Fulvestrant,Palbociclib,Abemaciclib,Passage))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                   total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                   total_median = median(total, na.rm = T),
                                                   AA = mean(AA), PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results_24wells = bind_rows(prediction_results_24wells, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name, "_posterior-predictions-passaging-24wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results_24wells, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_24wells_obs, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 3)) 
  
  print(ggplot(prediction_results_24wells %>% filter(Experiment == "Control"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_24wells_obs %>% filter(Experiment == "Control"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
  
}


simulate_validation_schedules_24well_differentIC50 = function(validation_fn, model_data, name,
                                                params,
                                                num_iters = 2250, dox = "Dox"){
  valid1_24wells = read.csv(validation_fn)
  valid1_24wells_tx = valid1_24wells %>% filter(Measurement == "Treatment", Dox == dox) %>% 
    arrange(Experiment, Day) %>% mutate(Fulvestrant = Fulvestrant * 1e6, Palbociclib = Palbociclib * 1e6, Abemaciclib = Abemaciclib*1e6)
  valid1_24wells_obs = valid1_24wells %>% filter(Measurement == "Observation", Dox == dox) %>% arrange(Experiment, Day)
  
  init_conditions = as.data.frame(model_data, pars = c("y0_ct_p")) 
  colnames(init_conditions) = c("G0G1","S","G2M")
  init_conditions$AA = 0
  init_conditions$PP = 0
  init_conditions$FF = 0
  
  init_conditions_24wells = init_conditions %>% 
    mutate(total = G0G1 + S + G2M, G0G1 = G0G1/total, S = S/total, G2M = G2M/total) %>% 
    mutate(total = mean((valid1_24wells_obs %>% filter(Day == 0, Experiment == "Control"))$total_count)) %>% 
    mutate(G0G1 = G0G1 * total,
           S = S * total,
           G2M = G2M * total) %>%
    select(-total)
  
  #params_for_prediction = as.data.frame(model_data, pars = params) 
  #names(params_for_prediction) = str_replace(names(params_for_prediction), "gamma","Gamma")
  #params_for_prediction_24wells = params_for_prediction %>% mutate(K = 6)
  
  prediction_results_24wells = data.frame()
  for (exp in unique(valid1_24wells_tx$Experiment)){
    pred = lapply(1:num_iters, function(i){
      ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper(params = params[i,], 
                                                                                      init = init_conditions_24wells[i,], 
                                                                                      tmax = 28, 
                                                                                      dose_times = (valid1_24wells_tx %>% filter(Experiment == exp))$Day,
                                                                                      drug_amount = valid1_24wells_tx %>% filter(Experiment == exp) %>% 
                                                                                        select (Fulvestrant,Palbociclib,Abemaciclib,Passage))
    })
    pred_summary = bind_rows(pred, .id = "rep") %>% group_by(time) %>% 
      mutate(total = G0G1 + S + G2M) %>% summarize(total_mean = mean(total, na.rm = T), 
                                                   total_low = quantile(total, 0.025, na.rm = T), total_high = quantile(total, 0.975, na.rm = T), 
                                                   total_median = median(total, na.rm = T),
                                                   AA = mean(AA), PP = mean(PP), FF = mean(FF))  %>% mutate(Experiment = exp)
    prediction_results_24wells = bind_rows(prediction_results_24wells, pred_summary)
  }
  
  pdf(paste("from_cluster_include_longterm_data/plots/posterior_predictions_",name, "_posterior-predictions-passaging-24wells.pdf", sep = "", collapse = ""), useDingbats = F)
  print(ggplot(prediction_results_24wells, aes(x = time, y = total_median)) + geom_line() + 
          geom_point(data = valid1_24wells_obs, aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25) + 
          facet_wrap(Experiment ~. , ncol = 3)) 
  
  print(ggplot(prediction_results_24wells %>% filter(Experiment == "Control"), aes(x = time, y = total_median)) + 
          geom_line() + 
          geom_point(data = valid1_24wells_obs %>% filter(Experiment == "Control"), aes(x = Day, y = total_count))  + 
          geom_ribbon(aes(ymin = total_low, ymax = total_high, linetype = NULL), alpha = 0.25))
  dev.off()
  
  return(prediction_results_24wells)
}


# ## plot beta as a function of palbocilib
# params_for_prediction_ND %>% summarize(b = median(b), c_P = median(c_P), b_P = median(b_P))
# 
# f_palbo <- function(x) 2.544512*(1/(1 + (x/0.02920118)^0.9529028))
# 
# curve(f_palbo, 0 , 1)
# 
# params_for_prediction_ND %>% summarize(b = median(b), c_F = median(c_F), b_F = median(b_F))
# f_fulv <- function(x) 2.544512*(1/(1 + (x/0.005342141)^0.753689))
# curve(f_fulv, 0, 1)

