library(tidyverse)
library(deSolve)

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper <- function(params, init, tmax, dose_times, drug_amount){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0,0,1), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), AA = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["AA"] = drug_amount[i,3]
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    init["G0G1"] = drug_amount[i,4] * init["G0G1"]
    init["S"] = drug_amount[i,4] * init["S"]
    init["G2M"] = drug_amount[i,4] * init["G2M"]
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1","S","G2M","AA","PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_mod <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                  passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] + init["G2M_5"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
      init["G2M_5"] = (init["G2M_5"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_mod,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1","S_1","G2M_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_2", "S_3", "S_4", "S_5", "G2M_2", "G2M_3", "G2M_4", "G2M_5", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_multiple <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                  passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] + init["G2M_5"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
      init["G2M_5"] = (init["G2M_5"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_multiple,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1","S_1","G2M_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_2", "S_3", "S_4", "S_5", "G2M_2", "G2M_3", "G2M_4", "G2M_5", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_IC <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                  passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] + init["G2M_5"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
      init["G2M_5"] = (init["G2M_5"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_IC,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_1", "S_2", "S_3", "S_4", "S_5", "G2M_1", "G2M_2", "G2M_3", "G2M_4", "G2M_5", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_15_even <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                  passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] + init["G2M_5"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
      init["G2M_5"] = (init["G2M_5"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_15_even,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1","S_1","G2M_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_2", "S_3", "S_4", "S_5", "G2M_2", "G2M_3", "G2M_4", "G2M_5", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_15_even_ND <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                          passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] + init["G2M_5"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
      init["G2M_5"] = (init["G2M_5"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_15_even_ND,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1","S_1","G2M_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_2", "S_3", "S_4", "S_5", "G2M_2", "G2M_3", "G2M_4", "G2M_5", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_15_even_ful_G2 <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                      passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] + init["G2M_5"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
      init["G2M_5"] = (init["G2M_5"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_15_even_ful_G2,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1","S_1","G2M_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_2", "S_3", "S_4", "S_5", "G2M_2", "G2M_3", "G2M_4", "G2M_5", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_15_even_ful_S <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                                 passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] + init["G2M_5"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
      init["G2M_5"] = (init["G2M_5"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_15_even_ful_S,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1","S_1","G2M_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_2", "S_3", "S_4", "S_5", "G2M_2", "G2M_3", "G2M_4", "G2M_5", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2 <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                        passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"] + init["G2M_1"] + init["G2M_2"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number

      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number

      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number

    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "G2M_1", "G2M_2", "PP", "FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2_palbo <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                            passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"]  + init["G2M_1"] + init["G2M_2"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number

      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number

      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number

    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2_palbo,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "G2M_1", "G2M_2",  "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2_ful <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                              passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"]  + init["G2M_1"] + init["G2M_2"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number
      
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number
      
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2_ful,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "G2M_1", "G2M_2",  "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_16_G2_ND <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                            passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["G0G1_7"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"] + init["S_7"]  + init["G2M_1"] + init["G2M_2"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number
      init["G0G1_7"] = (init["G0G1_7"]/total_count)*passage_number
      
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number
      init["S_7"] = (init["S_7"]/total_count)*passage_number
      
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_16_G2_ND,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "G0G1_7", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "S_7", "G2M_1", "G2M_2",  "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_16_G2_ful <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                            passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["G0G1_7"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"] + init["S_7"]  + init["G2M_1"] + init["G2M_2"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number
      init["G0G1_7"] = (init["G0G1_7"]/total_count)*passage_number
      
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number
      init["S_7"] = (init["S_7"]/total_count)*passage_number
      
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_16_G2_ful,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "G0G1_7", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "S_7", "G2M_1", "G2M_2",  "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2_both_G2_nonsyn <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                              passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"]  + init["G2M_1"] + init["G2M_2"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number
      
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number
      
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2_both_G2_nonsyn,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "G2M_1", "G2M_2",  "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_16_G2_both_G2_nonsyn_ND <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                                       passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["G0G1_7"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"] + init["S_7"] + init["G2M_1"] + init["G2M_2"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number
      init["G0G1_7"] = (init["G0G1_7"]/total_count)*passage_number
      
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number
      init["S_7"] = (init["S_7"]/total_count)*passage_number
      
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_16_G2_both_G2_nonsyn_ND,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "G0G1_7", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6",  "S_7", "G2M_1", "G2M_2",  "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_16_G2_both_G2_nonsyn <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                                       passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["G0G1_7"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"] + init["S_7"] + init["G2M_1"] + init["G2M_2"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number
      init["G0G1_7"] = (init["G0G1_7"]/total_count)*passage_number
      
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number
      init["S_7"] = (init["S_7"]/total_count)*passage_number
      
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_16_G2_both_G2_nonsyn,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "G0G1_7", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6",  "S_7", "G2M_1", "G2M_2",  "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_14_G2_both_G2_syn <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                                       passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"]  + init["G2M_1"] + init["G2M_2"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number
      
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number
      
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2_both_G2_syn,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "G2M_1", "G2M_2",  "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}
ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                  passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }

    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] + init["G2M_5"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
      init["G2M_5"] = (init["G2M_5"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1","S_1","G2M_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_2", "S_3", "S_4", "S_5", "G2M_2", "G2M_3", "G2M_4", "G2M_5", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_ful <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                  passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] + init["G2M_5"]
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
      init["G2M_5"] = (init["G2M_5"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_15_even_ful_G2,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1","S_1","G2M_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "S_2", "S_3", "S_4", "S_5", "G2M_2", "G2M_3", "G2M_4", "G2M_5", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}



ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_20_G2 <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                  passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["G0G1_7"] + init["G0G1_8"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"] + init["S_7"] + init["S_8"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number
      init["G0G1_7"] = (init["G0G1_7"]/total_count)*passage_number
      init["G0G1_8"] = (init["G0G1_8"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number
      init["S_7"] = (init["S_7"]/total_count)*passage_number
      init["S_8"] = (init["S_8"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_20_G2,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "G0G1_7", "G0G1_8", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "S_7", "S_8", "G2M_1", "G2M_2", "G2M_3", "G2M_4", "PP", "FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_validation_wrapper_20_G2_ful <- function(params, init, tmax, dose_times, drug_amount, 
                                                                                                      passage_day, passage_number){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    if(t_start == passage_day){
      total_count = init["G0G1_1"] + init["G0G1_2"] + init["G0G1_3"] + init["G0G1_4"] + init["G0G1_5"] + init["G0G1_6"] + init["G0G1_7"] + init["G0G1_8"] + init["S_1"]+ init["S_2"] + init["S_3"] + init["S_4"] + init["S_5"] + init["S_6"] + init["S_7"] + init["S_8"] + init["G2M_1"] + init["G2M_2"] + init["G2M_3"] + init["G2M_4"] 
      init["G0G1_1"] = (init["G0G1_1"]/total_count)*passage_number
      init["G0G1_2"] = (init["G0G1_2"]/total_count)*passage_number
      init["G0G1_3"] = (init["G0G1_3"]/total_count)*passage_number
      init["G0G1_4"] = (init["G0G1_4"]/total_count)*passage_number
      init["G0G1_5"] = (init["G0G1_5"]/total_count)*passage_number
      init["G0G1_6"] = (init["G0G1_6"]/total_count)*passage_number
      init["G0G1_7"] = (init["G0G1_7"]/total_count)*passage_number
      init["G0G1_8"] = (init["G0G1_8"]/total_count)*passage_number
      init["S_1"] = (init["S_1"]/total_count)*passage_number
      init["S_2"] = (init["S_2"]/total_count)*passage_number
      init["S_3"] = (init["S_3"]/total_count)*passage_number
      init["S_4"] = (init["S_4"]/total_count)*passage_number
      init["S_5"] = (init["S_5"]/total_count)*passage_number
      init["S_6"] = (init["S_6"]/total_count)*passage_number
      init["S_7"] = (init["S_7"]/total_count)*passage_number
      init["S_8"] = (init["S_8"]/total_count)*passage_number
      init["G2M_1"] = (init["G2M_1"]/total_count)*passage_number
      init["G2M_2"] = (init["G2M_2"]/total_count)*passage_number
      init["G2M_3"] = (init["G2M_3"]/total_count)*passage_number
      init["G2M_4"] = (init["G2M_4"]/total_count)*passage_number
    }
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_20_G2_ful,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1_1", "G0G1_2", "G0G1_3", "G0G1_4", "G0G1_5", "G0G1_6", "G0G1_7", "G0G1_8", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "S_7", "S_8", "G2M_1", "G2M_2", "G2M_3", "G2M_4", "PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}


ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_10day_no_passaging_validation_wrapper <- function(params, init, tmax, dose_times, drug_amount){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    init = out[out[,"time"] == t_end, c("G0G1","S","G2M","PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_bmin_palboOnly_validation_wrapper <- function(params, init, tmax, dose_times, drug_amount){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0,1), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    init["G0G1"] = drug_amount[i,3] * init["G0G1"]
    init["S"] = drug_amount[i,3] * init["S"]
    init["G2M"] = drug_amount[i,3] * init["G2M"]
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_bmin_palboOnly,
              y = init_num,
              times = seq(t_start, t_end, by = 0.05))
    
    d
    
    init = out[out[,"time"] == t_end, c("G0G1","S","G2M","PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_validation_wrapper_no_passaging <- function(params, init, tmax, dose_times, drug_amount){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0,0,1), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), AA = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["AA"] = drug_amount[i,3]
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    init["G0G1"] =  init["G0G1"]
    init["S"] = init["S"]
    init["G2M"] =  init["G2M"]
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro,
              y = init_num,
              times = seq(t_start, t_end, by = 0.1))
    init = out[out[,"time"] == t_end, c("G0G1","S","G2M","AA","PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_carrying_capacity_differentIC50_validation_wrapper <- function(params, init, tmax, dose_times, drug_amount){
  if (dose_times[1] > 0 ){
    dose_times = c(0, dose_times)
    drug_amount = rbind(c(0,0,0,1), drug_amount)
  }
  
  results = data.frame(time = numeric(), G0G1 = numeric(), S = numeric(),
                       G2M = numeric(), AA = numeric(), FF = numeric(),
                       PP = numeric())
  
  drug_idx = 1:length(dose_times)
  for ( i in drug_idx){
    t_start = dose_times[i]
    if (i < length(drug_idx)){
      t_end = dose_times[i+1]
    } else{
      t_end = tmax
    }
    
    init["AA"] = drug_amount[i,3]
    init["PP"] = drug_amount[i,2]
    init["FF"] = drug_amount[i,1]
    init["G0G1"] = drug_amount[i,4] * init["G0G1"]
    init["S"] = drug_amount[i,4] * init["S"]
    init["G2M"] = drug_amount[i,4] * init["G2M"]
    
    init_num = as.numeric(init)
    names(init_num) = names(init)
    
    out = ode(parms = params,
              func = ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_differentIC50,
              y = init_num,
              times = seq(t_start, t_end, by = 0.1))
    init = out[out[,"time"] == t_end, c("G0G1","S","G2M","AA","PP","FF")]
    outdf = data.frame(out)
    results = rbind(results, outdf)
  }
  return(results)
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_differentIC50 <- function(t, y, params){
  alpha = params$alpha
  b_A = params$b_A
  b_P = params$b_P
  b_F = params$b_F
  c_A = params$c_A
  c_P = params$c_P
  c_F = params$c_F
  c_P_eff = params$c_P_eff
  c_A_eff = params$c_A_eff
  a_FP = params$a_FP
  a_FA = params$a_FA
  K = params$K*1e5
  Gamma = params$Gamma
  delta = params$delta
  
  b = params$b
  
  logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + a_FP * ((y["PP"]/c_P_eff)/(1 + (y["PP"]/c_P_eff))))^(-1) * (1 + a_FA * ((y["AA"]/c_A_eff)/(1 + (y["AA"]/c_A_eff))))^(-1);
  beta = b*(1/(1 + (y["AA"]/c_A)^b_A))*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  
  dG0G1 = 2*alpha*y["G2M"]*logistic_term - beta*y["G0G1"]*logistic_term - delta*y["G0G1"]; # G0G1 compartment
  dS= beta*y["G0G1"]*logistic_term - Gamma*y["S"]*logistic_term - delta*y["S"]; # S compartment
  dG2M = Gamma*y["S"]*logistic_term - alpha*y["G2M"]*logistic_term - delta*y["G2M"];   # G2M compartment
  
  
  #print(c(dG0G1, dS, dG2M))
  dA = 0 ## abema
  dP = 0 ## palbo 
  dF = 0 ## fulv
  
  return(list(c(dG0G1,dS,dG2M, dA, dP, dF)))
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity <- function(t, y, params){
  alpha = params$alpha
  b_A = params$b_A
  b_P = params$b_P
  b_F = params$b_F
  c_A = params$c_A
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  a_FA = params$a_FA
  K = params$K*1e5
  Gamma = params$Gamma
  delta = params$delta
  
  b = params$b
  
  logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) * (1 + a_FA * ((y["AA"]/c_A)/(1 + (y["AA"]/c_A))))^(-1);
  beta = b*(1/(1 + (y["AA"]/c_A)^b_A))*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  
  dG0G1 = 2*alpha*y["G2M"]*logistic_term - beta*y["G0G1"]*logistic_term - delta*y["G0G1"]; # G0G1 compartment
  dS= beta*y["G0G1"]*logistic_term - Gamma*y["S"]*logistic_term - delta*y["S"]; # S compartment
  dG2M = Gamma*y["S"]*logistic_term - alpha*y["G2M"]*logistic_term - delta*y["G2M"];   # G2M compartment
  
  
  #print(c(dG0G1, dS, dG2M))
  dA = 0 ## abema
  dP = 0 ## palbo 
  dF = 0 ## fulv
  
  return(list(c(dG0G1,dS,dG2M, dA, dP, dF)))
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_mod <- function(t, y, params){
  alpha = 1.2*params$alpha #try increasing the rate
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/(5*c_F))^b_F)); #try less fulvestrant effect
  #beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2*10))^b_F_2))
  
  dG0G1_1 = 2*alpha*y["G2M_5"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - beta*y["G0G1_5"]; 
  dS_1= beta*y["G0G1_5"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dG2M_1 = alpha*y["S_5"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];
  dG2M_3 = alpha*y["G2M_2"] - alpha*y["G2M_3"];
  dG2M_4 = alpha*y["G2M_3"] - alpha*y["G2M_4"];
  dG2M_5 = alpha*y["G2M_4"] - alpha*y["G2M_5"];
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dS_1, dG2M_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dS_2, dS_3, dS_4, dS_5, dG2M_2, dG2M_3, dG2M_4, dG2M_5, dP, dF)))
}


ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_multiple <- function(t, y, params){
  alpha = params$alpha #try increasing the rate
  b = params$Beta
  c = params$Gamma
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/(c_F))^b_F)); #try less fulvestrant effect
  #beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2*10))^b_F_2))
  
  dG0G1_1 = 2*c*y["G2M_5"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - beta*y["G0G1_5"]; 
  dS_1= beta*y["G0G1_5"] - b*y["S_1"]; # S compartment
  dS_2= b*y["S_1"] - b*y["S_2"];
  dS_3= b*y["S_2"] - b*y["S_3"];
  dS_4= b*y["S_3"] - b*y["S_4"];
  dS_5= b*y["S_4"] - b*y["S_5"];
  dG2M_1 = b*y["S_5"] - c*y["G2M_1"];   # G2M compartment
  dG2M_2 = c*y["G2M_1"] - c*y["G2M_2"];
  dG2M_3 = c*y["G2M_2"] - c*y["G2M_3"];
  dG2M_4 = c*y["G2M_3"] - c*y["G2M_4"];
  dG2M_5 = c*y["G2M_4"] - c*y["G2M_5"];
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dS_1, dG2M_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dS_2, dS_3, dS_4, dS_5, dG2M_2, dG2M_3, dG2M_4, dG2M_5, dP, dF)))
}


ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_IC <- function(t, y, params){
  alpha = params$alpha #try increasing the rate
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/(c_F))^b_F)); #try less fulvestrant effect
  #beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2*10))^b_F_2))
  
  dG0G1_1 = 2*alpha*y["G2M_5"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - beta*y["G0G1_5"]; 
  dS_1= beta*y["G0G1_5"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dG2M_1 = alpha*y["S_5"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];
  dG2M_3 = alpha*y["G2M_2"] - alpha*y["G2M_3"];
  dG2M_4 = alpha*y["G2M_3"] - alpha*y["G2M_4"];
  dG2M_5 = alpha*y["G2M_4"] - alpha*y["G2M_5"];
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dS_1, dS_2, dS_3, dS_4, dS_5, dG2M_1, dG2M_2, dG2M_3, dG2M_4, dG2M_5, dP, dF)))
}


ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_15_even <- function(t, y, params){
  alpha = params$alpha #try increasing the rate
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/(c_F))^b_F)); #try less fulvestrant effect
  #beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2*10))^b_F_2))
  
  dG0G1_1 = 2*alpha*y["G2M_5"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - beta*y["G0G1_5"]; 
  dS_1= beta*y["G0G1_5"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dG2M_1 = alpha*y["S_5"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];
  dG2M_3 = alpha*y["G2M_2"] - alpha*y["G2M_3"];
  dG2M_4 = alpha*y["G2M_3"] - alpha*y["G2M_4"];
  dG2M_5 = alpha*y["G2M_4"] - alpha*y["G2M_5"];
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dS_1, dG2M_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dS_2, dS_3, dS_4, dS_5, dG2M_2, dG2M_3, dG2M_4, dG2M_5, dP, dF)))
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_15_even_ND <- function(t, y, params){
  alpha = params$alpha #try increasing the rate
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = 0.1*alpha_max + (alpha - 0.1*alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/(c_F))^(0.1*b_F))); #try less fulvestrant effect
  #beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2*10))^b_F_2))
  
  dG0G1_1 = 2*alpha*y["G2M_5"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - beta*y["G0G1_5"]; 
  dS_1= beta*y["G0G1_5"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dG2M_1 = alpha*y["S_5"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];
  dG2M_3 = alpha*y["G2M_2"] - alpha*y["G2M_3"];
  dG2M_4 = alpha*y["G2M_3"] - alpha*y["G2M_4"];
  dG2M_5 = alpha*y["G2M_4"] - alpha*y["G2M_5"];
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dS_1, dG2M_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dS_2, dS_3, dS_4, dS_5, dG2M_2, dG2M_3, dG2M_4, dG2M_5, dP, dF)))
}


ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_15_even_ful_G2 <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  b_F_2 = params$b_F_2
  c_F_2 = params$c_F_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2))^b_F_2))
  
  dG0G1_1 = 2*beta_2*y["G2M_5"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - beta*y["G0G1_5"]; 
  dS_1= beta*y["G0G1_5"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dG2M_1 = alpha*y["S_5"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];
  dG2M_3 = alpha*y["G2M_2"] - alpha*y["G2M_3"];
  dG2M_4 = alpha*y["G2M_3"] - alpha*y["G2M_4"];
  dG2M_5 = alpha*y["G2M_4"] - beta_2*y["G2M_5"];
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dS_1, dG2M_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dS_2, dS_3, dS_4, dS_5, dG2M_2, dG2M_3, dG2M_4, dG2M_5, dP, dF)))
}


ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_15_even_ful_S <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  b_F_2 = params$b_F_2
  c_F_2 = params$c_F_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2))^b_F_2))
  
  dG0G1_1 = 2*alpha*y["G2M_5"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - beta*y["G0G1_5"]; 
  dS_1= beta*y["G0G1_5"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - beta_2*y["S_5"];
  dG2M_1 = beta_2*y["S_5"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];
  dG2M_3 = alpha*y["G2M_2"] - alpha*y["G2M_3"];
  dG2M_4 = alpha*y["G2M_3"] - alpha*y["G2M_4"];
  dG2M_5 = alpha*y["G2M_4"] - alpha*y["G2M_5"];
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dS_1, dG2M_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dS_2, dS_3, dS_4, dS_5, dG2M_2, dG2M_3, dG2M_4, dG2M_5, dP, dF)))
}


ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2 <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  
  dG0G1_1 = 2*alpha*y["G2M_2"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - beta*y["G0G1_6"]; 

  dS_1= beta*y["G0G1_6"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];

  dG2M_1 = alpha*y["S_6"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];

  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dG2M_1, dG2M_2,  dP, dF)))
}


ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2_palbo <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  b_P_2 = params$b_P_2
  c_P_2 = params$c_P_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["PP"] /(c_P_2))^b_P_2))
  
  dG0G1_1 = 2*beta_2*y["G2M_2"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - beta*y["G0G1_6"]; 

  dS_1= beta*y["G0G1_6"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];

  dG2M_1 = alpha*y["S_6"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - beta_2*y["G2M_2"];

  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dG2M_1, dG2M_2,  dP, dF)))
  
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2_ful <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  b_F_2 = params$b_F_2
  c_F_2 = params$c_F_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2))^b_F_2))
  
  dG0G1_1 = 2*beta_2*y["G2M_2"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - beta*y["G0G1_6"]; 
  
  dS_1= beta*y["G0G1_6"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];
  
  dG2M_1 = alpha*y["S_6"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - beta_2*y["G2M_2"];
  
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dG2M_1, dG2M_2,  dP, dF)))
  
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_16_G2_ND <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  

  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = 0.1*alpha_max + (alpha - 0.1*alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^(0.1*b_F)));

  
  dG0G1_1 = 2*alpha*y["G2M_2"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - alpha*y["G0G1_6"]; 
  dG0G1_7 = alpha*y["G0G1_6"] - beta*y["G0G1_7"];
  
  
  dS_1= beta*y["G0G1_7"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];
  dS_7= alpha*y["S_6"] - alpha*y["S_7"];
  
  dG2M_1 = alpha*y["S_7"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];
  
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dG0G1_7, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dS_7, dG2M_1, dG2M_2,  dP, dF)))
  
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_16_G2_ful <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  b_F_2 = params$b_F_2
  c_F_2 = params$c_F_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2))^b_F_2))
  
  dG0G1_1 = 2*beta_2*y["G2M_2"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - alpha*y["G0G1_6"]; 
  dG0G1_7 = alpha*y["G0G1_6"] - beta*y["G0G1_7"];
  
  
  dS_1= beta*y["G0G1_7"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];
  dS_7= alpha*y["S_6"] - alpha*y["S_7"];
  
  dG2M_1 = alpha*y["S_7"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - beta_2*y["G2M_2"];
  
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dG0G1_7, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dS_7, dG2M_1, dG2M_2,  dP, dF)))
  
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2_both_G2_nonsyn <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F

  alpha_max = params$alpha_max
  
  b_P_2 = params$b_P_2
  b_F_2 = params$b_F_2
  c_P_2 = params$c_P_2
  c_F_2 = params$c_F_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  #drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (y["FF"]/c_F)^b_F));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["PP"]/c_P_2)^b_P_2))*(1/(1 + (y["FF"]/c_F_2)^b_F_2));
  
  dG0G1_1 = 2*beta_2*y["G2M_2"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - beta*y["G0G1_6"]; 
  
  dS_1= beta*y["G0G1_6"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];
  
  dG2M_1 = alpha*y["S_6"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - beta_2*y["G2M_2"];
  
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dG2M_1, dG2M_2,  dP, dF)))
  
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_16_G2_both_G2_nonsyn <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  
  alpha_max = params$alpha_max
  
  b_P_2 = params$b_P_2
  b_F_2 = params$b_F_2
  c_P_2 = params$c_P_2
  c_F_2 = params$c_F_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  #drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (y["FF"]/c_F)^b_F));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["PP"]/c_P_2)^b_P_2))*(1/(1 + (y["FF"]/c_F_2)^b_F_2));
  
  dG0G1_1 = 2*beta_2*y["G2M_2"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - alpha*y["G0G1_6"]; 
  dG0G1_7 = alpha*y["G0G1_6"] - beta*y["G0G1_7"]; 
  
  dS_1= beta*y["G0G1_7"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];
  dS_7= alpha*y["S_6"] - alpha*y["S_7"];
  
  dG2M_1 = alpha*y["S_7"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - beta_2*y["G2M_2"];
  
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dG0G1_7, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dS_7, dG2M_1, dG2M_2,  dP, dF)))
  
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_16_G2_both_G2_nonsyn_ND <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  
  alpha_max = params$alpha_max
  
  b_P_2 = params$b_P_2
  b_F_2 = params$b_F_2
  c_P_2 = params$c_P_2
  c_F_2 = params$c_F_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  #drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = 0.1*alpha_max + (alpha - 0.1*alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (y["FF"]/c_F)^(0.1*b_F)));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["PP"]/c_P_2)^b_P_2))*(1/(1 + (y["FF"]/c_F_2)^b_F_2));
  
  dG0G1_1 = 2*beta_2*y["G2M_2"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - alpha*y["G0G1_6"]; 
  dG0G1_7 = alpha*y["G0G1_6"] - beta*y["G0G1_7"]; 
  
  dS_1= beta*y["G0G1_7"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];
  dS_7= alpha*y["S_6"] - alpha*y["S_7"];
  
  dG2M_1 = alpha*y["S_7"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - beta_2*y["G2M_2"];
  
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dG0G1_7, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dS_7, dG2M_1, dG2M_2,  dP, dF)))
  
}


ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_14_G2_both_G2_syn <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP =  params$a_FP
  a_FP_2 =  params$a_FP_2
  alpha_max = params$alpha_max
  
  b_P_2 = params$b_P_2
  b_F_2 = params$b_F_2
  c_P_2 = params$c_P_2
  c_F_2 = params$c_F_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  drugF_eff_2 = y["FF"] * (1 + 0.1*a_FP_2 * ((y["PP"]/c_P_2)/(1 + (y["PP"]/c_P_2))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + ( drugF_eff /c_F)^b_F));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["PP"]/c_P_2)^b_P_2))*(1/(1 + (drugF_eff_2/c_F_2)^b_F_2));
  
  dG0G1_1 = 2*beta_2*y["G2M_2"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - beta*y["G0G1_6"]; 
  
  dS_1= beta*y["G0G1_6"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];
  
  dG2M_1 = alpha*y["S_6"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - beta_2*y["G2M_2"];
  
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dG2M_1, dG2M_2,  dP, dF)))
  
}


ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_20_G2 <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));

  dG0G1_1 = 2*alpha*y["G2M_4"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - alpha*y["G0G1_6"]; 
  dG0G1_7 = alpha*y["G0G1_6"] - alpha*y["G0G1_7"]; 
  dG0G1_8 = alpha*y["G0G1_7"] - beta*y["G0G1_8"]; 
  dS_1= beta*y["G0G1_8"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];
  dS_7= alpha*y["S_6"] - alpha*y["S_7"];
  dS_8= alpha*y["S_7"] - alpha*y["S_8"];
  dG2M_1 = alpha*y["S_8"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];
  dG2M_3 = alpha*y["G2M_2"] - alpha*y["G2M_3"];
  dG2M_4 = alpha*y["G2M_3"] - alpha*y["G2M_4"];
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dG0G1_7, dG0G1_8, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dS_7, dS_8, dG2M_1, dG2M_2, dG2M_3, dG2M_4, dP, dF)))
}



ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_palboOnly_20_G2_ful <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  alpha_max = params$alpha_max
  
  b_F_2 = params$b_F_2
  c_F_2 = params$c_F_2
  alpha_max_2 = params$alpha_max_2
  #K = params$K*1e5
  #Gamma = params$Gamma
  #delta = params$delta
  #b = params$b
  
  #logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + 0.1*a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) ;
  beta = alpha_max + (alpha - alpha_max)*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  beta_2 = alpha_max_2 + (alpha - alpha_max_2)*(1/(1 + (y["FF"] /(c_F_2*10))^b_F_2))
  
  dG0G1_1 = 2*beta_2*y["G2M_4"] - alpha*y["G0G1_1"]; # G0G1 compartment
  dG0G1_2 = alpha*y["G0G1_1"] - alpha*y["G0G1_2"]; 
  dG0G1_3 = alpha*y["G0G1_2"] - alpha*y["G0G1_3"]; 
  dG0G1_4 = alpha*y["G0G1_3"] - alpha*y["G0G1_4"]; 
  dG0G1_5 = alpha*y["G0G1_4"] - alpha*y["G0G1_5"]; 
  dG0G1_6 = alpha*y["G0G1_5"] - alpha*y["G0G1_6"]; 
  dG0G1_7 = alpha*y["G0G1_6"] - alpha*y["G0G1_7"]; 
  dG0G1_8 = alpha*y["G0G1_7"] - beta*y["G0G1_8"]; 
  dS_1= beta*y["G0G1_8"] - alpha*y["S_1"]; # S compartment
  dS_2= alpha*y["S_1"] - alpha*y["S_2"];
  dS_3= alpha*y["S_2"] - alpha*y["S_3"];
  dS_4= alpha*y["S_3"] - alpha*y["S_4"];
  dS_5= alpha*y["S_4"] - alpha*y["S_5"];
  dS_6= alpha*y["S_5"] - alpha*y["S_6"];
  dS_7= alpha*y["S_6"] - alpha*y["S_7"];
  dS_8= alpha*y["S_7"] - alpha*y["S_8"];
  dG2M_1 = alpha*y["S_8"] - alpha*y["G2M_1"];   # G2M compartment
  dG2M_2 = alpha*y["G2M_1"] - alpha*y["G2M_2"];
  dG2M_3 = alpha*y["G2M_2"] - alpha*y["G2M_3"];
  dG2M_4 = alpha*y["G2M_3"] - beta_2*y["G2M_4"];
  dP = 0
  dF = 0
  
  return(list(c(dG0G1_1, dG0G1_2, dG0G1_3, dG0G1_4, dG0G1_5, dG0G1_6, dG0G1_7, dG0G1_8, dS_1, dS_2, dS_3, dS_4, dS_5, dS_6, dS_7, dS_8, dG2M_1, dG2M_2, dG2M_3, dG2M_4, dP, dF)))
  
}



ode_sigmoid_effectiveDose_g1s_transition_1_invitro_carrying_capacity_bmin_palboOnly <- function(t, y, params){
  alpha = params$alpha
  b_P = params$b_P
  b_F = params$b_F
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  K = params$K*1e5
  Gamma = params$Gamma
  delta = params$delta
  
  b_max = params$b_max
  b_min = params$b_max
  
  logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1)
  beta = b_min*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F))
  
  dG0G1 = 2*alpha*y["G2M"]*logistic_term - beta*y["G0G1"]*logistic_term - delta*y["G0G1"]; # G0G1 compartment
  dS= beta*y["G0G1"]*logistic_term - Gamma*y["S"]*logistic_term - delta*y["S"]; # S compartment
  dG2M = Gamma*y["S"]*logistic_term - alpha*y["G2M"]*logistic_term - delta*y["G2M"];   # G2M compartment
  
  #print(c(dG0G1, dS, dG2M))
  dP = 0 ## palbo 
  dF = 0 ## fulv
  
  return(list(c(dG0G1,dS,dG2M, dP, dF)))
}

ode_sigmoid_effectiveDose_g1s_transition_1_invitro <- function(t, y, params){
  alpha = params$alpha
  b_A = params$b_A
  b_P = params$b_P
  b_F = params$b_F
  c_A = params$c_A
  c_P = params$c_P
  c_F = params$c_F
  a_FP = params$a_FP
  a_FA = params$a_FA
  K = params$K*1e5
  Gamma = params$Gamma
  delta = params$delta
  
  b = params$b
  
  logistic_term = (1 - (y["G0G1"] + y["S"] + y["G2M"])/K)
  
  #print(logistic_term)
  
  drugF_eff = y["FF"] * (1 + a_FP * ((y["PP"]/c_P)/(1 + (y["PP"]/c_P))))^(-1) * (1 + a_FA * ((y["AA"]/c_A)/(1 + (y["AA"]/c_A))))^(-1);
  beta = b*(1/(1 + (y["AA"]/c_A)^b_A))*(1/(1 + (y["PP"]/c_P)^b_P))*(1/(1 + (drugF_eff/c_F)^b_F));
  
  dG0G1 = 2*alpha*y["G2M"] - beta*y["G0G1"] - delta*y["G0G1"]; # G0G1 compartment
  dS= beta*y["G0G1"] - Gamma*y["S"] - delta*y["S"]; # S compartment
  dG2M = Gamma*y["S"] - alpha*y["G2M"] - delta*y["G2M"];   # G2M compartment
  
  
  #print(c(dG0G1, dS, dG2M))
  dA = 0 ## abema
  dP = 0 ## palbo 
  dF = 0 ## fulv
  
  return(list(c(dG0G1,dS,dG2M, dA, dP, dF)))
}
