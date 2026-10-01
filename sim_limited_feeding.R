################################################################################
# Code for "The effect of warming on the bioenergetics of juvenile common sole (Solea solea) in the Gironde estuary: implications for abundance and life-history traits in an energy-limited nursery ground"
# Author : Tristan Halna du Fretay

# Simulations for scenario of limited individual feeding

library(dplyr)
library(lubridate)
library(data.table)
library(ggplot2)
library(viridis)
library(deSolve)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(abind)
library(progress)
library(tidyverse)
library(lme4)
library(zoo)
library(patchwork)

################################################################################
# =========================== Data importation =============================== #

# f value

estim_f = read.table("data/estim_F_recalc.txt", dec = ",") %>% dplyr::select(dpf, f, cohort)

for (i in 0:6){                                          # 7 cohorts : 2010 - 2016
  nom_data_frame = paste("estim_f", 2010 + i, sep = "")  # Creating a data frame for a given cohort
  assign(nom_data_frame, estim_f %>%                     
           filter(cohort == 2010 + i) %>%                # Selecting the cohort of interest
           slice(1:940) %>%                              # Limiting to the juvenile phase (smallest amount of time a cohort has spent in the nursery)
           rename(time = dpf) %>%                        # Formatting data frame to the model
           dplyr::select(-cohort) %>%
           mutate(f = as.numeric(f))
  )
}

f_s <- list(estim_f2010, estim_f2011, estim_f2012, estim_f2013, estim_f2014, estim_f2015, estim_f2016)
f_mean <- Reduce("+", f_s)/ length(f_s) # Mean f value on the reference period

ggplot(f_mean) +
  theme_classic2() +
  geom_line(aes(x = time, y = f)) 

# Temperature data

temperature = read.table("data/Temperatures.txt")                   # Importation of temperature from MARS3D
temperature$date = as.Date(temperature$date, format = "%Y-%m-%d")

for (i in 0:6){  # 7 cohorts : 2010 - 2016
  nom_data_frame = paste("Temp_cohort_", 2010 + i, sep = "") 
  assign(nom_data_frame, temperature %>%                     
           filter(date >= paste(2010 + i,"-02-01", sep = "")) %>%   # Birth date on February 1st    
           slice(1:940) %>%    
           mutate(dpf = 1:940) %>%  
           dplyr::select(mean_daily_temp, dpf) %>%                  # Formatting the data for the model
           rename( TempC = mean_daily_temp, dpf = dpf)%>%
           'rownames<-'(1:940) %>%
           relocate(dpf)
  )
}


# Loading prediction temperature

load("data/predicted_rcp85_lag.RData")
load("data/predicted_rcp26_lag.RData")

# ------ Temperature under RCP 2.6 (exact same protocol)

for (i in 0:6){  
  nom_data_frame = paste("Temp_cohort_mdl_", "rcp_2_6_", 2090 + i, sep = "") 
  assign(nom_data_frame, predicted_rcp_26_with_lag %>%                     
           filter(date >= paste(2090 + i,"-02-01", sep = "")) %>%            # Birth date on February 1st 
           slice(1:940) %>%    
           mutate(dpf = 1:940) %>%  
           mutate(temperature = temperature) %>%
           dplyr::select(dpf, temperature) %>%  
           rename(dpf = dpf, TempC = temperature)
  )
}

# ------ Temperature under RCP 8.5

for (i in 0:6){ 
  nom_data_frame = paste("Temp_cohort_mdl_", "rcp_8_5_", 2090 + i, sep = "") 
  assign(nom_data_frame, predicted_rcp_85_with_lag %>%                    
           filter(date >= paste(2090 + i,"-02-01", sep = "")) %>%            # Birth date on February 1st 
           slice(1:940) %>%    
           mutate(dpf = 1:940) %>%   
           mutate(temperature = temperature) %>%
           dplyr::select(dpf, temperature) %>%  
           rename(dpf = dpf, TempC = temperature)
  )
}

# ---- Extraction of temperature for a birth in January
# ---- Reference temperature

for (i in 0:6){  
  nom_data_frame = paste("birth_jan_", 2010 + i, sep = "") 
  assign(nom_data_frame, temperature %>%                    
           filter(date >= paste(2010 + i,"-01-01", sep = "")) %>%    # Birth date on January 1st 
           slice(1:940) %>%    
           mutate(dpf = 1:940) %>%   
           dplyr::select(mean_daily_temp, dpf) %>%  
           rename( TempC = mean_daily_temp, dpf = dpf)%>%
           'rownames<-'(1:940) %>%
           relocate(dpf)
  )
}

# ---- RCP2.6 temperature born in January

for (i in 0:6){  
  nom_data_frame = paste("birth_jan_mdl_rcp_2_6_", 2090 + i, sep = "") 
  assign(nom_data_frame, predicted_rcp_26_with_lag %>%                    
           filter(date >= paste(2090 + i,"-01-01", sep = "")) %>%    # Birth date on January 1st 
           slice(1:940) %>%   
           mutate(dpf = 1:940) %>%   
           dplyr::select(temperature, dpf) %>%  
           rename(TempC = temperature, dpf = dpf)%>%
           'rownames<-'(1:940) %>%
           relocate(dpf)
  )
}

# ---- RCP8.5 temperature born in January

for (i in 0:6){  
  nom_data_frame = paste("birth_jan_mdl_rcp_8_5_", 2090 + i, sep = "") 
  assign(nom_data_frame, predicted_rcp_85_with_lag %>%                     
           filter(date >= paste(2090 + i,"-01-01", sep = "")) %>%    # Birth date on January 1st 
           slice(1:940) %>%    
           mutate(dpf = 1:940) %>%   
           dplyr::select(temperature, dpf) %>% 
           rename(TempC = temperature, dpf = dpf)%>%
           'rownames<-'(1:940) %>%
           relocate(dpf)
  )
}

# ---- Extraction of temperature for a birth in March
# ---- Reference temperature

for (i in 0:6){  #7 tours de boucles pour 7 cohort : 2010 - 2016
  nom_data_frame = paste("birth_mars_", 2010 + i, sep = "") 
  assign(nom_data_frame, temperature %>%                    
           filter(date >= paste(2010 + i,"-03-01", sep = "")) %>%    # Birth date on March 1st 
           slice(1:940) %>%    
           mutate(dpf = 1:940) %>%   
           dplyr::select(mean_daily_temp, dpf) %>%  
           rename( TempC = mean_daily_temp, dpf = dpf)%>%
           'rownames<-'(1:940) %>%
           relocate(dpf)
  )
}

# ---- RCP2.6 temperature born in Mars

for (i in 0:6){  
  nom_data_frame = paste("birth_mars_mdl_rcp_2_6_", 2090 + i, sep = "") 
  assign(nom_data_frame, predicted_rcp_26_with_lag %>%                    
           filter(date >= paste(2090 + i,"-03-01", sep = "")) %>%    # Birth date on March 1st
           slice(1:940) %>%    
           mutate(dpf = 1:940) %>%  
           dplyr::select(temperature, dpf) %>%  
           rename(TempC = temperature, dpf = dpf)%>%
           'rownames<-'(1:940) %>%
           relocate(dpf)
  )
}

# ---- RCP8.5 temperature born in Mars

for (i in 0:6){ 
  nom_data_frame = paste("birth_mars_mdl_rcp_8_5_", 2090 + i, sep = "") 
  assign(nom_data_frame, predicted_rcp_85_with_lag %>%                     
           filter(date >= paste(2090 + i,"-03-01", sep = "")) %>%    # Birth date on March 1st
           slice(1:940) %>%    
           mutate(dpf = 1:940) %>%   
           dplyr::select(temperature, dpf) %>%  
           rename(TempC = temperature, dpf = dpf)%>%
           'rownames<-'(1:940) %>%
           relocate(dpf)
  )
}

# ---- Vizualisation

ggplot() + theme_classic2() +
  geom_line(data = birth_jan_2010, mapping = aes(x = dpf, y = TempC), color = "blue") +
  geom_line(data = Temp_cohort_2010, mapping = aes(x = dpf, y = TempC), color = "green") +
  geom_line(data = birth_mars_2010, mapping = aes(x = dpf, y = TempC), color = "red") +
  geom_line(data = birth_jan_mdl_rcp_2_6_2090, mapping = aes(x = dpf, y = TempC), color = "lightblue2") +
  geom_line(data = Temp_cohort_mdl_rcp_2_6_2090, mapping = aes(x = dpf, y = TempC), color = "green2") +
  geom_line(data = birth_mars_mdl_rcp_2_6_2090, mapping = aes(x = dpf, y = TempC), color = "red2") +
  geom_line(data = birth_jan_mdl_rcp_8_5_2090, mapping = aes(x = dpf, y = TempC), color = "lightblue4") +
  geom_line(data = Temp_cohort_mdl_rcp_8_5_2090, mapping = aes(x = dpf, y = TempC), color = "green4") +
  geom_line(data = birth_mars_mdl_rcp_8_5_2090, mapping = aes(x = dpf, y = TempC), color = "red4") 

# ---- Creating f datasets for January and March

# Birth in january

estim_fjan = read.table("data/estim_optim_females_solea_January.txt", dec = ",") %>% dplyr::select(January.dpf, January.f, January.cohort) # Importation of back calculated f data

for (i in 0:6){ 
  nom_data_frame = paste("estim_fjan", 2010 + i, sep = "")  
  assign(nom_data_frame, estim_fjan %>%                     
           filter(January.cohort == 2010 + i) %>%  
           slice(1:940) %>%  
           rename(time = January.dpf, f = January.f) %>%  # Formatting the data for the model
           dplyr::select(-January.cohort) %>%
           mutate(f = as.numeric(f))
  )
}

f_s <- list(estim_fjan2010, estim_fjan2011, estim_fjan2012, estim_fjan2013, estim_fjan2014, estim_fjan2015, estim_fjan2016)
f_meanjan <- Reduce("+", f_s)/ length(f_s) # Mean f value for the reference period and a birth in January

# Birth in March

estim_fmars = read.table("data/estim_optim_females_solea_Mars.txt", dec = ",") %>% dplyr::select(Mars.dpf, Mars.f, Mars.cohort) # Importation of back calculated f data

for (i in 0:6){  
  nom_data_frame = paste("estim_fmars", 2010 + i, sep = "") 
  assign(nom_data_frame, estim_fmars %>%
           filter(Mars.cohort == 2010 + i) %>%
           slice(1:940) %>%  
           rename(time = Mars.dpf, f = Mars.f) %>%  # Formatting the data for the model
           dplyr::select(-Mars.cohort) %>%
           mutate(f = as.numeric(f))
  )
}

f_s <- list(estim_fmars2010, estim_fmars2011, estim_fmars2012, estim_fmars2013, estim_fmars2014, estim_fmars2016)
f_meanmars <- Reduce("+", f_s)/ length(f_s) # Mean f value for the reference period and a birth in March

# ---- Vizualisation

ggplot(f_mean) +
  theme_classic2() +
  geom_line(f_meanjan, mapping = aes(x = time, y = f, color = "January"), linewidth = 1.3) +
  geom_line(aes(x = time, y = f, color = "February"), linewidth = 1.3) +
  geom_line(f_meanmars, mapping = aes(x = time, y = f, color = "Mars"), linewidth = 1.3) +
  scale_color_manual(values = c("skyblue", "olivedrab", "goldenrod3"), labels = c("1st January", "1st February", "1st Mars")) +
  labs(color = "Birthdate") +
  xlab("Time (in dpf)") +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 1000)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 2)) +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 18))  

################################################################################
# ----------------------------- Simulations ---------------------------------- #

# ---- Individual growth under limited individual feeding scenario

time = 1:939                  # We set the simulation time - 939 days

PA_tot = data.frame()         # We keep reference simulations for future comparisons

for (date_birth in 1:3){   # 3 Birth dates - January, February, March
  
  J <- c("", "jan", "mars")[date_birth]
  B <- c("Temp_cohort", "birth_jan", "birth_mars")[date_birth]
  
  for(i in 0:6){    # 7 cohorts / birth date
    
    if(i == 5 & B == "birth_mars"){    # Excluding one cohort as we don't have a long enough f 
    } else {
      
      T = get(paste(B, "_201", i, sep = "")) # Temperature for a given cohort
      f = get(paste("f_mean", J, sep = ""))  # Mean f value for the birth date
      
      deb_model = DEBfunction(time, parameters = param_cont, T, f)  # Simulation of one cohort
      
      PA_tot = rbind(PA_tot, data.frame(dpf = deb_model$dpf,    
                                        p_A = deb_model$p_A,   # Assimilation flux, used for energy limitation
                                        L = deb_model$L,       # Structural length, used for energy limitation
                                        reserve_density = deb_model$reserve_density,   # Reserve density
                                        TC = deb_model$TC,    # Temperature correction, used for energy limitation
                                        SM = deb_model$s_M,   # s_M acceleration factor, used for energy limitation
                                        cohort = i,
                                        birth = date_birth))  
      
    }
  }
}

pa_tot_mean <- PA_tot  %>%  # Keeping all that is necessary for calculation of energy limitation
  select(-cohort)      %>%
  group_by(dpf, birth) %>%
  summarize(p_A = mean(p_A),
            L = mean(L),
            TC = mean(TC),
            SM = mean(SM))

len_tot_mean <- PA_tot %>%   # Keeping length and reserve density for calculation of anomalies
  select(-cohort)    %>%
  group_by(dpf)      %>%
  summarize(L = mean(L),
            reserve_density = mean(reserve_density))


test_birth_scenario_total_limited = data.frame()    # Defining the working data frame for simulations

for (temp in 1:3){     # 3 Temperature conditions
  
  temp_scen <- c("201", "mdl_rcp_2_6_209", "mdl_rcp_8_5_209")[temp]
  for (date_birth in 1:3){    # 3 Birth dates - January, February, March
    
    J <- c("", "jan", "mars")[date_birth]
    B <- c("Temp_cohort", "birth_jan", "birth_mars")[date_birth]
    
    for (i in 0:6){    # 7 Cohorts
      
      if (B == "birth_mars" & i == 5) {
      } else {
        
        T = get(paste(B, "_", temp_scen, i, sep = "")) # Temperature for a given cohort
        f = get(paste("f_mean", J, sep = ""))          # Mean f value for the birth date
        
        pa_mdl = pa_tot_mean%>%filter(birth == date_birth)%>%select(dpf, p_A, L, TC, SM)    # Reference simulations for a given birth date - used for energy limitation calculation

        deb_model = DEBfunction_PA_CORRECTED_TEMP(time, parameters = param_cont, T, f, ref = pa_mdl, scen_temp = temp) # Simulation (DEB model integrating energy limitation, see DEBfunction_limited_feeding.R)
        deb_model = energy_flux(deb_model, param_cont) # Adding the calculation of all energy fluxes (see Parameters.R)
        
        test_birth_scenario_total_limited <- rbind(test_birth_scenario_total_limited, data.frame(deb_model, birth = date_birth,      # Keeping the model and the birth dates
                                                                                                 scenario_temp = temp, cohort = i,   # The temperature condition and the cohort
                                                                                                 temperature = mean(T$TempC), scenario = 2,   # The mean temperature and it is scenario 2 for limited feeding
                                                                                                 delta_len = (deb_model$L - len_tot_mean$L)/len_tot_mean$L * 100,   # Length anomaly
                                                                                                 delta_E = (deb_model$reserve_density - len_tot_mean$reserve_density)/len_tot_mean$reserve_density * 100))   # Reserve density anomaly
        
        
      }
    }
  }
}

reserve_ref <- test_birth_scenario_total_limited %>%
  filter(scenario_temp == 1) %>%
  group_by(dpf) %>%
  summarize(reserve_ref = mean(reserve_density), .groups = "drop")

test_birth_scenario_total_limited_2 <- test_birth_scenario_total_limited %>%
  left_join(reserve_ref, by = "dpf") %>%
  group_by(scenario_temp) %>%
  mutate(diff_reserve = (reserve_density - reserve_ref) / reserve_ref * 100) %>%
  ungroup()

test_birth_mean_limited <- test_birth_scenario_total_limited_2%>%    # We aggregate all cohort to have a mean individual for each temperature condition
  group_by(dpf, scenario_temp) %>%
  summarize(length = mean(estim_L_cont),     # Mean physical length
            lower_range = min(estim_L_cont), # Lower boundary
            upper_range = max(estim_L_cont), # Upper boundary
            scaled_E = mean(reserve_density), # Mean reserve density
            lower_scaled_E = min(reserve_density),
            upper_scaled_E = max(reserve_density),
            age_mat = mean(tpcont),    # Mean age at maturity
            low_age_mat = min(tpcont),
            up_age_mat = max(tpcont),
            len_mat = mean(Lpcont),    # Mean length at maturity
            low_len_mat = min(Lpcont),
            up_len_mat = max(Lpcont),
            diff_len = mean(delta_len),   # Mean physical length anomaly
            low_diff_len = min(delta_len),
            up_diff_len = max(delta_len),
            diff_E = mean(diff_reserve),   # Mean reserve density anomaly
            low_diff_E = min(diff_reserve),
            up_diff_E = max(diff_reserve),
            p_R = mean(pR),       # Mean energy fluxes
            low_p_R = min(pR),
            up_p_R = max(pR),
            p_G = mean(pG_tot),
            low_p_G = min(pG_tot),
            up_p_G = max(pG_tot),
            p_C = mean(pC),
            low_p_C = min(pC),
            up_p_C = max(pC),
            p_M = mean(pM),
            p_J = mean(pJ),
            temperature = first(temperature),
            low_p_A = min(p_A),
            up_p_A = max(p_A),
            p_A = mean(p_A),
            L = mean(L),    # Mean structural length
            low_F = min(F),
            up_F = max(F),
            F = mean(F))  # Mean fecundity

save(test_birth_scenario_total_limited, pa_tot_mean, test_birth_mean_limited , file = "model_scenario2.RData")    # Save it for graphical representation (see graphical_representation.R)
