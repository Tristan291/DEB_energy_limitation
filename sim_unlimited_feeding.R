################################################################################
# Code for "The effect of warming on the bioenergetics of juvenile common sole (Solea solea) in the Gironde estuary: implications for abundance and life-history traits in an energy-limited nursery ground"
# Author : Tristan Halna du Fretay

# Simulations for scenario of unlimited individual feeding

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

# ------------------------ Observation on density ---------------------------- #
# ---- Data importation

# Densities in the estuary - Sturat campaign

densities <- read.csv2("data/Data_All_Sturat_2024.csv")

densities = densities %>%
  filter(espece_id == 182) %>%                          # We select Solea solea
  mutate(year = as.Date(Date_fmt, format = "%Y")) %>%
  group_by(year) %>%
  summarise(dens = mean(Densite),                       # We keep density and total number of indiv/year
            total = sum(total_number),
            traits = n())

densities$year = substr(densities$year, 1, 4)

# ---- Vizualisation

ggplot(densities) +
  geom_point(aes(x = year, y = dens)) +
  theme_classic2() 

# Age repartition in the estuary

ages <- read.table("data/Lengths_Solea_solea.txt")

ages <- ages %>%
  group_by(year, age) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(total = sum(count),
         proportion = count/total)


# Combination of tables

densities_ages <- merge(ages, densities, by = "year") # Merging tables to have both density and age repartition for each year

# ---- Vizualisation

ggplot(densities_ages) +
  geom_bar(aes(x = year, y = dens * proportion, fill = age),
           stat = "identity", position = "dodge") +
  theme_classic2()


################################################################################
# ----------------------------- Simulations ---------------------------------- #

# ---- Individual growth under unlimited individual feeding scenario

time = 1:939                  # We set the simulation time - 939 days

PA_tot = data.frame()         # We keep reference simulations for future comparisons

for (date_birth in seq(1:3)){  # 3 Birth dates - January, February, March
  
  J <- c("", "jan", "mars")[date_birth] 
  B <- c("Temp_cohort", "birth_jan", "birth_mars")[date_birth]
  
  for(i in 0:6){ # 7 cohorts / birth date
    
    if(i == 5 & B == "birth_mars"){ # We exclude one cohort as we don't have a long enough f 
    } else {
      
      T = get(paste(B, "_201", i, sep = "")) # Temperature for a given cohort
      
      f = get(paste("f_mean", J, sep = ""))  # Mean f value for the birth date
      
      deb_model = DEBfunction(time, parameters = param_cont, T, f)  # Simulation of one cohort
      
      PA_tot = rbind(PA_tot, data.frame(dpf = deb_model$dpf,                         # Time in dpf
                                        PA = deb_model$p_A,                          # Assimilation flux p_A
                                        L = deb_model$L,                             # Structural size L
                                        reserve_density = deb_model$reserve_density, # Reserve density
                                        tpcont = deb_model$tpcont,                   # Age at maturity
                                        cohort = i,
                                        birth = date_birth))  
      
    }
  }
}

# We keep the mean individual under reference period

pa_tot_mean <- PA_tot %>%
  select(-cohort)%>%
  group_by(dpf)%>%
  summarize(PA = mean(PA), L = mean(L), reserve_density = mean(reserve_density))

test_birth_scenario_total = data.frame() # Defining the working data frame for simulations

for (temp in seq(1:3)){    # 3 Temperature conditions
  
  temp_scen <- c("201", "mdl_rcp_2_6_209", "mdl_rcp_8_5_209")[temp]
  for (date_birth in 1:3){     # 3 Birth dates - January, February, March
    
    J <- c("", "jan", "mars")[date_birth]
    B <- c("Temp_cohort", "birth_jan", "birth_mars")[date_birth]
    
    for (i in 0:6){      # 7 Cohorts
      
      if (B == "birth_mars" & i == 5) {
      } else {
        
        T = get(paste(B, "_", temp_scen, i, sep = "")) #  Temperature for this cohort
        f = get(paste("f_mean", J, sep = ""))  # Mean f value (always the mean to allow comparison with future temperature conditions)
        
        deb_model = DEBfunction(time, parameters = param_cont, T, f)   # Simulation
        deb_model = energy_flux(deb_model, param_cont)   # Adding the calculation of all energy fluxes (see Parameters.R)
        
        test_birth_scenario_total <- rbind(test_birth_scenario_total, data.frame(deb_model, birth = date_birth,     # We keep the model and the birth dates
                                                                                 scenario_temp = temp, cohort = i,  # The temperature condition and the cohort
                                                                                 temperature = mean(T$TempC), scenario = 1,  # The mean temperature and it is scenario 1 for unlimited feeding
                                                                                 delta_len = (deb_model$L - pa_tot_mean$L)/pa_tot_mean$L * 100,  # Length anomaly
                                                                                 delta_E = (deb_model$reserve_density - pa_tot_mean$reserve_density)/pa_tot_mean$reserve_density * 100))  # Reserve density anomaly
        
        
      }
    }
  }
}

test_birth_scenario_total <- test_birth_scenario_total %>% mutate(cond_factor = 100 * (estim_W_cont / estim_L_cont^3))   # We calculate Fulton's K condition factor

test_birth_mean <- test_birth_scenario_total%>%   # We aggregate all cohort to have a mean individual for each temperature condition
  group_by(dpf, scenario_temp) %>%
  summarize(length = mean(estim_L_cont),  # Mean physical length
            lower_range = min(estim_L_cont),   # Lower boundary
            upper_range = max(estim_L_cont),   # Upper boundary
            scaled_E = mean(reserve_density),   # Mean reserve density
            lower_scaled_E = min(reserve_density),
            upper_scaled_E = max(reserve_density),
            age_mat = mean(tpcont),   # Mean age at maturity
            low_age_mat = min(tpcont),
            up_age_mat = max(tpcont),
            len_mat = mean(Lpcont),   # Mean length at maturity
            low_len_mat = min(Lpcont),
            up_len_mat = max(Lpcont),
            diff_len = mean(delta_len),   # Mean physical length anomaly
            low_diff_len = min(delta_len),
            up_diff_len = max(delta_len),
            diff_E = mean(delta_E),   # Mean reserve density anomaly
            low_diff_E = min(delta_E),
            up_diff_E = max(delta_E),
            p_R = mean(pR),   # Mean energy fluxes
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
            L = mean(L),   # Mean structural length
            low_F = min(F),
            up_F = max(F),
            F = mean(F),   # Mean fecundity
            low_cond = min(cond_factor),
            up_cond = max(cond_factor),
            cond = mean(cond_factor))

save(test_birth_scenario_total, test_birth_mean, file = "model_scenario1.RData")   # Save it for graphical representation (see graphical_representation.R)

# ------------------------------- Density ------------------------------------ #

################################################################################
# ====================== Trophic Carrying capacity  ========================== #

# Calculation of mean assimilation flux p_A 

simu_dens_actual_var = data.frame() 

for (temp in seq(1:3)){   # 3 Temperature conditions
  
  temp_scen <- c("201", "mdl_rcp_2_6_209", "mdl_rcp_8_5_209")[temp]
  for (date_birth in seq(1:3)){   # 3 Birth dates
    
    J <- c("", "jan", "mars")[date_birth]
    B <- c("Temp_cohort", "birth_jan", "birth_mars")[date_birth]
    
    for (i in 0:6){   # 7 Cohorts
      
      if (B == "birth_mars" & i == 5) {   # Cannort integrate one cohort due to lack of data
      } else {
        
        T = get(paste(B, "_", temp_scen, i, sep = ""))   #  Temperature for this cohort
        f = get(paste("f_mean", J, sep = ""))   # Mean f value 
    
        
        deb_model = DEBfunction(time, parameters = param_cont, T, f)   # Simulation 
        deb_model = energy_flux(deb_model, param_cont)
        
        
        mean_pa <- c(mean(deb_model$p_A[which(deb_model$dpf <= 365)]),    # Mean p_A of year-0
                     mean(deb_model$p_A[which(deb_model$dpf > 365 & deb_model$dpf <= 365*2)]),  # Mean p_A of year-1
                     mean(deb_model$p_A[which(deb_model$dpf > 365 * 2)]))   # Mean p_A of year-2
        
        mean_w  <- c(mean(deb_model$estim_W_cont[which(deb_model$dpf <= 365)]),   # Mean wet weight of year-0 (for biomass)
                     mean(deb_model$estim_W_cont[which(deb_model$dpf > 365 & deb_model$dpf <= 365*2)]),   # Mean wet weight of year-1
                     mean(deb_model$estim_W_cont[which(deb_model$dpf > 365 * 2)]))   # Mean wet weight of year-2
        
        if(temp == 1){
          year <- c(2010 + i, 2011 + i, 2012 + i)
        } else {
          year <- c(2090 + i, 2091 + i, 2092 + i)
        }
        
        age <- c("G0", "G1", "G2")
        
        simu_dens_actual_var = rbind(simu_dens_actual_var,  
                                     data.frame(mean_pa,   # Assimilation flux
                                                mean_w,    # Wet weigth
                                                year,     
                                                age,
                                                period = c("ref", "rcp2.6", "rcp8.5")[temp],
                                                temperature = mean(T$TempC),   # Mean temperature
                                                birth = date_birth))  
      }
    }
  }
}

# Mean trophic carrying capacity of the estuary on reference period

pa_tot <- merge(densities_ages, simu_dens_actual_var, by = c("year", "age"), all = FALSE)  # Merging density data and assimilation flux
pa_tot <- simu_dens_actual_var %>%
  filter(year != 2010,    # Not enough data for all ages
         year != 2011,    # Not enough data for all ages
         year != 2017,    # Not enough data for all ages
         year != 2018,    # Not enough data for all ages
         year != 2014,    # Abnormal year - 3 times bigger than every other
         year != 2090,    # Not enough data for all ages
         year != 2091,    # Not enough data for all ages
         year != 2097,    # Not enough data for all ages
         year != 2098)    # Not enough data for all ages

mean_prop   <- densities_ages%>%group_by(age)%>%summarize(mean_prop = mean(proportion)) # Mean proportion of each age class in a year
pa_tot <- pa_tot%>%left_join(mean_prop, by = "age")

pa_tot_ref <- pa_tot %>%    # Keeping only the reference period
  full_join(densities_ages, by = c("year", "age")) %>%
  filter(period == "ref") %>%
  mutate(pa_cohort = mean_pa * mean_prop) %>%   # Assimilation per age class for each cohort
  group_by(year, birth)   %>%
  summarise(pa_mean_indiv = mean(pa_cohort),  # Mean assimilation for a given year
            dens_estim = mean(dens),   # Mean density
            period = first(period),
            temp = mean(temperature))  # Mean temperature

mean_pa_tot   <- mean(pa_tot_ref$pa_mean_indiv)   # Mean carrying capcity of the estuary
mean_dens_tot <- mean(pa_tot_ref$dens_estim)      # Mean density on the reference period

# Estimation of future densities

dens_estimated <- pa_tot    %>%
  filter(period != "ref")   %>%  # Only keeping future simulations
  mutate(pa_cohort = mean_pa * mean_prop,   # Assimilation per age class for each cohort
         weight_cohort = mean_w * mean_prop) %>%   # Weight per age class for each cohort
  group_by(year, period, birth)    %>%
  summarize(pa_mean_indiv = mean(pa_cohort),   # Mean assimilation for a given year
            w_mean_indiv = mean(weight_cohort),   # Mean weight for a given year
            temp = mean(temperature)) %>%
  mutate(dens_estim = mean_pa_tot / pa_mean_indiv,   # Estimation of density
         mean_biomass = w_mean_indiv * dens_estim)   # Estimation of biomass

dens_estimated <- rbind(dens_estimated, pa_tot_ref)   # Associating reference and future simulations
dens_estimated$period <- factor(dens_estimated$period, levels = c("ref", "rcp2.6", "rcp8.5"))

A = dens_estimated %>% filter(period == "ref") %>% select(dens_estim)  # Only observed reference density
dens_moy_ref <- mean(A$dens_estim)

dens_estimated$diff_relative <- (dens_estimated$dens_estim - dens_moy_ref) / dens_moy_ref * 100   # Relative difference in density

save(dens_estimated, file = "unlimited_feed_dens.RData")   # Save it for graphical representation (see graphical_representation.R)