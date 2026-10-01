################################################################################
# Code for "The effect of warming on the bioenergetics of juvenile common sole (Solea solea) in the Gironde estuary: implications for abundance and life-history traits in an energy-limited nursery ground"
# Author : Tristan Halna du Fretay

# DEB and DEB_ode for unlimited feeding scenario


DEBfunction <- function(dpf, parameters, temperatures, f){
  
  # temperatures=TempC_obs
  # f=food
  # parameters=param_cont
  
  ### --- Provides E0
  E0=E0 # Value of the reserve at the beginning (in the egg)
  
  # Organics composition for oxygen flux calculation
  
  mu_X = 525000 # Chemical potential of food
  mu_V = 500000 # Chemical potential of structure
  mu_E = 550000 # Chemical potential of reserve
  mu_P = 480000 # Chemical potential of faeces
  
  n_X <- c(1, 1.8, 0.5, 0.15)  # Chemicals indices for food
  n_E <- c(1, 1.8, 0.5, 0.15)  # Chemicals indices for reserve
  n_V <- c(1, 1.8, 0.5, 0.15)  # Chemicals indices for structure
  n_P <- c(1, 1.8, 0.5, 0.15)  # Chemicals indices for faeces
  
  n_O <- matrix(cbind(n_X, n_E, n_V, n_P), nrow = 4)
  
  n_M <- matrix(c(2,1,2,0), nrow = 4) # Stoichiometric repartition of oxygen
  
  n_M_inv <- c(-1, -1/4, 1/2, 3/4)
  
  JM_JO <- -1 * n_M_inv %*% n_O
  
  # Calculations for flux of food, reserve, structure and faeces
  
  CHON <- c(12, 1, 16, 14)
  w_0 <- CHON %*% n_O # molar mass of structure
  w_V <- w_0[2]       # C-moles structure per volume
  M_V <- parameters$d_v / w_V
  
  y_EX <- parameters$kap_X * mu_X / mu_E # yield of reserve on food
  y_XE <- 1 / y_EX # yield of food on reserve
  y_VE <- mu_E * M_V / parameters$EG  # yield of structure on reserve
  y_PX <- parameters$kap_P * mu_X / mu_P # yield of faeces on food
  y_PE <- y_PX / y_EX # yield of faeces on reserve
  
  eta0 <- matrix(c(y_XE / mu_E * -1, 0, 1 / mu_E, y_PE / mu_E, 0, 0, -1 / mu_E, 0, 0, y_VE / mu_E, -1 / mu_E, 0), nrow = 4)
  
  parameters$eta0 <- eta0
  parameters$JM_JO <- JM_JO
  
  ### --- Create LEH (vector of state variables, and Lb and Lj)  
  LEH = numeric(8) 
  
  LEH[1] = 0.0001     # L, length
  LEH[2] = E0         # E, reserve
  LEH[3] = 0          # H, maturity
  LEH[4] = 0          # E_R, reproduction
  LEH[5] = 0          # Lb, Will be determined by the ode (when L reaches E_Hb)
  LEH[6] = 0          # Lj, Will be determined by the ode (when L reaches E_Hj)
  LEH[7] = 0          # q, ageing acceleration
  LEH[8] = 0          # h, hazard rate
  
  # Add the food conditions to parameters
  parameters$f = f
  
  # Add the temperature conditions to parameters
  parameters$TempC = temperatures
  
  ### --- Solve the ode :
  LEHovertime_cont = ode(y = LEH,
                         func = debODE_ABJ_Marion,
                         times = dpf,
                         parms = parameters,
                         method = "ode23")
  
  ### --- Create the dataframe, add the column names
  LEHovertime_cont = as.data.frame(LEHovertime_cont) 
  colnames(LEHovertime_cont) = c("dpf", "L", "E", "H", "E_R", "Lb", "Lj", "q", "h")
  ### --- estimated weights 
  
  # Estimated length : 
  # L is corrected by the shape coefficient, which is different before and after metamorphosis
  
  # Shape coefficient
  LEHovertime_cont$del_M=0
  LEHovertime_cont[which(LEHovertime_cont[,"H"]<param_cont$E_Hj),]$del_M = 0.166
  
  indices <- which(LEHovertime_cont[,"H"]>=param_cont$E_Hj)
  
  if (length(indices) > 0){
    LEHovertime_cont[indices, ]$del_M <- 0.171
  }
  
  #LEHovertime_cont[which(LEHovertime_cont[,"H"]>=param_cont$E_Hj),]$del_M = 0.171
  
  # Correction of L
  LEHovertime_cont$estim_L_cont = LEHovertime_cont[,"L"]/LEHovertime_cont[,"del_M"]
  
  # Estimated weights : 
  W_V = LEHovertime_cont[,"L"]^3*d_v   # Dry weight of structure
  W_E = LEHovertime_cont[,"E"] / q_E # Dry Weight of reserve 
  # W_R = Weight of Reproduction compartment,  null here as we consider a juvenile fish
  
  W_dry = W_V + W_E # dry weight
  
  k_D = 100*W_dry / LEHovertime_cont$estim_L_cont^3 # to calculate D (Fonds et al., 1989)
  LEHovertime_cont$D = 40.68*k_D^0.364  #  D = % of dry content
  
  LEHovertime_cont$estim_W_cont = W_dry / (LEHovertime_cont$D/100) # Wet weight
  
  # ---- Calculate tb, tj and tp (=age, days post fertilization. At birth, metamorphosis and puberty)
  # and Lb, Lj and Lp
  
  LEHovertime_cont$tbcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E_Hb) == min(abs(LEHovertime_cont[,"H"] - param_cont$E_Hb))),"dpf"]
  
  LEHovertime_cont$tjcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E_Hj) == min(abs(LEHovertime_cont[,"H"] - param_cont$E_Hj))),"dpf"]
  
  LEHovertime_cont$tpcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E_Hp) == min(abs(LEHovertime_cont[,"H"] - param_cont$E_Hp))),"dpf"][1]
  
  LEHovertime_cont$Lpcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E_Hp) == min(abs(LEHovertime_cont[,"H"] - param_cont$E_Hp))),"estim_L_cont"][1]
  
  LEHovertime_cont$Lbcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E_Hb) == min(abs(LEHovertime_cont[,"H"] - param_cont$E_Hb))),"L"]
  
  LEHovertime_cont$Ljcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E_Hj) == min(abs(LEHovertime_cont[,"H"] - param_cont$E_Hj))),"L"]
  
  
  # And add the conditions of temperature
  
  
  if(length(temperatures) > 1){
    LEHovertime_cont$TempC = temperatures[-dim(temperatures)[1],]$TempC 
  } else {
    LEHovertime_cont$TempC = temperatures
  }
  
  #######---------- CALCULATE THE FOOD ASSIMILATION
  
  # CREATE THE COLLUMN S_M IN ESTIM_RES_CONT :
  LEHovertime_cont$s_M=0
  # before birth, no acceleration
  LEHovertime_cont[which(LEHovertime_cont$H < E_Hb),]$s_M = 1 
  #between birth and metamorphosis, s_M = L/Lb
  LEHovertime_cont[which( LEHovertime_cont$H > E_Hb & LEHovertime_cont$H < E_Hj),]$s_M = 
    LEHovertime_cont[which( LEHovertime_cont$H > E_Hb & LEHovertime_cont$H < E_Hj),]$L / LEHovertime_cont[which( LEHovertime_cont$H > E_Hb & LEHovertime_cont$H < E_Hj),]$Lb
  #After metamorphosis, s_M = Lj/Lb
  LEHovertime_cont[which(LEHovertime_cont$H >= E_Hj),]$s_M = 
    LEHovertime_cont[which(LEHovertime_cont$H >= E_Hj),]$Lj / LEHovertime_cont[which(LEHovertime_cont$H >= E_Hj),]$Lb # after the metamorphosis
  
  
  # CREATE THE COLUMN TC IN ESTIM_RES_CONT
  
  LEHovertime_cont$TC=0 # creation of the column
  LEHovertime_cont$TC <- exp(((parameters$T_A)/(parameters$T_ref))-((parameters$T_A)/(LEHovertime_cont$TempC+273.15)))  # If we only had T_A and T_ref 
  LEHovertime_cont$TC <- LEHovertime_cont$TC/(1 + exp(parameters$T_AL/(LEHovertime_cont$TempC+273.15)-T_AL/T_L) + exp(parameters$T_AH/parameters$T_H -parameters$T_AH/(LEHovertime_cont$TempC+273.15)))
  
  # CALCULATE  p_A, p_D and p_G
  LEHovertime_cont$p_A = 0 # Creation of the column
  LEHovertime_cont$p_D = 0 # Creation of the column
  LEHovertime_cont$p_G = 0 # Creation of the column
  
  if (length(f)>1) {   # If f is variable over time
    LEHovertime_cont$p_A = f[-dim(f)[1],]$f*pAm*LEHovertime_cont$s_M*LEHovertime_cont$TC*(LEHovertime_cont$L)^2
  } else {
    LEHovertime_cont$p_A = f*pAm*LEHovertime_cont$s_M*LEHovertime_cont$TC*(LEHovertime_cont$L)^2
  } 
  
  LEHovertime_cont$pC = (parameters$EG * LEHovertime_cont$s_M * (parameters$v*LEHovertime_cont$TC) * ((LEHovertime_cont$L^3)^(-1/3)) + (parameters$pM*LEHovertime_cont$TC)) / ((parameters$EG / LEHovertime_cont$E) + (parameters$kap / (LEHovertime_cont$L)^3))
  LEHovertime_cont$pM = (parameters$pM*LEHovertime_cont$TC)*(LEHovertime_cont$L^3)  # Flux of energy to somatic maintenance (joules/day)
  LEHovertime_cont$p_G = (parameters$kap*LEHovertime_cont$pC)-LEHovertime_cont$pM   # Flux of energy to the soma (joules/day)
  
  LEHovertime_cont$pJ = (parameters$k_J*LEHovertime_cont$TC)*LEHovertime_cont$H
  LEHovertime_cont$pR = (1 - parameters$kap)*LEHovertime_cont$pC - LEHovertime_cont$pJ # Flux of energy that goes to maturity
  
  LEHovertime_cont$p_D = LEHovertime_cont$pM + LEHovertime_cont$pJ + LEHovertime_cont$pR
  
  #######---------- CALCULATE THE MAINTENANCE
  
  # SOMATIC MAINTENANCE
  
  LEHovertime_cont$som_maint = 0 # create the column
  LEHovertime_cont$som_maint = pM*LEHovertime_cont$TC*LEHovertime_cont$L^3
  
  # MATURITY MAINTENANCE 
  
  LEHovertime_cont$mat_maint = 0 # create the column
  LEHovertime_cont$mat_maint = k_J*LEHovertime_cont$TC*LEHovertime_cont$H
  
  # TOTAL MAINTENANCE
  LEHovertime_cont$total_maint = LEHovertime_cont$mat_maint + LEHovertime_cont$som_maint
  
  #######---------- CALCULATE THE FECUNDITY
  
  LEHovertime_cont$F = param_cont$kap_R * LEHovertime_cont[,'E_R'] / E0
  
  ######----------- MASS FLUX
  
  LEHovertime_cont$JOJx <- LEHovertime_cont$p_A * eta0[1,1] + LEHovertime_cont$p_D * eta0[1,2] + LEHovertime_cont$p_G * eta0[1,3] # molar flux of food (mol/time step)
  LEHovertime_cont$JOJv <- LEHovertime_cont$p_A * eta0[2,1] + LEHovertime_cont$p_D * eta0[2,2] + LEHovertime_cont$p_G * eta0[2,3] # molar flux of reserve (mol/time step)
  LEHovertime_cont$JOJe <- LEHovertime_cont$p_A * eta0[3,1] + LEHovertime_cont$p_D * eta0[3,2] + LEHovertime_cont$p_G * eta0[3,3] # molar flux of structure (mol/time step)
  LEHovertime_cont$JOJp <- LEHovertime_cont$p_A * eta0[4,1] + LEHovertime_cont$p_D * eta0[4,2] + LEHovertime_cont$p_G * eta0[4,3] # molar flux of faeces (mol/time step)
  
  LEHovertime_cont$JMO2 <- LEHovertime_cont$JOJx * JM_JO[1,1] + LEHovertime_cont$JOJv * JM_JO[1,2] + LEHovertime_cont$JOJe * JM_JO[1,3] + LEHovertime_cont$JOJp * JM_JO[1,4] # molar flux of O2 (mol/time step)
  
  cont_pA <- - (eta0[1,1] * JM_JO[1,1] + eta0[2,1] * JM_JO[1,2] + eta0[3,1] * JM_JO[1,3] + eta0[4,1] * JM_JO[1,4]) # This is how much of pA is taken to calculate the molar flux of oxygen
  cont_pD <- - (eta0[1,2] * JM_JO[1,1] + eta0[2,2] * JM_JO[1,2] + eta0[3,2] * JM_JO[1,3] + eta0[4,2] * JM_JO[1,4])
  cont_pG <- - (eta0[1,3] * JM_JO[1,1] + eta0[2,3] * JM_JO[1,2] + eta0[3,3] * JM_JO[1,3] + eta0[4,3] * JM_JO[1,4])
  
  #PV=nRT
  #T=273.15 #K
  #R=0.082058 #L*atm/mol*K
  #n=1 #mole
  #P=1 #atm
  #V=nRT/P=1*0.082058*273.15=22.41414
  #T=293.15
  #V=nRT/P=1*0.082058*293.15/1=24.0553
  T_REF <- 293.15
  P_atm <- 1
  R_const <- 0.082058
  gas_cor <- R_const * T_REF / P_atm * (LEHovertime_cont$TempC + 273.15) / T_REF * 1000 # 1 mole to ml/time at Tb and atmospheric pressure
  LEHovertime_cont$O2ML <- -1 * LEHovertime_cont$JMO2 * gas_cor # mlO2/time, temperature corrected (including SDA)
  
  return(LEHovertime_cont)
}

#############################
###  DEB ODE-ABJ MODEL   ####
#############################                 

# Marion Lefebvre Stage M2 2022
# DEB_abj model (acceleration of metabolism between birth and metamorphosis)


# Session

# This script is the function used in the ode(). 

debODE_ABJ_Marion <- function(t, LEH, parms){  # create the function
  
  with(as.list(parms), {
    
    
    dLEH=NULL     # LEH = List of the variables E,L,H,E_R, Lb and Lj
    
    # 1) DEB state variables :
    L       =   LEH[1]
    E       =   LEH[2]    # Reserve
    H       =   LEH[3]    # Maturity
    E.R     =   LEH[4]    # Reproduction
    Lb      =   LEH[5]    # Length at birth
    Lj      =   LEH[6]    # Length at puberty
    q       =   LEH[7]    # Ageing
    h       =   LEH[8]    # Hazard rate
    
    # 2) f and Temperature through time (interpolation if needed)
    
    # ---- If f varies over time
    
    if (length(f)>1) {   # If f is variable over time
      fdt = approx(c(floor(t),floor(t)+1),
                   c(f[which(f$time==floor(t)), 2], f[f$time==(floor(t)+1), 2]), xout=t)$y  # Interpolate
      fdt = max(fdt, 0)
    } else {
      fdt = f
    }     # else, f is always the same
    
    
    # ---- If TempC (Temperature, celsius) varies over time
    
    if (length(TempC)>1) {
      TempCdt = approx(c(floor(t),floor(t)+1),
                       c(TempC[TempC$dpf==floor(t), 2], TempC[TempC$dpf==(floor(t)+1), 2]), xout=t)$y
      
    } else {TempCdt = TempC}
    
    
    # 3) calculation of the shape coefficient (acceleration between birth and metamorphosis)
    
    if (H < E_Hb){  # before birth, no acceleration
      s_M = 1
    } else { 
      if (H < E_Hj){ # between birth and metamorphosis                              
        s_M = L/ Lb 
      } else{  # After metamorphosis, no acceleration anymore
        s_M = Lj/ Lb
      }
    }
    if ("sM" == FALSE){
      s_M=1
    }
    
    
    
    # 4) Correction of the parameters by the shape coefficient and thermic coefficient
    # --- Correction of p_Am and v by s_M 
    # Because acceleration of assimilation and mobilisation of energy
    
    v = v*s_M # For correction of Pc (mobilisation)
    pAm = pAm*s_M  #For correction of PA (assimilation)
    
    # ---- Correction of pM, kJ, p_Am and v by TC :
    
    # Calculate the Thermic coefficient 
    
    TC <- exp(((T_A)/(T_ref))-((T_A)/(TempCdt+273.15)))  # If we only had T_A and T_ref 
    
    TC <- TC/(1 + exp(T_AL/(TempCdt+273.15)-T_AL/T_L) + exp(T_AH/T_H -T_AH/(TempCdt+273.15)))
    
    # tempcdt = t, corrected to have kelvin degrees
    
    v = v*TC               
    pAm = pAm*TC
    pM = pM*TC               
    k_J = k_J*TC
    
    L_m = (kap * pAm / pM) # Maximum structural length
    
    V   = L^3
    V_m = L_m^3
    E_m = pAm  * v
    
    e  = E/E_m
    g = EG/(kap * E_m)
    
    #print(c(V_m, E_m))
    
    # 5)  Fluxes of assimilation (pA) and mobilization (pC)
    
    # Assimilation of energy, pA
    if (H>E_Hb){          # Assimilation only after birth
      pA=fdt*pAm*L^2
    } else {
      pA=0
    }
    
    # Mobilization rate (Out of reserves), pC  (J/cm^3)
    pC = E/(L^3) * (EG * v / L + pM) / (kap * E / (L^3) + EG)      
    
    # 6) Growth rate
    
    if (kap * pC < pM){ # if no available energy for the maintenance
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + EG * kap_G / kap)     # negative growth rate
    } else { 
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + EG / kap)
      #r2 = v * (e/L - (1 + L_T/L)/Lm) / e + g
    }
    # 7)  Generate dE, dH and dL (Evolution of reserve, length and maturity)
    
    #  Maturity(H) : Only before puberty
    if (H < E_Hp) { 
      dH = (1-kap) * pC * L^3 - k_J * H
    } else { 
      dH = 0
    }
    
    # For reproduction (E_R): Only after puberty
    if (H >= E_Hp) { 
      dE_R = kap_R * ((1-kap) * pC * L^3 - k_J * E_Hp)  # = kap_R * pR
    } else { 
      dE_R = 0
    }
    
    # Change in energy in reserve
    dE = pA - pC * (L^3) # reserve = pA - pC*V
    
    # Change in structural length
    dL = L * r / 3 
    
    # Lb and Lj depend on environment and animal's life history -> extracted dynamically
    
    # Lb increases until birth
    if (H <= E_Hb) { 
      dLb = max(0, dL)
    } else dLb=0
    
    # Lj increases until metamorphosis
    if (H <= E_Hj) { 
      dLj = max(0, dL)
    } else dLj=0
    
    # Ageing and hazard rate
    
    dq <- 0
    
    #dq <- (q * (V/V_m) * s_G + h_a) * e * ((v / L) - r) - r * q
    
    dh <- 0
    
    #dh <- q - r * h
    #print(c(t, q))
    # ---- Return dE, dH, dE_R and dL (states variables), but also Lb and Lj as they are dynamic too
    
    dLEH[1] <- dL
    dLEH[2] <- dE
    dLEH[3] <- dH
    dLEH[4] <- dE_R
    dLEH[5] <- dLb
    dLEH[6] <- dLj
    dLEH[7] <- dq
    dLEH[8] <- dh
    list(dLEH)
  })
}

