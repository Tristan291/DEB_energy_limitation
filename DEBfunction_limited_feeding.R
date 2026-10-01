################################################################################
# Code for "The effect of warming on the bioenergetics of juvenile common sole (Solea solea) in the Gironde estuary: implications for abundance and life-history traits in an energy-limited nursery ground"
# Author : Tristan Halna du Fretay

# DEB and DEB_ode for limited feeding scenario

DEBfunction_PA_CORRECTED_TEMP <- function(dpf, parameters, temperatures, f, ref){
  
  # temperatures=TempC_obs
  # f=food
  # PA = Energy consumed by one individual in natural conditions
  # parameters=param_cont
  
  ### --- Provides E0
  E0=E0 # Value of the reserve at the beginning (in the egg)
  
  ### --- Create LEH (vector of state variables, and Lb and Lj)  
  LEH = numeric(6) 
  
  LEH[1] = 0.0001     # L, length
  LEH[2] = E0         # E, reserve
  LEH[3] = 0          # H, maturity
  LEH[4] = 0          # E_R, reproduction
  LEH[5] = 0          # Lb, Will be determined by the ode (when L reaches E_Hb)
  LEH[6] = 0          # Lj, Will be determined by the ode (when L reaches E_Hj)
  
  # Add the food conditions to parameters
  parameters$f = f
  
  # Add the temperature conditions to parameters
  parameters$TempC = temperatures
  
  parameters$ref = ref
  
  ### --- Solve the ode :
  LEHovertime_cont = ode(y = LEH,
                         func = debODE_ABJ_PA_CORRECTED_TEMP,
                         times = dpf,
                         parms = parameters,
                         method = "ode23")
  
  ### --- Create the dataframe, add the column names
  LEHovertime_cont = as.data.frame(LEHovertime_cont) 
  colnames(LEHovertime_cont) = c("dpf", "L", "E", "H", "E_R", "Lb", "Lj")
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
  
  # cREATE THE COLLUMN S_M IN ESTIM_RES_CONT :
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
  
  # CALCULATE  PA
  
  LEHovertime_cont$p_A_ref = ref$p_A
  
  LEHovertime_cont$p_A = 0 # Creation of the column
  LEHovertime_cont$p_A = pmin(LEHovertime_cont$p_A_ref / ((ref$L**2) * ref$SM * ref$TC) * (LEHovertime_cont$L**2) * LEHovertime_cont$s_M * LEHovertime_cont$TC,
                             LEHovertime_cont$p_A_ref)
  
  
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
  
  LEHovertime_cont$f_estim = LEHovertime_cont$p_A / (param_cont$pAm*LEHovertime_cont$TC*LEHovertime_cont$s_M * LEHovertime_cont$L^2)
  
  return(LEHovertime_cont)
}

# This script is the function used in the ode(). 

debODE_ABJ_PA_CORRECTED_TEMP <- function(t, LEH, parms){  # create the function
  
  with(as.list(parms), {
    
    dLEH=NULL     # LEH = List of the variables E,L,H,E_R, Lb and Lj
    
    # 1) DEB state variables :
    L      =   LEH[1]
    E      =   LEH[2]     # Reserve
    H      =   LEH[3]     # Maturity
    E.R    =   LEH[4]     # Reproduction
    Lb     =   LEH[5]     # Length at birth
    Lj     =   LEH[6]     # Length at puberty
    
    # 2) Temperature through time (interpolation if needed)
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
    
    #fdt = pA/pAm*L^2
    
    # 5)  Fluxes of assimilation (pA) and mobilization (pC)
    
    # Assimilation of energy, pA
    
    if (H<E_Hb) {
      
      pA = 0
      
    } else {
      
      pA = min(ref[which(ref$dpf == floor(t)),]$p_A / ((ref[which(ref$dpf == floor(t)),]$L**2) * ref[which(ref$dpf == floor(t)),]$SM * ref[which(ref$dpf == floor(t)),]$TC) * (L**2) * s_M * TC ,
               ref[which(ref$dpf == floor(t)),]$p_A)
      
    }
    
    #print(fdt)
    # Mobilization rate (Out of reserves), pC  (J/cm^3)
    pC = E/(L^3) * (EG * v / L + pM) / (kap * E / (L^3) + EG)
    
    # 6) Growth rate
    if (kap * pC < pM){ # if no available energy for the maintenance
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + EG * kap_G / kap)     # negative growth rate
    } else { 
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + EG / kap)
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
    
    
    # ---- Return dE, dH, dE_R and dL (states variables), but also Lb and Lj as they are dynamic too
    
    dLEH[1] <- dL
    dLEH[2] <- dE
    dLEH[3] <- dH
    dLEH[4] <- dE_R
    dLEH[5] <- dLb
    dLEH[6] <- dLj
    list(dLEH)
  })
}

