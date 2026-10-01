######################
### DEB PARAMETERS ###
######################

# The parameters of the DEB_abj model, for Solea solea.
# Provided by Mounier et al. (2020)


v = 0.0724/2.7612    # energy conductance
kap = 0.7682   # fraction allocated to the soma

pM = 39.18     # Somatic maintenance rate
EG = 5430      # Costs of structure
k_J = 0.002    # Maturity maintenance rate coefficient

kap_R = 0.95   # Fraction allocated to gametes
kap_X = 0.8    # Digestion efficiency
kap_P = 0.1    # Faecation efficiency
E_Hb = 0.285   # Threshold of birth
E_Hj = 6.039   # Threshold of metamorphosis
E_Hp = 150600  # Threshold of puberty
del_l = 0.166  # shape coefficient before metamorphosis
del_m = 0.171  # shape coefficient after metamorphosis
d_v = 0.2      # Density of structure
q_E = 23012    # Energy density of reserve
pAm = 710/2.7612     # Max assimilation rate
D = 0.2          # % of dry content
w_V = 23.9
mu_V = 500000
E0 = 2.6
h_a = 2.529e-09  # Weibull aging acceleration
s_G = 0.0001     # Gompertz stress coefficient

M_V = d_v/ w_V     # mol/cm^3, volume-specific mass of structure
kap_G = mu_V * M_V / EG 

# ----- TEMPERATURE PARAMETERS
T_L = 276      # Lower bound. of tolerance range
T_H = 303      # Upper bound. of tolerance range
T_A = 5119     # Arrhenius T
T_AL = 50000   # rate of decrease at lower bound.
T_AH = 100000  # rate of decrease at upper bound
T_ref = 293.5  # Demander a bastien

# Time step
dt = dt  

# create the list of parameters
param_cont = as.list(c(v=v, kap=kap,kap_X=kap_X ,kap_P=kap_P , pM=pM, EG=EG, k_J=k_J, 
                       kap_R=kap_R, E_Hb=E_Hb,E_Hj=E_Hj, E_Hp=E_Hp, del_l=del_l,
                       del_m = del_m, d_v=d_v,q_E=q_E, kap_G=kap_G, T_L=T_L, T_H=T_H,
                       T_A=T_A, T_AL=T_AL,T_AH=T_AH, T_ref=T_ref, dt=dt, D=D, pAm=pAm,
                       w_V=w_V, mu_V,E0 = E0, h_a = h_a, s_G = s_G))


# Calculation of energy fluxes

energy_flux <- function(debmodel, params) {
  
  debmodel$pC = (params$EG * debmodel$s_M * (params$v*debmodel$TC) * ((debmodel$L^3)^(-1/3)) + (params$pM*debmodel$TC)) / ((params$EG / debmodel$E) + (params$kap / (debmodel$L)^3))
  debmodel$pM = (params$pM*debmodel$TC)*(debmodel$L^3)  # Flux of energy to somatic maintenance (J/day)
  debmodel$pG_tot = (params$kap*debmodel$pC)-debmodel$pM   # Flux of energy to the soma (J/day)
  debmodel$pG_growth = debmodel$pG_tot*params$kap_G  
  debmodel$G_sur_M <- debmodel$pG_tot/debmodel$pM # Ratio between growth & maintenance
  debmodel$K <- debmodel$pG_growth / (debmodel$pC/0.8) # growth efficiency coefficient
  debmodel$pJ = (params$k_J*debmodel$TC)*debmodel$H
  debmodel$pR = (1 - params$kap)*debmodel$pC - debmodel$pJ # Flux of energy that goes to maturity
  
  return(debmodel)
}