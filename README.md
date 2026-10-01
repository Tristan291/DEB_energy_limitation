This is the code used for the simulations and graphical representation of "The effect of warming on the bioenergetics of juvenile common sole (Solea solea) in the Gironde estuary: implications for abundance and life-history traits in an energy-limited nursery ground".

DEBfunction_limited_feeding.R: Function and ode used for the simulations under the limited individual feeding. This model includes a limitation in assimilation related to the assimilation under reference conditions.

DEBfunction_unlimited_feeding.R: Function and ode used for the simulations under the unlimited individual feeding. This model is a classical abj DEB model.

parameters.R: Set of parameters used in the DEB model (from Mounier et al. 2020) + a function used to calculate all energy fluxes.

sim_limited_feeding.R: Simulations under limited feeding scenario made under 3 birth dates and 3 temperature conditions.

sim_unlimited_feeding.R: Simulations under unlimited feeding scenario made under 3 birth dates and 3 temperature conditions.

graphical_representation.R: All graphical representation in the manuscript and the supplementary materials.
