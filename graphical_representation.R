################################################################################
# -------- Code for all graphical representations in the manuscript ---------- #
# -------- Tristan Halna du Fretay - PhD

# Importation

library(dplyr)
library(lubridate)
library(data.table)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(progress)
library(tidyverse)
library(patchwork)

load("model_scenario1.RData")
load("model_scenario2.RData")
load("unlimited_feed_dens.RData")
load("temp_annual.RData")
load("predicted_rcp26_lag.RData")
load("predicted_rcp85_lag.RData")

# Creation of maturity THV dataset

test_birth_maturity_s1 <- test_birth_scenario_total %>%   # Unlimited feeding scenario
  group_by(scenario_temp, birth, cohort) %>%
  slice(1) %>%
  select(temperature, cohort, tpcont, Lpcont, scenario, scenario_temp, birth)

test_birth_maturity_s2 <- test_birth_scenario_total_limited %>% # Limited feeding scenario
  filter(scenario_temp > 1) %>%
  group_by(scenario_temp, birth, cohort) %>%
  slice(1) %>%
  select(temperature, cohort, tpcont, Lpcont, scenario, scenario_temp, birth)

test_birth_scenario_maturity_total <- rbind(test_birth_maturity_s1, test_birth_maturity_s2) # Selecting only life history traits at maturity

# Creation of daily temperature dataset between 2090-2098

rcp26_2090 <- predicted_rcp_26_with_lag%>%filter(date>="2090-01-01")
rcp85_2090 <- predicted_rcp_85_with_lag%>%filter(date>="2090-01-01")

# Graphs of the main manuscript


# Fig Temperature (now in Supplementary)
ggplot(data = graph, aes(x = year, y = temperature_moyenne)) +
  theme_classic2() +
  geom_point(aes(color = scenario)) +
  geom_line(aes(color = scenario)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(12, 18)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(2000, 2100)) +
  labs(color = "Scenario") +
  scale_color_manual(values = c("#234a12", "#2986cc"), labels = c("RCP 2.6", "RCP 8.5")) +
  xlab("") +
  ylab("Mean annual temperature (°C)") 

ggsave("Graphiques/temp.png", width = 1220, height = 820, units = "px", scale = 2)

# Fig Temperature 2090-2100
ggplot() +
  theme_classic2() +
  geom_line(data = rcp26_2090, mapping = aes(x = date, y = temperature, color = "1")) +
  geom_line(data = rcp85_2090, mapping = aes(x = date, y = temperature, color = "2")) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 30),
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  labs(color = "Scenario") +
  scale_color_manual(values = c("#234a12", "#2986cc"), labels = c("RCP 2.6", "RCP 8.5")) +
  scale_x_date(limits = as.Date(c("2090-01-01", "2100-01-01")),
               expand = c(0,0),
               date_breaks = "2 year",
               date_labels = "%Y",
               ) +
  xlab("") +
  ylab("Annual temperature (°C)") 

ggsave("Graphiques/new_temp.png", width = 1220, height = 820, units = "px", scale = 2)

# Fig Length diff
a = ggplot(test_birth_mean) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = low_diff_len, ymax = up_diff_len, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = diff_len, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_point(data = test_birth_mean[which(test_birth_mean$dpf == round(test_birth_mean$age_mat)),], 
             mapping = aes(x = age_mat, y = diff_len, fill = as.factor(scenario_temp)),
             shape = 23,
             size = 3) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  scale_y_continuous(expand = c(0,0),limits = c(-100,100))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill = "Climatic conditions")    +
  ylab("Relative difference in length (%)") +
  xlab("Time (dpf)")  +
  guides(color = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 45,
           y = 100, label = "A", size = 8, vjust = 2
  )

b = ggplot(test_birth_mean_limited) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = low_diff_len, ymax = up_diff_len, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = diff_len, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_point(data = test_birth_mean_limited[which(test_birth_mean_limited$dpf == round(test_birth_mean_limited$age_mat)),], 
             mapping = aes(x = age_mat, y = diff_len, fill = as.factor(scenario_temp)),
             shape = 23,
             size = 3) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  scale_y_continuous(expand = c(0,0),limits = c(-100,100))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill = "Climatic conditions")    +
  ylab("Relative difference in length (%)") +
  xlab("Time (dpf)")  +
  guides(color = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 45,
           y = 100, label = "B", size = 8, vjust = 2)

a + b + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Graphiques/length dynamic.png", width = 1220, height = 820, units = "px", scale = 2)

#Fig density diff
letters_signif_dens = data.frame(period = c("ref", "rcp2.6", "rcp8.5"), y = c(60, 33, 7), letters = c("a", "a", "b"))

ggplot(dens_estimated, (aes(y = diff_relative, x = period, fill = period))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.5) +
  stat_summary(fun.y = mean, geom = "point", shape = 18, size = 4) +
  theme_classic2() +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"),
                    labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(-100, 100)) +
  scale_x_discrete(labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  ylab("Relative difference in estimated density (%)") +
  xlab("") +
  guides(fill = FALSE) +
  geom_text(data = letters_signif_dens, mapping = aes(x = period, y = y, label = letters)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12)) 

ggsave("Graphiques/dens.png", width = 1220, height = 820, units = "px", scale = 2)

#Fig LHT maturity

test_birth_scenario_maturity_total$fill_graph <- with(
  test_birth_scenario_maturity_total,
  ifelse(
    scenario_temp == 1,
    "Reference",
    as.character(scenario)
  )
)

test_birth_scenario_maturity_total$fill_graph <- factor(test_birth_scenario_maturity_total$fill_graph,
                                                        levels = c("Reference", "1", "2"))

lettre_age_mat = data.frame(x = c(1, 1.82, 2.18, 2.82, 3.18),
                            y =c(700, 660, 700, 630, 700),
                            lettre = c("a", "bc", "a", "c", "ab"))

A = ggplot(test_birth_scenario_maturity_total) +
  theme_classic2() +
  geom_boxplot(aes(x = as.factor(scenario_temp), y = tpcont, fill = as.factor(fill_graph))) +
  scale_fill_manual(values = c("#fce7a7" , "darkolivegreen4", "darkorange3"), labels = c("Reference", "Unlimited", "Limited")) +
  scale_x_discrete(labels = c("1" = "Reference", "2" = "RCP 2.6", "3" = "RCP 8.5")) +
  ylab("Age at maturity (dpf)") +
  xlab("") +
  labs(fill = "Food limitation scenario") +
  scale_y_continuous(limits = c(400, 800),
                     expand = c(0,0)) +
  annotate("text", x = 1,
           y = 800, label = "A", size = 6, vjust = 2, hjust = 1.5) +
  geom_text(data = lettre_age_mat, mapping = aes(x = x, y = y, label = lettre)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))


lettre_len_mat = data.frame(x = c(1, 1.82, 2.18, 2.82, 3.18),
                            y =c(22.7, 22.7, 22.6, 22.6, 22.5),
                            lettre = c("a", "ab", "bc", "c", "c")) 

B = ggplot(test_birth_scenario_maturity_total) +
  theme_classic2() +
  geom_boxplot(aes(x = as.factor(scenario_temp), y = Lpcont, fill = as.factor(fill_graph))) +
  scale_fill_manual(values = c("#fce7a7" , "darkolivegreen4", "darkorange3"), labels = c("Reference", "Unlimited", "Limited")) +
  scale_x_discrete(labels = c("1" = "Reference", "2" = "RCP 2.6", "3" = "RCP 8.5")) +
  ylab("Length at maturity (cm)") +
  xlab("") +
  labs(fill = "Food limitation scenario") +
  scale_y_continuous(limits = c(21, 24),
                     expand = c(0,0)) +
  annotate("text", x = 1,
           y = 24, label = "B", size = 6, vjust = 2, hjust = 1.5) +
  geom_text(data = lettre_len_mat, mapping = aes(x = x, y = y, label = lettre)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

A + B + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Graphiques/THV mat.png", width = 1220, height = 820, units = "px", scale = 2)

# Fig reserve energy
a = ggplot(test_birth_mean%>%filter(dpf > 50)) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = low_diff_E, ymax = up_diff_E, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = diff_E, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_point(data = test_birth_mean[which(test_birth_mean$dpf == round(test_birth_mean$age_mat)),], 
             mapping = aes(x = age_mat, y = diff_E, fill = as.factor(scenario_temp)),
             shape = 23,
             size = 3) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  scale_y_continuous(expand = c(0,0), limits = c(-100, 100))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill = "Climatic conditions")    +
  ylab("Relative change in reserve density (%)") +
  xlab("Time (dpf)")  +
  guides(color = FALSE, fill = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 80,
           y = 100, label = "A", size = 8, vjust = 2)+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

b = ggplot(test_birth_mean_limited%>%filter(dpf > 50)) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = low_diff_E, ymax = up_diff_E, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = diff_E, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_point(data = test_birth_mean_limited[which(test_birth_mean_limited$dpf == round(test_birth_mean_limited$age_mat)),], 
             mapping = aes(x = age_mat, y = diff_E, fill = as.factor(scenario_temp)),
             shape = 23,
             size = 3) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  scale_y_continuous(expand = c(0,0), limits = c(-100, 100))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill = "Climatic conditions")    +
  ylab("") +
  xlab("Time (dpf)")  +
  guides(color = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 85,
           y = 100, label = "B", size = 8, vjust = 2)+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

a + b + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Graphiques/reserve dens dynamic.png", width = 1220, height = 820, units = "px", scale = 2)

# SUPPLEMENTARY MATERIALS

# Fig Length dynamic
a = ggplot(test_birth_mean) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = lower_range, ymax = upper_range, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = length, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_point(data = test_birth_mean[which(test_birth_mean$dpf == round(test_birth_mean$age_mat)),], 
             mapping = aes(x = age_mat, y = length, fill = as.factor(scenario_temp)),
             shape = 23,
             size = 3) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  scale_y_continuous(expand = c(0,0),limits = c(-0,40))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill = "Climatic conditions")    +
  ylab("Length (cm)") +
  xlab("Time (dpf)")  +
  guides(color = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 45,
           y = 100, label = "A", size = 8, vjust = 2
  )

b = ggplot(test_birth_mean_limited) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = lower_range, ymax = upper_range, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = length, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_point(data = test_birth_mean_limited[which(test_birth_mean_limited$dpf == round(test_birth_mean_limited$age_mat)),], 
             mapping = aes(x = age_mat, y = length, fill = as.factor(scenario_temp)),
             shape = 23,
             size = 3) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  scale_y_continuous(expand = c(0,0),limits = c(-0,40))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill = "Climatic conditions")    +
  ylab("Length (cm)") +
  xlab("Time (dpf)")  +
  guides(color = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 45,
           y = 100, label = "B", size = 8, vjust = 2)

a + b + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Graphiques/supplementary/total length dynamic.png", width = 1220, height = 820, units = "px", scale = 2)


# Annex 4 - Dynamic of p_A

# Estimation of the difference between theoretical assimilation (p_A wihtout limitation) and limited p_A

ref = pa_tot_mean%>%rename(no_lim_ref_p_A = p_A, ref_L = L, ref_TC = TC, ref_SM = SM)

test_birth_scenario_total_limited <- test_birth_scenario_total_limited %>%
  group_by(cohort, scenario_temp, birth) %>%
  left_join(ref, by = c("birth", "dpf")) %>%
  mutate(theoretical_p_A = p_A_ref / (ref_L**2 * ref_SM * ref_TC) * L**2 * s_M * TC,
         limit_p_A = (p_A - theoretical_p_A)/theoretical_p_A * 100)

test_birth_mean_limited <- test_birth_scenario_total_limited%>%
  group_by(dpf, scenario_temp) %>%
  summarize(length = mean(estim_L_cont),
            lower_range = min(estim_L_cont),
            upper_range = max(estim_L_cont),
            weight = mean(estim_W_cont),
            lower_weight = min(estim_W_cont),
            upper_weight = max(estim_W_cont),
            scaled_E = mean(E / (L^3)),
            lower_scaled_E = min(E / (L^3)),
            upper_scaled_E = max(E / (L^3)),
            age_mat = mean(tpcont),
            low_age_mat = min(tpcont),
            up_age_mat = max(tpcont),
            len_mat = mean(Lpcont),
            low_len_mat = min(Lpcont),
            up_len_mat = max(Lpcont),
            diff_len = mean(delta_len),
            low_diff_len = min(delta_len),
            up_diff_len = max(delta_len),
            diff_E = mean(delta_E),
            low_diff_E = min(delta_E),
            up_diff_E = max(delta_E),
            p_R = mean(pR),
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
            f_estim = mean(f_estim),
            temperature = first(temperature),
            low_p_A = min(p_A),
            up_p_A = max(p_A),
            p_A = mean(p_A),
            low_p_A = min(p_A),
            up_p_A = max(p_A),
            p_A_ref = mean(p_A_ref),
            theoretical_p_A = mean(theoretical_p_A),
            limit_p_A = mean(limit_p_A),
            L = mean(L))

# Figure S6

ggplot(test_birth_mean_limited %>% filter(scenario_temp != 1)) + theme_classic2() +
  geom_line(aes(x = dpf, y = p_A, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_line(aes(x = dpf, y = p_A_ref, color = "Reference"), linewidth = 0.75) +
  scale_color_manual(values = c("Reference" = "black", "2" = "#234a12", "3" = "#2986cc"), 
                     labels = c("Reference" = "Reference" ,"2" = "RCP 2.6","3" = "RCP 8.5"),
                     breaks = c("Reference", "2", "3")) +
  scale_y_continuous(expand = c(0,0))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(color = "Climatic scenario")    +
  ylab("Assimilation flux p_A (in J/day)") +
  xlab("Time (in dpf)")  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12))

ggsave("Graphiques/supplementary/p_A_dynamic_lim_feeding.png", width = 1020, height = 820, units = "px", scale = 2)

# Figure S7

ggplot(test_birth_mean_limited%>%filter(scenario_temp != 1)) + theme_classic2() +
  geom_line(aes(x = dpf, y = limit_p_A, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_hline(aes(yintercept = 0, color = "0"), linewidth = 2) +
  scale_color_manual(values = c("0" = "black", "2" = "#234a12", "3" = "#2986cc"), labels = c("0" = "Reference", "2" = "RCP2.6","3" = "RCP8.5"),
                     breaks = c("0", "2", "3")) +
  geom_hline(aes(yintercept = 0)) +
  scale_y_continuous(expand = c(0,0), limits = c(-50, 10))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(color = "Climatic scenario")    +
  ylab("Limitation of assimilation flux (%)") +
  xlab("Time (dpf)")  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12))

ggsave("Graphiques/supplementary/limitation energy.png", width = 1220, height = 820, units = "px", scale = 2)

# Annex 5 - f values for varying birthdates

load("f_graph.RData")

ggplot(f_mean) +
  theme_classic2() +
  geom_line(f_meanjan, mapping = aes(x = time, y = f, color = "January"), linewidth = 1.3) +
  geom_line(aes(x = time, y = f, color = "February"), linewidth = 1.3) +
  geom_line(f_meanmars, mapping = aes(x = time, y = f, color = "Mars"), linewidth = 1.3) +
  scale_color_manual(values = c("skyblue", "olivedrab", "goldenrod3"), labels = c("1st January", "1st February", "1st March")) +
  labs(color = "Birthdate") +
  xlab("Time (dpf)") +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 1000)) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 2)) +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.text = element_text(size = 18))  

ggsave("Graphiques/supplementary/fbackcalc.png", width = 1120, height = 820, units = "px", scale = 2)

# Annex 6 - Reserve density dynamics

a = ggplot(test_birth_mean%>%filter(dpf > 50)) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = lower_scaled_E, ymax = upper_scaled_E, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = scaled_E, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_point(data = test_birth_mean[which(test_birth_mean$dpf == round(test_birth_mean$age_mat)),], 
             mapping = aes(x = age_mat, y = scaled_E, fill = as.factor(scenario_temp)),
             shape = 23,
             size = 3) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  scale_y_continuous(expand = c(0,0))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill = "Climatic conditions")    +
  ylab("Reserve density (J/cm3)") +
  xlab("Time (dpf)")  +
  guides(color = FALSE, fill = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 80,
           y = 100, label = "A", size = 8, vjust = 2)

b = ggplot(test_birth_mean_limited%>%filter(dpf > 50)) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = lower_scaled_E, ymax = upper_scaled_E, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = scaled_E, color = as.factor(scenario_temp)), linewidth = 0.75) +
  geom_point(data = test_birth_mean_limited[which(test_birth_mean_limited$dpf == round(test_birth_mean_limited$age_mat)),], 
             mapping = aes(x = age_mat, y = scaled_E, fill = as.factor(scenario_temp)),
             shape = 23,
             size = 3) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  scale_y_continuous(expand = c(0,0))  +
  scale_x_continuous(expand = c(0,0)) +
  labs(fill = "Climatic conditions")    +
  ylab("") +
  xlab("Time (dpf)")  +
  guides(color = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 85,
           y = 100, label = "B", size = 8, vjust = 2)

a + b + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Graphiques/supplementary/reserve density dynamic.png", width = 1220, height = 820, units = "px", scale = 2)


# Annex 7 - Metabolic fluxes

# Reference

ggplot(test_birth_mean%>%filter(scenario_temp == 1)) + theme_classic2() +
  geom_line(aes(x = dpf, y = p_G, color = "p_G", linetype = "Unlimited feeding"), linewidth = 0.75) +
  geom_line(aes(x = dpf, y = p_R, color = "p_R"), linewidth = 0.75) +
  geom_line(aes(x = dpf, y = p_C, color = "p_C"), linewidth = 0.75) +
  geom_line(data = test_birth_mean_limited%>%filter(scenario_temp == 1),
            mapping = aes(x = dpf, y = p_G, color = "p_G", linetype = "Limited feeding"), linewidth = 0.75) +
  geom_line(data = test_birth_mean_limited%>%filter(scenario_temp == 1),
            mapping = aes(x = dpf, y = p_R, color = "p_R", linetype = "Limited feeding"), linewidth = 0.75) +
  geom_line(data = test_birth_mean_limited%>%filter(scenario_temp == 1),
            mapping = aes(x = dpf, y = p_C, color = "p_C", linetype = "Limited feeding"), linewidth = 0.75) +
  scale_y_continuous(expand = c(0,0), limits = c(-500, 15000))  +
  scale_x_continuous(expand = c(0,0)) +
  scale_linetype_manual(values = c("Unlimited feeding" = "solid", "Limited feeding" = "dashed")) +
  labs(color = "Flux", linetype = "Scenario")    +
  ylab("Energy repartition (J/day)") +
  xlab("Time (dpf)")  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) 

ggsave("Graphiques/supplementary/flux_comparison_ref.png", width = 1020, height = 820, units = "px", scale = 2)

# RCP2.6

ggplot(test_birth_mean%>%filter(scenario_temp == 2)) + theme_classic2() +
  geom_line(aes(x = dpf, y = p_G, color = "p_G", linetype = "Unlimited feeding"), linewidth = 0.75) +
  geom_line(aes(x = dpf, y = p_R, color = "p_R"), linewidth = 0.75) +
  geom_line(aes(x = dpf, y = p_C, color = "p_C"), linewidth = 0.75) +
  geom_line(data = test_birth_mean_limited%>%filter(scenario_temp == 2),
            mapping = aes(x = dpf, y = p_G, color = "p_G", linetype = "Limited feeding"), linewidth = 0.75) +
  geom_line(data = test_birth_mean_limited%>%filter(scenario_temp == 2),
            mapping = aes(x = dpf, y = p_R, color = "p_R", linetype = "Limited feeding"), linewidth = 0.75) +
  geom_line(data = test_birth_mean_limited%>%filter(scenario_temp == 2),
            mapping = aes(x = dpf, y = p_C, color = "p_C", linetype = "Limited feeding"), linewidth = 0.75) +
  scale_y_continuous(expand = c(0,0), limits = c(-500, 15000))  +
  scale_x_continuous(expand = c(0,0)) +
  scale_linetype_manual(values = c("Unlimited feeding" = "solid", "Limited feeding" = "dashed")) +
  labs(color = "Flux", linetype = "Scenario")    +
  ylab("Energy repartition (J/day)") +
  xlab("Time (dpf)")  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) 

ggsave("Graphiques/supplementary/flux_comparison_rcp2.6.png", width = 1020, height = 820, units = "px", scale = 2)

# RCP8.5

ggplot(test_birth_mean%>%filter(scenario_temp == 3)) + theme_classic2() +
  geom_line(aes(x = dpf, y = p_G, color = "p_G", linetype = "Unlimited feeding"), linewidth = 0.75) +
  geom_line(aes(x = dpf, y = p_R, color = "p_R"), linewidth = 0.75) +
  geom_line(aes(x = dpf, y = p_C, color = "p_C"), linewidth = 0.75) +
  geom_line(data = test_birth_mean_limited%>%filter(scenario_temp == 3),
            mapping = aes(x = dpf, y = p_G, color = "p_G", linetype = "Limited feeding"), linewidth = 0.75) +
  geom_line(data = test_birth_mean_limited%>%filter(scenario_temp == 3),
            mapping = aes(x = dpf, y = p_R, color = "p_R", linetype = "Limited feeding"), linewidth = 0.75) +
  geom_line(data = test_birth_mean_limited%>%filter(scenario_temp == 3),
            mapping = aes(x = dpf, y = p_C, color = "p_C", linetype = "Limited feeding"), linewidth = 0.75) +
  scale_y_continuous(expand = c(0,0), limits = c(-500, 15000))  +
  scale_x_continuous(expand = c(0,0)) +
  scale_linetype_manual(values = c("Unlimited feeding" = "solid", "Limited feeding" = "dashed")) +
  labs(color = "Flux", linetype = "Scenario")    +
  ylab("Energy repartition (J/day)") +
  xlab("Time (dpf)")  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) 

ggsave("Graphiques/supplementary/flux_comparison_rcp8.5.png", width = 1020, height = 820, units = "px", scale = 2)


# Annex 8 - Biomass

ggplot(dens_estimated, mapping = aes(y = (mean_biomass - mean_biomass_tot) / mean_biomass_tot * 100, x = period, fill = period)) +
  geom_boxplot(width = 0.5) +
  theme_classic2() +
  scale_y_continuous(limits = c(-50, 50), expand = c(0,0)) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP2.6", "RCP8.5")) +
  scale_x_discrete(labels = c("Reference", "RCP 2.6", "RCP 8.5")) +
  ylab("Relative change in biomass (%)") +
  xlab("") +
  labs(fill = "Climatic conditions") +
  geom_signif(comparisons = list(c("ref", "rcp2.6"),
                                 c("ref", "rcp8.5"),
                                 c("rcp2.6", "rcp8.5")),
              map_signif_level = function(p) {
                ifelse(p > 0.05, "ns",   
                       ifelse(p > 0.01, "*",
                              ifelse(p > 0.001, "**", "***")))
              },
              y_position = c(15, 22, 10))

ggsave("Graphiques/supplementary/biomass_relative.png", width = 1020, height = 820, units = "px", scale = 2)

data.frame(moyenne = aggregate((dens_estimated$mean_biomass - mean_biomass_tot) / mean_biomass_tot * 100, by = list(dens_estimated$period), FUN = mean),
           sd = aggregate((dens_estimated$mean_biomass - mean_biomass_tot) / mean_biomass_tot * 100, by = list(dens_estimated$period), FUN = sd))


# Annex 9 - Fecundity

fecundity_unlimited_feed <- test_birth_scenario_total%>%
  group_by(dpf, scenario_temp) %>%
  summarize(fecundity = mean(F),
            lower_fecundity = min(F),
            upper_fecundity = max(F),
  )

fecundity_limited_feed <- test_birth_scenario_total_limited%>%
  group_by(dpf, scenario_temp) %>%
  summarize(fecundity = mean(F),
            lower_fecundity = min(F),
            upper_fecundity = max(F),
  )

A = ggplot(fecundity_unlimited_feed) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = lower_fecundity, ymax = upper_fecundity, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = fecundity, color = as.factor(scenario_temp)), linewidth = 0.75) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP2.6", "RCP8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP2.6", "RCP8.5")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 250000))  +
  scale_x_continuous(expand = c(0,0), limits = c(600, 940)) +
  labs(fill = "Climatic conditions")    +
  ylab("Fecundity (number of eggs)") +
  xlab("Time (dpf)")  +
  guides(fill = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 620,
           y = 250000, label = "A", size = 8, vjust = 2) +
  labs(color = "Climatic conditions")

B = ggplot(fecundity_limited_feed) + theme_classic2() +
  geom_ribbon(aes(x = dpf, ymin = lower_fecundity, ymax = upper_fecundity, fill = as.factor(scenario_temp), color = as.factor(scenario_temp)),
              alpha = 0.25,
              linetype = "dashed",
              linewidth = 0.25) +
  geom_line(aes(x = dpf, y = fecundity, color = as.factor(scenario_temp)), linewidth = 0.75) +
  scale_fill_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP2.6", "RCP8.5"))  +
  scale_color_manual(values = c("#fce7a7", "#234a12", "#2986cc"), labels = c("Reference", "RCP2.6", "RCP8.5")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 250000))  +
  scale_x_continuous(expand = c(0,0), limits = c(600, 940)) +
  labs(fill = "Climatic conditions")    +
  ylab("") +
  xlab("Time (dpf)")  +
  guides(fill = FALSE)  +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  annotate("text", x = 620,
           y = 250000, label = "B", size = 8, vjust = 2)+
  labs(color = "Climatic conditions")

A + B + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Graphiques/supplementary/fecundity dynamic.png", width = 1220, height = 820, units = "px", scale = 2)

