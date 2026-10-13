## Validation of GEE design effects for partial clustered trials using blocked randomisation
## Kylie Lange (kylie.lange@sahmri.com)

# Figures of theoretical design effect functions in manuscript

library(tidyverse)
library(patchwork)

source("R/deff functions.R")

### Binary DEFFs for individual and opposite simple randomisation
simpledeff_logit_indiv = function(gee_corr, rho, gamma_s, pC, pI) {
  pi_term = sqrt(pI*pC*(1-pI)*(1-pC)) / (pI*(1-pI) + pC*(1-pC))
  if (gee_corr == "independence") {
    deff = 1 + rho*(1-gamma_s)*(0.5 - pi_term)
  }
  else if (gee_corr == "exchangeable") {
    deff = ((1 - rho^2)*(1 - gamma_s*rho^2 - (1-gamma_s)*rho*(0.5 + pi_term))) /
      ((1 - gamma_s*rho^2)*(1 - gamma_s*rho^2 - (1-gamma_s)*rho))
  }
  deff
}
simpledeff_logit_opp = function(gee_corr, rho, gamma_s, pC, pI) {
  pi_term = sqrt(pI*pC*(1-pI)*(1-pC)) / (pI*(1-pI) + pC*(1-pC))
  if (gee_corr == "independence") {
    deff = 1 - rho*(1-gamma_s)*2*pi_term
  }
  else if (gee_corr == "exchangeable") {
    deff = (1 - rho)*(1 - gamma_s*rho^2 - (1-gamma_s)*rho*2*pi_term) / 
      ((1 - gamma_s*rho)*(1 - gamma_s*rho^2 - (1-gamma_s)*rho))
  }
  deff
}
simpledeff_log_indiv = function(gee_corr, rho, gamma_s, pC, pI) {
  pi_term = sqrt(pI*pC*(1-pI)*(1-pC)) / (pI*(1-pC) + pC*(1-pI))
  if (gee_corr == "independence") {
    deff = 1 + rho*(1-gamma_s)*(0.5 - pi_term)
  }
  else if (gee_corr == "exchangeable") {
    deff = ((1 - rho^2)*(1 - gamma_s*rho^2 - (1-gamma_s)*rho*(0.5 + pi_term))) /
      ((1 - gamma_s*rho^2)*(1 - gamma_s*rho^2 - (1-gamma_s)*rho))
  }
  deff
}
simpledeff_log_opp = function(gee_corr, rho, gamma_s, pC, pI) {
  pi_term = sqrt(pI*pC*(1-pI)*(1-pC)) / (pI*(1-pC) + pC*(1-pI))
  if (gee_corr == "independence") {
    deff = 1 - rho*(1-gamma_s)*2*pi_term
  }
  else if (gee_corr == "exchangeable") {
    deff = (1 - rho)*(1 - gamma_s*rho^2 - (1-gamma_s)*rho*2*pi_term) / 
      ((1 - gamma_s*rho)*(1 - gamma_s*rho^2 - (1-gamma_s)*rho))
  }
  deff
}

## Figure 1: Continuous DEFFs with and without stratification - ICC = 0.8
b = c(2,4,6,8,10,12,14)
gee_corr = c("independence", "exchangeable")
stratified = c(TRUE, FALSE)
gamma_s = 2/3
rho = 0.8
data1 = expand.grid(b = b, gee_corr = gee_corr, stratified = stratified) %>%
  mutate(gamma_s, rho,
         mod = case_when(
           ((stratified == TRUE) & (gee_corr == "independence")) ~ "GEE-ind with stratification",
           ((stratified == FALSE) & (gee_corr == "independence")) ~ "GEE-ind without stratification",
           ((stratified == TRUE) & (gee_corr == "exchangeable")) ~ "GEE-exch with stratification",
           ((stratified == FALSE) & (gee_corr == "exchangeable")) ~ "GEE-exch without stratification")) %>% 
  rowwise() %>% 
  mutate(deff = blockdeff_cont(gee_corr, rho, gamma_s, b, stratified))
data1$mod = factor(data1$mod, levels = c("GEE-ind without stratification", "GEE-ind with stratification",
                                            "GEE-exch without stratification", "GEE-exch with stratification"))

simple_ind_indiv1  = 1
simple_exch_indiv1 = (1 - rho^2) / (1 - (rho^2)*gamma_s)
simple_ind_opp1    = 1 - rho*(1-gamma_s)
simple_exch_opp1   = (1 - rho) / (1 - rho*gamma_s)

plot1 = ggplot(data = data1, aes(x = b, y = deff, group = mod, shape = mod, colour = mod)) +
  geom_point(aes(x = b, y = deff, colour = mod), size = 4, stroke = 1.5) +
  geom_hline(yintercept = simple_ind_indiv1, linetype = "dashed", linewidth = 1, colour = "#FF0000") +
  geom_hline(yintercept = simple_exch_indiv1, linetype = "dashed", linewidth = 1, colour = "#66B2FF") +
  geom_hline(yintercept = simple_ind_opp1, linetype = "dotted", linewidth = 1.5, colour = "#FF0000") +
  geom_hline(yintercept = simple_exch_opp1, linetype = "dotted", linewidth = 1.5, colour = "#66B2FF") +
  scale_shape_manual(name = "", labels = c("GEE-ind without stratification", "GEE-ind with stratification",
                                           "GEE-exch without stratification", "GEE-exch with stratification"),
                     values = c(1, 0, 1, 0)) +
  scale_color_manual(name = "", labels = c("GEE-ind without stratification", "GEE-ind with stratification",
                                           "GEE-exch without stratification", "GEE-exch with stratification"),
                     values = c("#FF0000", "#FF0000", "#66B2FF", "#66B2FF")) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14)) +
  labs(x = "Block size", y = "DEFF") +
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size= 16))

## Figure S1: Continuous DEFFs with and without stratification - ICC = 0.2
b = c(2,4,6,8,10,12,14)
gee_corr = c("independence", "exchangeable")
stratified = c(TRUE, FALSE)
gamma_s = 2/3
rho = 0.2
dataS1 = expand.grid(b = b, gee_corr = gee_corr, stratified = stratified) %>%
  mutate(gamma_s, rho,
         mod = case_when(
           ((stratified == TRUE) & (gee_corr == "independence")) ~ "GEE-ind with stratification",
           ((stratified == FALSE) & (gee_corr == "independence")) ~ "GEE-ind without stratification",
           ((stratified == TRUE) & (gee_corr == "exchangeable")) ~ "GEE-exch with stratification",
           ((stratified == FALSE) & (gee_corr == "exchangeable")) ~ "GEE-exch without stratification")) %>% 
  rowwise() %>% 
  mutate(deff = blockdeff_cont(gee_corr, rho, gamma_s, b, stratified))
dataS1$mod = factor(dataS1$mod, levels = c("GEE-ind without stratification", "GEE-ind with stratification",
                                           "GEE-exch without stratification", "GEE-exch with stratification"))

simple_ind_indivS1  = 1
simple_exch_indivS1 = (1 - rho^2) / (1 - (rho^2)*gamma_s)
simple_ind_oppS1    = 1 - rho*(1-gamma_s)
simple_exch_oppS1   = (1 - rho) / (1 - rho*gamma_s)

plotS1 = ggplot(data = dataS1, aes(x = b, y = deff, group = mod, shape = mod, colour = mod)) +
  geom_point(aes(x = b, y = deff, colour = mod), size = 4, stroke = 1.5) +
  geom_hline(yintercept = simple_ind_indivS1, linetype = "dashed", linewidth = 1, colour = "#FF0000") +
  geom_hline(yintercept = simple_exch_indivS1, linetype = "dashed", linewidth = 1, colour = "#66B2FF") +
  geom_hline(yintercept = simple_ind_oppS1, linetype = "dotted", linewidth = 1.5, colour = "#FF0000") +
  geom_hline(yintercept = simple_exch_oppS1, linetype = "dotted", linewidth = 1.5, colour = "#66B2FF") +
  scale_shape_manual(name = "", labels = c("GEE-ind without stratification", "GEE-ind with stratification",
                                           "GEE-exch without stratification", "GEE-exch with stratification"),
                     values = c(1, 0, 1, 0)) +
  scale_color_manual(name = "", labels = c("GEE-ind without stratification", "GEE-ind with stratification",
                                           "GEE-exch without stratification", "GEE-exch with stratification"),
                     values = c("#FF0000", "#FF0000", "#66B2FF", "#66B2FF")) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14)) +
  labs(x = "Block size", y = "DEFF") +
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size= 16))
