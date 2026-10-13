## Validation of GEE design effects for partial clustered trials using blocked randomisation
## Kylie Lange (kylie.lange@sahmri.com)

##### HELPER FUNCTIONS #####

# Calculate gamma_s, the proportion of observations that belong to clusters of size 1, from the number of clusters of each size

calc_gamma_s = function(num_clus1, num_clus2) {
  N = num_clus1 + 2*num_clus2
  gamma_s = num_clus1 / N
  gamma_s
}

##### DESIGN EFFECTS - CONTINUOUS #####

# DEFFs under simple randomisation
# From Yelland et al 2017 (DOI 10.1002/sim.7201)
# individual randomisation, continuous outcome

# Inputs: gee working correlation structure, rho, proportion of obs that belong to a cluster of size 1
# Output: design effect

simpledeff_cont = function(gee_corr, rho, gamma_s) {
  if (gee_corr == "independence") {
    deff = 1
  }
  else if (gee_corr == "exchangeable") {
    deff = (1 - rho^2) / (1 - (rho^2)*gamma_s)
  }
  deff
}

# DEFFs under blocked randomisation
# continuous outcomes

blockdeff_cont = function(gee_corr, rho, gamma_s, b, stratified = TRUE) {
  
  # if using variable block sizes, convert the string of block sizes to a vector
  if (is.character(b)) {
    b = as.numeric(unlist(strsplit(b, ",")))
  } else {
    b = b
  }
  # set equal probabilities for each block size
  p = rep(1/length(b), times = length(b))
  # weightings for each block size
  w = p*b / sum(unlist(p*b))
  
  if (stratified == TRUE) {
    psi_II = sum(unlist(w*(b - 2) / (4*b - 4)))
    psi_CC = sum(unlist(w*(b - 2) / (4*b - 4)))
    psi_IC = 1 - 2*sum(unlist(w*(b - 2) / (4*b - 4)))
  } else if (stratified == FALSE) {
    psi_II = sum(unlist(w*(b - 1) / (4*b)))
    psi_CC = sum(unlist(w*(b - 1) / (4*b)))
    psi_IC = 1 - 2*sum(unlist(w*(b - 1) / (4*b))) }
  
  if (gee_corr == "independence") {
    deff = 1 + rho*(1 - gamma_s)*(psi_II + psi_CC - psi_IC)
  }
  else if (gee_corr == "exchangeable") {
    deff = (1 - rho^2) / (1 - (rho^2)*gamma_s - rho*(1 - gamma_s)*(psi_II + psi_CC - psi_IC))
  }
  deff
}

# Expected DEFF and power under blocked randomisation
# Inputs: trial params including rho, number of clusters of each size, treatment effect, block size
# Outputs: expected power (power of a Wald test for independent data for the specified treatment effect
#          and effective sample size), effective sample size (number of observations / expected deff under blocked randomisation), for both gee correlation structures

calculate_blockdeffs_cont = function(params) {
  gamma_s = calc_gamma_s(params$num_clus1, params$num_clus2)
  N       = params$num_clus1 + 2*params$num_clus2
  
  # independence gee
  ind_deff      = blockdeff_cont(gee_corr = "independence", rho = params$ICC, gamma_s = gamma_s, b = params$block_size,
                                 stratified = params$stratified)
  ind_N_eff     = N / ind_deff
  ind_power_exp = pwrss.t.2means(mu1 = params$b0, mu2 = params$b1, sd1 = 1, sd2 = 1, kappa = 1,
                                 n2 = ind_N_eff/2, alpha = 0.05, verbose=FALSE)$power
  ind_res       = data.frame(deff = ind_deff, N_eff = ind_N_eff, power_exp = ind_power_exp, model = "gee-ind")
  
  # exchangeable gee
  exch_deff      = blockdeff_cont(gee_corr = "exchangeable", rho = params$ICC, gamma_s = gamma_s, b = params$block_size,
                                  stratified = params$stratified)
  exch_N_eff     = N / exch_deff
  exch_power_exp = pwrss.t.2means(mu1 = params$b0, mu2 = params$b1, sd1 = 1, sd2 = 1, kappa = 1,
                                  n2 = exch_N_eff/2, alpha = 0.05, verbose=FALSE)$power
  exch_res       = data.frame(deff = exch_deff, N_eff = exch_N_eff, power_exp = exch_power_exp, model = "gee-exch")
  
  res = bind_rows(ind_res, exch_res)
  res
}

##### DESIGN EFFECTS - LOGISTIC #####

# DEFFs under simple randomisation
# From Yelland et al 2017 (DOI 10.1002/sim.7201)
# individual randomisation, binary outcome with logit link

# Inputs: gee working correlation structure, rho, proportion of obs that belong to a cluster of size 1, prevalence
# Output: design effect

simpledeff_logit = function(gee_corr, rho, gamma_s, pC, pI) {
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

# DEFFs under blocked randomisation
# binary outcome with logit link

blockdeff_logit = function(gee_corr, rho, gamma_s, b, pC, pI, stratified = TRUE) {
  
  # if using variable block sizes, convert the string of block sizes to a vector
  if (is.character(b)) {
    b = as.numeric(unlist(strsplit(b, ",")))
  } else {
    b = b
  }
  # set equal probabilities for each block size
  p = rep(1/length(b), times = length(b))
  # weightings for each block size
  w = p*b / sum(unlist(p*b))
  
  if (stratified == TRUE) {
    psi_II = sum(unlist(w*(b - 2) / (4*b - 4)))
    psi_CC = sum(unlist(w*(b - 2) / (4*b - 4)))
    psi_IC = 1 - 2*sum(unlist(w*(b - 2) / (4*b - 4)))
  } else if (stratified == FALSE) {
    psi_II = sum(unlist(w*(b - 1) / (4*b)))
    psi_CC = sum(unlist(w*(b - 1) / (4*b)))
    psi_IC = 1 - 2*sum(unlist(w*(b - 1) / (4*b))) }
  
  if (gee_corr == "independence") {
    deff = 1 + 2*rho*(1 - gamma_s) / (pI*(1-pI)+pC*(1-pC)) *
      (psi_II*pC*(1-pC) + psi_CC*pI*(1-pI) - psi_IC*sqrt(pI*pC*(1-pI)*(1-pC)))
  }
  else if (gee_corr == "exchangeable") {
    deff = (1 - rho^2) * 
      (1 - (rho^2)*gamma_s - rho*(1 - gamma_s)*(psi_II + psi_CC + psi_IC*(2*sqrt(pI*pC*(1-pI)*(1-pC)) / 
                                                                           (pI*(1-pI) + pC*(1-pC))))) /
      ((1 - (rho^2)*gamma_s - rho*(1 - gamma_s)*(psi_II + psi_CC - psi_IC)) *
      (1 - (rho^2)*gamma_s - rho*(1 - gamma_s)))
  }
  deff
}

# Expected DEFF and power under blocked randomisation
# Inputs: trial params including rho, number of clusters of each size, outcome prevalence, block size
# Outputs: expected power (power of a Wald test for independent data for the specified treatment effect
#          and effective sample size), effective sample size, for both gee correlation structures

calculate_blockdeffs_logit = function(params) {
  gamma_s = calc_gamma_s(params$num_clus1, params$num_clus2)
  N       = params$num_clus1 + 2*params$num_clus2
  
  # independence gee
  ind_deff      = blockdeff_logit(gee_corr = "independence", rho = params$ICC, gamma_s = gamma_s, b = params$block_size,
                                 pI = params$pI, pC = params$pC, stratified = params$stratified)
  ind_N_eff     = N / ind_deff
  ind_power_exp = pwrss.z.logreg(p0=params$pC, p1=params$pI, n = ind_N_eff, alpha = 0.05,
                                 distribution="bernoulli", verbose=FALSE)$power
  ind_res       = data.frame(deff = ind_deff, N_eff = ind_N_eff, power_exp = ind_power_exp, model = "gee-ind")
  
  # exchangeable gee
  exch_deff      = blockdeff_logit(gee_corr = "exchangeable", rho = params$ICC, gamma_s = gamma_s, b = params$block_size,
                                   pI = params$pI, pC = params$pC, stratified = params$stratified)
  exch_N_eff     = N / exch_deff
  exch_power_exp = pwrss.z.logreg(p0=params$pC, p1=params$pI, n = exch_N_eff, alpha = 0.05,
                                  distribution="bernoulli", verbose=FALSE)$power
  exch_res       = data.frame(deff = exch_deff, N_eff = exch_N_eff, power_exp = exch_power_exp, model = "gee-exch")
  
  res = bind_rows(ind_res, exch_res)
  res
}

##### DESIGN EFFECTS - LOG BINOMIAL #####

# DEFFs under simple randomisation
# From Yelland et al 2017 (DOI 10.1002/sim.7201)
# individual randomisation, binary outcome with log link

# Inputs: gee working correlation structure, rho, proportion of obs that belong to a cluster of size 1, prevalence
# Output: design effect

simpledeff_log = function(gee_corr, rho, gamma_s, pC, pI) {
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

# DEFFs under blocked randomisation
# binary outcome with log link

blockdeff_log = function(gee_corr, rho, gamma_s, b, pC, pI, stratified = TRUE) {
  
  # if using variable block sizes, convert the string of block sizes to a vector
  if (is.character(b)) {
    b = as.numeric(unlist(strsplit(b, ",")))
  } else {
    b = b
  }
  # set equal probabilities for each block size
  p = rep(1/length(b), times = length(b))
  # weightings for each block size
  w = p*b / sum(unlist(p*b))
  
  if (stratified == TRUE) {
    psi_II = sum(unlist(w*(b - 2) / (4*b - 4)))
    psi_CC = sum(unlist(w*(b - 2) / (4*b - 4)))
    psi_IC = 1 - 2*sum(unlist(w*(b - 2) / (4*b - 4)))
  } else if (stratified == FALSE) {
    psi_II = sum(unlist(w*(b - 1) / (4*b)))
    psi_CC = sum(unlist(w*(b - 1) / (4*b)))
    psi_IC = 1 - 2*sum(unlist(w*(b - 1) / (4*b))) }
  
  if (gee_corr == "independence") {
    deff = 1 + 2*rho*(1 - gamma_s) / (pI*(1-pC)+pC*(1-pI)) *
      (psi_II*pC*(1-pI) + psi_CC*pI*(1-pC) - psi_IC*sqrt(pI*pC*(1-pI)*(1-pC)))
  }
  else if (gee_corr == "exchangeable") {
    deff = (1 - rho^2) * 
      (1 - (rho^2)*gamma_s - rho*(1 - gamma_s)*(psi_II + psi_CC + psi_IC*(2*sqrt(pI*pC*(1-pI)*(1-pC)) / 
                                                                            (pC*(1-pI) + pI*(1-pC))))) /
      ((1 - (rho^2)*gamma_s - rho*(1 - gamma_s)*(psi_II + psi_CC - psi_IC)) *
         (1 - (rho^2)*gamma_s - rho*(1 - gamma_s)))
  }
  deff
}

# Expected DEFF and power under blocked randomisation
# Inputs: trial params including gee working correlation structure, rho, number of clusters of each size,
#         treatment effect, block size
# Outputs: expected power (power of a Wald test for independent data for the specified treatment effect
#          and effective sample size), effective sample size, for both gee correlation structures

calculate_blockdeffs_log = function(params) {
  gamma_s = calc_gamma_s(params$num_clus1, params$num_clus2)
  N       = params$num_clus1 + 2*params$num_clus2
  
  # independence gee
  ind_deff      = blockdeff_log(gee_corr = "independence", rho = params$ICC, gamma_s = gamma_s, b = params$block_size,
                                  pC = params$pC, pI=params$pI, stratified = params$stratified)
  ind_N_eff     = N / ind_deff
  ind_power_exp = pwrss.z.2props(p1=params$pC, p2=params$pI, kappa=1, n2 = ind_N_eff/2, alpha = 0.05, verbose=FALSE)$power
  ind_res       = data.frame(deff = ind_deff, N_eff = ind_N_eff, power_exp = ind_power_exp, model = "gee-ind")
  
  # exchangeable gee
  exch_deff      = blockdeff_log(gee_corr = "exchangeable", rho = params$ICC, gamma_s = gamma_s, b = params$block_size,
                                 pI = params$pI, pC = params$pC, stratified = params$stratified)
  exch_N_eff     = N / exch_deff
  exch_power_exp = pwrss.z.2props(p1=params$pC, p2=params$pI, kappa = 1, n2 = exch_N_eff/2, alpha = 0.05, verbose=FALSE)$power
  exch_res       = data.frame(deff = exch_deff, N_eff = exch_N_eff, power_exp = exch_power_exp, model = "gee-exch")
  
  res = bind_rows(ind_res, exch_res)
  res
}
