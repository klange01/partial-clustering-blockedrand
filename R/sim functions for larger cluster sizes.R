## Validation of GEE design effects for partial clustered trials using blocked randomisation
## Kylie Lange (kylie.lange@sahmri.com)

# Functions for substudy in settings with larger cluster sizes
# Consider maximum cluster size k = 4


##### DESIGN EFFECTS #####

# Calculate the proportions of observations in clusters of each size
calc_gamma = function(num_clus1, num_clus2, num_clus3, num_clus4){
  M_k   = c(num_clus1, num_clus2, num_clus3, num_clus4)
  N     = sum(unlist(M_k %>% imap(~ .y*.x)))
  gamma = M_k %>% imap(~ .y*.x/N)   # gamma_k = k*M_k/N
  return(gamma)
}

# Design effects for partially clustered designs with simple individual randomisation and cluster sizes > 2
# Formulas derived in: Lange KM et al. Stats Med 2025;44(15-17):e70172. doi: 10.1002/sim.70172

simpledeff_cont_sizek = function(gee_corr, rho, gamma) {
  if (gee_corr == "independence") {
    deff = 1
  }
  else if (gee_corr == "exchangeable") {
    sum1 = gamma %>% imap(~ (1 + (.y-2)*rho)/((1 - rho)*(1 + (.y-1)*rho))*.x)
    deff = 1 / sum(unlist(sum1))
  }
  deff
}

# Expected DEFF and power under simple randomisation
# expected power = power of a Wald test for independent data for the specified treatment effect and effective sample size
# effective sample size = number of observations / expected deff under simple randomisation

calculate_simpledeffs_sizek = function(params) {
  gamma = calc_gamma(params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  M_k   = c(params$num_clus1, params$num_clus2, params$num_clus3, params$num_clus4)
  N     = sum(unlist(M_k %>% imap(~ .y*.x)))

  # independence gee
  ind_deff      = simpledeff_cont_sizek(gee_corr = "independence", rho = params$ICC, gamma = gamma)
  ind_N_eff     = N / ind_deff
  ind_power_exp = pwrss.t.2means(mu1 = params$b0, mu2 = params$b1, sd1 = 1, sd2 = 1, kappa = 1,
                                 n2 = ind_N_eff/2, alpha = 0.05, verbose=FALSE)$power
  ind_res       = data.frame(deff = ind_deff, N_eff = ind_N_eff, power_exp = ind_power_exp, model = "gee-ind")

  # exchangeable gee
  exch_deff      = simpledeff_cont_sizek(gee_corr = "exchangeable", rho = params$ICC, gamma = gamma)
  exch_N_eff     = N / exch_deff
  exch_power_exp = pwrss.t.2means(mu1 = params$b0, mu2 = params$b1, sd1 = 1, sd2 = 1, kappa = 1,
                                  n2 = exch_N_eff/2, alpha = 0.05, verbose=FALSE)$power
  exch_res       = data.frame(deff = exch_deff, N_eff = exch_N_eff, power_exp = exch_power_exp, model = "gee-exch")

  res = bind_rows(ind_res, exch_res)
  res
}

##### GENERATING DATA #####

generate_seed = function() {
  my_seed = as.integer(Sys.time()) %% 100000
  return(my_seed)
}

create_trial_nostrat_fixed_k4 = function(block_size, num_clus1, num_clus2, num_clus3, num_clus4) {
  
  k        = c(1, 2, 3, 4)                                  # cluster sizes
  M_k      = c(num_clus1, num_clus2, num_clus3, num_clus4)  # numbers of clusters of each size
  num_clus = sum(unlist(M_k))                               # total number of clusters
  N_k      = M_k %>% imap(~ .y*.x)                          # number of obs within each cluster size
  N        = sum(unlist(N_k))                               # total number of obs 
  b        = ceiling(N / block_size)                        # number of blocks needed (rounded up to to a full block)
  
  # create the clusters
  clusters = rep(k, times = M_k)
  dat = data.frame(cluster_id = seq(1:num_clus), clus_size = clusters, block_size = block_size) %>% 
    slice_sample(n = num_clus)        # randomly order the clusters
  
  # create the treatment indicator variable
  # randomisation is via randomly permuted blocks, but without stratification by cluster size
  tmt = unlist(replicate(n = b, expr={sample(rep(c(0,1), length = block_size))}, simplify = FALSE))
  # truncate the tmt list to N, to handle cases where N is not a multiple of block size
  tmt = tmt[1:N]
  
  # final dataframe
  dat = dat %>% uncount(clus_size, .remove = FALSE) %>%    # expand the clusters into observation level data
    mutate(tmt = tmt) %>%
    group_by(cluster_id) %>%
    mutate(withinclus_id = row_number(cluster_id),
           t = sum(tmt)) %>%
    ungroup() %>% 
    mutate(obs_id = seq(1:nrow(.)),
           block_size = block_size) %>% 
    select(obs_id, cluster_id, withinclus_id, clus_size, block_size, t, tmt)
  return(dat)
}

create_trial_nostrat_var_k4 = function(block_size, num_clus1, num_clus2, num_clus3, num_clus4) {
  
  k        = c(1, 2, 3, 4)                                  # cluster sizes
  M_k      = c(num_clus1, num_clus2, num_clus3, num_clus4)  # numbers of clusters of each size
  num_clus = sum(unlist(M_k))                               # total number of clusters
  N_k      = M_k %>% imap(~ .y*.x)                          # number of obs within each cluster size
  N        = sum(unlist(N_k))                               # total number of obs
  
  # create the clusters
  clusters = rep(k, times = M_k)
  dat = data.frame(cluster_id = seq(1:num_clus), clus_size = clusters) %>% 
    slice_sample(n = num_clus)        # randomly order the clusters
  # expand the clusters into observation level data
  dat = dat %>% uncount(clus_size, .remove = FALSE)
  
  # maximum number of blocks needed of each size.
  # by setting this the same for each block size, each block size will have equal probability of selection
  max_total = ceiling(N/min(unlist(block_size)))
  max_b = ceiling(max_total/length(block_size))
  
  # create a tmt indicator variable by looping through each block size
  tmt = data.frame()
  for (i in (1:length(block_size))) {
    b = block_size[i]
    
    t = unlist(replicate(n = max_b, expr={sample(rep(c(0,1), length = b))}, simplify = FALSE))
    block_id = rep(1:(length(t)/b), each = b) + max_b*(i-1)
    
    res = data.frame(block_id, b, t)
    tmt = bind_rows(tmt, res)
  }
  
  # randomly order the blocks
  tmt_blocks = split(tmt, tmt$block_id)
  tmt_shuffled = do.call(rbind, sample(tmt_blocks))
  rownames(tmt_shuffled) = NULL
  
  # truncate the tmt list to N, the required number of observations
  tmt_shuffled = tmt_shuffled[1:N,]
  
  # final dataframe
  dat = bind_cols(dat, tmt_shuffled) %>% 
    rename(tmt = t,
           blk_size = b) %>%
    group_by(cluster_id) %>%
    mutate(withinclus_id = row_number(cluster_id),
           t = sum(tmt)) %>%
    ungroup() %>% 
    mutate(obs_id = seq(1:nrow(.))) %>% 
    select(obs_id, cluster_id, withinclus_id, block_id, blk_size, clus_size, t, tmt)
  return(dat)
}

# Create a dataframe for a specified trial design
# Inputs: block size, number of clusters of sizes 1-4
# Outputs: a dataframe of observation and cluster ids and treatment allocation
create_trial_k4 = function(params) {
  
  # if using variable block sizes, convert the string of block sizes to a vector
    if (is.character(params$block_size)) {
    b = as.numeric(unlist(strsplit(params$block_size, ",")))
  } else {
    b = params$block_size
  }

  # set the create_trial function depending on block type (fixed or varying)
    if (length(b) == 1) {
      create_trial = function() { create_trial_nostrat_fixed_k4(block_size = b, num_clus1 = params$num_clus1,
                                                                num_clus2 = params$num_clus2,
                                                                num_clus3 = params$num_clus3,
                                                                num_clus4 = params$num_clus4) }
    } else if (length(b) > 1) {
      create_trial = function() { create_trial_nostrat_var_k4(block_size = b, num_clus1 = params$num_clus1,
                                                              num_clus2 = params$num_clus2,
                                                              num_clus3 = params$num_clus3,
                                                              num_clus4 = params$num_clus4) }
    }
  
  # create the dataset
  dat = create_trial()
  return(dat)
}

# Generate a single simulated dataset with a continuous outcome
# Inputs: simulation parameters
# Output: a dataframe
generate_cont_data_k4 = function(params) {
  
  # create dataframe of observations and the treatment variable
  dat = create_trial_k4(params)
  num_clus   = n_distinct(dat$cluster_id)
  clus_sizes = dat %>% summarise(clus_sizes=first(clus_size), .by=cluster_id) %>% select(clus_sizes)
  
  # random variables
  sd_alpha = sqrt(params$ICC)
  sd_eps   = sqrt(1 - params$ICC)
  alpha_j  = rnorm(num_clus, 0, sd_alpha)
  alpha_j  = rep(alpha_j, times = unlist(clus_sizes))
  eps_ij   = rnorm(nrow(dat), 0, sd_eps)
  dat      = data.frame(dat, alpha_j, eps_ij)
  
  # response variable
  dat = dat %>% mutate(Y = params$b0 + params$b1*tmt + alpha_j + eps_ij)
  
  # attach sim params
  dat = dat %>% 
    mutate(block_size = params$block_size,
           ICC        = params$ICC,
           trueb0     = params$b0,
           trueb1     = params$b1,
           num_clus1  = params$num_clus1,
           num_clus2  = params$num_clus2,
           num_clus3  = params$num_clus3,
           num_clus4  = params$num_clus4)
  return(dat)
}

# Generate multiple simulated datasets for a single scenario
# Inputs: simulation parameters for a single scenario
# Outputs: a list of datasets created by generate_data() and a matrix of random seeds for each dataset
generate_multiple_datasets_k4 = function(sims, params, outcome) {
  # initialise
  nsims = 0
  simulated_data = list()
  states = matrix(ncol=626, nrow=sims)
  # simulation loop
  while (nsims < sims) {
    states[nsims + 1, ] = .Random.seed
    simulated_data[[nsims + 1]] = generate_cont_data_k4(params)
    nsims = nsims + 1
  }
  simulated_data = Map(cbind, simulated_data, index=seq_along(simulated_data))
  list(data=simulated_data, states=states)
}

# Run simulation for one scenario
# Inputs: a dataframe of simulation parameters; number of simulated datasets to create per scenario; type of outcome
# Outputs: a list of lists of datasets created by generate_data()
# Data and random seeds are saved to file
simulate_k4 = function(params, nsims, outcome) {
  print('Simulating data with params:')
  print(params)
  # simulate the datasets
  result = generate_multiple_datasets_k4(nsims, params, outcome)
  data = result$data
  data = Map(cbind, scenario_id=params$scenario_id, data)
  # save the results
  saveRDS(data, file=paste0("simdata/data_", params$scenario_id, ".rds"))
  saveRDS(result$states, file=paste0("simdata/states_", params$scenario_id, ".rds"))
  return(data)
  #  }
}

##### ANALYSING DATA #####

# Fit models for partially clustered data with a continuous outcome to a dataset
# Inputs: a dataset as produced by generate_data()
# Output: a dataframe of simulation parameters and model estimates, one row per analysis
# Analyses: geeglm (ind and exch correlation structures)
fit_cont_models_k4 = function(dat) {
  # logging
  if ((dat$index[1] == 1) | (dat$index[1] %% 50 == 0)) {
    print(dat[1,c("scenario_id", "index")]) 
    print(Sys.time())
  }
  # sim params
  d = data.frame(scenario_id = dat$scenario_id[1],
                 dataset     = dat$index[1],
                 block_size  = dat$block_size[1],
                 N           = nrow(dat),
                 M1          = dat$num_clus1[1],
                 M2          = dat$num_clus2[1],
                 M3          = dat$num_clus3[1],
                 M4          = dat$num_clus4[1],
                 ICC         = dat$ICC[1],
                 trueb0      = dat$trueb0[1],
                 trueb1      = dat$trueb1[1],
                 stringsAsFactors = FALSE, row.names = NULL)
  
  # independence gee via geepack::geeglm
  m1 = geeglm(Y ~ tmt, id=cluster_id, family=gaussian, corstr="independence", data=dat)
  s1 = summary(m1)
  d1 = data.frame(model  = "geeglm-ind",
                  b0     = s1$coefficients["(Intercept)", "Estimate"],
                  se_b0  = s1$coefficients["(Intercept)", "Std.err"],
                  b1     = s1$coefficients["tmt", "Estimate"],
                  se_b1  = s1$coefficients["tmt", "Std.err"],
                  var_b1 = m1$geese$vbeta[2,2],
                  p_b1   = s1$coefficients["tmt", "Pr(>|W|)"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d1 = bind_cols(d, d1)
  
  # exchangeable gee via geepack::geeglm
  m2 = geeglm(Y ~ tmt, id=cluster_id, family=gaussian, corstr="exchangeable", data=dat)
  s2 = summary(m2)
  d2 = data.frame(model    = "geeglm-exch",
                  b0       = s2$coefficients["(Intercept)", "Estimate"],
                  se_b0    = s2$coefficients["(Intercept)", "Std.err"],
                  b1       = s2$coefficients["tmt", "Estimate"],
                  se_b1    = s2$coefficients["tmt", "Std.err"],
                  var_b1   = m2$geese$vbeta[2,2],
                  p_b1     = s2$coefficients["tmt", "Pr(>|W|)"],
                  alpha    = s2$corr["alpha", "Estimate"],
                  alpha_se = s2$corr["alpha", "Std.err"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d2 = bind_cols(d, d2)
  
  # independence gee with effect for cluster size via geepack::geeglm
  m3 = geeglm(Y ~ tmt + clus_size, id=cluster_id, family=gaussian, corstr="independence", data=dat)
  s3 = summary(m3)
  d3 = data.frame(model  = "geeglm-ind-adj",
                  b0     = s3$coefficients["(Intercept)", "Estimate"],
                  se_b0  = s3$coefficients["(Intercept)", "Std.err"],
                  b1     = s3$coefficients["tmt", "Estimate"],
                  se_b1  = s3$coefficients["tmt", "Std.err"],
                  var_b1 = m3$geese$vbeta[2,2],
                  p_b1   = s3$coefficients["tmt", "Pr(>|W|)"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d3 = bind_cols(d, d3)

  # exchangeable gee with effect for cluster size via geepack::geeglm
  m4 = geeglm(Y ~ tmt + clus_size, id=cluster_id, family=gaussian, corstr="exchangeable", data=dat)
  s4 = summary(m4)
  d4 = data.frame(model  = "geeglm-exch-adj",
                  b0     = s4$coefficients["(Intercept)", "Estimate"],
                  se_b0  = s4$coefficients["(Intercept)", "Std.err"],
                  b1     = s4$coefficients["tmt", "Estimate"],
                  se_b1  = s4$coefficients["tmt", "Std.err"],
                  var_b1 = m4$geese$vbeta[2,2],
                  p_b1   = s4$coefficients["tmt", "Pr(>|W|)"],
                  alpha    = s4$corr["alpha", "Estimate"],
                  alpha_se = s4$corr["alpha", "Std.err"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d4 = bind_cols(d, d4)
  
  # linear regression assuming independence
  m5 = lm(Y ~ tmt, data=dat)
  s5 = summary(m5)
  d5 = data.frame(model    = "linreg",
                  b0       = s5$coefficients["(Intercept)", "Estimate"],
                  se_b0    = s5$coefficients["(Intercept)", "Std. Error"],
                  b1       = s5$coefficients["tmt", "Estimate"],
                  se_b1    = s5$coefficients["tmt", "Std. Error"],
                  var_b1   = vcov(m5)["tmt", "tmt"],
                  sd_resid = sd(s5$residuals),
                  p_b1     = s5$coefficients["tmt", "Pr(>|t|)"],
                  stringsAsFactors = FALSE, row.names = NULL)
  d5 = bind_cols(d, d5)
  
  result = bind_rows(d1, d2, d3, d4, d5)
  return(result)
}

##### REPORTING RESULTS #####

# Selects the results of the models adjusted or unadjusted for cluster size
# Only include partially (not fully) clustered datasets for adjusted analyses;
# ie, a non-zero number of clusters of both sizes

select_adjustment = function(dat, adjusted = FALSE) {
  if (adjusted == TRUE) {
    dat = dat %>% filter(adjusted == "adj" | adjusted == "linreg" | adjusted == "glm") %>% 
      filter((M1 > 0) & (M2 > 0))
  }
  else if (adjusted == FALSE) {
    dat = dat %>% filter(adjusted == "unadj" | adjusted == "linreg" | adjusted == "glm")
  }
  return(dat)
}

# Summarise the model performance outcomes for each scenario
# Inputs: model_results dataframe
# Outputs: table of non-convergence and non-positive definite rates and estimated ICC for each scenario
create_modelling_result_table = function (model_results) {
  
  # calculate performance measures
  model_res = model_results %>%
    filter(model == "geeglm-ind" | model == "geeglm-exch") %>%
    group_by(scenario_id, model) %>%
    summarise(conv    = sum(.dropbig)/sims*100,
              npd     = sum(abs(alpha) > 1)/sims*100,
              estICC  = median(alpha),
              .groups = "drop")
  
  # restructure for presentation table
  model_res_wide = model_res %>%
    pivot_wider(names_from = model, values_from = c(conv, npd, estICC)) %>%
    rename(nonconv_ind  = "conv_geeglm-ind",
           nonconv_exch = "conv_geeglm-exch",
           npd_exch  = "npd_geeglm-exch",
           estICC    = "estICC_geeglm-exch") %>%
    select(-c("npd_geeglm-ind", "estICC_geeglm-ind")) %>%
    left_join(parameters, by="scenario_id") %>% 
    relocate(block_size, .after = scenario_id)

  model_table = flextable(model_res_wide[,1:6]) %>%
    colformat_double(j = c("estICC"), digits = 2) %>%
    fontsize(size=9)
  return(model_table)
}

# Restructure simulation results for final presentation table
# Input: summary table of simulation results as created by create_summary_table()
# Output: dataframe of 1 row per scenario of deff and power results
convert_wide = function(res){
  res_wide = res %>%
    pivot_wider(names_from = model, values_from = c(deff_exp, deff_obs_med, deff_reldiff_med,
                                                    pwr_exp, pwr_obs, pwr_diff)) %>% 
    rename(deff_obs_ind      = "deff_obs_med_geeglm-ind",
           deff_exp_ind      = "deff_exp_geeglm-ind",
           deff_reldiff_ind  = "deff_reldiff_med_geeglm-ind",
           deff_obs_exch     = "deff_obs_med_geeglm-exch",
           deff_exp_exch     = "deff_exp_geeglm-exch",
           deff_reldiff_exch = "deff_reldiff_med_geeglm-exch",
           pwr_obs_ind       = "pwr_obs_geeglm-ind",
           pwr_exp_ind       = "pwr_exp_geeglm-ind",
           pwr_diff_ind      = "pwr_diff_geeglm-ind",
           pwr_obs_exch      = "pwr_obs_geeglm-exch",
           pwr_exp_exch      = "pwr_exp_geeglm-exch",
           pwr_diff_exch     = "pwr_diff_geeglm-exch",) %>% 
    arrange(block_size, N, scenario_id) %>% 
    relocate(block_size, N, ICC,
             deff_obs_ind, deff_exp_ind, deff_reldiff_ind,
             deff_obs_exch, deff_exp_exch, deff_reldiff_exch,
             pwr_obs_ind, pwr_exp_ind, pwr_diff_ind,
             pwr_obs_exch, pwr_exp_exch, pwr_diff_exch)
  return(res_wide)
}

# Creates a table of observed and expected design effect results
deff_table_k4 = function(dat, numsim = numsim, caption = caption) {
  
  dat = dat %>%
    mutate(ppairs = round(1 - calc_gamma_s(num_clus1 = M1, num_clus2 = M2), digits =2)) %>% 
    select(c(ppairs, block_size, N, ICC, deff_obs_ind, deff_exp_ind, deff_reldiff_ind,
             deff_obs_exch, deff_exp_exch, deff_reldiff_exch)) %>% 
    arrange(ppairs, block_size, N)
  
  tab = flextable(dat) %>% 
    set_caption(caption) %>%
    set_header_labels(ppairs            = "Propn of paired obs",
                      block_size        = "Block size",
                      N                 = "Sample size",
                      ICC               = "ICC",
                      deff_obs_ind      = "Observed DEFF*",
                      deff_exp_ind      = "Expected DEFF^",
                      deff_reldiff_ind  = "Relative difference (%)*",
                      deff_obs_exch     = "Observed DEFF*",
                      deff_exp_exch     = "Expected DEFF^",
                      deff_reldiff_exch = "Relative difference (%)*") %>% 
    add_header_row(colwidths = c(4, 3, 3), values = c("", "GEE independence", "GEE exchangeable")) %>% 
    align(align = "center", part = "header") %>%
    vline(j = c(4, 7)) %>%
    line_spacing(space=0.5) %>% 
    add_footer_lines(paste0("*Median value across ", numsim, " simulated datasets")) %>% 
    add_footer_lines(paste0("^Based on the formula for simple randomisation")) %>% 
    add_footer_lines("GEE = generalised estimating equation, DEFF = design effect, ICC = intracluster correlation coefficient") %>% 
    fontsize(size=9) %>%
    fontsize(part = "header", size = 10) %>% 
    fontsize(part="footer", size=9) %>% 
    colformat_double(j = c("deff_obs_ind", "deff_exp_ind", "deff_reldiff_ind",
                           "deff_obs_exch", "deff_exp_exch", "deff_reldiff_exch"), digits = 2)
  return(tab)
}

# Creates a table of observed and expected power results
power_table_k4 = function(dat, caption = caption){
  
  dat = dat %>%
    mutate(ppairs = round(1 - calc_gamma_s(num_clus1 = M1, num_clus2 = M2), digits =2)) %>% 
    select(c(ppairs, block_size, N, ICC, pwr_obs_ind, pwr_exp_ind, pwr_diff_ind,
             pwr_obs_exch, pwr_exp_exch, pwr_diff_exch)) %>% 
    arrange(ppairs, block_size, N)
  
  tab = flextable(dat) %>% 
    set_caption(caption) %>%
    set_header_labels(ppairs        = "Propn of paired obs",
                      block_size    = "Block size",
                      nobs          = "Sample size",
                      ICC           = "ICC",
                      pwr_obs_ind   = "Observed power",
                      pwr_exp_ind   = "Expected power^",
                      pwr_diff_ind  = "Difference",
                      pwr_obs_exch  = "Observed power",
                      pwr_exp_exch  = "Expected power^",
                      pwr_diff_exch = "Difference") %>% 
    add_header_row(colwidths = c(4, 3, 3), values = c("", "GEE independence", "GEE exchangeable")) %>% 
    align(align = "center", part = "header") %>%
    vline(j = c(4, 7)) %>%
    line_spacing(space=0.5) %>% 
    add_footer_lines(paste0("^Based on the formula for simple randomisation")) %>% 
    add_footer_lines("GEE = generalised estimating equation, ICC = intracluster correlation coefficient") %>% 
    fontsize(size=9) %>%
    fontsize(part = "header", size = 10) %>% 
    fontsize(part="footer", size=9) %>% 
    colformat_double(j = c("pwr_obs_ind", "pwr_exp_ind", "pwr_diff_ind",
                           "pwr_obs_exch", "pwr_exp_exch", "pwr_diff_exch"), digits = 2)
  return(tab)
}
