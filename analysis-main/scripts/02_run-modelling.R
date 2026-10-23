############################
# ---- Modelling Stage ----
############################

run_modelling <- function(cfg) {

  # ----------------------------
  # 00 Space: MPL, Flags, Drift, Samples
  # ----------------------------
  # These four steps form one space because they are chained: the
  # Holt-Laury (MPL) model yields the posterior draws of r_i; the
  # expected-utility optimal stake a*(r) computed from those draws
  # defines the participant flags (the confirmatory sample); the
  # block-drift tests are run on that sample; and the analysis samples
  # collect, in one table, who is in which analysis. Each later step is
  # stored inside the folder of the MPL model it was built from, so a
  # change in the MPL fit invalidates flags, drift decisions and samples
  # together.
  
  

  # primary preregistered model specification
  source(here::here("scripts", "modelling", "01_mpl-stan-bernoulli.R"))
  mpl_stan_bernoulli(cfg)   
  
  # alternative: hierarchical (partially pooled) Holt-Laury model, robustness for the full samples
  source(here::here("scripts", "modelling", "02_mpl-stan-hierarchical.R"))
  mpl_stan_hierarchical(cfg)
  
  # draw a* for the selected model, produce participant flags
  source(here::here("scripts", "modelling", "07_a-star-draws.R"))
  a_star_draws(cfg)             
  
  # drift decisions on the full, consistent and confirmatory samples
  source(here::here("scripts", "modelling", "08_drift-check.R"))
  drift_check(cfg)              

  # analysis samples: who is in which analysis (one table, read by all models)
  source(here::here("scripts", "modelling", "09_analysis-sample.R"))
  analysis_sample(cfg)


  # ----------------------------
  # 10 Space: RQ1
  # ----------------------------
 
  # primary preregistered model of betting, full and confirmatory samples
  source(here::here("scripts", "modelling", "11_rq1-stan-bernoulli.R"))
  rq1_stan_bernoulli(cfg)   
            

  # ----------------------------
  # 20 Space: RQ2 
  # ----------------------------
  
  # primary preregistered model of staking, Gaussian model of standardized stake deviations
  source(here::here("scripts", "modelling", "21_rq2-stan-gaussian.R"))
  rq2_stan_gaussian(cfg)        
  
  # preregistered alternative, one-inflated Beta-Binomial
  source(here::here("scripts", "modelling", "22_rq2-stan-one-inflated.R"))
  rq2_stan_one_inflated(cfg)    
  
  source(here::here("scripts", "modelling", "23_rq2-stan-one-inflated-seq.R"))
  rq2_stan_one_inflated_seq(cfg)

   
   # ----------------------------
   # 30 Space: RQ3 
   # ----------------------------
   source(here::here("scripts", "modelling", "31_rq3-stan-hurdle-gamma.R"))
   rq3_stan_hurdle_gamma(cfg)
   
   source(here::here("scripts", "modelling", "32_rq3-stan-gamma.R"))
   rq3_stan_gamma(cfg)
   
   source(here::here("scripts", "modelling", "33_rq3-stan-gaussian.R"))
   rq3_stan_gaussian(cfg)

   # sensitivity of RQ3 to the wealth floor (refits the selected RQ3 model at each floor)
   source(here::here("scripts", "modelling", "35_rq3-stan-floor-robustness.R"))
   rq3_stan_floor_robustness(cfg)

   
   # ----------------------------
   # 40 Space: RQ4
   # ----------------------------
   source(here::here("scripts", "modelling", "41_rq4-stan-bernoulli.R"))
   rq4_stan_bernoulli(cfg)

   
   # ----------------------------
   # 50 Space: EX1
   # ----------------------------

   source(here::here("scripts", "modelling", "51_ex1-chi-draws.R"))
   ex1_chi_draws(cfg)
   
   source(here::here("scripts", "modelling", "54_ex1-stan-regression.R"))
   ex1_stan_regression(cfg)

   
   # ----------------------------
   # 60 Space: EX2
   # ----------------------------
   source(here::here("scripts", "modelling", "61_ex2-stan-associations.R"))
   ex2_stan_associations(cfg)

   # ----------------------------
   # 90 Space: Synthesis (joint model of betting and side choice)
   # ----------------------------
   source(here::here("scripts", "modelling", "91_synthesis-stan-joint.R"))
   synthesis_stan_joint(cfg)

   



  
  msg("\nMODELLING completed.\n")
}