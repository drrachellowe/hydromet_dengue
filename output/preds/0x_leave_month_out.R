# !/usr/bin/env Rscript

yyear = as.integer(commandArgs(trailingOnly=TRUE)[1])
mmonth = as.integer(commandArgs(trailingOnly=TRUE)[2])
print(c(yyear, mmonth))
inla.setOption(num.threads = "4:1")


# Script to run INLA models in cross validation prediction mode

# Step 0: load packages pre-processed data and functions
# Step 1: rerun the selected model (fitted with config = TRUE for sampling) 
# Step 2: produce cross-validated posterior predictive samples leaving out one year and one month at a time

# Step 0: load packages and pre-processed data
source("00_load_packages_data.R")

# Step 1: rerun the selected model (fitted with config = TRUE for sampling) 

# fit full final model to use with control.mode = list(result = model1, restart = TRUE)
formula <- Y ~ 1 + 
  f(T1, replicate = S2, model = "rw1", scale.model = TRUE, cyclic = TRUE, 
    constr = TRUE, hyper = precision.prior) +
  f(S1, model = "bym2", replicate = T2, graph = "output/map.graph", 
    scale.model = TRUE, hyper = precision.prior) +
  basis_tmin + basis_pdsi + urban_basis1_pdsi + Vu


if (FALSE) {
    model1 <- inla(formula, data = df, family = "nbinomial", offset = log(E), 
                   control.inla = list(strategy = 'adaptive'), 
                   control.compute = list(dic = TRUE, config = TRUE, 
                                          cpo = TRUE, return.marginals = FALSE),
                   control.fixed = list(correlation.matrix = TRUE, 
                                        prec.intercept = 1, prec = 1),
                   control.predictor = list(link = 1, compute = TRUE), 
                   verbose = FALSE)
    model1 <- inla.rerun(model1)
    save(model1, file = "output/model1_config.RData")
} else {
    load(file = "output/model1_config.RData")
    model1$misc$configs <- NULL
    gc()
}

# Step 2: produce cross-validated posterior predictive samples leaving out one year and one month at a time

# define number of samples
s <- 1000

if (TRUE) {
    # replace dengue data in testing period with NA for out of sample prediction
    casestopred <- data$dengue_cases # response variable
    idx.pred <- which(data$year_index == yyear & data$month == mmonth)
    casestopred[idx.pred] <- NA # replace cases in year and month of interest to NA
    mpred <- length(idx.pred)
    
    # set response variable and year indicator
    df$Y <- casestopred
    
    # final model
    formula <- Y ~ 1 + 
        f(T1, replicate = S2, model = "rw1", cyclic = TRUE, constr = TRUE, 
          scale.model = TRUE, hyper = precision.prior) +
        f(S1, model = "bym2", replicate = T2, graph = "output/map.graph", 
          scale.model = TRUE, hyper = precision.prior) +
        basis_tmin + basis_pdsi + urban_basis1_pdsi + Vu
    
    if (TRUE) {
        model <- inla(formula, data = df, family = "nbinomial", offset = log(E), 
                      control.inla = list(strategy = 'adaptive'), 
                      control.compute = list(dic = TRUE, config = TRUE, 
                                             cpo = TRUE, return.marginals = FALSE),
                      control.fixed = list(correlation.matrix = TRUE, 
                                           prec.intercept = 1, prec = 1),
                      control.predictor = list(link = 1, compute = TRUE), 
                      control.mode = list(result = model1, restart = TRUE),
                      verbose = TRUE)
        model <- inla.rerun(model)
    } else {
        model <- model1
    }
    
    xx <- inla.posterior.sample(s, model)
    xx.s <- inla.posterior.sample.eval(function(...) c(theta[1], Predictor[idx.pred]), xx)
    y.pred <- matrix(NA, mpred, s)
    for(s.idx in 1:s) {
        xx.sample <- xx.s[, s.idx]
        y.pred[, s.idx] <- rnbinom(mpred, mu = exp(xx.sample[-1]), size = xx.sample[1])
    }
    preds <- list(year = 2000 + yyear, month = mmonth, idx.pred = idx.pred, 
                  mean = apply(y.pred, 1, mean), median = apply(y.pred, 1, median),
                  lci = apply(y.pred, 1, quantile, probs = c(0.025)),
                  uci = apply(y.pred, 1, quantile, probs = c(0.975)))
    save(preds, file = paste0("output/preds_",2000 + yyear, "_", mmonth, ".RData"))
}    
