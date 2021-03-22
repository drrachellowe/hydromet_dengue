# Combined effects of hydrometeorological hazards and urbanisation on dengue risk in Brazil: a spatiotemporal modelling study

# Rachel Lowe (2021)
# https://github.com/drrachellowe/hydromet_dengue

# R script to run INLA models of increasing complexity
# WARNING: the script may take over a day to run

# Step 0: load packages and pre-processed data 
# Step 1: formulate a baseline model including spatiotemporal random effects and test different combinations of DLNM climate indicators
# Step 2: using best fitting model in Step 1 to test interactions between PDSI-DLNM and 
# 1) % urban population from highly urbanised to more rural
# 2) frequency of reported water supply shortages from high to low frequency
# centred at high (upper quartile), medium (median) and low (lower quartile)

# Step 0: load packages and pre-processed data
source("00_load_packages_data.R")

# run models of increasing complexity in INLA

# Step 1: fit a baseline model including spatiotemporal random effects

# formulate a base model including: 
# state-specific monthly random effects to account for variation in seasonality between states (random walk cyclic prior)
# year-specific spatial random effects to account for interannual variation in spatial overdisperson and dependency structures (modified Besag-York-Mollie prior bym2)

# baseline model
baseformula <- Y ~ 1 + f(T1, replicate = S2, model = "rw1", cyclic = TRUE, constr = TRUE,
                     scale.model = TRUE,  hyper = precision.prior) +
  f(S1, model = "bym2", replicate = T2, graph = "output/map.graph", 
    scale.model = TRUE, hyper = precision.prior) 

# test baseline model with Poisson distribution
# model <- mymodel(baseformula, family = "poisson")
# model$dic$dic

# define formulas by updating the baseline formula with different combinations of Tmin, Tmax and PDSI cross-basis functions
formula0.1 <- update.formula(baseformula, ~. + basis_tmin)
formula0.2 <- update.formula(baseformula, ~. + basis_tmax)
formula0.3 <- update.formula(baseformula, ~. + basis_pdsi)
formula0.4 <- update.formula(baseformula, ~. + basis_tmin + basis_pdsi)
formula0.5 <- update.formula(baseformula, ~. + basis_tmax + basis_pdsi)

# create a list of formulas
formulas <- list(baseformula, formula0.1, formula0.2, formula0.3, formula0.4, formula0.5)
# create model label string
lab <- c("basemodel", "model0.1", "model0.2", "model0.3", "model0.4", "model0.5")

# create a function to run a model for each formula in the list and save the model output to file
# WARNING: this may take a long time to run
models <- lapply(1:length(formulas), 
              function(i) {
                model <- mymodel(formulas[[i]], df)
                save(model, file = paste0("output/", lab[i],".RData"))})

# create table to store DIC and select best model 
table0 <- data.table(Model  = c("base", "tmin", "tmax", "pdsi", "tmin + pdsi", "tmax + pdsi"), 
                     DIC = NA)

for(i in 1:length(formulas))
  {
  load(paste0("output/",lab[i],".RData"))
  table0$DIC[i] <- round(model$dic$dic, 0)
}

# view table
table0

# define position of best fitting model
best.fit <- which.min(table0$DIC)

# Step 3: use best fitting model in Step 2 to test interactions betweeen 
# PDSI-DLNM and % population living in urban areas and PDSI-DLNM and water supply shortage frequency
# centred on high (model1.1), medium (model1.2) and low (model1.3) levels of urbanisation and
# centred on high (model1.4), medium (model1.5) and low (model1.6) frequency of water supply shortages

# assign formula for best fitting model to the new baseformula
# baseformula <- formulas[[best.fit]]

# redefine baseformula as best fitting model from above
baseformula <- Y ~ 1 + f(T1, replicate = S2, model = "rw1", cyclic = TRUE, constr = TRUE,
                         scale.model = TRUE,  hyper = precision.prior) +
  f(S1, model = "bym2", replicate = T2, graph = "output/map.graph", 
    scale.model = TRUE, hyper = precision.prior) + basis_tmin + basis_pdsi

# define formulas by updating the best fitting model0 formula with interactions between pdsi cross-basis and socio-economic indicators

# best fitting model0 + interaction between pdsi cross-basis and urbanisation
formula1.1 <- update.formula(baseformula, ~. + urban_basis1_pdsi + Vu)
formula1.2 <- update.formula(baseformula, ~. + urban_basis2_pdsi + Vu)
formula1.3 <- update.formula(baseformula, ~. + urban_basis3_pdsi + Vu)

# best fitting model0 + interaction between pdsi cross-basis and water supply shortage
formula1.4 <- update.formula(baseformula, ~. + water_basis1_pdsi + Vw)
formula1.5 <- update.formula(baseformula, ~. + water_basis2_pdsi + Vw)
formula1.6 <- update.formula(baseformula, ~. + water_basis3_pdsi + Vw)

# create a list of formulas
formulas <- list(formula1.1, formula1.2, formula1.3, formula1.4, formula1.5, formula1.6)

# create model label string
lab <- c("model1.1_urban", "model1.2_urban", "model1.3_urban",
         "model1.1_water", "model1.2_water", "model1.3_water")

# create a function to run a model for each formula in the list and save the model output to file
# WARNING: this may take a long time to run
models <- lapply(1:length(formulas), 
                 function(i) {
                   model <- mymodel(formulas[[i]], df)
                   save(model, file = paste0("output/", lab[i],".RData"))})

# create table to store DIC
table1 <- data.frame(Model  = c("high urban", "intermediate urban", "low urban",
                                "high water shortage", "intermediate water shortage", "low water shortage"), 
                     DIC = NA)

for(i in 1:length(formulas))
{
  load(paste0("output/",lab[i],".RData"))
  table1$DIC[i] <- round(model$dic$dic, 0)
}

# view table 
# note: high, intermediate and low urban (water) models are different parametrisations of the same urban (water) interaction model 
table1




