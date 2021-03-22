# Combined effects of hydrometeorological hazards and urbanisation on dengue risk in Brazil: a spatiotemporal modelling study

# Rachel Lowe (2021)
# https://github.com/drrachellowe/hydromet_dengue

# R script to prepare data and lagged variables for INLA-DLNM modelling

# install INLA
# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

# load INLA
library(INLA)

#  select other packages
packages <- c("data.table", "tidyverse", "sf", "sp", "spdep",
              "dlnm", "tsModel", "hydroGOF","RColorBrewer", 
              "geofacet", "ggpubr", "ggthemes")

# install.packages
# lapply(packages, install.packages, character.only = TRUE)

# load packages
lapply(packages, library, character.only = TRUE)

# load shape file for Brazil
map <- read_sf("data/shape_brazil.shp")
# dim(map)

# Create adjacency matrix
nb.map <- poly2nb(as_Spatial(map$geometry))
g.file <- "output/map.graph"
if (!file.exists(g.file)) nb2INLA(g.file, nb.map)

# load pre-defined grid of Brazilian states for geofacet plots
# note: could use pre-loaded grid = "br_states_grid1" in geofacet package, would need match state names
grid <- read.csv("data/br_states_grid.csv")
# head(grid)

# load data
# note dengue data available from Jan 2001 
# Climate data included for 2000 to obtained lagged values prior to 2001
data <- fread("data/data_2000_2019.csv", header = T)
# head(data)

# summary climate variables
summary(data$tmin)
summary(data$tmax)
summary(data$pdsi)

# Note, Palmer Drought Severity Index data not available for island Fernando de Noronha (code 26019)

# Create lagged variables
# define matrices of lagged terms for monthly mean climate variables

# set maximum lag
nlag = 6

# Minimum temperature (Tmin)
lag_tmin <- tsModel::Lag(data$tmin, group = data$micro_code, k = 0:nlag)
# Maximum temperature (Tmax)
lag_tmax <- tsModel::Lag(data$tmax, group = data$micro_code, k = 0:nlag)
# Palmer drought severity index (PDSI)
lag_pdsi <- tsModel::Lag(data$pdsi, group = data$micro_code, k = 0:nlag)

# Remove year 2000 from lagged climate variables
lag_tmin <- lag_tmin[data$year > 2000,]
lag_tmax <- lag_tmax[data$year > 2000,]
lag_pdsi <- lag_pdsi[data$year > 2000,]

# remove year 2000 from dengue dataframe
data <- data[data$year > 2000,]

# define dimensions

# re-define time indicator to set 1 to Jan 2001
data$time <- data$time - 12
# total number of months
ntime <- length(unique(data$time))
# total number of years
nyear <- length(unique(data$year))
# total number of microregions
nmicro <- length(unique(data$micro_code))
# total number of states
nstate <- length(unique(data$state_code))

# define cross-basis matrix (combining nonlinear exposure and lag functions)
# set lag knots
lagknot = equalknots(0:nlag, 2)

# Tmin
var <- lag_tmin
basis_tmin <- crossbasis(var,
                    argvar = list(fun = "ns", knots = equalknots(data$tmin, 2)),
                    arglag = list(fun = "ns", knots = nlag/2))

# Tmax
var <- lag_tmax
basis_tmax <- crossbasis(var,
                         argvar = list(fun = "ns", knots = equalknots(data$tmax, 2)),
                         arglag = list(fun = "ns", knots = nlag/2))

# PDSI
var <- lag_pdsi
basis_pdsi <- crossbasis(var,
                         argvar = list(fun="ns", knots = equalknots(data$pdsi, 2)),
                         arglag = list(fun="ns", knots = lagknot))

# test linear interaction with % residents living in urban areas
# centre the urban variable at different levels of urbanisation (25th, 50th and 75th percentiles)
# from highly urbanised to more rural

summary(data$urban)
# set indicator to zero at point of interest (centring point)
# re-parameterise model to extract different predictions
urban_ind1 <- data$urban - quantile(data$urban, p = 0.75) # highly urbanised 
urban_ind2 <- data$urban - quantile(data$urban, p = 0.5) # intermediate
urban_ind3 <- data$urban - quantile(data$urban, p = 0.25) # more rural

# test linear interaction with frequency of water shortages between 2000-2016
# centre the water shortage variable at different levels of frequency of shortages (25th, 50th and 75th percentiles)
# from high to low frequency

summary(data$water_shortage)
# set indicator to zero at point of interest (centring point)
# re-parameterise model to extract different predictions
water_ind1 <- data$water_shortage - quantile(data$water_shortage, p = 0.75) # high frequency shortages
water_ind2 <- data$water_shortage - quantile(data$water_shortage, p = 0.5) # intermediate
water_ind3 <- data$water_shortage - quantile(data$water_shortage, p = 0.25) # low frequency shortages

# Multiply each cross-basis variable by the linear terms (see Gasparrini et al. EHP 2015)
# note: exploit the column by column product 

# multiply the PDSI cross-basis variables by the urban linear terms
urban_basis1_pdsi <- basis_pdsi*urban_ind1
urban_basis2_pdsi <- basis_pdsi*urban_ind2
urban_basis3_pdsi <- basis_pdsi*urban_ind3

# multiply the PDSI cross-basis variables by the water shortage linear terms
water_basis1_pdsi <- basis_pdsi*water_ind1
water_basis2_pdsi <- basis_pdsi*water_ind2
water_basis3_pdsi <- basis_pdsi*water_ind3

# assign unique column names to cross-basis matrix for inla() model
# note: not necessary for glm(), gam() or glm.nb() models
colnames(basis_tmin) = paste0("basis_tmin.", colnames(basis_tmin))
colnames(basis_tmax) = paste0("basis_tmax.", colnames(basis_tmax))
colnames(basis_pdsi) = paste0("basis_pdsi.", colnames(basis_pdsi))

colnames(urban_basis1_pdsi) = paste0("urban_basis1_pdsi.", colnames(urban_basis1_pdsi))
colnames(urban_basis2_pdsi) = paste0("urban_basis2_pdsi.", colnames(urban_basis2_pdsi))
colnames(urban_basis3_pdsi) = paste0("urban_basis3_pdsi.", colnames(urban_basis3_pdsi))

colnames(water_basis1_pdsi) = paste0("water_basis1_pdsi.", colnames(water_basis1_pdsi))
colnames(water_basis2_pdsi) = paste0("water_basis2_pdsi.", colnames(water_basis2_pdsi))
colnames(water_basis3_pdsi) = paste0("water_basis3_pdsi.", colnames(water_basis3_pdsi))

# create indices for INLA models
# note: for INLA models an index should start with 1 and with the max value equal to the length of unique values

# create microregion index 
data$micro_index <- rep(1:nmicro, ntime)

# create state index
# state length
k <- unique(data$state_code)

for (j in 1:nstate)
{
  data$state_index[data$state_code == k[j]] <- j 
}

# create year index
# set first year (in this case 2001) to 1
data$year_index <- data$year - 2000 

# set up data and priors for INLA model

# set data for models
Y  <- data$dengue_cases # response variable
N  <- length(Y) # total number of data points
E  <- data$population/10^5 # model offset so that response is equivalent to an incidence rate per 100,000 people
T1 <- data$month # for random effect to account for annual cycle (seasonality)
T2 <- data$year_index # for random effect to account for inter-annual variability
S1 <- data$micro_index # for microregion spatial random effect
S2 <- data$state_index # for state interaction with month random effect
Vu <- data$urban # include level of urbanisation (% pop living in urban areas) variable along with linear urban interaction
Vw <- data$water_shortage # include frequency of water shortages along with linear water shortage interaction

# create dataframe for model testing
df <- data.frame(Y, E, T1, T2, S1, S2, Vu, Vw)

# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

# inla model function

# include formula and set defaults for data, family (to allow other prob dist models e.g. Poisson) and config (to allow for sampling)
mymodel <- function(formula, data = df, family = "nbinomial", config = FALSE)

  {
  model <- inla(formula = formula, data = data, family = family, offset = log(E),
       control.inla = list(strategy = 'adaptive'), 
       control.compute = list(dic = TRUE, config = config, 
                              cpo = TRUE, return.marginals = FALSE),
       control.fixed = list(correlation.matrix = TRUE, 
                            prec.intercept = 1, prec = 1),
       control.predictor = list(link = 1, compute = TRUE), 
       verbose = FALSE)
  model <- inla.rerun(model)
  return(model)
}
