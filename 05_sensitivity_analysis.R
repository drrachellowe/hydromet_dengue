# Combined effects of hydrometeorological hazards and urbanisation on dengue risk in Brazil: a spatiotemporal modelling study

# Rachel Lowe (2021)
# https://github.com/drrachellowe/hydromet_dengue

# R script to perform sensitivity analysis, testing exposure-lag-response associations using different variables
# including Tmax instead of Tmin and PDSI DLNM interaction with continuous water shortage frequency instead of urbanisation variable.

# Step 0: load packages and pre-processed data
# Step 1: load models and plot overall response across all lags for Tmin and Tmax 
# Step 2: load models and plot exposure lag response associations for PDSI at high and low frequency of water supply shortages
# Step 3: plot lag response associations for extreme values of PDSI given high and low frequency of water supply shortages

# Step 0: load packages and pre-processed data
source("00_load_packages_data.R")

# Step 1: load models and plot overall response across all lags for Tmin and Tmax 

# plot overall Tmin and Tmax response across all lags (Appendix Fig A11)
pdf("figs/fig_S11_Tmin_Tmax_overall.pdf", width = 12, height = 6)
par(mfrow = c(1,2))

# Tmin
# load model Tmin DLNM + PDSI DLNM
load("output/model0.4.RData")

# extract full coef and vcov and create indicators for each term
coef <- model$summary.fixed$mean
vcov <- model$misc$lincomb.derived.covariance.matrix

# create indicators for the terms associated with Tmin crossbasis
indt <- grep("basis_tmin", model$names.fixed)

# extract predictions from the Tmin DLNM centred on overall mean Tmin (19 deg C)
predt <- crosspred(basis_tmin, coef = coef[indt], vcov=vcov[indt,indt],
                   model.link = "log", bylag = 0.25, cen = round(mean(data$tmin), 0))

# define x values (lag, by lag)
lagbylag <- seq(0, nlag, 0.25)

# plot overall response across all lags 
plot(predt,"overall", xlab = expression(paste("Temperature (",degree,"C)")), 
     ylab = "Relative risk", main = "", 
     ylim = c(range(predt$allRRlow,predt$allRRhigh)))
box()
mtext(side = 2,text="a",cex = 1.2, las = 2, 
      at = max(predt$allRRhigh)*1.2,line = 2.5)

# Tmax
# load model Tmax DLNM + PDSI DLNM
load("output/model0.5.RData")

# extract full coef and vcov and create indicators for each term
coef <- model$summary.fixed$mean
vcov <- model$misc$lincomb.derived.covariance.matrix

# select the position of the terms associated with Tmax crossbasis
indt <- grep("basis_tmax", model$names.fixed)

# extract predictions from the Tmax DLNM centred on overall mean Tmax (30 deg C)
predtmx <- crosspred(basis_tmax, coef = coef[indt], vcov=vcov[indt,indt],
                   model.link = "log", bylag = 0.25, cen = round(mean(data$tmax), 0)) 

# plot overall response across all lags 
plot(predtmx,"overall",xlab = expression(paste("Temperature (",degree,"C)")), 
     ylab = "Relative risk", main = "", 
     ylim = c(range(predt$allRRlow,predt$allRRhigh)))
box()
mtext(side = 2,text="b",cex = 1.2, las = 2, 
      at = max(predt$allRRhigh)*1.2,line = 2.5)

dev.off()

# find Tmin/Tmax value at max RR
predt$predvar[which.max(predt$allRRhigh)]
predtmx$predvar[which.max(predtmx$allRRhigh)]

# Step 2: plot exposure lag response associations for PDSI interacted with water supply shortage frequency

# load best fitting model with climate DLNMs but no interactions
# Tmin DLNM + PDSI DLNM
load("output/model0.4.RData")
model0 <- model
# load urban interaction models
load("output/model1.1_water.RData")
model1.1 <- model
load("output/model1.2_water.RData")
model1.2 <- model
load("output/model1.3_water.RData")
model1.3 <- model

# Step 2: load models and plot exposure lag response associations for PDSI at high and low frequency of water supply shortages

# create model name and label strings
mod.name <- c("model0", "model1.1", "model1.2", "model1.3")
lab <- c("a", "a", "b", "b") #  this avoids plotting the overall and intermediate scenario, to view all change second a to b and final c to d

# make table to save relative risks for wet and dry scenarios
table1 <- as.data.frame(matrix(NA, 4, 11))
colnames(table1) <- c("Setting", 
                      "wet_var", "wet_lag", "wet_rr", "wet_lci","wet_uci",
                      "dry_var", "dry_lag", "dry_rr", "dry_lci","dry_uci")
table1[,1] <- c("Overall", "High","Intermediate", "Low")

for (j in 1:length(mod.name))
{
  
  model <- eval(parse(text = as.name(mod.name[j]))) 
  
  # extract coefficients and variance-covariance matrix
  coef <- model$summary.fixed$mean
  vcov <- model$misc$lincomb.derived.covariance.matrix
  
  # create indicators for terms associated with PDSI cross basis
  indp <- grep("basis_pdsi",model$names.fixed)
  
  # extract predictions from the PDSI DLNM centred on zero (normal conditions)
  predp <- crosspred(basis_pdsi, coef = coef[indp], vcov = vcov[indp,indp],
                     model.link = "log", bylag = 0.25, cen = 0)
  
  # define range limits
  minlim <- min(lag_pdsi, na.rm = T)
  maxlim <- max(lag_pdsi, na.rm = T)
  
  # contour plot of exposure-lag-response associations (Appendix Fig A12a and A12b)
  pdf(paste0("figs/fig_S12",lab[j],"_pdsi_surface.pdf"), width = 6.5, height = 6)
  
  y <- predp$predvar
  x <- seq(0, nlag, 0.25)
  z <- t(predp$matRRfit)
  
  pal <- rev(brewer.pal(11, "PRGn"))
  levels <- pretty(c(0, z, 2), 20)
  col1 <- colorRampPalette(pal[1:6])
  col2 <- colorRampPalette(pal[6:11])
  cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))
  
  filled.contour(x, y, z,
                 xlab = "Lag", ylab = "Drought severity index", main = "",
                 ylim = c(minlim, maxlim), col = cols,levels = levels,
                 plot.axes = { axis(1, at = 0:nlag, seq(0,nlag,1)) 
                   axis(2)})
  mtext(side = 2, at = max(y)*1.2, text = lab[j], las = 2, cex = 1.2, line = 2)
  
  dev.off()
  
  # Table results
  
  # get exposures values
  vars<-predp$predvar
  
  # make dataframe with relative risk for each variable and lag
  extremes <- data.frame(rr = as.vector(predp$matRRfit),
                         rr.lci = as.vector(predp$matRRlow),
                         rr.uci = as.vector(predp$matRRhigh),
                         var = rep(vars, length(lagbylag)),
                         lag = rep(lagbylag, each = length(vars)))
  
  dry <- extremes[extremes$var <= 0,]
  wet <- extremes[extremes$var > 0,]
  
  dry_ind <- which.max(dry$rr)
  wet_ind <- which.max(wet$rr)
  
  table1$wet_var[j] <- wet$var[wet_ind]
  table1$wet_lag[j] <- round(wet$lag[wet_ind],0)
  table1$wet_rr[j]  <- round(wet$rr[wet_ind], 2)
  table1$wet_lci[j] <- round(wet$rr.lci[wet_ind], 2)
  table1$wet_uci[j] <- round(wet$rr.uci[wet_ind], 2)
  table1$dry_var[j] <- dry$var[dry_ind]
  table1$dry_lag[j] <- round(dry$lag[dry_ind], 0)
  table1$dry_rr[j]  <- round(dry$rr[dry_ind], 2)
  table1$dry_lci[j] <- round(dry$rr.lci[dry_ind], 2)
  table1$dry_uci[j] <- round(dry$rr.uci[dry_ind], 2)
}

# save relative risk results (Appendix Table A3)
write.csv(table1, file = "figs/table_S03.csv", quote = FALSE, row.names = FALSE)

# Step 3: plot lag response associations for extreme values of PDSI given high and low frequency of water supply shortages


# plot lag response at extreme PDSI values in areas with: 
# a) high frequency water supply shortages (model1.1 parameterisation) and 
# b) low frequency water supply shortages (model1.3 parameterisation) (Appendix Fig S13)

pdf("figs/fig_S13_pdsi_scenario.pdf", width = 12, height = 6)
par(mfrow = c(1,2))

lab <- c("", "a", "", "b")

for (j in c(2,4))
{
  
  model <- eval(parse(text = as.name(mod.name[j]))) 
  
  # extract coeficients and variance-covariance matrix
  coef <- model$summary.fixed$mean
  vcov <- model$misc$lincomb.derived.covariance.matrix
  
  # create indicators for terms associated with PDSI cross basis
  indp <- grep("basis_pdsi",model$names.fixed)
  
  # extract predictions from the PDSI DLNM centred on zero (normal conditions)
  predp <- crosspred(basis_pdsi, coef = coef[indp], vcov = vcov[indp,indp],
                     model.link = "log", bylag = 0.25, cen = 0)
  
  # get exposures values
  vars<-predp$predvar
  
  # obtain RR fit and upper and lower confidence limits for all exposure variables
  rr <-predp$matRRfit
  rr.lci <- predp$matRRlow
  rr.uci <- predp$matRRhigh
  
  # set relative risk range
  r1 <- min(range(rr, rr.lci, rr.uci))
  r2 <- max(range(rr, rr.lci, rr.uci))
  
  # selected values
  mn <- which(round(vars, 2) == -7)
  mx <- which(round(vars, 2) == 7)
  
  # define colors
  col1 <- brewer.pal(11, "BrBG")[2]
  tcol1 <- do.call(rgb, c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))
  
  col2 <- brewer.pal(11, "BrBG")[10]
  tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/4, max = 255))
  
  # define x values (lag, by lag)
  lagbylag <- seq(0, nlag, 0.25)
  
  # dry
  plot(lagbylag, rr[mn,],
       col = col1, type = "l", lwd = 1, 
       xlab = "Lag", ylab = "Relative risk", main = "", 
       ylim = c(r1, max(r2, 2.7)), frame.plot = T, axes = F)
  axis(1, at = 0:nlag, labels = 0:nlag)
  axis(2)
  # points(lagbylag, rr[mn,], col = col1, pch = 20)
  xx <- c(lagbylag, rev(lagbylag))
  yy <- c(rr.lci[mn,], rev(rr.uci[mn,]))
  polygon(xx, yy, col = tcol1, border = tcol1)
  # wet
  lines(lagbylag, rr[mx,], col = col2, lwd = 1)
  xx<-c(lagbylag, rev(lagbylag))
  yy<-c(rr.lci[mx,], rev(rr.uci[mx,]))
  polygon(xx, yy, col = tcol2, border = tcol2)
  abline(h = 1, lty = 3)
  # legend("top", legend = c(paste0("PDSI = ",vars[mn]) , paste0("PDSI = ",vars[mx])), col = c(col1, col2), pch=20, lwd = 2, lty = 1, bty = "n", y.intersp = 1.5)
  legend("top", legend = c("Extreme drought","Extremely wet"), col = c(col1, col2), lwd = 2, lty = 1, bty = "n", y.intersp = 1.5)
  mtext(side = 2, text = lab[j], cex = 1.2, las = 2, at = max(r2, 2.7)*1.1, line = 2)
}

dev.off() 


