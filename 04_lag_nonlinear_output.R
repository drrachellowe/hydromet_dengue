# Combined effects of hydrometeorological hazards and urbanisation on dengue risk in Brazil: a spatiotemporal modelling study

# Rachel Lowe (2021)
# https://github.com/drrachellowe/hydromet_dengue

# R script to explore exposure-lag-response relationships

# Step 0: load packages, pre-processed data and models
# Step 1: plot Tmin output 
# Step 2: plot exposure lag response associations for PDSI overall and at high and low levels of urbanisation
# Step 3: plot lag response associations for extreme values of PDSI given high and low levels of urbanisation

# Step 0: load packages, pre-processed data and models

# load packages and pre-processed data
source("00_load_packages_data.R")

# load model output

# load best fitting model with climate DLNMs but no interactions
# Tmin DLNM + PDSI DLNM
load("output/model0.4.RData")
model0 <- model
# load urban interaction models
load("output/model1.1_urban.RData")
model1.1 <- model
load("output/model1.2_urban.RData")
model1.2 <- model
load("output/model1.3_urban.RData")
model1.3 <- model

# Step 1: plot Tmin output 

# Select the hydrometeorological  model (without interactions) 
model <- model0

# extract full coef and vcov and create indicators for each term
coef <- model$summary.fixed$mean
vcov <- model$misc$lincomb.derived.covariance.matrix

# find position of the terms associated with Tmin crossbasis
indt <- grep("basis_tmin", model$names.fixed)

# extract predictions from the Tmin DLNM centred on overall mean Tmin (19 deg C)
predt <- crosspred(basis_tmin, coef = coef[indt], vcov=vcov[indt,indt],
                   model.link = "log", bylag = 0.25, cen = round(mean(data$tmin), 0)) 

# contour and scenario plots for Tmin (Main text Fig 3)

# contour plot of exposure-lag-response associations (Main text Fig 3a)
pdf("figs/fig_03a_Tmin_contour.pdf", width = 6.5, height = 6)

y <- predt$predvar
x <- seq(0, nlag, 0.25)
z <- t(predt$matRRfit)

pal <- rev(brewer.pal(11, "PRGn"))
levels <- pretty(z, 20)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 1)), col2(sum(levels > 1)))

filled.contour(x,y,z,
               xlab = "Lag", ylab = expression(paste("Temperature (",degree,"C)")), main = "",
               col = cols,levels = levels,
               plot.axes = { axis(1, at = 0:nlag, c(0:nlag)) 
                 axis(2)})
mtext(side = 2, at = max(y)*1.1, text = "a", las = 2, cex = 1.2, line = 2)

dev.off()

# lag response for different Tmin scenarios (Main text Fig 3b)
pdf("figs/fig_03b_Tmin_scenario.pdf", width = 6, height = 6)

# get exposures values
vars <- predt$predvar

# obtain relative risk (RR) fit and upper and lower confidence limits for all exposure variables
rr <- predt$matRRfit
rr.lci <- predt$matRRlow
rr.uci <- predt$matRRhigh

# set relative risk range 
r1 <- min(range(rr, rr.lci, rr.uci))
r2 <- max(range(rr, rr.lci, rr.uci))

# get selected exposure variable positions
mn <- which(round(vars, 2) == 15)
mx <- which(round(vars, 2) == 21)
mx2 <- which(round(vars, 2) == 25)

# define colours
col1 <- brewer.pal(11, "RdBu")[9]
tcol1 <- do.call(rgb, c(as.list(col2rgb(col1)), alpha = 255/4, max = 255))

col2 <- brewer.pal(11, "RdBu")[3]
tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/4, max = 255))

col3 <- brewer.pal(11, "RdBu")[1]
tcol3 <- do.call(rgb, c(as.list(col2rgb(col3)), alpha = 255/4, max = 255))

# define x values (lag, by lag)
lagbylag <- seq(0, nlag, 0.25)

# cool
plot(lagbylag, rr[mn,], col = col1, type = "l", lwd = 1, 
     xlab = "Lag", ylab = "Relative risk", main = "", 
     ylim = range(r1, r2*1.11), frame.plot = T, axes = F)
axis(1, at = 0:nlag, labels = 0:nlag)
axis(2)
xx <- c(lagbylag, rev(lagbylag))
yy <- c(rr.lci[mn,], rev(rr.uci[mn,]))
polygon(xx, yy, col = tcol1, border = tcol1)
# warm
lines(lagbylag, rr[mx,], col = col2, lwd = 1)
xx <- c(lagbylag, rev(lagbylag))
yy <- c(rr.lci[mx,], rev(rr.uci[mx,]))
polygon(xx, yy, col = tcol2, border = tcol2)
abline(h = 1, lty = 3)
# warmest
lines(lagbylag, rr[mx2,], col = col3, lwd = 1)
xx <- c(lagbylag, rev(lagbylag))
yy <- c(rr.lci[mx2,],rev(rr.uci[mx2,]))
polygon(xx, yy, col = tcol3, border = tcol3)
abline(h = 1, lty = 3)

legend("topleft",
       legend = c(paste0("Tmin = ",vars[mn]," deg C"),
                            paste0("Tmin = ", vars[mx]," deg C"),
                            paste0("Tmin = ", vars[mx2]," deg C")),
       col = c(col1, col2, col3), 
       lwd = 2, lty = 1, bty = "n", 
       y.intersp = 1.5, horiz = F)

mtext(side = 2, at = r2*1.25, text = "b", las = 2, cex = 1.2, line = 2)

dev.off()

# Step 2: explore exposure lag response associations overall and in urban and rural areas

# create model name and label strings
mod.name <- c("model0", "model1.1", "model1.2", "model1.3")
lab <- c("a", "b", "c", "c") #  this avoid plotting the intermediate scenario, to view all change final c to d

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
  
  # create indicators for terms associated with PDSI cross-basis
  indp <- grep("basis_pdsi",model$names.fixed)
  
  # extract predictions from the PDSI DLNM centred on zero (normal conditions)
  predp <- crosspred(basis_pdsi, coef = coef[indp], vcov = vcov[indp,indp],
                     model.link = "log", bylag = 0.25, cen = 0)
  
  # define range limits
  minlim <- min(lag_pdsi, na.rm = T)
  maxlim <- max(lag_pdsi, na.rm = T)
  
  # contour plot of exposure-lag-response associations (Main text Fig 4a, 4b and 4c)
  pdf(paste0("figs/fig_04",lab[j],"_pdsi_surface.pdf"), width = 6.5, height = 6)
  
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

# save relative risk results (Main text Table)
write.csv(table1, file = "figs/table_01.csv", quote = FALSE, row.names = FALSE)

# Step 3: plot lag response associations for extreme values of PDSI in urban and rural areas

# plot lag response at extreme PDSI values in areas with:
# a) high levels of urbanisation (model1.1 parameterisation) and 
# b) low levels of urbanisation (model1.3 parameterisation) (Main text Fig 5)

pdf("figs/fig_05_pdsi_scenario.pdf", width = 12, height = 6)
par(mfrow = c(1,2))

lab <- c("", "a", "", "b")

for (j in c(2,4))
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

