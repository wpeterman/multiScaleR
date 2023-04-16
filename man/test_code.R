
# USE PACKAGE CODE --------------------------------------------------------
# code_files <- list.files("R/", full.names = T)
# sapply(code_files, source)
library(multiScaleR)
# library(terra)
# library(sf)
# library(fields)

data('count_data')
hab <- rast(system.file('data/hab.tif', package = 'multiScaleR'))
data('pts')

pts <- pts
count_data <- count_data

mod <- glm(y ~ hab,
           family = 'poisson',
           data = count_data)

summary(mod)

kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = hab,
                             max_D = 300,
                             kernel = 'gaussian')

opt_pois <- multiScale_optim(fitted_mod = mod,
                             # par = 100/250,
                             kernel_inputs = kernel_inputs)
summary(opt_pois)
plot(opt_pois)

# source("man/kernel_scale_funcs--20221222.R")
source("man/sim_functions.R")
# code_files <- list.files("R/", full.names = T)
# sapply(code_files, source)
library(multiScaleR)
# library(terra)
# library(sf)
#  *** Test Run *** -------------------------------------------------------

# Simulated surfaces ------------------------------------------------------
r_stack <- rast_sim(dim = 300)
raster_stack <- subset(r_stack, c(1,2))

# raster_stack <- subset(r_stack, 1)

plot(raster_stack)
# plot(kernel_scale.raster(r_stack, sigma = 25))
# plot(kernel_scale.raster(r_stack, c(75,250, 25, 50)))


# Simulate kernel-scale data ----------------------------------------------
alpha <- -1.5
sigma <- c(100, 40)
shape <- c(2,5)
beta <- c(1., -0.75)
n_points <- 100
max_D <- 400
kernel <- 'gaussian'
# unit_conv <- 1e5
StDev <- 1

count_dat <- sim_dat.kernel(alpha = alpha,
                            beta = beta,
                            n_points = n_points,
                            raster_stack = raster_stack,
                            kernel = kernel,
                            sigma = sigma,
                            shape = shape,
                            max_D = max_D,
                            # unit_conv = unit_conv,
                            type = 'count',
                            user_seed = 321)

plot(count_dat$df$y ~ count_dat$df$bin1)
plot(count_dat$df$y ~ count_dat$df$bin2)

# length(unique(cellFromXY(raster_stack[[1]], count_dat$pts)))
# plot(raster_stack[[1]]); plot(count_dat$pts, add = T, pch = 20)
# corr <- ncf::correlog(count_dat$pts@coords[,1],
#                       count_dat$pts@coords[,2],
#                       count_dat$obs,
#                       increment = 10)
# plot(corr)

## For creating variogram: https://gsp.humboldt.edu/olm/R/04_01_Variograms.html
## http://rstudio-pubs-static.s3.amazonaws.com/78780_9354d81680f74fa7808be320c4769039.html
# vario <- gstat::variogram(obs ~ 1, data = count_dat$pts)
# plot(vario)
# vario_mod <- gstat::vgm(psill = 8, model = "Exp", nugget = 0.01, range = 100)
# plot(vario, model = vario_mod)
# fit_vario <- gstat::fit.variogram(vario, model = vario_mod)
# plot(vario, model = fit_vario)

count_dat.nb <- sim_dat.kernel(alpha = alpha,
                               beta = beta,
                               n_points = n_points,
                               raster_stack = raster_stack,
                               # kernel = 'gaussian',
                               sigma = sigma,
                               shape = shape,
                               max_D = max_D,
                               # unit_conv = unit_conv,
                               type = 'count_nb',
                               StDev = StDev,
                               user_seed = 112182)
plot(count_dat.nb$df$y ~ count_dat.nb$df$bin1)

occ_dat <- sim_dat.kernel(alpha = alpha,
                          beta = beta,
                          n_points = n_points,
                          raster_stack = raster_stack,
                          kernel = 'gaussian',
                          sigma = sigma,
                          max_D = max_D,
                          # unit_conv = unit_conv,
                          type = 'occ')
mean(occ_dat$obs)
plot(occ_dat$df$y ~ occ_dat$df$bin1)

gaus_dat <- sim_dat.kernel(alpha = alpha,
                           beta = beta,
                           n_points = n_points,
                           raster_stack = raster_stack,
                           kernel = 'gaussian',
                           sigma = sigma,
                           max_D = max_D,
                           # unit_conv = unit_conv,
                           StDev = StDev,
                           type = 'gauss')
plot(gaus_dat$df$y ~ gaus_dat$df$bin1)
abline(lm(gaus_dat$df$y ~ gaus_dat$df$bin1))


summary((cnt_mod <- glm(y ~ .,
                        family = 'poisson',
                        data = count_dat$df)))

acf(resid(cnt_mod))

summary(cnt_mod.nb <- MASS::glm.nb(y ~ .,
                                   data = count_dat.nb$df))
acf(resid(cnt_mod.nb))


summary(occ_mod <- glm(y ~ .,
                       family = 'binomial',
                       data = occ_dat$df))
acf(resid(occ_mod))

summary(gaus_mod <- lm(y ~ .,
                       data = gaus_dat$df))
acf(resid(gaus_mod))


# Optimize kernel --------------------------------------------------------

# >> Poisson --------------------------------------------------------------

fit_pois1 <- glm(y ~ bin2 + bin1,
                 family = 'poisson',
                 data = count_dat$df)

fit_pois2 <- glm(y ~ .,
                 family = 'poisson',
                 data = count_dat$df)
summary(fit_pois1)
# summary(fit_pois2)


# r_stack2 <- c(r_stack$cont1,r_stack$bin1)

opt_input.pois <- kernel_prep(pts = count_dat$pts,
                              kernel = 'expow',
                              max_D = max_D,
                              raster_stack = r_stack)

# opt_input.pois2 <- kernel_prep(pts = count_dat$pts,
#                                kernel = 'gaussian',
#                                max_D = max_D,
#                                raster_stack = r_stack)

opt_pois1 <- multiScale_optim(fitted_mod = fit_pois1,
                              kernel_inputs = opt_input.pois,
                              # par = sigma/unit_conv,
                              opt_parallel = T,
                              n_cores = 6)
# opt_pois2 <- multiScale_optim(fitted_mod = fit_pois2,
                              # kernel_inputs = opt_input.pois,
                              # method = ""
                              # par = c(75,75)/unit_conv,
                              # opt_parallel = T,
                              # n_cores = 6)


summary(opt_pois1)
# summary(opt_pois2)

plot(opt_pois1)
smooth_rast <- kernel_scale.raster(raster_stack = r_stack,
                                   scale_opt = opt_pois1)
plot(smooth_rast)


data.frame(true_mod = logLik(cnt_mod),
           opt_mod = logLik(opt_pois$opt_mod))


# scalescape --------------------------------------------------------------
#
# library(scalescape)
# l_mat <- landscape_matrix(raster(r_stack),
#                           sites = count_dat$pts,
#                           25)
# dist_weight(mod0 = fit_pois2,
#             landscape.vars = l_mat,
#             landscape.formula = ~bin2 + bin1,
#             data=count_dat$df)

# >> Occupancy --------------------------------------------------------------

fit_occ <- glm(y ~ .,
               family = 'binomial',
               data = occ_dat$df)

opt_input.occ <- kernel_prep(pts = occ_dat$pts,
                             max_D = max_D,
                             kernel = 'gaussian',
                             unit_conv = unit_conv,
                             raster_stack = raster_stack)

opt_occ <- multiScale_optim(fitted_mod = fit_occ,
                            kernel_inputs = opt_input.occ,
                            n_cores = 6)

summary(opt_occ)

data.frame(true_mod = logLik(occ_mod),
           opt_mod = logLik(opt_occ$opt_mod))

# >> Gaussian --------------------------------------------------------------

fit_gaus <- glm(y ~ .,
                # family = 'binomial',
                data = gaus_dat$df)

opt_input.gaus <- kernel_prep(pts = gaus_dat$pts,
                              kernel = 'gaussian',
                              max_D = max_D,
                              unit_conv = unit_conv,
                              raster_stack = raster_stack)

opt_gaus <- multiScale_optim(fitted_mod = fit_gaus,
                             kernel_inputs = opt_input.gaus,
                             n_cores = 6)

opt_gaus
summary(opt_gaus)

data.frame(true_mod = logLik(gaus_mod),
           opt_mod = logLik(opt_gaus$opt_mod))



# Simulations to conduct: -------------------------------------------------

## 1) Effect of sample size on precision of scale estimate
## 2) Optimization time as raster size increases
## 3) Precision of scale estimate and correct sign of covariates for count, occ, and gaus data
## 4) XXX
## 5) XXX
