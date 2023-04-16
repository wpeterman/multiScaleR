# load("C:/Users/peterman.73/OneDrive - The Ohio State University/Teaching/ENR5374_LANDECO/Labs/Lab_Dev/Lab5/sim_dat.RData")
# library(terra)
# library(sf)
#
# pts <- count_dat$pts
# count_data <- count_dat$df
# hab <- rast("C:/Users/peterman.73/OneDrive - The Ohio State University/Teaching/ENR5374_LANDECO/Labs/Lab_Dev/Lab5/hab.tif")
#
# usethis::use_data(pts, count_data, hab)
#

## 3/26/2023 Update
# source('man/sim_functions.R')
library(multiScaleR)

n_points <- 100
alpha <- 0.25
sigma <- c(250, # Spatial
           450, # Local
           500, # Spatial
           300) # Not included
kernel <- 'gaussian'#'expow'
# kernel <- 'expow'
shape = 3.5
# beta <- c(0.75,-0.5,0.5,0)
# beta <- c(0.3,-0.6,0.5,0)
beta <- c(-0.5,0.3,0.7,0)
max_D <- 2000
StDev <- 1
cnt <- 1
r_dim <- 200

r_stack <- sim_rast(dim = r_dim,
                    resolution = 20,
                    user_seed = cnt * 555)
# raster_stack <- subset(r_stack, c(1))
raster_stack <- r_stack

## Poisson
count_dat <- sim_dat(alpha = alpha,
                     beta = beta,
                     n_points = n_points,
                     raster_stack = raster_stack,
                     sigma = sigma,
                     # kernel = kernel,
                     shape = shape,
                     max_D = max_D,
                     # type = 'count',
                     user_seed = cnt * 543210)

plot(count_dat$df$y ~ count_dat$df$bin1)
plot(count_dat$df$y ~ count_dat$df$bin2)
plot(count_dat$df$y ~ count_dat$df$cont1)

plot(raster_stack);
plot(kernel_scale.raster(raster_stack,kernel = kernel, sigma = sigma))
plot(raster_stack$bin1);plot(st_geometry(count_dat$pts), add=T, pch=19)
plot(raster_stack$bin2);plot(st_geometry(count_dat$pts), add=T, pch=19)
plot(raster_stack$cont1);plot(st_geometry(count_dat$pts), add=T, pch=19)


# fit_mod <- glm(y ~ .,
#                family = 'poisson',
#                data = count_dat$df)
# summary(fit_mod)
#
# opt_input <- kernel_prep(pts = count_dat$pts,
#                          raster_stack = raster_stack,
#                          kernel = kernel,
#                          max_D = max_D)
#
# opt <- try(multiScale_optim(fitted_mod = fit_mod,
#                             kernel_inputs = opt_input,
#                             par = sigma/max_D,
#                             n_cores = 6))
# summary(opt)
# plot(opt)
#
# ## Opt True
max_D <- 2000
fit_mod2 <- glm(y ~ bin1 + bin2 + cont1,
                family = 'poisson',
                data = count_dat$df)
summary(fit_mod2)

opt_input2 <- kernel_prep(pts = count_dat$pts,
                          raster_stack = subset(raster_stack,c(1,3,4)),
                          # kernel = kernel,
                          max_D = max_D)

system.time(opt2 <- try(multiScale_optim(fitted_mod = fit_mod2,
                                         kernel_inputs = opt_input2,
                                         par = sigma[c(1,3)]/max_D,
                                         n_cores = 6))) ## <100 seconds
summary(opt2)
plot(opt2)


# Create package data -----------------------------------------------------

count_data <- count_dat$df[,c(1,3)]
names(count_data) <- c('counts', 'site')

surv_pts <- count_dat$pts

landscape_rast <- subset(raster_stack, c(1,3,4))
names(landscape_rast) <- c('land1', 'land2', 'land3')

## Create buffers to get initial summary of landscapes around points
pt_buff <- buffer(vect(surv_pts), 500)

# land_vars <- extract(landscape, pt_buff, 'mean', ID = FALSE)
land_vars <- extract(landscape_rast, vect(surv_pts), 'mean', ID = FALSE)
landscape_counts <- data.frame(count_data, (land_vars))

kernel_inputs <- kernel_prep(pts = surv_pts,
                             raster_stack = landscape_rast,
                             max_D = 2000,
                             # kernel = 'gaussian',
                             projected = T)

mod <- glm(counts ~ site + land1 + land2,
           data = landscape_counts,
           family = 'poisson')
summary(mod)

system.time(opt <- try(multiScale_optim(fitted_mod = mod,
                                       kernel_inputs = kernel_inputs,
                                       # par = c(175,350)/max_D,
                                       n_cores = 6)))

summary(opt)
plot(opt)

# sigma <- c(175, # Spatial (land1 -- 0.75)
#            450, # Local (site -- -0.5)
#            350, # Spatial (land2 -- 0.5)
#            300) # Not included (land3 -- 0)
# beta <- c(0.75,-0.5,0.5,0)





# usethis::use_data(surv_pts, landscape_counts, overwrite = TRUE)
# writeRaster(landscape_rast, "data/landscape.tif")

