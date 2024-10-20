

# List of packages required:
packages <- c("tidyverse", "PNWColors", "janitor", "broom", "DHARMa", "performance", "sdmTMB", "sf")

# Load packages into session
lapply(packages, require, character.only = TRUE)
rm(packages)

# Ensure functions with duplicate names are from the correct package
select <- dplyr::select
map <- purrr::map
summarize <- dplyr::summarize
clean_names <- janitor::clean_names
margin <- ggplot2::margin

set.seed(123) #Set seed for pseudo-random number generator, for reproducibility

mytheme <- theme_light()+ #define custom theme for ggplots
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, l = 0)),
        text=element_text(size=15))

crabs <- read_csv("./data/crab_data.csv", show_col_types = FALSE)

glimpse(crabs)

crabs <- crabs %>% filter(source == "NEFSC", latitude > 39.5) %>%
  select(
    zone,
    longitude,
    latitude,
    year,
    station_id,
    crab_id,
    ch,
    cw) %>% mutate(zone = as.factor(zone))


sdm_df <- add_utm_columns(crabs, units="km") %>% 
  mutate(loc = paste0("(", X, ", ", Y, ")"))


mesh <- make_mesh(sdm_df, c("X", "Y"), cutoff = 7)
plot(mesh)



coast <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", country = "united states of america") %>% 
  st_crop(c(xmin = -75, ymin = 10, xmax = -50, ymax = 50)) %>% 
  st_transform(crs=32619)

males_sf <- st_as_sf(crabs, coords=c("longitude", "latitude"), crs=4326)

ggplot(data=males_sf)+
  facet_wrap(~year)+
  geom_sf(aes(color=zone))+
  mytheme+
  geom_sf(data=coast)+labs(color="Zone")

library(marmap)
bat <- getNOAA.bathy(-65, -75, 45, 39.5, resolution = 4)

# Creating color palettes
blues <- c("lightsteelblue4", "lightsteelblue3",
           "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

plot(bat, image = TRUE, land = TRUE, lwd = 0.03,
     bpal = list(c(0, max(bat), greys),
                 c(min(bat), 0, blues)))

sdm_coords <- crabs %>% select(longitude, latitude)

depth1 <- get.depth(bat, sdm_coords, locator = FALSE)

plot(bat, image = TRUE, land = TRUE, n=1,
     bpal = list(c(0, max(bat), greys),
                 c(min(bat), 0, blues)))

points(depth1$lon, depth1$lat, pch = 21, col = "black",
       bg = "yellow", cex = 1.3)

text(-70, 40, "Gulf of \nMaine", col = "white", font = 3)

# Add coastline
plot(bat, n = 1, lwd = 0.4, add = TRUE)

plot(bat, image=TRUE, deep=-6000, shallow=0, step=1000, add=TRUE)
summary(bat)

# url_base <- "https://coastwatch.pfeg.noaa.gov/erddap/"
# parameter <- "sea_floor_depth"
# lats <- c(39.5, 45)
# lons <- c(-73, -66)
# dat_info <- rerddap::info("erdEtopo22SeafloorGradient_Atlantic", url = url_base)
# #depth_data <- griddap(dat_info, longitude = lons, latitude = lats, fields = parameter, url = url_base)
# 
# depth_data <- rxtracto_3D(dat_info, parameter = parameter, xcoord = lons, ycoord = lats)
# 
# plotBBox(depth_data, maxpixels = 1000)
# 
# depth_tidy <- tidy_grid(depth_data)
# 
# sdm_df %>% left_join(depth_tidy) %>% view()

sdm_df2 <- males_sf %>% 
  bind_cols(depth1) %>% add_utm_columns(units="km",
                                        ll_names = c("lon", "lat")) %>% 
  mutate(depth = abs(depth))
  
sdm_df2 %>% 
  ggplot()+geom_sf(aes(color = depth))+
  mytheme+
  geom_sf(data=coast)+labs(color="Depth (m)") + scale_color_viridis_c()




lm1 <- lm(cw ~ depth + lat + lon, data = sdm_df2)

check_residuals(lm1)

summary(lm1)
hist(log(sdm_df2$cw))

sdm_df3 <- st_drop_geometry(sdm_df2)

sdm_df3$scaled_year <- (sdm_df3$year - mean(sdm_df3$year))

m <- sdmTMB(
  data = sdm_df3,
  formula = cw ~ s(depth),
  mesh = mesh, # can be omitted for a non-spatial model
  spatial = "on",
  spatial_varying = ~ 0 + scaled_year,
  #family = Gamma()
  family = tweedie(link = "log")
)
m


# residuals ---------------------------------------------------------------

qqnorm(residuals(m, type = "mle-mvn"))
qqline(residuals(m, type = "mle-mvn"))
hist(residuals(m, type = "mle-mvn"))
# set.seed(123)
# samps <- sdmTMBextra::predict_mle_mcmc(m, mcmc_warmup = 300, mcmc_iter = 1000)
# 
# r <- residuals(m, "mle-mcmc", mcmc_samples = samps)
# qqnorm(r)
# qqline(r)

s_mod <- simulate(m, nsim = 500, type = "mle-mvn")
dharma_residuals(s_mod, m)

r_mod <- dharma_residuals(s_mod, m, return_DHARMa = TRUE)
plot(r_mod)

library(mgcv)
library(gratia)


gam1 <- gam(cw ~ s(depth) + year + lat + lon, data = sdm_df4)

appraise(gam1)

hist(sdm_df3$cw)

# multivariate outlier identification using Mahalanobis distance and
# Minimum Covariance Determinant (MCD), a robust version
outliers_m <- check_outliers(
  x = sdm_df3 %>% select(cw, ch, crab_id),
  method = c("mahalanobis", "mcd"),
  ID = "crab_id"
) %>%
  as_tibble() %>%
  select(crab_id, Outlier)

sdm_df4 <- sdm_df3 %>% left_join(outliers_m) %>% 
  filter(Outlier < 1)

dat <- cbind(sdm_df4$cw, sdm_df4$ch) # create array for input to mclust function

# Perform Gaussian mixture model-based clustering
jonah_class_raw <- Mclust(data = dat, G = 2, modelNames = "EVV")

# Add cluster-based classifications to the data frame
jonah_class <- sdm_df4 %>%
  mutate(pred_class = jonah_class_raw$classification)

# which cluster number corresponds to the adult males?
adult_num <- jonah_class %>%
  slice_max(ch) %>%
  pull(pred_class)

# standardize adults to be cluster 1 and juveniles to be cluster 0
jonah_class <- jonah_class %>%
  mutate(
    pred_class = if_else(pred_class == adult_num, as.factor(1), as.factor(0)),
    pred_class_num = if_else(pred_class == adult_num, as.numeric(1), as.numeric(0))
  ) %>%
  mutate(uncertainty = jonah_class_raw$uncertainty)

# fit linear equation to juvenile males
lin_eq <- lm(ch ~ cw, data = jonah_class %>% filter(pred_class == 0))

# keep only the juvenile males and the males classified as adults
# that lie below the juvenile line
res_remove_orig <- jonah_class %>%
  filter(!((ch > coef(lin_eq)[1] + coef(lin_eq)[2] * cw) & pred_class == 1))

res_remain <- jonah_class %>%
  anti_join(res_remove_orig) %>%
  filter(cw > 70)

# Visualize points to be evaluated for outlier removal
ggplot(data = res_remove_orig) +
  geom_point(aes(x = cw, y = ch, color = pred_class)) +
  mytheme +
  labs(
    x = "Carapace width (mm)", y = "Claw height (mm)",
    title = "Points to be evaluated for outlier removal",
    color = "Predicted cluster"
  )

max_rem <- 40
mse_vec <- rep(NA, max_rem)

res_remove_temp <- res_remove_orig

for (i in c(1:max_rem)) {
  mod <- lm(data = res_remove_temp, log(ch) ~ log(cw))
  mse_vec[i] <- mse(mod)
  most_neg <- which.min(residuals(mod))
  res_remove_temp <- res_remove_temp %>% slice(-most_neg)
}

mse_diffs <- diff(mse_vec)

ggplot() +
  geom_smooth(aes(x = c(1:(max_rem - 1)), y = abs(mse_diffs)),
              se = FALSE, method = "gam", color = "#A5AAAF"
  ) +
  geom_point(aes(x = c(1:(max_rem - 1)), y = abs(mse_diffs)), size = 2) +
  mytheme +
  labs(x = "Points removed", y = "Change in MSE")

res_remove <- res_remove_orig %>%
  mutate(res = residuals(lm(data = res_remove_orig, log(ch)~ log(cw)))) %>%
  slice_max(res, n = -25)

jonah <- res_remove %>%
  select(-res) %>%
  bind_rows(res_remain) %>%
  mutate(pred_class = if_else((ch < coef(lin_eq)[1] + coef(lin_eq)[2] * cw) &
                                pred_class == 1, as.factor(0), pred_class))


ggplot(data = jonah) +
  geom_point(aes(x = cw, y = ch, color = pred_class)) +
  mytheme +
  labs(
    x = "Carapace width (mm)", y = "Claw height (mm)",
    color = "Predicted cluster"
  )

mesh_class <- make_mesh(jonah, c("X", "Y"), cutoff = 8)
plot(mesh_class)

m2 <- sdmTMB(
  data = jonah,
  formula = cw ~ depth + year + pred_class,
  mesh = mesh_class, # can be omitted for a non-spatial model
  spatial = "on",
)
m2


# residuals classification ------------------------------------------------


qqnorm(residuals(m2, type = "mle-mvn"))
qqline(residuals(m2, type = "mle-mvn"))
hist(residuals(m2, type = "mle-mvn"))
# set.seed(123)
# samps <- sdmTMBextra::predict_mle_mcmc(m, mcmc_warmup = 300, mcmc_iter = 1000)
# 
# r <- residuals(m, "mle-mcmc", mcmc_samples = samps)
# qqnorm(r)
# qqline(r)

s_mod <- simulate(m2, nsim = 500, type = "mle-mvn")
dharma_residuals(s_mod, m2)

r_mod <- dharma_residuals(s_mod, m2, return_DHARMa = TRUE)
plot(r_mod)

testResiduals(r_mod)
testDispersion(r_mod)
