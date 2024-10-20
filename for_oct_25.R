#---------------------------------------------------#
# Exploring spatial variation in Jonah crab size
# Ruby Krasnow
# Last updated Oct 20, 2024
#---------------------------------------------------#


# 1: Initial setup -----------------------------------------------------------

# List of packages required:
packages <- c("tidyverse", "PNWColors", "janitor", "broom", "DHARMa", "performance", "sdmTMB", "sf", "marmap")

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


# Read in crab data -------------------------------------------------------

crabs_raw <- read_csv("./data/crab_data.csv", show_col_types = FALSE)

crabs_initial <- crabs_raw %>% filter(source == "NEFSC", latitude > 39.5) %>%
  select(longitude, latitude, year, station_id, crab_id, ch, cw, zone) %>% 
  mutate(zone = as.factor(zone))

# Outlier removal - round 1 ---------------------------------------------------------

# multivariate outlier identification using Mahalanobis distance and
# Minimum Covariance Determinant (MCD), a robust version
outliers_m <- check_outliers(
  x = crabs_initial %>% select(cw, ch, crab_id),
  method = c("mahalanobis", "mcd"),
  ID = "crab_id"
) %>%
  as_tibble() %>%
  select(crab_id, Outlier)

crabs <- crabs_initial %>% left_join(outliers_m) %>%
  filter(Outlier < 1)


# Clustering --------------------------------------------------------------

dat <- cbind(crabs$cw, crabs$ch) # create array for input to mclust function

# Perform Gaussian mixture model-based clustering
crabs_class_raw <- Mclust(data = dat, G = 2, modelNames = "EVV")

# Add cluster-based classifications to the data frame
crabs_class <- crabs %>%
  mutate(pred_class = crabs_class_raw$classification)

# which cluster number corresponds to the adult males?
adult_num <- crabs_class %>%
  slice_max(ch) %>%
  pull(pred_class)

# standardize adults to be cluster 1 and juveniles to be cluster 0
crabs_class <- crabs_class %>%
  mutate(
    pred_class = if_else(pred_class == adult_num, as.factor(1), as.factor(0)),
    pred_class_num = if_else(pred_class == adult_num, as.numeric(1), as.numeric(0))
  )



# Outlier removal - from clustering ---------------------------------------

# fit linear equation to juvenile males
lin_eq <- lm(ch ~ cw, data = crabs_class %>% filter(pred_class == 0))

# keep only the juvenile males and the males classified as adults
# that lie below the juvenile line
res_remove_orig <- crabs_class %>%
  filter(!((ch > coef(lin_eq)[1] + coef(lin_eq)[2] * cw) & pred_class == 1))

res_remain <- crabs_class %>%
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
  labs(x = "Carapace width (mm)", y = "Claw height (mm)",
    color = "Predicted cluster")

# Set up bathymetry -------------------------------------------------------

bat <- getNOAA.bathy(-65, -75, 45, 39.5, resolution = 4)

sdm_coords <- jonah %>% select(longitude, latitude)

depth_df <- get.depth(bat, sdm_coords, locator = FALSE)


# Set-up for sdmTMB model -------------------------------------------------

jonah <- jonah %>% bind_cols(depth_df) %>% add_utm_columns(units="km") %>% 
  mutate(depth = abs(depth))
  
mesh <- make_mesh(jonah, c("X", "Y"), cutoff = 8)
plot(mesh)


no_spat <- sdmTMB(
  data = jonah,
  formula = cw ~ depth + pred_class,
  mesh = mesh, # can be omitted for a non-spatial model
  spatial = "off")

no_year <- sdmTMB(
  data = jonah,
  formula = cw ~ depth + pred_class,
  mesh = mesh,
  spatial = "on")

lin_mod <- sdmTMB(
  data = jonah,
  formula = cw ~ depth + year + pred_class,
  mesh = mesh,
  spatial = "on")

gam_mod <- sdmTMB(
  data = jonah,
  formula = cw ~ s(depth) + pred_class,
  mesh = mesh, # can be omitted for a non-spatial model
  spatial = "on")

AIC(no_spat, no_year, lin_mod, gam_mod)

summary(lin_mod)

jonah$scaled_year <- (jonah$year - mean(jonah$year))

scaled_mod <- sdmTMB(
  data = jonah,
  formula = cw ~ depth + scaled_year + pred_class,
  mesh = mesh,
  spatial = "on")

summary(scaled_mod)
tidy(scaled_mod, conf.int = TRUE)

summary(no_year)
tidy(no_year, conf.int = TRUE)

mod <- no_year

# Residual checking -------------------------------------------------------

qqnorm(residuals(mod, type = "mle-mvn"))
qqline(residuals(mod, type = "mle-mvn"))
hist(residuals(mod, type = "mle-mvn"))

s_mod <- simulate(mod, nsim = 500, type = "mle-mvn")
dharma_residuals(s_mod, mod)

r_mod <- dharma_residuals(s_mod, mod, return_DHARMa = TRUE)
plot(r_mod)
testDispersion(r_mod)

# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html 
# unconcerned about underdispersion bc it will just bias p-values to be on 
# the conservative side


# Visualization -----------------------------------------------------------

fdepth<- ggeffect(mod, "depth")
plot(fdepth)

jonah_sf <- st_as_sf(jonah, coords=c("longitude", "latitude"), crs=4326)

coast <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", country = "united states of america") %>% 
  st_crop(c(xmin = -75, ymin = 10, xmax = -50, ymax = 50)) %>% 
  st_transform(crs=32619)

ggplot(data=jonah_sf)+
  facet_wrap(~year)+
  geom_sf(aes(color=zone))+
  mytheme+
  geom_sf(data=coast)+labs(color="Zone")

grid0 <- st_make_grid(jonah_sf, n=50, what = "centers") # a square grid over the region with observed pts

points_reg <- jonah_sf %>%  # an outline of the region with observed points (polygon)
  st_union() %>% 
  st_convex_hull() 

cont <- st_contains(points_reg, grid0)[[1]] # which points from the square grid are inside the non-square (convex hull) polygon with observed points

grid1 <- grid0[cont] %>% # take only the points in the observed polygon, 
  st_transform(crs=32619) %>% 
  st_coordinates() %>% 
  as_tibble()

plotting_coords <- grid0[cont] %>% 
  st_coordinates() %>% 
  as_tibble()

depth_df <- get.depth(bat,
                      x = plotting_coords$X,
                      y = plotting_coords$Y,
                      locator = FALSE) 


plotting_coords <- grid1 %>% bind_cols(depth_df) %>% 
  as_tibble() %>% 
  mutate(X=X/1000, Y=Y/1000) #convert from m to km

grid_classes <- replicate_df(plotting_coords, "pred_class", unique(jonah$pred_class))

predictions <- predict(mod, newdata=grid_classes, type="response") %>% as.data.frame()

preds_plot <- ggplot() + 
  geom_tile(data=predictions, aes(x = X*1000, y = Y*1000, color = est), linewidth=4) +
  geom_sf(data=coast) +
  theme_light() +
  labs(color = "CW (mm)") +
  labs(x = NULL, y = NULL)+
  scale_color_viridis_c()+
  theme(text = element_text(size=14))

preds_plot
