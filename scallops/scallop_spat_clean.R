#---------------------------------------------------#
# Modeling scallop spat count data
# Ruby Krasnow
# Last updated Nov 8, 2024
#---------------------------------------------------#

# Setup ----------------------------------------------------------

## List of packages required:
packages <- c(
  "tidyverse",
  "ggthemes",
  "gt",
  "janitor",
  "glmmTMB",
  "performance",
  "DHARMa",
  "rstatix",
  "broom.mixed",
  "marginaleffects",
  "emmeans",
  "patchwork",
  "clipr"
)

# Load packages into session
lapply(packages, require, character.only = TRUE)
rm(packages)

set.seed(123)

mytheme <- theme_light()+ # define custom theme for ggplots
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, l = 0)),
        text = element_text(size = 13))

spat2 <- read.csv("./spat_filled.csv") %>%
  clean_names() %>% # clean up column names
  rename(count = count_per_bag) %>%
  select(-c(bag_line, status_open_closed)) %>% # remove unnecessary columns
  complete(line, year, bag_position) %>% # fill in gaps where there should be bags but they were not recovered
  # create new line variable that is not duplicated between years, since line 7 in year X is not necessarily in the same area as line 7 in year Y
  mutate(line2 = paste(year, line, sep = "_")) %>% 
  filter(!is.na(bag_position)) %>% # it also added bag_position = NA for each year and line; get rid of those
  mutate(area = case_when((line < 7) ~ "south", 
                          (line > 6 & line < 10) ~ "north", .default = "closure"))  %>% # add areas
  mutate( # convert variables to factors; don't really have enough years to use as continuous predictor
    bag_position = as_factor(bag_position),
    line2 = as_factor(line2),
    year = as_factor(year),
    area = as_factor(area))

# Calculate percentage of bags recovered by closure status (open/closed)
perc_rec <- spat2 %>% mutate(status = if_else(area == "closure", "closed", "open")) %>%
  filter(!is.na(count)) %>%
  count(year, status) %>%
  rename(bags = n) %>%
  mutate(deployed = if_else(status == "open", 45, 15)) %>%
  mutate(perc = bags / deployed)

# Calculate percentage of bags recovered by area (north, closure, south)
perc_rec_area <- spat2 %>% mutate(status = if_else(area == "closure", "closed", "open")) %>%
  filter(!is.na(count)) %>%
  count(year, area) %>%
  rename(bags = n) %>%
  add_row(year = as.factor(2015),
          area = "south",
          bags = 0) %>%
  mutate(deployed = if_else(area == "south", 30, 15)) %>%
  mutate(perc = bags / deployed)

# Figures -----------------------------------------------------------------

# Percent of bags recovered by year and closure status
ggplot(perc_rec) + 
  geom_bar(aes(x = year, y = perc*100,
               fill = str_to_sentence(status)), 
           stat = "identity",
           position = "dodge", 
           alpha = 0.8, color = NA,
           linewidth = 1, width = 0.8) + 
  mytheme + 
  scale_color_tableau() +
  scale_fill_tableau() +
  labs(x = "Year", y = "% of bags recovered", fill = NULL, color = NULL)

perc_rec %>% mutate(across(where(is.numeric), \(x) round(x, digits = 2))) %>% write_clip()

# Percent of bags recovered by year and closure status
ggplot(perc_rec_area) + 
  geom_bar(aes(x = year, y = perc*100, 
               fill = str_to_sentence(area)),
           color = NA,
           stat = "identity",
           position = "dodge", 
           alpha = 0.8, 
           linewidth = 1, width = 0.8) + 
  mytheme + 
  scale_color_tableau() +
  scale_fill_tableau() +
  labs(x = "Year", y = "% of bags recovered", fill = NULL, color = NULL)

perc_rec_area %>% mutate(across(where(is.numeric), \(x) round(x, digits = 2))) %>% write_clip()

# Number of spat per recovered bag by year and area
spat2 %>% filter(!is.na(count)) %>% 
  ggplot(aes(x = year, y = count, color = str_to_sentence(area), fill = str_to_sentence(area))) +
  geom_boxplot(alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0)) +
  mytheme +
  scale_fill_tableau() + scale_color_tableau() +
  labs(x = "Year", y = "# spat per recovered bag", fill = NULL, color = NULL) +
  scale_y_continuous(labels = scales::comma) +
  theme(panel.grid.minor = element_blank())
  # legend.position = "inside", legend.position.inside = c(0.8, 0.8),


# Modeling ----------------------------------------------------------------

# Fixed effect of interest could be status (open vs closed), but using area is
# more balanced and also allows for a greater understanding of how the closed
# area might affecting larval dispersal in relation to the Western Maine Coastal
# Current (according to the manuscript draft?). We are including area as fixed
# effect because there are arguably too few areas to adequately serve as a
# random effect (say, in a random intercept like (1|area)), but also because we
# care about the differences between the north and south open areas

# Also, glmmTMB cannot estimate all parameters when area and status are in the
# same model, giving a warning that the model is "rank-deficient". Since only 1
# area is closed, we cannot independently estimate the effects of area and
# closure status. The model cannot contain interaction terms with year and area
# because not all combinations of year and area represented due to gear loss (no
# lines were recovered from south of closed area in 2015, despite half of all
# lines being set in the south area)
  
########## "rank-deficient" formulas:
# count ~ year + status + area
# count ~ area + status
# count ~ status + area + bag_position + year
# count ~ year * area (no south lines in 2015)

########## Convergence issue formula
# (model too complex for such a small sample size)
# modD = count ~ status + year + (1|bag_position)

# Define list of formulas to test
zi_formula <- list(
  modA1 = count ~ year,
  modA2 = count ~ area,
  modA3 = count ~ year + area,
  modB1 = count ~ area + (1 | line2),
  modB2 = count ~ year + (1 | line2),
  modB3 = count ~ area + year + (1 | line2),
  modC1 = count ~ area + (1 | bag_position),
  modC2 = count ~ year + (1 | bag_position),
  modC3 = count ~ year + area + (1 | bag_position),
  modD1 = count ~ area + (1 | bag_position) + (1 | line2),
  modD2 = count ~ year + (1 | bag_position) + (1 | line2),
  modD3 = count ~ year + area + (1 | bag_position) + (1 | line2)
)

# fit a new glmmTMB model for each formula in the list using a zero-inflated
# negative binomial distribution note that it is not including the zeros from
# missing bags - those are coded as NA and dropped when fitting the model
zi_models <- tibble(zi_formula, models = purrr::map(
  zi_formula,
  \(x)
  glmmTMB(
    data = spat2,
    formula = x,
    family = nbinom2(),
    ziformula = ~ 1
  )
)) %>%
  add_column(id = names(zi_formula), .before = 1)

# Note: the nbinom2() parameterization of the negative binomial distribution
# lowered AIC values and improved residual normality compared to nbinom1()

zi_models_tidy <- zi_models %>%
  mutate(tidy_model = map(models, tidy),
         AIC = map(models, AIC)) %>% 
  unnest(cols = c(AIC)) %>% 
  mutate(logLik = map(models, ~ (logLik(.x)[1]))) %>% 
  unnest(cols = c(logLik)) %>% 
  arrange(AIC)

tab1 <- bbmle::ICtab(zi_models_tidy$models, mnames = zi_models_tidy$zi_formula, type = c("AIC"), logLik = T)
as.data.frame(tab1) %>% mutate(across(where(is.numeric), \(x) round(x, digits = 2))) %>% write_clip()

bbmle::ICtab(zi_models_tidy$models, mnames = zi_models_tidy$zi_formula, type = c("BIC"), logLik = T)

# refit model with lowest AIC and highest log-likelihood
final_mod <- glmmTMB(count ~ year + area + (1 | bag_position) + (1 | line2),
                   data = spat2,
                   family = nbinom2(),
                   zi =  ~ 1)

# Model summaries
summary(final_mod)
tidy(final_mod) # %>% mutate(across(where(is.numeric), \(x) round(x, digits = 4))) %>% write_clip()
Anova(final_mod)	# %>% mutate(across(where(is.numeric), \(x) round(x, digits = 4))) %>% write_clip()
# Anova(final_mod) is equivalent to emmeans::joint_tests(final_mod)

# Extensive diagnostic checks -----------------------------------------------------------------------

# Residuals look good
testResiduals(final_mod) # same as testUniformity(testmod)
plot(check_residuals(final_mod))+labs(subtitle = NULL)

# Posterior predictive checks look good
check_predictions(final_mod, type = "discrete_interval", check_range = TRUE)

testZeroInflation(final_mod) # No remaining issues with zero-inflation

testOutliers(final_mod) # outlier test looks good

# No within-group deviations from uniformity
# Levene Test for homogeneity of variance looks good
final_mod_simres <- simulateResiduals(final_mod)
testCategorical(final_mod_simres, catPred = na.omit(spat2)$area)
testCategorical(final_mod_simres, catPred = na.omit(spat2)$year)

# Efron's pseudo r-squared
r2_efron(final_mod)


# Post-hoc contrasts (tables) --------------------------------------------------------

names <- c(
  "Year" = "year",
  "Area" = "area",
  "SE" = "std.error",
  "P-value adj." = "p.value",
  "P-value adj." = "adj.p.value",
  "P-value adj." = "p.value.adj",
  "Odds ratio" = "odds.ratio",
  "Odds ratio" = "ratio"
)

create_supp_gt_avg_comps <- function(mod, variable, by = TRUE) {
  marginaleffects::avg_comparisons(
    mod,
    variables = variable,
    by = by,
    p_adjust = "bonferroni",
    re.form = NA,
  ) %>%
    select(any_of(c("term", "contrast", "year", "area", "estimate", "std.error",
                    "statistic", "p.value", "conf.low", "conf.high"))) %>%
    dplyr::rename(any_of(names)) %>%
    gt() %>%
    fmt_number(decimals = 4) %>%
    sub_small_vals(threshold = 0.001) %>%
    tab_style(
      locations = cells_column_labels(columns = any_of(
        c("term", "contrast", "year", "area", "estimate", "statistic",
          "conf.low", "conf.high"))),
      style = cell_text(transform = "capitalize")) %>%
    tab_style(
      locations = cells_body(columns = "term"),
      style = cell_text(transform = "capitalize"))
}


avg_predictions(final_mod, by = "area", re.form = NA)

p1 <- plot_predictions(
  final_mod,
  by = "year",
  re.form = NA,
  vcov = TRUE) + mytheme + labs(x = "Year", y = "Predicted # spat per bag")

avg_predictions(final_mod, by = "year", re.form = NA)

p2 <- plot_predictions(
  final_mod,
  by = "area",
  re.form = NA,
  vcov = TRUE) + mytheme + labs(x = "Area", y = "Predicted # spat per bag") +
  scale_x_discrete(limits = c("closure", "south", "north"), labels = c("south" = "South", "closure" = "Closure", "north" = "North"))

(p1 + p2) & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 16, hjust = -1, vjust = -1),
        plot.margin = margin(20,30,5,10)) & plot_annotation(tag_levels = "A")


# Area, across all years
s1 <- create_supp_gt_avg_comps(final_mod, variable = list(area = "pairwise"))

# Year, across all areas
s2 <- create_supp_gt_avg_comps(final_mod, variable = list(year = "pairwise"))

# Area contrasts by year
s3 <- create_supp_gt_avg_comps(final_mod, variable = list("area" = "pairwise"),
                               by = "year")

# Year contrasts by area
s4 <- create_supp_gt_avg_comps(final_mod, variable = list("year" = "pairwise"),
                               by = "area")

supp_gts <- gt_group(s1, s2, s3, s4)
supp_gts

