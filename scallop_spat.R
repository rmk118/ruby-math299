#---------------------------------------------------#
# Modeling scallop spat count data
# Ruby Krasnow
# Last updated Nov 9, 2024
#---------------------------------------------------#

## List of packages required:
packages <- c("tidyverse", "ggthemes", "gt", "janitor", "glmmTMB", "performance", "DHARMa", "rstatix", "broom.mixed", "marginaleffects", "emmeans")

# Load packages into session
lapply(packages, require, character.only = TRUE)
rm(packages)

set.seed(123)

mytheme <- theme_light()+ #define custom theme for ggplots
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, l = 0)),
        text=element_text(size=15))

# Bags recovered ----------------------------------------------------------

## Setup -------------------------------------------------------------------

perc_rec <- tribble(
  ~year,	~bags,	~lines, ~status,
  2013,	25,	5,	"open",	
  2014,	20,	4,	"open",
  2015,	4,	1,	"open",
  2016,	14,	3,	"open",	
  2017,	8,	1.5,	"open",
  2013,	4,	1,	"closed",
  2014,	5,	1,	"closed",
  2015,	14,	3,	"closed",
  2016,	6,	1.25,	"closed",
  2017,	5,	1,	"closed"
)

perc_rec <- perc_rec %>% 
  mutate(deployed = if_else(status == "open", 45, 15)) %>% 
  mutate(perc = bags/deployed,
         inv_dep = 1/deployed)

spat_orig <- read.csv("./spat_orig.csv") %>% clean_names() %>% rename(status = status_open_closed, line = bag_line, count = count_per_bag) %>% select(-site) %>% filter(year != 2018) %>% mutate(bag_pos = case_when(
  bag_position == 2 ~ "b",
  bag_position == 3 ~ "c",
  bag_position == 4 ~ "d",
  bag_position == 5 ~ "e",
  .default = bag_position
)) %>% mutate(
  bag_position = as_factor(bag_position),
  line = as_factor(line),
  year = as_factor(year),
  area = as_factor(area),
  status = as_factor(status)
)

# spat_orig %>% count(year, area, line)
# spat_orig %>% select(year, area, line) %>% unique() %>% count(year, area)

# full_grid <- expand_grid(
#   line = c(1:12),
#   bag_pos = c("a", "b", "c", "d", "e"),
#   year = c(2013:2017)) %>% 
#   mutate(area = case_when((line <= 3) ~ "n",
#                           (line > 3 & line < 7) ~ "c",
#                           .default = "s")) %>% 
#   mutate(line = paste0(line, area, year))

expand_grid(area = c(1, 2, 3, 4), line = c(1:3)) %>% 
  expand_grid(bag_pos = c(-0.2, -0.1, 0, 0.1, 0.2)) %>% ggplot() + 
  geom_tile(aes(x = line,y = area), fill = "white", color = "black") + 
  geom_point(aes(x = line, y = area + bag_pos, color = as.factor(bag_pos)))

## Modeling -------------------------------------------------------------------

perc_lm <- lm(data = perc_rec, perc ~ status)
# perc_lm_area <- lm(data = perc_rec_area, perc ~ area)

check_model(perc_lm)
# check_residuals(perc_lm_area)

# shapiro.test(perc_lm_area$residuals)

weighted1 <- glmmTMB(data = perc_rec, bags ~ status, weights = deployed)
weighted2 <- glmmTMB(data = perc_rec, bags ~ status, weights = inv_dep)
no_weight <- glmmTMB(data = perc_rec, bags ~ status)
null_mod <- glmmTMB(data = perc_rec, bags ~ 1)

# weighted1_area <- glmmTMB(data = perc_rec_area, bags ~ area, weights = deployed)
# no_weight_area  <- glmmTMB(data = perc_rec_area, bags ~ area)
# null_mod_area  <- glmmTMB(data = perc_rec_area, bags ~ 1)

# counts_weighted_area <- glmmTMB(data = perc_rec_area,
#                             perc ~ area,
#                             weights = deployed,
#                             family = binomial())
# 
# bbmle::AICtab(weighted1_area,
#               no_weight_area,
#               null_mod_area,
#               counts_weighted_area)

r2(counts_weighted)
check_residuals(counts_weighted)
shapiro.test(residuals(counts_weighted))
summary(counts_weighted)
Anova(counts_weighted)

(em1 <- emmeans(counts_weighted, "status", type="response"))
plot(em1) + theme_light() + labs(y = "Status", x = "Estimated marginal means for % bags recovered")

m1 <- DHARMa::simulateResiduals(counts_weighted)
plot(m1)

perc_rec %>% wilcox_test(perc ~ status)


# Spat abundance ----------------------------------------------------------

## Setup -------------------------------------------------------------------

spat <- read.csv("./spat_clean.csv") %>%
  janitor::clean_names() %>%
  select(-c(x, site)) %>% filter(year != 2018) %>% mutate(
    bag_position = as_factor(bag_position),
    line = as_factor(line),
    year = as_factor(year),
    area = as_factor(area),
    status = as_factor(status)
  )

## Modeling -------------------------------------------------------------------


# Fixed effect of interest is status (open vs closed) "rank-deficient" when area
# and status in same model. Only 1 area is closed, cannot independently estimate
# effects. Too few areas to add as random effect, so included as fixed effect
# (area) rather than random intercept (1|area) or slope. Cannot assess
# interaction terms with year and area because not all combinations of year and
# area represented due to gear loss (no lines recovered from south of closed
# area in 2015, despite half of all lines being set in that area)

########## "rank-deficient" formulas:
# full = count ~ year + status + area
# modB4 = count ~ area + status
# modE = count ~ status + area + bag_position + year
# modB2 = count ~ year * area --> no south lines in 2015

########## Convergence issue formulas:
# modD = count ~ status + year + (1|bag_position)

########## Other formulas including status:
# modA3 = count ~ status
# modB3 = count ~ year * status
# modB4 = count ~ bag_position * status
# modC3 = count ~ year + status
# modC4 = count ~ bag_position + status
# modC6 = count ~ area + status
# modD2 = count ~ status + (1|year),
# modF = count ~ status + year + (1|area))

########## Formulas that include status rather than area - 
# count ~ status
# count ~ year * status
# count ~ bag_position * status
# count ~ year + status
# count ~ bag_position + status
# count ~ area + status
# count ~ status + (1|year),
# count ~ status + year + (1|area))

nb_formula <- list(
  modA1 = count ~ year,
  modA2 = count ~ area, 
  modA3 = count ~ bag_position,
  modB1 = count ~ year * bag_position,
  modB2 = count ~ bag_position * area,
  modC1 = count ~ year + bag_position,
  modC2 = count ~ year + area,
  modC3 = count ~ bag_position + area,
  modD1 = count ~ area + (1|year))

# https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf

nb1_models <- tibble(zi_formula, models = purrr::map(zi_formula, \(x)
                                                    glmmTMB(
                                                      #data = spat,
                                                      data = zero_df,
                                                      formula = x,
                                                      family = nbinom2(),
                                                      #family = poisson,
                                                       ziformula = ~1
                                                    ))) %>%
  add_column(id = names(zi_formula), .before = 1)

nb1_models_tidy <- nb1_models %>%
  mutate(tidy_model = map(models, tidy),
         AIC = map(models, AIC),
         resids = map(models, residuals)) %>% 
  unnest(cols = c(AIC)) %>% 
  mutate(resids = map(resids, as.numeric()), 
         normal_p = map(resids, ~ shapiro.test(.x)$p.value)) %>% 
  unnest(cols = c(normal_p)) %>% 
  mutate(logLik = map(models, ~ (logLik(.x)[1]))) %>% 
  unnest(cols = c(logLik)) %>% 
  arrange(AIC)

nb2_models <- tibble(nb_formula, models = purrr::map(nb_formula, \(x)
                                                     glmmTMB(
                                                       data = spat,
                                                       formula = x,
                                                       family = nbinom2(),
                                                       # ziformula = ~1
                                                     ))) %>%
  add_column(id = names(nb_formula), .before = 1)

nb2_models_tidy <- nb2_models %>%
  mutate(tidy_model = map(models, tidy),
         AIC = map(models, AIC),
         resids = map(models, residuals)) %>% 
  unnest(cols = c(AIC)) %>% 
  mutate(resids = map(resids, as.numeric()), 
         normal_p = map(resids, ~ shapiro.test(.x)$p.value)) %>% 
  unnest(cols = c(normal_p)) %>% 
  mutate(logLik = map(models, ~ (logLik(.x)[1]))) %>% 
  unnest(cols = c(logLik)) %>% 
  arrange(AIC)

bbmle::AICtab(nb2_models_tidy$models, mnames = nb2_models_tidy$id)

# do diagnostic checking on model w/ lowest AIC

c2 <- (nb2_models_tidy %>% filter(id =="modC2") %>% pull(models))[[1]]
testResiduals(c2) # equivalent to check_residuals(c2)
testDispersion(c2) # equivalent to check_overdispersion(c2)
check_predictions(c2, type = "discrete_interval", check_range = TRUE, plot = FALSE)
c2_simres <- simulate_residuals(c2)
testCategorical(c2_simres, catPred = na.omit(spat)$area)
testCategorical(c2_simres, catPred = spat$year)
plot(simulateResiduals(c2, n = 1000))

c2b <- glmmTMB(count ~ year + area, data = spat, family = nbinom2(), zi=~1)
testResiduals(c2b)
testDispersion(c2)
check_predictions(c2b, type = "discrete_interval", check_range = TRUE, plot = FALSE)
c2b_simres <- simulate_residuals(c2b)
testCategorical(c2b_simres, catPred = na.omit(spat)$area)
testCategorical(c2b_simres, catPred = spat$year)
plot(simulateResiduals(c2b, n = 1000))
DHARMa::testQuantiles(c2b) # slight quantile deviation detected but not statistically significant

AIC(c2, c2b)

tidy(c2b)
# then use emmeans, ggeffects, or marginaleffects to isolate impact of area

# https://strengejacke.github.io/ggeffects/articles/introduction_randomeffects.html
# https://marginaleffects.com/
# https://www.andrewheiss.com/blog/2022/11/29/conditional-marginal-marginaleffects/
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/model_evaluation.pdf

# AIC orders of magnitude higher when using poisson


spat2 <- read.csv("./spat_filled.csv") %>%
  clean_names() %>%
  rename(count = count_per_bag) %>%
  select(-c(bag_line, status_open_closed)) %>%
  complete(line, year, bag_position) %>%
  mutate(line2 = paste(year, line, sep = "_")) %>%
  filter(!is.na(bag_position)) %>%
  mutate(area = case_when((line < 7) ~ "south",
                          (line > 6 & line < 10) ~ "north", .default = "closure"))  %>%
  mutate(
    bag_position = as_factor(bag_position),
    line2 = as_factor(line2),
    year = as_factor(year),
    area = as_factor(area) )


spat2 %>% mutate(status = if_else(area == "closure", "closed", "open")) %>% 
  filter(!is.na(count)) %>% 
  count(year, status) %>% 
  # mutate(lines = n/5) %>% 
  rename(bags = n) %>% 
  arrange(status)

zero_df <- spat2 %>% mutate( # create version of data set where not-recovered bags are 0
  true_zero = if_else(count!=0 | is.na(count), FALSE, TRUE),
  count = replace_na(count, 0))

zi_formula <- list(
  modA1 = count ~ year,
  modA2 = count ~ area, 
  modA3 = count ~ bag_position,
  modB1 = count ~ year * bag_position,
  modB2 = count ~ area * bag_position,
  modC1 = count ~ year + bag_position,
  modC2 = count ~ year + area,
  modC3 = count ~ bag_position + area,
  modD1 = count ~ area + (1|line2),
  modD2 = count ~ year + (1|line2),
  modD3 = count ~ bag_position + (1|line2))

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

zi_models <- tibble(zi_formula, models = purrr::map(zi_formula, \(x)
                                                     glmmTMB(
                                                       data = spat2,
                                                       formula = x,
                                                       family = nbinom2(),
                                                      # ziformula = ~1
                                                     ))) %>%
  add_column(id = names(zi_formula), .before = 1)

zi_models_tidy <- zi_models %>%
  mutate(tidy_model = map(models, tidy),
         AIC = map(models, AIC),
         resids = map(models, residuals)) %>% 
  unnest(cols = c(AIC)) %>% 
  mutate(resids = map(resids, as.numeric()), 
         normal_p = map(resids, ~ shapiro.test(.x)$p.value)) %>% 
  unnest(cols = c(normal_p)) %>% 
  mutate(logLik = map(models, ~ (logLik(.x)[1]))) %>% 
  unnest(cols = c(logLik)) %>% 
  arrange(AIC)

testmod <- glmmTMB(count ~ year + area + (1|line2) + (1|bag_position), data = spat2, family = nbinom2()) #ziformula = ~ 1)
testmod2 <- glmmTMB(count ~ year + area + (1|line2) + (1|bag_position), 
                    data = spat2, 
                    #data = zero_df, 
                    family = nbinom2(), ziformula = ~ 1)
summary(testmod)
tidy(testmod)

testResiduals(testmod) # same as testUniformity(testmod)
testDispersion(testmod)
check_model(testmod)
check_predictions(testmod, type = "discrete_interval", check_range = TRUE)
testmod_simres <- simulateResiduals(testmod)
testZeroInflation(testmod)
testCategorical(testmod_simres, 
                #catPred = na.omit(spat2)$year)
                catPred = zero_df$year)
testCategorical(testmod_simres, 
                #catPred = na.omit(spat2)$area)
                catPred = zero_df$area)
plot(simulateResiduals(testmod, n = 1000))
plot(simulateResiduals(testmod, n = 1000), quantreg = F)
DHARMa::testQuantiles(testmod)

# Efron's pseudo r-squared
r2_efron(testmod)

Anova(testmod)
# same as
emmeans::joint_tests(testmod) %>%
  select(`model term`, F.ratio, Chisq, p.value) %>%
  dplyr::rename("Term" = "model term",
                "F-ratio" = "F.ratio",
                "P-value" = "p.value")

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

avg_predictions(testmod, variables = "area", re.form = NA)

plot_predictions(
  testmod,
  condition = "year",
  type = "response",
  re.form = NA,
  vcov = TRUE) + mytheme

avg_predictions(testmod, variables = "year", re.form = NA)

plot_predictions(
  testmod,
  condition = "area",
  type = "response",
  re.form = NA,
  vcov = TRUE) + mytheme


# Area, across all years
s1 <- create_supp_gt_avg_comps(testmod, variable = list(area = "pairwise"))

# Year, across all areas
s2 <- create_supp_gt_avg_comps(testmod, variable = list(year = "pairwise")) 

# Area contrasts by year
s3 <- create_supp_gt_avg_comps(testmod, variable = list("area" = "pairwise"), by = "year") 
emmeans(testmod, ~ area | year, type = "response")

# Year contrasts by area
s4 <- create_supp_gt_avg_comps(testmod,variable=list("year"="pairwise"), by="area") 

supp_gts <- gt_group(s1, s2, s3, s4)
supp_gts


## Visualization ---------------------------------------------------------------

ggplot(perc_rec) + geom_bar(aes(x = year, y = perc, color = str_to_sentence(status), fill = str_to_sentence(status)), stat = "identity", position = "dodge", alpha = 0.7, linewidth = 1) + mytheme + scale_color_tableau() + scale_fill_tableau() + labs(x = "Year", y = "% of bags recovered", fill = NULL, color = NULL)

ggplot(spat) + geom_boxplot(aes(x = year, y = count, fill = str_to_sentence(area)), alpha = 0.7) + mytheme + scale_fill_tableau() + labs(x = "Year", y = "# spat per recovered bag", fill = NULL) + scale_y_continuous(labels = scales::comma)


# Hurdle models -----------------------------------------------------------


h_models <- tibble(zi_formula, models = purrr::map(
  zi_formula,
  \(x)
  glmmTMB(
    data = zero_df,
    formula = x,
    family = truncated_nbinom2(),
    ziformula = ~ 1
  )
)) %>%
  add_column(id = names(zi_formula), .before = 1)

h_models_tidy <- h_models %>%
  mutate(tidy_model = map(models, tidy),
         AIC = map(models, AIC),
         resids = map(models, residuals)) %>% 
  unnest(cols = c(AIC)) %>% 
  mutate(resids = map(resids, as.numeric()), 
         normal_p = map(resids, ~ shapiro.test(.x)$p.value)) %>% 
  unnest(cols = c(normal_p)) %>% 
  mutate(logLik = map(models, ~ (logLik(.x)[1]))) %>% 
  unnest(cols = c(logLik)) %>% 
  arrange(AIC)

library(bbmle)
AICtab(testmod, final_mod, logLik = T)
