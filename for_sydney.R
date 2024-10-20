

library(tidyverse)
library(janitor)
library(mgcv)
library(broom)

dat <- iris %>% clean_names()

mod1 <- lm(data = dat, sepal_length ~ sepal_width)
mod2 <- lm(data = dat, sepal_length ~ sepal_width + petal_length)
mod3 <- lm(data = dat, sepal_length ~ sepal_width * petal_length)

bind_rows(list(
  DEI = tidy(mod1),
  DMC = tidy(mod2),
  MSF = tidy(mod3)
), .id = "Hatchery")

# example with GAMs
mod2 <- gam(data = dat, sepal_length ~ s(sepal_width))
mod3 <- gam(data = dat, sepal_length ~ s(sepal_width) + s(petal_length))

tidy(mod1) %>%
  bind_rows(tidy(mod2), .id = "model")
