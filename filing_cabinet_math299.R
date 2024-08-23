

# *A)* The approximate percentage of variation in gambling expenditure explained by these predictors is `r paste0(round(modsum$r.squared*100, 1), "%")`.

tibble(x2.5 = confint(full_mod)[,1], x97.5=confint(full_mod)[,2]) %>% 
  mutate(across(where(is.numeric), ~round(.x, 3))) %>% 
  add_column(term = tidy(full_mod)$term, .before="x2.5")