
library(FITfileR)
library(tidyverse)
library(lubridate)

nats_data <- readFitFile("~/Downloads/nats.fit")
nats_data

nats_records <- records(nats_data)
nats_records

nats_allrecords <- records(nats_data) %>% 
  bind_rows() %>% 
  arrange(timestamp) 

nats2 <- nats_allrecords %>% 
  select(timestamp, distance, heart_rate, power) %>% 
  mutate(time = as.character(timestamp)) %>% 
  mutate(time = substr(time, 15, 19)) %>% 
  mutate(time = ms(time))


ggplot(nats2, aes(y = time, x = distance)) +
  geom_line()

nats2 %>% filter(round(distance, -1) == 1000)
nats2 %>% filter(round(distance, -1) == 2000)
nats2 %>% filter(round(distance, -1) == 3000)
nats2 %>% filter(round(distance, -1) == 4000)
nats2 %>% filter(round(distance, -1) == 5000)
nats2 %>% filter(round(distance, -1) == 6000) 


# library(tabulapdf)
# 
# tables_list <- paste0("~/Downloads/MYiClubOnline", c(0:7))
# tables_list <- paste0(tables_list, ".pdf")
# 
# f1 <- extract_tables(tables_list[1])[[1]]
# f2 <- extract_tables(
#   tables_list[2],
#   pages = c(1, 2),
#   guess = FALSE,
#   area = list(
#     c(361.33333, 19.33333, 766.66667, 598),
#     c(14.66667, 24.66667, 481.33333, 595.33333)
#   )
# )[[1]]
# 
# f3 <- extract_tables(
#   tables_list[3],
#   pages = c(1, 2),
#   guess = FALSE,
#   area = list(
#     c(401.3333, 22, 780, 595.3333),
#     c(12, 24.66667, 638.66667, 595.33333)
#   )
# )[[1]]
# 
# f4 <- extract_tables(
#   tables_list[4],
#   pages = c(1, 2),
#   guess = FALSE,
#   area = list(c(332, 22, 774.6667, 590), c(12, 22, 350.6667, 598))
# )[[1]]
# 
# f5 <- extract_tables(
#   tables_list[5],
#   pages = c(1, 2),
#   guess = FALSE,
#   area = list(
#     c(334.6667, 19.3333, 774.6667, 590),
#     c(9.333333, 24.66667, 57.33333, 595.33333)
#   )
# )[[1]]
# 
# f6 <- extract_tables(tables_list[6])[[1]]
# f7 <- extract_tables(tables_list[7])[[1]]
# f8 <- extract_tables(tables_list[8], guess = FALSE,
#                      area = list(c(316, 22, 577.3333, 600.6667)))[[1]]
# 
# all_tables <- bind_rows(f1, f2, f3, f4, f5, f6, f7, f8)
# 
# write_csv(all_tables, "planetfitness.csv")

all_tables <- read.csv("./planetfitness.csv") %>% janitor::clean_names()

pf <- all_tables %>% 
  select(-home_club) %>% 
  mutate(date = mdy(date)) %>% 
  arrange(date)

unique(pf$location)
