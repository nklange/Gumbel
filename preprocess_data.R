library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(readr)
library(tibble)
library(ggplot2)

loaddata <- function(expt){
  # Get a list of all csvs in rawdata directory
  csvlist <- list.files(path=paste0(getwd(), "/SpantonBerry2021Data/rawdata/", expt, "/"))
  # Import all data files into a list without prints
  bind_rows(map(csvlist, ~ read_csv(paste0(getwd(), "/SpantonBerry2021Data/rawdata/", expt, "/", .x)) ) )
}

experimentfile <- bind_rows(loaddata("Experiment1") %>% mutate(exp = "SB2021_e1"),
                     loaddata("Experiment2") %>% mutate(exp = "SB2021_e2")) %>%
  mutate(id = paste0(exp,"_",subject_nr))

selecttest <- experimentfile %>% filter(phase == "test") %>%
  select(exp,id,condition,conditionorders,test_word,oldnew,response)

demographics <- selecttest %>% group_by(exp,id) %>%
  dplyr::slice_tail(n=3) %>% select(-oldnew)

testphase <- selecttest %>% group_by(exp,id) %>% dplyr::slice(1:480) %>% mutate(trial = c(1:480)) %>%
  mutate(oldnew = ifelse(oldnew == 0,"New","Old"))

# A: High Strength, high EV
# B: High Strength, low EV
# C: Low strength, high EV
# D: Low Strength, low EV

