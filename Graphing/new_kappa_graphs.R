setwd("C://Users//insan_000//Desktop")
library(tidyverse)

dfr = read_csv("all_averages_for_SNP_sim_data.csv")
dfr$kappa = as.factor(dfr$kappa)
dfr$`GC%` = as.factor(dfr$`GC%`)
dfr$phi = as.factor(dfr$phi)

ggplot(dfr,aes(x=`Mutations on Each Strain`,y=`Average CMs`,color=`GC%`)) + 
  geom_point() + 
  theme_classic() + 
  ggtitle("CMs vs total mutations for variant Kappa") + 
  facet_wrap(~`kappa`) +
  labs(x="Muts.",y="conv. muts.")
