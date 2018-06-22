#this script imports the large CSV and makes graphs.
setwd("C://Users//insan_000//Documents//summer 2018//simulation data//bigg data")
library(tidyverse)

dfr = read_csv("all_averages_for_SNP_sim_1000.csv")
dfr$Kappa = as.factor(dfr$Kappa)
dfr$`GC%` = as.factor(dfr$`GC%`)

dfrset = subset(dfr, dfr$`GC%` == 0.5)

ggplot(dfrset,aes(x=`Mutations on Each Strain`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic()
