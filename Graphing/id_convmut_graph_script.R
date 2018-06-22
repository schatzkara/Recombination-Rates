#this script imports the large CSV and makes graphs.
setwd("C://Users//insan_000//Documents//summer 2018//simulation data//bigg data")
library(tidyverse)

dfrid = read_csv("all_averages_for_ID_sim_1000.csv")
dfrid$Kappa = as.factor(dfrid$Kappa)
dfrid$`GC%` = as.factor(dfrid$`GC%`)

dfrsub = subset(dfrid, dfrid$`GC%` == 0.5)

ggplot(dfrsub,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic()
