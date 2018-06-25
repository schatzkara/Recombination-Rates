setwd("C://Users//insan_000//Documents//summer 2018//simulation data//bigg data")
library(tidyverse)

dfr = read_csv("all_averages_for_ID_sim_1000.csv")
dfr$Kappa = as.factor(dfr$Kappa)
dfr$`GC%` = as.factor(dfr$`GC%`)

ggplot(dfr,aes(x=`Average ID%`,y=`Average CMs`,color=`GC%`)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for variant Kappa") + facet_wrap(~`Kappa`)
