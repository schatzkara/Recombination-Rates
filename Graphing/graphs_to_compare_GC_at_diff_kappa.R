setwd("C://Users//insan_000//Documents//summer 2018//simulation data//bigg data")
library(tidyverse)

dfr = read_csv("all_averages_for_ID_sim_1000.csv")
dfr$Kappa = as.factor(dfr$Kappa)
dfr$`GC%` = as.factor(dfr$`GC%`)

dfr.k1 <- subset(dfr,dfr$Kappa == 1)

ggplot(dfr,aes(x=`Average ID%`,y=`Average CMs`,color=`GC%`)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for variant Kappa") + facet_wrap(~`Kappa`)

aovmod <- aov(dfr.k1$`Average CMs`~dfr.k1$`GC%`,data=dfr.k1)

summary(aovmod)
TukeyHSD(aovmod)
pairwise.t.test(dfr.k1$`Average CMs`,dfr.k1$`GC%`,p.adjust.method="b")

ggplot(dfr,aes(x=`GC%`,y=`Average CMs`)) + geom_boxplot()
