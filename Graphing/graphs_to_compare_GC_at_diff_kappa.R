setwd("C://Users//insan_000//Documents//summer 2018//simulation data//bigg data")
library(tidyverse)

dfr = read_csv("all_averages_for_SNP_sim_1000.csv")
dfr$Kappa = as.factor(dfr$Kappa)
dfr$`GC%` = as.factor(dfr$`GC%`)

dfr.k1 <- subset(dfr,dfr$Kappa == 1)
dfr.k2 <- subset(dfr,dfr$Kappa == 2)
dfr.k3 <- subset(dfr,dfr$Kappa == 3)

ggplot(dfr,aes(x=`Mutations on Each Strain`,y=`Average CMs`,color=`GC%`)) + geom_smooth() + theme_classic() + ggtitle("CMs vs total mutations for variant Kappa") + facet_wrap(~`Kappa`)

aovmodk1 <- aov(dfr.k1$`Average CMs`~dfr.k1$`GC%`,data=dfr.k1)
aovmodk2 <- aov(dfr.k2$`Average CMs`~dfr.k2$`GC%`,data=dfr.k2)
aovmodk3 <- aov(dfr.k3$`Average CMs`~dfr.k3$`GC%`,data=dfr.k3)

summary(aovmodk1)
summary(aovmodk2)
summary(aovmodk3)

TukeyHSD(aovmodk1)
TukeyHSD(aovmodk2)
TukeyHSD(aovmodk3)

pairwise.t.test(dfr.k1$`Average CMs`,dfr.k1$`GC%`,p.adjust.method="b")
pairwise.t.test(dfr.k2$`Average CMs`,dfr.k2$`GC%`,p.adjust.method="b")
pairwise.t.test(dfr.k3$`Average CMs`,dfr.k3$`GC%`,p.adjust.method="b")

ggplot(dfr,aes(x=`GC%`,y=`Average CMs`)) + geom_boxplot()

kappaaov <- aov(dfr$`Average CMs`~dfr$Kappa,data=dfr)
summary(kappaaov)
TukeyHSD(kappaaov)
pairwise.t.test(dfr$`Average CMs`,dfr$`Kappa`,p.adjust.method = "b")
