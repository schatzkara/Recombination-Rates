setwd("C://Users//insan_000//Desktop")
library(tidyverse)

dfr = read_csv("all_averages_for_id_sim_1000.csv")
dfr$Kappa = as.factor(dfr$Kappa)
dfr$`GC%` = as.factor(dfr$`GC%`)

dfr.k1 <- subset(dfr,dfr$Kappa == 1)
dfr.k2 <- subset(dfr,dfr$Kappa == 2)
dfr.k3 <- subset(dfr,dfr$Kappa == 3)

ggplot(dfr,aes(x=`Mutations on Each Strain`,y=`Average CMs`,color=`GC%`)) + 
  geom_smooth() + 
  theme_classic() + 
  ggtitle("CMs vs total mutations for variant Kappa") + 
  facet_wrap(~`Kappa`) +
  labs(x="Muts.",y="conv. muts.")
ggplot(dfr,aes(x=`Mutations on Each Strain`,y=`Average CMs`,color=`Kappa`)) + 
    geom_smooth(method="loess") + 
    theme_classic() + 
    facet_wrap(~`GC%`)

aovmodk1 <- aov(dfr.k1$`Average CMs`~dfr.k1$`GC%`,data=dfr.k1)
aovmodk2 <- aov(dfr.k2$`Average CMs`~dfr.k2$`GC%`,data=dfr.k2)
aovmodk3 <- aov(dfr.k3$`Average CMs`~dfr.k3$`GC%`,data=dfr.k3)

aovmod.allk <- aov(dfr$`Average CMs`~dfr$`GC%`,data=dfr)

summary(aovmodk1)
summary(aovmodk2)
summary(aovmodk3)

summary(aovmod.allk)
TukeyHSD(aovmod.allk)

TukeyHSD(aovmodk1)
TukeyHSD(aovmodk2)
TukeyHSD(aovmodk3)

pairwise.t.test(dfr.k1$`Average CMs`,dfr.k1$`GC%`,p.adjust.method="b")
pairwise.t.test(dfr.k2$`Average CMs`,dfr.k2$`GC%`,p.adjust.method="b")
pairwise.t.test(dfr.k3$`Average CMs`,dfr.k3$`GC%`,p.adjust.method="b")

ggplot(data=dfr.k3) +
  geom_smooth(aes(x=`Mutations on Each Strain`,y=`Average CMs`,color=`GC%`))

ggplot(dfr,aes(x=`GC%`,y=`Average CMs`)) + geom_boxplot()

pairwise.t.test(dfr$`Average ID%`,dfr$`GC%`,p.adjust="b")

kappaaov <- aov(dfr$`Average CMs`~dfr$Kappa,data=dfr)
summary(kappaaov)
TukeyHSD(kappaaov)
pairwise.t.test(dfr$`Average CMs`,dfr$`Kappa`,p.adjust.method = "b")

k3.g0.qm <- lm(
  formula=`Average CMs`~(`Mutations on Each Strain`+(`Mutations on Each Strain`)^2),
  data = dfr.k3.g0) 
summary(k3.g0.qm)

k3.g1.qm <- lm(
  formula=`Average CMs`~(`Mutations on Each Strain`+(`Mutations on Each Strain`)^2),
  data = dfr.k3.g1) 
summary(k3.g1.qm)

anova(k3.g0.qm,k3.g1.qm)

dfr.k3.g0 = subset(dfr.k3,dfr.k3$`GC%` == 0)
#dfr.k3.g1 = subset(dfr.k3,dfr.k3$`GC%` == 0.1)
#dfr.k3.g2 = subset(dfr.k3,dfr.k3$`GC%` == 0.2)

x.g0 <- dfr.k3.g0$`Mutations on Each Strain`
#x.g1 <- dfr.k3.g1$`Mutations on Each Strain`
#x.g2 <- dfr.k3.g2$`Mutations on Each Strain`
y.g0 <- dfr.k3.g0$`Average CMs`
#y.g1 <- dfr.k3.g1$`Average CMs`
#y.g2 <- dfr.k3.g2$`Average CMs`

ggplot() +
  geom_point(aes(x=x.g0,y=y.g0,color="GC = 0")) +
  geom_point(aes(x=x.g1,y=y.g1,color="GC = 10")) +
  geom_point(aes(x=x.g2,y=y.g2,color="GC = 20"))

t.test(y.g0,y.g1)
t.test(y.g0,y.g2)

setwd("D://NEW//SNP Sim")
dfr.k3.g0 <- read_csv("averages_for_SNP_sim_1000_300_0.0_3.0_0.5.csv")

x.g0 <- dfr.k3.g0$`Mutations on Each Strain`
y.g0 <- dfr.k3.g0$`Average CMs`

ggplot() +
  geom_point(aes(x=dfr.k3.g0$`Mutations on Each Strain`,y=dfr.k3.g0$`Average CMs`,color="0")) +
  geom_point(aes(x=dfr.k3$`Mutations on Each Strain`,y=dfr.k3$`Average CMs`,color=dfr.k3$`GC%`)) 
  
