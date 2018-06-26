#this script imports the large CSV and makes graphs.
setwd("C://Users//insan_000//Documents//summer 2018//simulation data//bigg data")
library(tidyverse)

dfr = read_csv("all_averages_for_SNP_sim_1000.csv")
dfr$Kappa = as.factor(dfr$Kappa)
dfr$`GC%` = as.factor(dfr$`GC%`)

dfrset0 = subset(dfr, dfr$`GC%` == 0)

dfrset0k1 = subset(dfrset0,dfrset0$`Kappa` == 1)
x <- dfrset0k1$`Mutations on Each Strain`
y <- dfrset0k1$`Average CMs`
x2 <- x^2

ythr <- ((x/1000)^2*(1.5/4))*1000
yq <- y/ythr
yr <- y-ythr

linear.model <- lm(y~x)
quadratic.model <- lm(y~(x2+x))

summary(linear.model)
summary(quadratic.model)

#ggplot(dfrset0k1,aes(x=x,y=y)) + geom_point(aes(color="data")) + stat_smooth(method="lm",formula= y~x+I(x^2),linetype=2,aes(color="quadratic model")) + stat_smooth(method="lm",linetype=2,se=F,aes(color="Linear model")) + theme_classic() + labs(title="Modeling convergent mutations vs. total mutations",x="Total mutations",y="Convergent mutations")
ggplot(dfrset0k1,aes(x=x,y=y)) + geom_point(aes(color="data")) + geom_smooth(aes(y=ythr,color="model"),linetype=2)
qplot(data=dfrset0k1,x=x,y=yr)

#kappa=2
dfrset0k2 = subset(dfrset0,dfrset0$`Kappa` == 2)
x <- dfrset0k2$`Mutations on Each Strain`
y <- dfrset0k2$`Average CMs`                                                                                                                                 
ymodel <- ((x/1000)^2*(1/2))*1000

ggplot(dfrset0k2,aes(x=x,y=y)) + geom_point(aes(color="data")) + geom_line(aes(y=ymodel,color="model"),linetype=2)
