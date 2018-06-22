#this script imports the large CSV and makes graphs.
setwd("C://Users//insan_000//Documents//summer 2018//simulation data//bigg data")
library(tidyverse)

dfr = read_csv("all_averages_for_ID_sim_1000.csv")
dfr$Kappa = as.factor(dfr$Kappa)
dfr$`GC%` = as.factor(dfr$`GC%`)

dfrset0 = subset(dfr, dfr$`GC%` == 0)

dfrset0k1 = subset(dfrset0,dfrset0$`Kappa` == 1)
x <- dfrset0k1$`Average ID%`
y <- dfrset0k1$`Average CMs`
x2 <- x^2

linear.model <- lm(y~x)
quadratic.model <- lm(y~(x2+x))

summary(linear.model)
summary(quadratic.model)

ggplot(dfrset0k1,aes(x=x,y=y)) + geom_point(aes(color="data")) + stat_smooth(method="lm",formula= y~x+I(x^2),linetype=2,aes(color="quadratic model")) + stat_smooth(method="lm",linetype=2,se=F,aes(color="Linear model")) + theme_classic() + labs(title="Modeling Convergent mutations vs. ID%",x="ID%",y="Convergent mutations")

ggsave("models_cms_v_id.png",device="png")
