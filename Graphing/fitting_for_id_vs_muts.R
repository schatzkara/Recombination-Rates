#this script imports the large CSV and makes graphs.
setwd("C://Users//insan_000//Documents//summer 2018//simulation data//bigg data")
library(tidyverse)

dfr = read_csv("all_averages_for_ID_sim_1000.csv")
dfr$Kappa = as.factor(dfr$Kappa)
dfr$`GC%` = as.factor(dfr$`GC%`)

dfrset0 = subset(dfr, dfr$`GC%` == 0.5)

dfrset0k1 = subset(dfrset0,dfrset0$`Kappa` == 1)
x <- dfrset0k1$`Average ID%`
y <- dfrset0k1$`Average CMs`
x2 <- x^2

quadratic.model.k1 <- lm(y~(x2+x))

ggplot(dfrset0k1,aes(x=x,y=y)) +
  geom_point(aes(color="data")) +
  stat_smooth(method="lm",formula= y~x+I(x^2),linetype=2,aes(color="quadratic model")) +
  theme_classic() + labs(title="Modeling Convergent mutations vs. ID%, kappa=1",x="ID%",y="Convergent mutations")

#ggsave("models_cms_v_id.png",device="png")

dfrset0k2 = subset(dfrset0,dfrset0$`Kappa` == 2)
x <- dfrset0k2$`Average ID%`
y <- dfrset0k2$`Average CMs`
x2 <- x^2

quadratic.model.k2 <- lm(y~(x2+x))

summary(quadratic.model)

ggplot(dfrset0k2,aes(x=x,y=y)) +
  geom_point(aes(color="data")) +
  stat_smooth(method="lm",formula= y~x+I(x^2),linetype=2,aes(color="quadratic model")) +
  theme_classic() + 
  labs(title="Modeling Convergent mutations vs. ID%, kappa=2",x="ID%",y="Convergent mutations")

#ggsave("models_cms_v_id.png",device="png")

dfrset0k3 = subset(dfrset0,dfrset0$`Kappa` == 3)
x <- dfrset0k3$`Average ID%`
y <- dfrset0k3$`Average CMs`
x2 <- x^2

quadratic.model.k3 <- lm(y~(x2+x))

ggplot(dfrset0k3,aes(x=x,y=y)) +
  geom_point(aes(color="data")) +
  stat_smooth(method="lm",formula= y~x+I(x^2),linetype=2,aes(color="quadratic model")) +
  theme_classic() +
  labs(title="Modeling Convergent mutations vs. ID%, kappa=3",x="ID%",y="Convergent mutations")

summary(quadratic.model.k1)
summary(quadratic.model.k2)
summary(quadratic.model.k3)
