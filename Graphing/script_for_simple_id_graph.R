setwd("C://Users//wzbillin//Documents//cm vs id graph")
library(tidyverse)

#dfr = read_csv("all_identity_sim_data_1000_allGCs.csv")
#dfr$kappa = as.factor(dfr$kappa)

dfr <- read_csv("all_identity_sim_data_1000_allGCs.csv", 
                col_types = cols(`GC%` = col_factor(levels = c("0", 
                                                               "0.1", "0.2", "0.3", "0.4", "0.5", 
                                                               "0.6", "0.7", "0.8", "0.9", "1")), 
                                 kappa = col_factor(levels = c("0.5", 
                                                               "1", "2", "3", "4", "5", "6"))))

dfr.gc5 = subset(dfr,dfr$`GC%` == 0.5)
dfr.gc5.khalf = subset(dfr.gc5,dfr.gc5$kappa == 0.5)

#dfr$`GC%` = transform(dfr,`GC%` = `GC%`*10)
#dfr$`GC%` = as.factor(dfr$`GC%`)

# dfrset0 = subset(dfr, dfr$`GC%` == 0)
# 
# dfrset0k1 = subset(dfrset0,dfrset0$`Kappa` == 1)

# x2 <- x^2

# linear.model <- lm(y~x)
# quadratic.model <- lm(y~(x2+x))
# 
# summary(linear.model)
# summary(quadratic.model)

ggplot(dfr.gc5.khalf,
       mapping = aes(x=`Average ID%`,y=`Average c`)) + 
  geom_smooth(mapping = aes(x=`Average ID%`,
                  y=`Average c`)) + 
  #stat_smooth(method="lm",
  #            formula= y~x+I(x^2),
  #            linetype=2,
  #            aes(color="quadratic model")) + 
  #stat_smooth(method="lm",
  #            linetype=2,
  #            se=F,
  #            aes(color="Linear model")) + 
  theme_classic() + 
  labs(title="",
      x = "",
      y = "")
  #facet_wrap(~`GC%`)

ggsave("simple_cm_vs_id_k_0.5_gc_0.5.png",device="png")
