setwd("C://Users//insan_000//Desktop//")

library(gdata)
library(grid)
library(tidyverse)

#set up dataframe
mutdfr <- read_csv("all_averages_for_id_sim_1000.csv")
mutdfr <- subset(mutdfr,mutdfr$`GC%` == 0.5)
mutdfr$Kappa <- as.factor(mutdfr$Kappa)
moddfr <- read_csv("all_model_c_vs_id_data.csv")
moddfr <- subset(moddfr,
                 moddfr$kappa == 1 | 
                 moddfr$kappa == 2 |
                 moddfr$kappa == 3)
moddfr$kappa <- as.factor(mutdfr$Kappa)

mutdfr.k1 <- subset(mutdfr,mutdfr$Kappa == 1)
mutdfr.k2 <- subset(mutdfr,mutdfr$Kappa == 2)
mutdfr.k3 <- subset(mutdfr,mutdfr$Kappa == 3)
moddfr.k1 <- subset(moddfr,moddfr$kappa == 1)
moddfr.k2 <- subset(moddfr,moddfr$kappa == 2)
moddfr.k3 <- subset(moddfr,moddfr$kappa == 3)

p.k1 <- ggplot() +
  geom_point(
    mapping = aes(x = mutdfr.k1$`Average ID%`,
                  y = mutdfr.k1$`Average CMs`,
                  color = "Simulated data")
    ) +
  geom_smooth(
    mapping = aes(x = moddfr.k1$`E[ID%]`,
                  y = moddfr.k1$`E[c]`,
                  color = "Predicted data"),
    linetype = 2,
    method = "loess"
  ) +
  labs(
    title = "Identity percentage vs. total mutations for Kappa = 1",
    x = "Identity percentage",
    y = "Number of convergent mutations",
    legend = "Data series"
  ) +
  theme_classic()

p.k2 <- ggplot() +
  geom_point(
    mapping = aes(x = mutdfr.k2$`Average ID%`,
                  y = mutdfr.k2$`Average CMs`,
                  color = "Simulated data")
  ) +
  geom_smooth(
    mapping = aes(x = moddfr.k2$`E[ID%]`,
                  y = moddfr.k2$`E[c]`,
                  color = "Predicted data"),
    linetype = 2,
    method = "loess"
  ) +
  labs(
    title = "Identity percentage vs. total mutations for Kappa = 2",
    x = "Identity percentage",
    y = "Number of convergent mutations",
    legend = "Data series"
  ) +
  theme_classic()

p.k3 <- ggplot() +
  geom_point(
    mapping = aes(x = mutdfr.k3$`Average ID%`,
                  y = mutdfr.k3$`Average CMs`,
                  color = "Simulated data")
  ) +
  geom_smooth(
    mapping = aes(x = moddfr.k3$`E[ID%]`,
                  y = moddfr.k3$`E[c]`,
                  color = "Predicted data"),
    linetype = 2,
    method = "loess"
  ) +
  labs(
    title = "Identity percentage vs. total mutations for Kappa = 3",
    x = "Identity percentage",
    y = "Number of convergent mutations",
    legend = "Data series"
  ) +
  theme_classic() 

print(p.k1)
print(p.k2)
print(p.k3)

x.k3 <- mutdfr$`Average ID%`
y.k3 <- mutdfr$`Average CMs`
x.k3.sq <- (x.k3)^2
regress.k3 <- lm(
  formula = y.k3~(x.k3+I(x.k3.sq))
)

summary(regress.k3)

# ggplot(data=mutdfr,
#        mapping = aes(
#          x = `Mutations on Each Strain`,
#          y = `Average ID%`
#        )) +
#   geom_point()
