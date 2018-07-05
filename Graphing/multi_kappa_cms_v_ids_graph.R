setwd("C://Users//wzbillin//Documents//cm vs id graph")
library(tidyverse)


#import simulated data and cut to only GC of 50%
sim.dfr <- read_csv("all_identity_sim_data_1000_allGCs.csv", 
                    col_types = cols(`GC%` = col_factor(levels = c("0", 
                                                                   "0.1", "0.2", "0.3", "0.4", "0.5", 
                                                                   "0.6", "0.7", "0.8", "0.9", "1")), 
                                     kappa = col_factor(levels = c("0.5", 
                                                                   "1", "2", "3", "4", "5", "6"))))
sim.dfr.gc5 = subset(sim.dfr,sim.dfr$`GC%` == 0.5)

#plot the data, separated by kappa
ggplot(
  data = sim.dfr.gc5,
  mapping = aes(x = `Average ID%`,
                y = `Average c`,
                color = `kappa`)
  ) +
  #adds smooth curve for each kappa
  geom_smooth(
    method = "loess",
    linetype = 1) +
  #adds axis labels and title
  labs(title = "",
       x = "this is a filler label",
       y = "this label is filler also") +
  #removes gridlines  
  theme_classic() +
  #edits theme elements e.g. axis mark sizes
  theme(
    axis.text.x = element_text(
      face = "bold",
      size = 15),
    axis.text.y = element_text (
      face = "bold",
      size = 15),
    text = element_text (
      size = 20),
    legend.position = "none")
  )