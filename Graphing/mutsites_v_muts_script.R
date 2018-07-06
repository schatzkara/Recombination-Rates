setwd("C://Users//wzbillin//Documents//mut sites")
library(tidyverse)


#import simulated data and cut to only GC of 50%
dfr <- read_csv("mutation_sites_sim_average_1000_0.001_300_00.5_00.5.csv", 
                                                                col_types = cols(kappa = col_factor(levels = c("0.5")), 
                                                                                 phi = col_factor(levels = c("0.5"))))
x = 1:300
y = 1:300

#plot the data, separated by kappa
ggplot(
  data = dfr,
  mapping = aes(x = Mutations,
                y = `Mutation Sites`)) +
  #adds smooth curve for each kappa
  geom_smooth(
    method = "loess",
    linetype = 1) +
  geom_line(
   mapping = aes(
     x = x,
     y = y),
   linetype = 2
  ) +
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