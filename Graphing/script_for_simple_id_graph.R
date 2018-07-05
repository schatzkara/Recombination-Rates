setwd("C://Users//wzbillin//Documents//cm vs id graph")
library(tidyverse)


#import simulated data
sim.dfr <- read_csv("all_identity_sim_data_1000_allGCs.csv", 
                col_types = cols(`GC%` = col_factor(levels = c("0", 
                                                               "0.1", "0.2", "0.3", "0.4", "0.5", 
                                                               "0.6", "0.7", "0.8", "0.9", "1")), 
                                 kappa = col_factor(levels = c("0.5", 
                                                               "1", "2", "3", "4", "5", "6"))))
sim.dfr.gc5 = subset(sim.dfr,sim.dfr$`GC%` == 0.5)
sim.dfr.gc5.khalf = subset(sim.dfr.gc5,sim.dfr.gc5$kappa == 0.5)


#import model data
mod.dfr <- read_csv("all_model_c_vs_id_data.csv", 
                    col_types = cols(kappa = col_factor(levels = c("0.5", 
                                                                   "1", "2", "3", "4", "5", "6"))))
mod.dfr.khalf = subset(mod.dfr,mod.dfr$kappa == 0.5)

#plot the data
ggplot(sim.dfr.gc5.khalf,
       mapping = aes(x = `Average ID%`,y = `Average c`)) + 
  # adds a curve for the model
  geom_smooth(mapping = aes(x = mod.dfr.khalf$`E[ID%]`,
                            y = mod.dfr.khalf$`E[c]`),
              method = "loess",
              color = "#ff5c5c", #change to #ffffff for simple graph
              linetype = 1) +
  # adds a curve for the sim data
  geom_smooth(mapping = aes(x = `Average ID%`,
                            y = `Average c`),
              method = "loess",
              color = "#376595",
              linetype = 1) + 
  #changes the labels on the plot
  labs(title = "",
      x = "this is a filler label",
      y = "this label is filler also") +
  #removes gridlines, makes everything white
  theme_classic() +
  #adjusts font sizes and can do other things.
  theme(
    axis.text.x = element_text(
      face = "bold",
      size = 15),
    axis.text.y = element_text (
      face = "bold",
      size = 15),
    text = element_text (
      size = 20
    )) 
