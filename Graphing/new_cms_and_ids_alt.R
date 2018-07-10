setwd("C://Users//insan_000//Desktop//")

library(gdata)
library(grid)
library(tidyverse)

#set up dataframe
mutdfr <- read_csv("all_model_c_vs_id_data.csv")

mutdfr.subk <- subset(mutdfr,mutdfr$kappa == 0.5)
mutdfr.subk.subp <- subset(mutdfr.subk,mutdfr.subk$phi == 0)

n <- 299

modeldfr <- read_csv("model_c_vs_id_data_1000_300_0.5_0.0.csv")
modeldfr <- modeldfr[-300,]

x <- 1:n
y <- 

loopindex <- 1:n

for (k in loopindex) {
  jerry = ym[,k]
  #ev = evs[[k]]
  
  plot <- ggplot() +
    geom_point(mapping = aes(x = x, y = jerry)) +
    labs(x = "Number of Simulations", 
         y = "Average convergent mutations", 
         title = paste("CMs v. ID for mu = ",mus[[k]],",kappa =",kappas[[k]],", phi = ",phis[[k]],", L = 1000.png",sep="") 
    ) +
    theme_classic()
  #annotate("text",x = 750, y = (ev*IQR(y)), label = paste("E[c] = ",ev,sep=""), size = 4) 
  
  print(plot)
  #fn = paste("MC_Sim_Graph_mu_",mus[[k]],"_kappa_",kappas[[k]],"_phi_",phis[[k]],"_L_1000.png",sep="")
  #ggsave(fn,device="png")
  
}
