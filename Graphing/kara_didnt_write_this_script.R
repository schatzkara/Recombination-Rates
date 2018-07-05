setwd("C://Users//wzbillin//Documents//monte carlo garph")
library(gdata)
library(grid)
library(tidyverse)

#set up dataframe
mutdfr <- read_csv("all_mutation_sim_data_1000.csv")
# mutdfr$mu <- as.factor(mutdfr$mu)
# mutdfr$kappa <- as.factor(mutdfr$kappa)
# mutdfr$phi <- as.factor(mutdfr$phi)

# rows = function(x) lapply(seq_len(nrow(x)), function(i) lapply(x,"[",i))
# tab <- rows(mutdfr)

evlist <- mutdfr$`Expected CMs`
kappalist <- mutdfr$kappa
mulist <- mutdfr$mu
philist <- mutdfr$phi

evs <- vector("numeric",1000)
kappas <- vector("numeric",1000)
mus <- vector("numeric",1000)
phis <- vector("numeric",1000)

i <- 1
j <- 1

while (i <= length(evlist)) {
  evs[[j]] = evlist[[i]]
  kappas[[j]] = kappalist[[i]]
  mus[[j]] = mulist[[i]]
  phis[[j]] = philist[[i]]
  i = i + 1000
  j = j + 1
}

i <- 1
j <- 1

x <- 1:1000
y <- mutdfr$`Average CMs`

ym <- matrix(y,nrow = 1000,byrow=F)

colsiny <- dim(ym)[[2]]
loopindex <- 1:colsiny

for (k in loopindex) {
  jerry = ym[,k]
  ev = evs[[k]]
  
  plot <- ggplot() +
    geom_point(mapping = aes(x = x, y = jerry)) +
    geom_hline(yintercept = ev, linetype = 2, color = "red", size = 1.25) +
    labs(x = "Number of Simulations", 
         y = "Average convergent mutations", 
         title = paste("MC_Sim_Graph_mu_",mus[[k]],"_kappa_",kappas[[k]],"_phi_",phis[[k]],"_L_1000.png",sep=""), 
         subtitle = paste("E[c] = ",ev,sep="")) +
    theme_classic()
  #annotate("text",x = 750, y = (ev*IQR(y)), label = paste("E[c] = ",ev,sep=""), size = 4) 
  
  print(plot)
  fn = paste("MC_Sim_Graph_mu_",mus[[k]],"_kappa_",kappas[[k]],"_phi_",phis[[k]],"_L_1000.png",sep="")
  #ggsave(fn,device="png")
  
}