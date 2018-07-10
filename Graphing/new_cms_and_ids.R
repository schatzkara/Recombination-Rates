setwd("C://Users//insan_000//Desktop//")

library(gdata)
library(grid)
library(tidyverse)

#set up dataframe
mutdfr <- read_csv("all_averages_for_id_sim_1000.csv")
mutdfr <- subset(mutdfr,mutdfr$`GC%` == 0.5)
moddfr <- read_csv("all_model_c_vs_id_data.csv")
# mutdfr$mu <- as.factor(mutdfr$mu)
# mutdfr$kappa <- as.factor(mutdfr$kappa)
# mutdfr$phi <- as.factor(mutdfr$phi)

# rows = function(x) lapply(seq_len(nrow(x)), function(i) lapply(x,"[",i))
# tab <- rows(mutdfr)

n <- 299

kappalist <- mutdfr$kappa
mulist <- mutdfr$mu
philist <- mutdfr$phi

evs <- vector("numeric",n)
kappas <- vector("numeric",n)
mus <- vector("numeric",n)
phis <- vector("numeric",n)

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

x <- 1:n
y <- mutdfr$`E[c]`

ym <- matrix(y,nrow = n,byrow=F)

colsiny <- dim(ym)[[2]]
loopindex <- 1:colsiny

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
