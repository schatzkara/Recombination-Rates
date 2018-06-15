setwd("S://simulation data")

library(ggplot2)
library(gdata)
library(grid)

#Load all .csv in directory into list
dataFiles <- list.files(pattern="*.csv")
N <- length(dataFiles)
dataList <- vector("list",N)
j <- 1
paramList = list()

for(i in dataFiles) {
  #read all of the csv files
  name = gsub("-",".",i)
  name = gsub(".csv","",name)
  params = unlist(strsplit(name,split="_"))[c(4,5)]
  paramList[[j]]=params
  i = paste(".\\",i,sep="")
  assign(name,read.csv(i, header=T))
  
  #add to dataList
  dataList[[j]] = assign(name,read.csv(i, header=T))
  j = j+1 
}

#reset j, initialize vector of lists
j <- 1
evList <- vector("list",N)
mkList <- vector("list",N)

#get out expected values
for (k in dataList) {
  #extract expected values
  evk = k[1001,1]
  name = paste("expected", j, sep = "")
  evList[[j]] = assign(name, evk)
  
  #extract mu^2/3
  muk = k[1002,1]
  vname = paste("mu squared over three", j, sep = "")
  mkList[[j]] = assign(vname, muk)
  
  k = k[-c(1002), ]
  k = k[-c(1001), ]
  dataList[[j]] = k
  j = j+1
}

#clear unnecessary variables
keep(dataList,evList,mkList,paramList,sure=T)
ls()
i = 1

#ggplot those graphs
for (k in dataList) {
  cp = ggplot(k, aes(x=k[[1]],y=k[[2]]))
  cp = cp+geom_point()
  cp = cp+ggtitle(paste("Monte Carlo Simulation with mu = ",paramList[[i]][[2]], "and L = ",paramList[[i]][[1]]))
  cp = cp+labs(x="Number of simulations",y="Expected number of convergent mutations")
  cp = cp+geom_hline(yintercept=evList[[i]], linetype="dashed", color = "red", size=1.25)
  cp = cp+theme_classic()
  str = paste("E[m] = ",round(evList[[i]],2))
  iqrk = IQR(k[ ,2])
  madk = mad(k[ ,2],center = median(k[ ,2]))
  sf = 3*(madk+iqrk)
  if (k[750,2]>=evList[[i]]) {
    cp = cp + annotate("text", x = 750, y = (k[750,2]+sf), label = str, size=4)
  } else {
    cp = cp + annotate("text", x = 750, y = (evList[[i]]+sf), label = str, size=4)
  }
  print(cp)
  fn = paste("sim_mu",paramList[[i]][[2]], "_L",paramList[[i]][[1]],".png")
  ggsave(fn,device = "png")
  i = i+1
}