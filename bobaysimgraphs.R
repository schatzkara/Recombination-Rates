setwd("S://simulationdata//kappaGCsims//GC = 0.5")

library(ggplot2)
library(gdata)
library(grid)

#Import K1 files
setwd("S://simulationdata//kappaGCsims//GC = 0.5//kappa = 1")

#Load all .csv in directory into list
k1Files <- list.files(pattern="*.csv")
Nk1 <- length(k1Files)
k1List <- vector("list",Nk1)
k1Params = list()

j <- 1

for(i in k1Files) {
  #read all of the csv files
  name = gsub("-",".",i)
  name = gsub(".csv","",name)
  params = unlist(strsplit(name,split="_"))[c(5,7)]
  k1Params[[j]]=params
  i = paste(".\\",i,sep="")
  assign(name,read.csv(i, header=T))
  
  #add to dataList
  k1List[[j]] = assign(name,read.csv(i, header=T))
  j = j+1 
}

#Import K2 files
setwd("S://simulationdata//kappaGCsims//GC = 0.5//kappa = 2")

#Load all .csv in directory into list
k2Files <- list.files(pattern="*.csv")
Nk2 <- length(k2Files)
k2List <- vector("list",Nk2)
j <- 1
k2Params = list()

for(i in k2Files) {
  #read all of the csv files
  name = gsub("-",".",i)
  name = gsub(".csv","",name)
  params = unlist(strsplit(name,split="_"))[c(4,5)]
  k2Params[[j]]=params
  i = paste(".\\",i,sep="")
  assign(name,read.csv(i, header=T))
  
  #add to dataList
  k2List[[j]] = assign(name,read.csv(i, header=T))
  j = j+1 
}

#Import K3 files
setwd("S://simulationdata//kappaGCsims//GC = 0.5//kappa = 3")

#Load all .csv in directory into list
k3Files <- list.files(pattern="*.csv")
Nk3 <- length(k3Files)
k3List <- vector("list",Nk3)
j <- 1
k3Params = list()

for(i in k3Files) {
  #read all of the csv files
  name = gsub("-",".",i)
  name = gsub(".csv","",name)
  params = unlist(strsplit(name,split="_"))[c(5,6)]
  k3Params[[j]]=params
  i = paste(".\\",i,sep="")
  assign(name,read.csv(i, header=T))
  
  #add to dataList
  k3List[[j]] = assign(name,read.csv(i, header=T))
  j = j+1 
}

#clear unnecessary variables
keep(k1List,k1Params,k2List,k2Params,k3List,k3Params,Nk1,Nk2,Nk3,sure=T)
ls()
i <- 1
fn = list()

for (n in 1:1000) {
  fn[[i]] <- paste("plot",i,sep = "_")
  i = i+1
}

i <- 1
avk1 = vector("numeric",length(k1List[[1]][[2]]))

snps = 1:300

for (k in snps) {
  avk1[[i]] = sum(k1List[[k]][[2]])/300
  i = i+1
}



k1Data = data.frame(snps,avk1)

sampleplot = ggplot(data=k1Data,aes(x=k1Data[[1]],y=k1Data[[2]]))
print(sampleplot + geom_smooth())

# #ggplot the data
# for (m in k1List) {
#   fn[[i]] = ggplot(data = m, aes(x=m[[1]]))
#   fn[[i]] = fn[[i]] + geom_smooth(aes(x=m[[1]],y=m[[2]],color = "K = 1"),method="loess")
#   i = i+1
# }
# 
# i <- 1
# 
# for (n in k2List) {
#   fn[[i]] = fn[[i]] + geom_smooth(aes(x=n[[1]],y=n[[2]],color = "K = 2"),method="loess")
#   i = i+1
# }
# 
# i <- 1
# 
# for (p in k3List) {
#   filename = paste("Simulation #",k1Params[[i]][[2]],"with %GC =",k1Params[[i]][[1]],sep=" ")
#   fn[[i]] = fn[[i]] + geom_smooth(aes(x=p[[1]],y=p[[2]],color = "K = 3"),method="loess")
#   fn[[i]] = fn[[i]] + labs(title=filename,x="Mutations",y="Conv. Mutations")
#   print(fn[[i]])
#   i = i+1
# }
# 
