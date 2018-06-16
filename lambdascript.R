#Creates a matrix of lambda values for contour plotting

setwd("C://Users//insan_000//Documents//summer 2018//simulation data//export")

#set up parameters
L = seq(100, 10000, by=100)
mu = seq(0.001, 0.1, by=0.001)
alphas = c(0.05,0.01,0.005,0.001)
lambdas = t(outer(L,mu,"*"))
write.csv(lambdas,'lambdamatrix.csv')
qvals = lambdas

#run a for loop to replace qval with qpois value
for (a in alphas) {
  for (l in lambdas) {
    for (q in qvals) {
      q = qpois(a,l)
    }
  }
  #export current qval as csv before it is overwritten
  fn = paste("qpois",a,"matrix.csv")
  write.csv(qvals,fn)
}