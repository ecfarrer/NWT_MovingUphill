
library(foreach)
library(doParallel)


getDoParWorkers()



#number of iterations in the loop
iters<-50
#vector for appending output
ls<-vector('list',length=iters)
#start time
strt<-Sys.time()
#loop
for(i in 1:iters){
    #counter
    cat(i,'\n')
    to.ls<-rnorm(1e6)
    to.ls<-summary(to.ls)
    #export
    ls[[i]]<-to.ls    
    }
#end time
print(Sys.time()-strt)
# 50 iterations takes 9 sec

#number of iterations
iters<-50
#setup parallel backend to use 4 processors
cl<-makeCluster(4)
registerDoParallel(cl)
#start time
strt<-Sys.time()
#loop
ls<-foreach(i=1:iters,.combine=rbind) %dopar% {
    to.ls<-rnorm(1e6)
    to.ls<-c(summary(to.ls))
    }
print(Sys.time()-strt)
stopCluster(cl)
ls
#50 iterations take 4 sec in parallel
#.packages=c("iterators")




#playing with the i=1:5 format
#setup parallel backend to use 4 processors
cl<-makeCluster(4)
registerDoParallel(cl)
#start time
strt<-Sys.time()
#loop
ls<-foreach(i=1:5,.combine=rbind) %dopar% {
    sqrt(i)
    }
print(Sys.time()-strt)
stopCluster(cl)
ls
#.packages=c("iterators")
