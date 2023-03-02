rm(list=ls())
sup <- seq(0,100,0.01)
shapes <- c(0.1,1,2,5,10) 
rates <- c(0.1,1,2,5,10)
try(dev.off())
x11(width=10,height=8)
layout(matrix(1:25,ncol=5))
par(mar=rep(2,4))
for(sh in shapes){
for(rt in rates){
pdf <- dgamma(sup,rate=rt,shape=sh)
plot(sup,pdf,xlim=c(0,10),ylim=c(0,1.5),lwd=3,type='l')
text(6,1,paste('shape=',sh,', rate=',rt,sep=''))
}
}
