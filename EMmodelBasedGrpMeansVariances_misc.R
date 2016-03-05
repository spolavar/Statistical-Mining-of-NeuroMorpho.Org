m <- max(abs(parms$mean[,1]))
c1 <- parms$mean[,1]/m
m <- max(abs(parms$mean[,2]))
c2 <- parms$mean[,2]/m
m <- max(abs(parms$mean[,3]))
c3 <- parms$mean[,3]/m
m <- max(abs(parms$mean[,4]))
c4 <- parms$mean[,4]/m
m <- max(abs(parms$mean[,5]))
c5 <- parms$mean[,5]/m
m <- max(abs(parms$mean[,6]))
c6 <- parms$mean[,6]/m
m <- max(abs(parms$mean[,7]))
c7 <- parms$mean[,7]/m


plot(c1,type="b",ylim=c(-1,1),xlim=c(1,32),xlab="principal components",ylab="normalized mean values")
lines(c2,type="b",lty=2,pch=0)#,col="darkred")
lines(c3,type="b",lty=3,pch=6)#,col="darkorange")
lines(c4,type="b",lty=4,pch=12,col="red")
lines(c5,type="b",lty=5,pch=7)#,col="darkgreen")
lines(c6,type="b",lty=6,pch=13)#,col="lightgreen")
lines(c7,type="b",lty=7,pch=2)#,col="red")
legend("bottomright",text.width = strwidth("1,00"), legend = c(1:params$variance$G), pch=c(1,0,6,12,7,13,2), cex=0.7, col=c(rep(1,3),2,rep(1,3)))
parms$variance$sigma[]
parms$variance$sigma[,,1]
parms$variance$sigma[1:2,1:2,1]
str(params$variance$scale)


metricmeans <- matrix(data=NA, nrow= 7,ncol=1)
rownames(metricmeans) <- c(1:7)
colnames(metricmeans) <- "mean"#colnames(sigmetricdf[,33:64])
for(c in 1:7){
  x <- subset(sigmetricdf[,c(33:64,72)],sigmetricdf$classification == c)
  dim(x)
  v <- aggregate(x[,1:32],by=list(x$classification),FUN=mean)
  dim(v)
  colnames(v)
  metricmeans[c] <- list(v[,2:33])
}

m <- max(unlist(metricmeans[1]))
c1 <- unlist(metricmeans[1])/m
m <- max(unlist(metricmeans[2]))
c2 <-unlist(metricmeans[2])/m
m <- max(unlist(metricmeans[3]))
c3 <- unlist(metricmeans[3])/m
m <- max(unlist(metricmeans[4]))
c4 <- unlist(metricmeans[4])/m
m <- max(unlist(metricmeans[5]))
c5 <- unlist(metricmeans[5])/m
m <-  max(unlist(metricmeans[6]))
c6 <-  unlist(metricmeans[6])/m
m <- max(unlist(metricmeans[7]))
c7 <- unlist(metricmeans[7])/m


parms$variance$sigma[]
parms$variance$sigma[,,1]
parms$variance$sigma[1:2,1:2,1]
parms$variance$sigma[1:2,1:2,2]
parms$variance$sigma[1:2,1:2,3]
parms$variance$sigma[1:2,1:2,4]
parms$variance$sigma[1:2,1:2,5]
parms$variance$sigma[1:2,1:2,6]
parms$variance$sigma[1:2,1:2,7]
str(params$variance$scale)
parms$variance$sigma[c(1,4),c(1,4),1]
parms$variance$sigma[c(1,4),c(1,4),2]
parms$variance$sigma[c(1,4),c(1,4),3]
parms$variance$sigma[c(1,4),c(1,4),4]
parms$variance$sigma[c(1,4),c(1,4),5]
parms$variance$sigma[c(1,4),c(1,4),6]
parms$variance$sigma[c(1,4),c(1,4),7]