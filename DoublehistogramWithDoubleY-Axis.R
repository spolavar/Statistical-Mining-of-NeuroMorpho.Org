#first distribution
d1 <- rnorm(1000,mean=0,sd=3)
#second distribution
d2 <- rnorm(40,mean=1,sd=1)

#scale factor 
scaleupby <- 10
str(scaleupby)
#compute the breaks and counts for the second distribution
newp <- hist(d2,plot=FALSE)
length(newp$breaks)
#scale the counts
newp$counts[newp$counts > 0] <- aggregate(d2,by=list(cut(d2,newp$breaks)),FUN=sum)$x*scaleupby
#the new ticks values
axpos <- seq(from=floor(min(newp$counts)),to=ceiling(max(newp$counts)),
             by=diff(range(newp$counts))/(length(newp$breaks)-1))
axpos <- seq(from=floor(min(newp$counts)))
hist(d1)
plot(newp,add=T,col='red',axes=FALSE)
axis(4)
axis(side=4,at=axpos,col='red',labels=round(axpos,0),cex=0.7)
hist(d2,add=T,col='blue')