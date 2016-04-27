###input signal###
rm(list=ls())
uLen<-1:10000
u<-sin(uLen/5)

# plot some of it
while( dev.cur() != 1 ) dev.off() # close all previous plots
dev.new()
plot(u[1:1000],type='l')
title('input signal')

###define the initial length and train length
initLen<-100
trainLen<-300

# generate the ESN reservoir 100
inSize = outSize = 1
resSize = 100

set.seed(100)
Win<-matrix(sample(c(-1,1),size=inSize*resSize,replace = TRUE,prob=c(0.5,0.5)),resSize)
W<-matrix(sample(c(0,0.4,-0.4),size=10000,replace=TRUE,prob=c(0.95,0.025,0.025)),resSize)
cat('Computing spectral radius...')
rhoW = abs(eigen(W,only.values=TRUE)$values[1])
print('done.')

# allocated memory for the design (collected states) matrix
X = matrix(0,resSize,trainLen-initLen)
# set the corresponding target matrix directly
Yt = 1/2*(sin(uLen/5))^7

# run the reservoir with the data and collect X
x = rep(0,resSize)
for (t in 1:trainLen){
  un = u[t]
  x = tanh( Win %*% un + W %*% x )
  if (t > initLen)
    X[,t-initLen] = x
}

# plot some signals
dev.new()
matplot( t(X[(1:6),(1:100)]), type='l' )
title(main=expression(paste('Some reservoir activations ', bold(x)(italic(n)))))

# regression
rownames(X)<-paste("x",c(1:100),sep="")
Y<-Yt[(initLen+1):(trainLen)]

X<-t(X)
dataRegress<-cbind(Y,X)
dataRegress<-as.data.frame(dataRegress)
fit<-lm(Y~.,data = dataRegress)


