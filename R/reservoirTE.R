rm(list=ls())
library(myCCM)
library(TransferEntropy)
data("cstr")
Ca<-cstr$Ca[-1000:-1]
Cb<-cstr$Cb[-1000:-1]
data<-cbind(as.vector(Ca),as.vector(Cb))
tevalue<-c(1:4000)
for(i in 1:4000){
  x<-Ca[i:(i+999)]
  y<-Cb[i:(i+999)]
  tevalue[i]<-computeTE(x,y,embedding = 10, k=3)$TE
}

# plot some of it
while( dev.cur() != 1 ) dev.off() # close all previous plots
dev.new()
plot(Ca[1000:2000],type='l')
title('input signal Ca')

# plot the value of cross correlation function
dev.new()
plot(tevalue,type='l',main = 'The TE value',xlab = 'Tag',
     ylab = 'TEValue')

dev.new()
plot(Cb[1000:2000],type='l')
title('input signal Cb')

# define the initial length and train length
initLen<-1000
trainLen<-2000
testLen<-2000

# generate the ESN reservoir 100
inSize = 2
outSize = 1
resSize = 1000
a = 0.3 # leaking rate

set.seed(42)
Win = matrix(runif(resSize*(1+inSize),-0.5,0.5),resSize)
W = matrix(runif(resSize*resSize,-0.5,0.5),resSize)
# Option 1 - direct scaling (quick&dirty, reservoir-specific):
#W = W * 0.135 
# Option 2 - normalizing and setting spectral radius (correct, slow):
cat('Computing spectral radius...')
rhoW = abs(eigen(W,only.values=TRUE)$values[1])
print('done.')
W = W * 1.25 / rhoW

# allocated memory for the design (collected states) matrix
X = matrix(0,1+inSize+resSize,trainLen-initLen)
# set the corresponding target matrix directly
Yt = matrix(tevalue[1:1000],1)

# run the reservoir with the data and collect X
x = rep(0,resSize)
for (t in 1:trainLen){
  u = as.matrix(data[t,])
  x = (1-a)*x + a*tanh( Win %*% rbind(1,u) + W %*% x )
  if (t > initLen)
    X[,t-initLen] = rbind(1,u,x)
}

# train the output
reg = 1e-8  # regularization coefficient
X_T = t(X)
Wout = Yt %*% X_T %*% solve( X %*% X_T + reg*diag(1+inSize+resSize) )

# run the trained ESN in a generative mode. no need to initialize here, 
# because x is initialized with training data and we continue from there.
Y = matrix(0,outSize,testLen)
u = as.matrix(data[trainLen+1,])
for (t in 1:testLen){
  x = (1-a)*x + a*tanh( Win %*% rbind(1,u) + W %*% x )
  y = Wout %*% rbind(1,u,x)
  Y[,t] = y
  #   # generative mode:
  #   u = y
  # this would be a predictive mode:
  u = as.matrix(data[trainLen+t+1,])
}

# compute MSE for the first errorLen time steps
errorLen = 500
mse = ( sum( (tevalue[1001:(1000+errorLen)] - Y[1,1:errorLen])^2 )
        / errorLen )
print( paste( 'MSE = ', mse ) )

# plot some signals
dev.new() 
plot( c(Y), type='l', col='green',lwd=3)
lines( tevalue[1001:3000], col='blue',lwd=3)
title(main=expression(paste('Target and generated signals ', bold(y)(italic(n)), 
                            ' starting at ', italic(n)==0 )))
legend('bottomleft',legend=c('Free-running predicted signal','Target signal'), 
       col=c('green','blue'),lty=1, bty='n' )
dev.new()
matplot( t(X[(1:20),(1:200)]), type='l' )
title(main=expression(paste('Some reservoir activations ', bold(x)(italic(n)))))

dev.new()
barplot( Wout )
title(main=expression(paste('Output weights ', bold(W)[out])))
