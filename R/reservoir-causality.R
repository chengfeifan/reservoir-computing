
library(myCCM)
data("cstr")
Ca<-cstr$Ca[-1000:-1]
Cb<-cstr$Cb[-1000:-1]

# plot some of it
while( dev.cur() != 1 ) dev.off() # close all previous plots
dev.new()
plot(Ca[1000:2000],type='l')
title('input signal Ca')

dev.new()
plot(Cb[1000:2000],type='l')
title('input signal Cb')

# define the initial length and train length
initLen<-100
trainLen<-300
testLen<-500

# generate the ESN reservoir 100
inSize = outSize = 1
resSize = 100

