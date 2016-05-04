#' The markey Glass model function
#' 
#' Generate the data sequence of markey glass model
#' @param initialValue The initial value of the sequence value
#' @param lengthx The sequence length you want
#' @param stepSize The stepSize of the integration
#' @param alpha The parameter of the markey glass differential equation
#' @param beta The parameter of the markey glass differential equation
#' @param gama The parameter of the markey glass differential equation
#' @param tal The time delay of the markey glass differential equation
#' 
#' @return xNeed The markey glass sequence
#'@examples 
#' data<-markeyGlass()
#' plot(data[1000:2000],type='l',xlab='Time',ylab='x')
#' title('The sequence of markey glass')
#' dev.new()
#' plot(data[1000:2000],data[1017:2017],type='l',xlab = expression('x'[Tau]),ylab = expression('x'[Tau+17]))
#'@export

markeyGlass<-function(initialValue=0.5,lengthx=3000,stepSize=0.1,alpha=0.2,beta=10,gama=0.1,tal=17,...){
  xlen<-(lengthx-1)/stepSize+1
  x<-c(1:xlen)
  x[1:(tal/stepSize+1)]<-initialValue
  for(i in (tal/stepSize+1):(xlen-1)){
    x[i+1]<-x[i]+stepSize*(alpha*x[i-(tal/stepSize)]/(1+x[i-(tal/stepSize)]^beta)-gama*x[i])
  }
  xNeed<-x[seq(1,xlen,1/stepSize)]
  return(xNeed)
}