# initial parameter
M<-10;T<-10;n<-100;
beta1<-beta2 <- pu <- sigma1 <- sigma2 <- numeric(M)
beta1[1]<-0.9;beta2[1]<-4;
pu[1]<-0.5;sigma1[1]<-1.5;sigma2[1]<-9.5;K<-50

# generate fixed effect X
X<-matrix(rnorm(1000,0,1),n,T)
X1<-X2<-matrix(0,n,T)
U <- sample(c(1,2),size = n,replace = T,prob = c(pu[1],1-pu[1]))
X1[U==1,] <- X[U==1,]
X2[U==2,] <- X[U==2,]
wi1 <- matrix(0,nrow=M,ncol=n)
wi1[1,U==1] <- 1




#generate Y
U_true <- sample(c(1,2),size = n,replace = T,prob = c(0.6,0.4))
X1_true <- X[U_true==1,]
X2_true <- X[U_true==2,]
Z1_true <- rnorm(length(X1_true),0,2)
Z2_true <- rnorm(length(X2_true),0,10)
Y1 <- ifelse(exp(X1_true+Z1_true)/(1+exp(X1_true+Z1_true))>0.5,1,0)
Y2 <- ifelse(exp(5*X2_true+Z2_true)/(1+exp(5*X2_true+Z2_true))>0.5,1,0)
Y<-matrix(0,nrow=100,ncol=10)
Y[U_true==1,]<-Y1
Y[U_true==2,]<-Y2



flag<-1
while(flag==1){
  
  
F1 <- function(pai,sigma,z,x,beta){
  return(pai*dnorm(z,0,sigma)*(exp(beta*x+z)/(1+exp(beta*x+z))))
}
  

MH<-function(x,beta,sigma,pai){
    z<-numeric(1000)
    z[1]<-rnorm(1,0,sigma)
    for(i in 2:1000){
      zt<-z[i-1]
      Z<-rnorm(1,zt,sigma)
      num<-F1(pai,sigma,Z,x,beta)
      den<-F1(pai,sigma,zt,x,beta)
      u<-runif(1)
      ifelse(u<(num/den),z[i]<-Z,z[i]<-zt)
    }
    z1<-sample(z[-c(1:100)],1)
    return(z1)
  }
# generate #500*100 z
Z1 <- matrix(0,nrow = K,ncol=n)
Z2 <- matrix(0,nrow = K,ncol=n)

Xmean <- apply(X, 1, mean)
for(k in 1:K){
 for(i in 1:n){

     Z1[k,i]<-MH(Xmean[i],beta1[m-1],sigma1[m-1],pu[m-1])
     Z2[k,i]<-MH(Xmean[i],beta2[m-1],sigma2[m-1],pu[m-1])
  
   
 }
}

#generate u(update pai)
num <- 0
den <- 0
for(k in 1:K){
 for(i in 1:n){
  for(j in 1:T){
      num <- num+exp(beta1[m-1]*X[i,j]+Z1[k,i])^Y[i,j]/(1+exp(beta1[m-1]*X[i,j]+Z1[k,i]))
      den <- den+exp(beta2[m-1]*X[i,j]+Z2[k,i])^Y[i,j]/(1+exp(beta2[m-1]*X[i,j]+Z2[k,i]))

  }
 }
}
Num <- num*pu[m-1] 
Den <- num*pu[m-1]+den*(1-pu[m-1])
pu[m] <- Num/Den
U <- sample(c(1,2),size = n,replace = T,prob = c(pu[m],1-pu[m]))
wi1[m,U==1] <- 1


#update sigma 1 & 2
num1 <- num2 <- 0
for(k in 1:K){
  for(i in 1:n){
      num1 <- num1+wi1[m-1,i]*Z1[k,i]^2
      num2 <- num2+(1-wi1[m-1,i])*Z2[k,i]^2
  }
}
sigma1[m] <- sqrt(1/K*num1/sum(wi1[m-1,]))
sigma2[m] <- sqrt(1/K*num2/sum(1-wi1[m-1,]))


#update beta 1&2
s1<-s2<-s3<-s4<-0
for(k in 1:K){
  for(i in 1:n){
    for(j in 1:T){
      s1 <- s1+wi1[m-1,i]*(Y[i,j]*X[i,j]-X[i,j]*exp(beta1[m-1]*X[i,j]+Z1[k,i])/(1+exp(beta1[m-1]*X[i,j]+Z1[k,i])))
      s2 <- s2-wi1[m-1,i]*(X[i,j]^2*exp(beta1[m-1]*X[i,j]+Z1[k,i])/(1+exp(beta1[m-1]*X[i,j]+Z1[k,i]))^2)
      s3 <- s3+(1-wi1[m-1,i])*(Y[i,j]*X[i,j]-X[i,j]*exp(beta2[m-1]*X[i,j]+Z2[k,i])/(1+exp(beta2[m-1]*X[i,j]+Z2[k,i])))
      s4 <- s4+(wi1[m-1,i]-1)*(X[i,j]^2*exp(beta2[m-1]*X[i,j]+Z2[k,i])/(1+exp(beta2[m-1]*X[i,j]+Z2[k,i]))^2)
    }
  }
}
beta1[m] <- beta1[m-1]-s1/s2
beta2[m] <- beta2[m-1]-s3/s4

#convergence rule
while(abs(beta1[m]-beta1[m-1])<=0.5 & abs(beta2[m]-beta2[m-1])<=0.5 & abs(sigma1[m]-sigma1[m-1])<=0.5 & abs(pu[m]-pu[m-1])<=0.5)
{flag<-0}
}