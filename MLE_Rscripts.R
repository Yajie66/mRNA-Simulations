

#original function,ssa
ori.f<-function(Mt,theta,S0){
  
  Sini<-S0
  
  
  S.can<-S0
  
  t<-0
  
  Reaction<-matrix(c(-1,0,0,+1,0,-1),byrow=TRUE,ncol=2)
  
  
  
  while(t[length(t)] <  Mt){
    
    
    if ( S.can[1] != 0  |   S.can[2] != 0 ){
      
      
      theta.can<-c(theta[1]*S.can[1], theta[2]*S.can[1], theta[3]*S.can[2])
      
      lambda<-sum(theta.can, na.rm = T) 
      
      wait<-rexp(1,rate=lambda)
      
      
      
      t[length(t)+1] <-t[length(t)]+wait
      
      
      if(t[length(t)]<= Mt){
        
        if ( S.can[1] > 0  ){
          
          
          prob<-theta.can/lambda   
          
          R.can<-sample(size =1,x=c(1,2,3), prob = prob)
          
          
          
          if(R.can==1){
            S.can<-S.can+Reaction[1,] 
            
            
          }else{
            
            if(R.can==3){
              #protein death process
              
              
              
              S.can<-S.can+Reaction[3,]
              
              
            }else{
              
              S.can<-S.can+Reaction[2,] #protein birth process
            }
          }
          
        }else{
          
          if(S.can[2] > 0){
            S.can<- S.can+Reaction[3,]  #protein death process
            
          }
          
        }
        
      }else{
        t[length(t)]=Mt
        
      }   
      
    }
    
    
    Sini<-rbind.data.frame(Sini,S.can)
    
  }
  
  total.d<- cbind.data.frame(t,Sini)
  
  colnames(total.d)<-c("t","S1.ori","S2.ori")
  
  return(total.d)
  
}



#Extract the Ei(event.type up or down) and Ti(integral of integral) from original Protein process(S2)

OriS2.ei.ti<-function(original.out){
  
  S2.ori<-original.out[,3]
  
  Ei<-rep(NA,l=length(S2.ori))
  Ti<-original.out[,1]
  
  ori.I2<-0
  
  
  for(i in 2:(length(S2.ori)-1)){
    
    if(S2.ori[i] == (S2.ori[i-1]+1)){
      Ei[i]<-1 #birth_of_protein
      
    }else{
      if(S2.ori[i] == (S2.ori[i-1]-1)){
        Ei[i]<-0  #death_of_protein
        
      }else{
        Ei[i]<-NA 
        Ti[i]<-NA
      }}
  }
  
  
  
  #the integral of protein
  
  index<-which(!is.na(Ti))
  
  
  for(i in 2:length(index)){
    
    ori.I2<-ori.I2+(Ti[index[i]]-Ti[index[i-1]])*S2.ori[index[i-1]]
  }
  
  S2.data.na<-cbind(Ti,S2.ori,Ei)
  
  S2.data<-S2.data.na[index,]
  
  list(total.data=S2.data[,1:2],Event.type=S2.data[,3],Integral.protein=ori.I2)
}



#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################

#Approximation

s_1.value<-function(theta,S0,t){
  S0[1]*exp(-theta[1]*t)
}#input whole vector

Next_birth_time<-function(ti,theta,S0){
  
  s<-rexp(1,rate=1)
  
  In.log<--(theta[1]*s)/(theta[2]*S0[1])+exp(-theta[1]*ti) 
  
  while(In.log <= 0){
    
    s<-rexp(1,rate=1)
    
    In.log<--(theta[1]*s)/(theta[2]*S0[1])+exp(-theta[1]*ti) 
  }
  
  Tb<-(-1/theta[1])*log(In.log)-ti
  return(Tb)
}
#input whole vector


#################################################################################
###################################################################################
#################################################################################




#The function of one simulation result 

# Secondly, setting up the approximation function

Approx.mrna<-function(Time.max,theta,S0){
  
  
  Sini<-S0 # initial value
  Sini.new<-S0 # change with the process of reaction
  Event.type<-NA #e_i as vector
  
  
  Integral.protein<-0
  
  
  t<-0 # time start from 0
  
  
  while(t[length(t)] < Time.max){
    
    
    if ( Sini.new[1] != 0  |   Sini.new[2] != 0 ){
      
      
      if(Sini.new[1] >0 ){
        #next_birth_time_start_from_t{i-1}
        tb<- Next_birth_time(ti=t[length(t)],theta=theta,S0=S0)
        
      }else{
        Event.type=c(Event.type,0) # Death process
        tb=Inf
      }
      
      death.rate=theta[3]*Sini.new[2]
      
      
      if(death.rate != 0){
        td<-rexp(1,rate=death.rate)
      }else{
        # birth process
        Event.type=c(Event.type,1)
        td=Inf}
      
      
      
      
      wait<-min(tb,td) # the minimal waiting time
      
      Direction<-which(wait==c(tb,td))#provide the direction of next jump
      
      t[length(t)+1] <-t[length(t)]+wait
      
      
      
      
      #control the final jump to the T.max
      if(t[length(t)] <=  Time.max){
        
        Integral.protein<- Integral.protein + Sini.new[2]* wait
        
        #calculation the I_2
        
        
        if(Direction==1){
          
          # the birth process
          
          Sini.new[2] <-Sini.new[2] +1 
          Event.type=c(Event.type,1)# birth process
          
        }else{
          #death process
          Sini.new[2] <-Sini.new[2] -1 
          Event.type=c(Event.type,0) # Death process
          
        } # the death process
        
        Sini.new[1]<-s_1.value(theta=theta,t=t[length(t)],S0=S0)
        
      }else{ 
        
        Integral.protein<- Integral.protein + Sini.new[2]*(Time.max-t[length(t)-1]) 
        
        t[length(t)]<-Time.max
        
        Sini.new <- c( s_1.value(theta=theta, t=t[length(t)], S0=S0),Sini.new[2]) 
        
        
        Event.type[length(Event.type)]<-NA
        
        #Event.type<-c(Event.type,NA)
        
        
      } 
      
    }
    
    Sini<-rbind.data.frame(Sini,Sini.new)
    
  }
  
  total.data<-cbind.data.frame(t,Sini,Event.type)
  
  colnames(total.data)<-c("times","mRNA","Protein","Event.type")
  
  
  list(total.data=total.data, Event.type=Event.type, Integral.protein=Integral.protein)
}




##Parameter Estimator

#input the (theta_1, log(theta_2) and log(theta_3))

Neg.loglike.fun<-function(theta,ei,Ti,S0.1,I.2){
  
  theta1<-theta[1]
  
  # ln.theta2<-log(theta[2])
  
  ln.theta2<-theta[2]
  
  theta2<-exp(theta[2])
  
  # ln.theta3<-log(theta[3])
  
  ln.theta3<-theta[3]
  theta3<-exp(theta[3])
  
  
  OUT<- -( sum(ei)*ln.theta2- theta1*sum(ei*Ti)+ sum(1-ei)*ln.theta3+ ((theta2*S0.1)/theta1)*(exp(-theta1*t.max)-1)-theta3*I.2)
  
  return(OUT)
}



MLE.es<-function(single.out,S0){
  
  Ti<-single.out$total.data[,1][-c(1,length(single.out$total.data[,1]))]
  
  ei<-single.out$Event.type[-c(1,length(single.out$Event.type))]
  
  S0.1<-S0[1]
  
  I.2<-single.out$Integral.protein
  
  #hat.mle<-optim(par=c(0.01,0.01,0.01),fn=Neg.loglike.fun,Ti=Ti,ei=ei,S0.1=S0.1,I.2=I.2,lower = Boundary.parameter[1], upper =Boundary.parameter[2],method="L-BFGS-B")$par
  
  
  hat.ln.mle<-optim(par=c(0.01,0.01,0.01),fn=Neg.loglike.fun,Ti=Ti,ei=ei,S0.1=S0.1,I.2=I.2)$par
  
  hat.mle<-c(hat.ln.mle[1],exp(hat.ln.mle[2:3]))                 
  return(hat.mle)
}


#MLE function of original data form 




##Iteration with MLE estimator with multiply trajectories to check the MlE behavior.

# Let us define the General MLE estimator function

###change the data to original form

# MLE estimator with multiply iterations 




iter.mle<-function(iter,t.m,theta,S0){
  
  theta_ori_mle<-matrix(rep(NA,iter*3),ncol=3)
  theta_approx_mle<-matrix(rep(NA,iter*3),ncol=3)
  
  colnames(theta_ori_mle)<-c("theta1.ori.mle","theta2.ori.mle","theta3.ori.mle")
  
  colnames(theta_approx_mle)<-c("theta1.approx.mle","theta2.approx.mle","theta3.approx.mle")
  
  for(i in 1:iter){
    
    single.ori.data<-ori.f(Mt=t.m,theta=theta,S0=S0)
    single.ei.ti.I2<-OriS2.ei.ti(original.out = single.ori.data)
    
    single.approx.out<- Approx.mrna(Time.max=t.m,theta = theta, S0=S0)
    
    theta_ori_mle[i,] <-MLE.es(single.out=single.ei.ti.I2,S0=S0)
    theta_approx_mle[i,] <- MLE.es(single.out=single.approx.out,S0=S0)
  }
  
  list(theta_ori_mle=theta_ori_mle,theta_approx_mle=theta_approx_mle)
}



N=1000
#Boundary.parameter<-c(zapsmall(0+(1e-5), digits = 6),zapsmall(1-(1e-5), digits = 6))

theta<-c(0.11,0.001*N,0.03) 

S0<-c(N,0)

t.max=100
iter=2000


print("start simulation with iteration")

start.time <- Sys.time()

mle.results<-iter.mle(iter=iter,t.m=t.max,theta=theta,S0=S0)


end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)


print(paste("run",iter,"simulation with N=",N,"takes",time.taken,"s"))

print("finished simulation with iteration")


mle.ori<-mle.results$theta_ori_mle
mle.approx<-mle.results$theta_approx_mle


#write.csv(mle.ori, "C:\\Users\\pmxyg3\\downloads\\mle-ori-n=10000,iter=1000.csv", row.names=TRUE)
#write.csv(mle.approx, "C:\\Users\\pmxyg3\\downloads\\mle-approx-n=10000,iter=1000.csv", row.names=TRUE)

write.csv(mle.ori, "rep3.mle-ori-n=1000,iter=2000.csv", row.names=TRUE)
write.csv(mle.approx, "rep3.mle-approx-n=1000,iter=2000.csv", row.names=TRUE)

print("finished.mle.data saved with file")
