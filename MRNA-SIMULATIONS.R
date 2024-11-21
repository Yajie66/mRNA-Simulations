
# The original stochastic model with three reactions

ssa.f<-function(Mt,theta,S0){
  
  Sini<-S0
  
  names(Sini)<-c('S1','S2')
  
  S.can<-S0
  
  t<-0
  
  Reaction<-matrix(c(-1,0,0,+1,0,-1),byrow=TRUE,ncol=2)
  
  while(t[length(t)] <= Mt){
    
    if ( S.can[1] != 0  |   S.can[2] != 0 ){
      
      theta.can<-c(theta[1]*S.can[1], theta[2]*S.can[1], theta[3]*S.can[2])
      
      lambda<-sum(theta.can, na.rm = T) 
      
      wait<-rexp(1,rate=lambda)
      
      
      t[length(t)+1] <-t[length(t)]+wait
      
      if ( S.can[1] > 0  ){
        
        
        prob<-theta.can/lambda   
        
        R.can<-sample(size =1,x=c(1,2,3), prob = prob)
        
        
        
        if(R.can==1){
          S.can<-S.can+Reaction[1,] } else{
            
            if(R.can==3){
              S.can<-S.can+Reaction[3,]
            }else{
              
              S.can<-S.can+Reaction[2,]
            }
          }
        
      }else{
        
        if(S.can[2] > 0){
          S.can<- S.can+Reaction[3,] }else{ return( cbind.data.frame(t,Sini) )  } 
        
      }
      
    }else{ return(cbind.data.frame(t,Sini) )}
    
    Sini<-rbind.data.frame(Sini,S.can)
    
  }
  cbind.data.frame(t,Sini)
}




#####################################################################
#####################################################################
#####################################################################
#the approximated model

library(deSolve)
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

# Secondly, setting up the approximation function

Approx.mrna<-function(Time.max,theta,S0){
  
  
  Sini<-S0 # initial value
  Sini.new<-S0 # change with the process of reaction
  #Event.type<-NA #e_i as vector
  
  
  #Integral.protein<-0
  
  
  t<-0 # time start from 0
  
  
  while(t[length(t)] < Time.max){
    
    
    if ( Sini.new[1] != 0  |   Sini.new[2] != 0 ){
      
      
      if(Sini.new[1] >0 ){
        #next_birth_time_start_from_t{i-1}
        tb<- Next_birth_time(ti=t[length(t)],theta=theta,S0=S0)
        
      }else{
       # Event.type=c(Event.type,0) # Death process
        tb=Inf
      }
      
      death.rate=theta[3]*Sini.new[2]
      
      
      if(death.rate != 0){
        td<-rexp(1,rate=death.rate)
      }else{
        # birth process
       # Event.type=c(Event.type,1)
        td=Inf}
      
      
      
      
      wait<-min(tb,td) # the minimal waiting time
      
      Direction<-which(wait==c(tb,td))#provide the direction of next jump
      
      t[length(t)+1] <-t[length(t)]+wait
      
      
      
      
      #control the final jump to the T.max
      if(t[length(t)] <=  Time.max){
        
       # Integral.protein<- Integral.protein + Sini.new[2]* wait
        
        #calculation the I_2
        
        
        if(Direction==1){
          
          # the birth process
          
          Sini.new[2] <-Sini.new[2] +1 
         # Event.type=c(Event.type,1)# birth process
          
        }else{
          #death process
          Sini.new[2] <-Sini.new[2] -1 
         # Event.type=c(Event.type,0) # Death process
          
        } # the death process
        
        Sini.new[1]<-s_1.value(theta=theta,t=t[length(t)],S0=S0)
        
      }else{ 
        
       # Integral.protein<- Integral.protein + Sini.new[2]*(Time.max-t[length(t)-1]) 
        
        t[length(t)]<-Time.max
        
        Sini.new <- c( s_1.value(theta=theta, t=t[length(t)], S0=S0),Sini.new[2]) 
        
        
      #  Event.type[length(Event.type)]<-NA
        
        #Event.type<-c(Event.type,NA)
        
        
      } 
      
    }
    
    Sini<-rbind.data.frame(Sini,Sini.new)
    
  }
  
 # total.data<-cbind.data.frame(t,Sini,Event.type)
  total.data<-cbind.data.frame(t,Sini)
 # colnames(total.data)<-c("times","mRNA","Protein","Event.type")
  colnames(total.data)<-c("times","mRNA","Protein")
  
 # list(total.data=total.data, Event.type=Event.type, Integral.protein=Integral.protein)
  list(total.data=total.data)
}





#################################################################

#################################################################
# Iteration with the Approximated model and Original model


#library(dplyr)
# The iteration of approximated model

iter.appro.out<-function(Mt,Time.max,scale,theta,S0,iter){
  
  times<-seq(0,Mt,by=Mt/scale) 
  
  # Mt is the maximum time length of the series,scale is the time division 
  length.t<- scale+1
  
  si.out<-matrix(NA,ncol=iter+2, nrow = length.t )
  
  si.out[,1]<-times #first col is the time
  
  s2.name<-as.character(paste("S2",c(1:iter),sep="."))
  
  colnames(si.out)<-c("times","S1",s2.name)
  
  
  si.out[1,2]<-S0[1]
  
  
  #As for S1, since each iteration the S1 value is the determined model, so simply give the s1 ode solution with specific time  
  for(j in 2:length.t ){
    t.j<-times[j]
    si.out[j,2]<-s_1.value(theta=theta,S0=S0,t=t.j)
  } 
  
  # As for S2  
  for(i in 1:iter){
    
    #t.S2<-na.omit( Approx.mrna(Time.max = Time.max,theta=theta,S0=S0)$total.data %>% select("times","Protein") ) #the reduce function results
    t.S2<-na.omit( Approx.mrna(Time.max = Time.max,theta=theta,S0=S0)$total.data) #the reduce function results
    
    index<-findInterval(times,t.S2[,1]) # select the interval index of S2 at specific time
    
    si.out[,i+2]=t.S2[index,2]
    
  }
  si.out
}

#We can control the iteration number, get the data from sigle cell time lapse with equal 





##########################################################################

##########################################################################

##########################################################################
#

##########################################################################
#
##########################################################################
# the iteration of original stochastic model
#simulation based, one time scale to plot the graphs


ssa.out<-function(Mt,scale,theta,S0,iter){
  
  times=seq(0,Mt,by=Mt/scale)
  
  si.out<-matrix(NA,ncol=2*iter+1,nrow = length(times))
  si.out[,1]<-times
  #get the name
  s1.name<- paste("S1.ssa",c(1:iter),sep = ".")
  #s1.name<-rep('S1.ssa',iter) 
  s2.name<-paste("S2.ssa",c(1:iter),sep = ".")
  #s2.name<-rep('S2.ssa',iter)
  colnames(si.out) <- c("times",s1.name,s2.name)
  
  for(i in 1:iter){
    
    result<-na.omit(ssa.f(Mt,theta,S0)) 
    index<-findInterval(times,result[,1])
    si.out[,i+1]<-result[index,2]
    si.out[,i+1+iter]<-result[index,3]
  }
  
  si.out
}



# Graphic simulation with stochastic model and its approximation
# The final results

#library(reshape)
#library(ggplot2)
#library(gghighlight)

#one iteration

N=1000
S0<-c(N,0)
theta<-c(0.11, 0.001*N, 0.03)
iter=1
Mt=10
scale=1000

print("start")
start<-Sys.time()
Ori<-ssa.out(Mt=Mt,scale=scale,theta=theta,S0=S0,iter=iter)
end<-Sys.time()
print("Original.finished")

(total.time<-end-start)








Approx<-iter.appro.out(Mt=Mt,Time.max=Mt,scale=scale,theta=theta,S0=S0,iter=iter) 

print("approx.finished")








write.csv(Ori, "ori.sim.data-n=1000.csv", row.names=TRUE)

write.csv(Approx, "approx.sim.data-n=1000.csv", row.names=TRUE)
print("finished.original.data")








