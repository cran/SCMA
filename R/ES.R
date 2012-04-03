ES<-function(design,ES,data=read.table(file.choose(new=FALSE))){

  if(design=="CRD"|design=="ATD"|design=="RBD"|design=="AB"|design=="ABA"){
    if(design=="CRD"|design=="ATD"|design=="RBD"|design=="AB"){
      A<-data[,2][data[,1]=="A"]
      B<-data[,2][data[,1]=="B"]
    }    
    if(design=="ABA"){
      A<-c(data[,2][data[,1]=="A1"],data[,2][data[,1]=="A2"])
      B<-data[,2][data[,1]=="B1"]
    }
    if(ES=="SMD"){
      es<-(mean(B)-mean(A))/sd(A) 
    }
    if(ES=="SMDpool"){
      es<-(mean(B)-mean(A))/(sqrt((((length(A)-1)*var(A))+((length(B)-1)*var(B)))/(length(A)+length(B)-2))) 
    }
    if(ES=="PND+"){
      es<-sum(B>max(A))/length(B)*100
    }
    if(ES=="PND-"){
      es<-sum(B<min(A))/length(B)*100
    }
    if(ES=="PEM+"){
      es<-sum(B>median(A))/length(B)*100
    }
    if(ES=="PEM-"){
      es<-sum(B<median(A))/length(B)*100
    }
  }
  
  if(design=="ABAB"){
    A1<-data[,2][data[,1]=="A1"]
    B1<-data[,2][data[,1]=="B1"]
    A2<-data[,2][data[,1]=="A2"]
    B2<-data[,2][data[,1]=="B2"]
  
    if(ES=="SMD"){
      es1<-(mean(B1)-mean(A1))/sd(A1)
      es2<-(mean(B2)-mean(A2))/sd(A2) 
    }
    if(ES=="SMDpool"){
      es1<-(mean(B1)-mean(A1))/(sqrt((((length(A1)-1)*var(A1))+((length(B1)-1)*var(B1)))/(length(A1)+length(B1)-2))) 
      es2<-(mean(B2)-mean(A2))/(sqrt((((length(A2)-1)*var(A2))+((length(B2)-1)*var(B2)))/(length(A2)+length(B2)-2))) 
    }
    if(ES=="PND+"){
      es1<-sum(B1>max(A1))/length(B1)*100
      es2<-sum(B2>max(A2))/length(B2)*100
    }
    if(ES=="PND-"){
      es1<-sum(B1<min(A1))/length(B1)*100
      es2<-sum(B2<min(A2))/length(B2)*100
    }
    if(ES=="PEM+"){
      es1<-sum(B1>median(A1))/length(B1)*100
      es2<-sum(B2>median(A2))/length(B2)*100
    }
    if(ES=="PEM-"){
      es1<-sum(B1<median(A1))/length(B1)*100
      es2<-sum(B2<median(A2))/length(B2)*100
    }
    es<-mean(c(es1,es2))
  }
  
  if(design=="MBD"){
    N<-ncol(data)/2
    effectsizes<-numeric(N)
    for(it in 1:N){
        A<-data[,it*2][data[,(it*2)-1]=="A"]	
        B<-data[,it*2][data[,(it*2)-1]=="B"]
      
      if(ES=="SMD"){
        effectsizes[it]<-(mean(B)-mean(A))/sd(A) 
      }
      if(ES=="SMDpool"){
        effectsizes[it]<-(mean(B)-mean(A))/(sqrt((((length(A)-1)*var(A))+((length(B)-1)*var(B)))/(length(A)+length(B)-2))) 
      }
      if(ES=="PND+"){
        effectsizes[it]<-sum(B>max(A))/length(B)*100
      }
      if(ES=="PND-"){
        effectsizes[it]<-sum(B<min(A))/length(B)*100
      }
      if(ES=="PEM+"){
        effectsizes[it]<-sum(B>median(A))/length(B)*100
      }
      if(ES=="PEM-"){
        effectsizes[it]<-sum(B<median(A))/length(B)*100
      }
    } 
    es<-mean(effectsizes)
  }
  
  return(es)

}