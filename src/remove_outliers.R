args<-commandArgs(TRUE)
#Collect input arguments
fin<-args[1]
iterations<-args[2]
#Print to stdout
print(paste("Removing outlier points from ",fin ,"...",sep=""))
#Create output file name
fout<-paste("fix.",fin,sep="")
# Reading the input file
data<-read.table(fin,sep=" ",skip = 1,header=F)
# Function to compute midpoint
midpoint_form<-function(rowi,rowj){
  x1<-rowi[1]
  y1<-rowi[2]
  z1<-rowi[3]
  x2<-rowj[1]
  y2<-rowj[2]
  z2<-rowj[3]  
  midpoint<-matrix(c((x1+x2)/2,(y1+y2)/2,(z1+z2)/2))
  return(midpoint)
}

# Function to compute the distance between consecutive points
distance_form<-function(rowi,rowj){
  x1<-rowi[1]
  y1<-rowi[2]
  z1<-rowi[3]
  x2<-rowj[1]
  y2<-rowj[2]
  z2<-rowj[3]  
  distance<-sqrt(((x2-x1)^2+(y2-y1)^2+(z2-z1)^2))
  return(distance)
}

# number of rows in the data
row_len<-length(data[,1])


data_old<-data
for( z in c(1:iterations)){
# Empty object to collect distances
distances<-c()

# looping to find distances
for(i in c(1:(row_len-1))){
  j=i+1
  d<-distance_form(data[i,],data[j,])
  distances<-c(distances,d)
}


out<-boxplot(distances[-length(distances)])
idx<-(1+which(distances%in%out$out))

midpoints<-matrix(nrow=0,ncol=3)
for( i in idx){
  up_i<-(i+1)
  down_i<-(i-1)
  m<-midpoint_form(data[up_i,],data[down_i,])
  midpoints<-rbind(midpoints,(t(m)))  
}
data<-data
data[idx,c(1,2,3)]<-data.frame(midpoints)
data<-format(data[,c(1,2,3)],digits=8)
data[,1]<-sub(" ","",data[,1])
data[,2]<-sub(" ","",data[,2])
data[,3]<-sub(" ","",data[,3])
data[,4]<-""
data[,1]<-as.numeric(as.character(data[,1]))
data[,2]<-as.numeric(as.character(data[,2]))
data[,3]<-as.numeric(as.character(data[,3]))
}
data_new<-data


lines<-data.frame(length(data_new[,1]))
write.table(lines,paste(fout),sep=" ",quote=F,col.names = F,row.names=F)
write.table(data_new,paste(fout),sep=" ",quote=F,col.names = F,row.names=F,append=TRUE)
print(paste("   ...Output written to:",fout))