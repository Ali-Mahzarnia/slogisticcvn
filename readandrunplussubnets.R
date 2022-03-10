
#run this right after running Step A matlab that would produce the connectivity data by reading and re-arranging
#install.packages("R.matlab")
library(R.matlab)
path="/Users/ali/Desktop/mar/mice/cvn/resultsconnectivity_all_ADDecode_Dipy.mat"
path2="/Users/ali/Desktop/mar/mice/cvn/resultsresponse_array.mat"

data=readMat(path)
connectivity=data$connectivity
temp=connectivity[,,1]
indexlower=lower.tri(temp, diag=FALSE)
indexlowertrue=which(indexlower==TRUE)
temp=temp[indexlower]
len=sum(indexlower)



data2=readMat(path2)
response=data2$response.array
riskfactors=matrix(NA,  dim(response)[1], (dim(response)[2]-1))
#sum(riskfactors[,2]==3)

subjnameofconnectivity=data$subjlist

for (i in 1:dim(riskfactors)[1]) {
  ind=which(response[i,1]==subjnameofconnectivity)
  if (i!=ind) cat("here", i)
  riskfactors[ind,]=response[ind, 2:(dim(response)[2])     ]
}


riskfactorind=riskfactors[,4]<2
riskfactorind=riskfactors>0
sum(riskfactorind)
riskfactors=riskfactors[riskfactorind,] # removing riskfactor 2,3


image=matrix(NA,  dim(connectivity)[3], len)

for (i in 1:dim(connectivity)[3]) {
  temp=connectivity[,,i]
  indexlower=lower.tri(temp, diag=FALSE)
  temp=temp[indexlower]
image[i,]=temp
}
image=image[riskfactorind,]

#recordzerocols
indd=0
for (i in 1:dim(image)[2]) if(sd(image[,i])==0 ) {indd=rbind(indd,i);  cat ( i , sd(image[,i]), "\n" );}
indd=indd[2:dim(indd)[1]]
image=image[,-indd]




#riskfactors=riskfactors[,c(1,3)]
#image=t(image);
#riskfactors=t(riskfactors);



#lets run
## Not run:
#install.packages("PMA")
library(PMA)



set.seed(3189) #for reproductivity

# Can run CCA with default settings, and can get e.g. 3 components
??glmnet
#out <- CCA(x=image,z=riskfactors,typex="standard",typez="standard")
#install.packages("glmnet")
library(glmnet)

#out=cv.glmnet(x=image,y=riskfactors, family = "multinomial")
out=cv.glmnet(x=image,y=riskfactors,  family = "binomial", type.measure = "class")
print(out)
plot(out)
coef=coef(out, s = "lambda.min")
coef=coef
sum(coef!=0)


#lambdas=seq(1:30)/30 
#lambdas=lambdas[2:(length(lambdas)-1)]
#perm.out <- CCA.permute(x=image,z=riskfactors,typex="standard",typez="standard",nperms=50, standardize=FALSE, penaltyxs=lambdas, penaltyzs=lambdas)
#perm.out <- CCA.permute(x=image,z=riskfactors,typex="standard",typez="standard",nperms=1000, standardize=FALSE, penaltyz=1)
# if you dont want to run the 1000 permuttion again just load the 1000 permutation R data and run from here:

#print(perm.out)
#plot(perm.out)
#out <- CCA(x=image,z=riskfactors,typex="standard",typez="standard",
#           penaltyx=perm.out$bestpenaltyx,penaltyz=1,
#           v=perm.out$v.init)

#print(out) # could do print(out,verbose=TRUE)
#print(image[out$u!=0]) 

u=coef
#v=out$v
sum(u==0)
#len=length(u)
sum(u!=0)

sum(u==0)/len #sparsity 

uout=matrix(NA, dim(u)[1]+length(indd),1 )
#put those zeros back
uout[indd]=0
uout[-indd]=coef
#make it square again



#take a look at zero rois through all subjects:

# ##### Example to reverse lower tri
# A=matrix(4*1:16,4,4)
# indexample=lower.tri(A, diag=FALSE)
# indexampletrue=which(indexample==TRUE)
# tempex=A[indexample]
# #position 2 of tempex is "#12" A(3,1)
# posindex=2
# tempex[posindex]
# indexampletrue[posindex]
# A[indexampletrue[posindex]]
# # multiple positions
# # position c(2,4) tempex of it are c(12,28)
# posindex=c(2,4)
# tempex[posindex]
# indexampletrue[posindex]
# A[indexampletrue[posindex]]
# # works just fine

indd
indexlowertrue=which(indexlower==TRUE)
temp[indd]
indexlowertrue[indd]

connectivityexample=connectivity[,,1]
connectivityexample[indexlowertrue[indd]] ##yes the're them
connectivityexample[indexlowertrue[indd]]="zeros" # lest make them known for a special word
indexofzeros=which(connectivityexample=="zeros", arr.ind=TRUE)

indexofzeros[,1]
indexofzeros

# #lest check really quick:
# for (j in 1:dim(indexofzeros)[1]) { 
# for (i in 1:dim(connectivity)[3]) {  cat(  "subject", i, "at position", indexofzeros[j,] , connectivity[matrix(c(indexofzeros[j,],i),1,3)] , "\n")
# }
# } ## yes theyre all odly zeros

#results of connectivities that matter:
nonzeroindex=which(uout!=0)
connectivityexample=connectivity[,,1]
connectivityexample[]=0
connectivitvals=connectivityexample
nonzerouout=uout[uout!=0]
for (i in 1:length(nonzeroindex)) {
  connectivityexample[indexlowertrue[nonzeroindex[i]]]=c("nonzero") # lest make them known for a special word
  connectivitvals[indexlowertrue[nonzeroindex[i]]]=nonzerouout[i] #store their coefitient values
}



library('igraph');
connectivitvalsones=connectivitvals
t=which(connectivitvalsones!=0, arr.ind=TRUE)
t <- cbind(t, connectivitvals[which(connectivitvals!=0,arr.ind=TRUE)]) 
t.graph=graph.data.frame(t,directed=F)
E(t.graph)$color <- ifelse(E(t.graph)$V3 > 0,'blue','red') 
#t.names <- colnames(cor.matrix)[as.numeric(V(t.graph)$name)]
minC <- rep(-Inf, vcount(t.graph))
maxC <- rep(Inf, vcount(t.graph))
minC[1] <- maxC[1] <- 0
l <- layout_with_fr(t.graph, minx=minC, maxx=maxC,
                    miny=minC, maxy=maxC)      

pathnames='/Users/ali/Desktop/mar/mice/mouse_anatomy.csv'
datanmes=read.csv(pathnames, header = TRUE, sep = ",", quote = "")
 datanmes$ROI

 par(mfrow=c(1,1))
 
#set.vertex.attribute(t.graph, "name", value=datanmes$ROI   )
plot(t.graph, layout=l, 
     rescale=T,
     asp=0,
     edge.arrow.size=0.1, 
     vertex.label.cex=0.8, 
     vertex.label.family="Helvetica",
     vertex.label.font=4,
     #vertex.label=t.names,
     vertex.shape="circle", 
     vertex.size=5, 
     vertex.color="deepskyblue2",
     vertex.label.color="black", 
     #edge.color=E(t.graph)$color, ##do not need this since E(t.graph)$color is already defined.
     edge.width=as.integer(cut(abs(E(t.graph)$V3), breaks = 5)))

connectivitvals=connectivitvals+t(connectivitvals) #symetric


nonzeroposition=which(connectivityexample=="nonzero", arr.ind=TRUE)
getwd()
filename=paste(getwd(), "/", "valandpos.mat", sep = "")
writeMat(filename, nonzeroposition = nonzeroposition, connectivitvals = connectivitvals , oddzeroposition=indexofzeros)




subnets=groups(components(t.graph))
subnetsresults=vector(mode = "list", length = length(subnets))
colsumabs=colSums(abs(connectivitvals))
colsum=colSums(connectivitvals)

for (i in 1:length(subnets)) {
  temp=subnets[[i]]
  temp=as.numeric(temp)
  net=matrix(NA,7,length(temp) )
  net[2,]=datanmes$ROI[temp]
  net[1,]=as.numeric(temp)
  net[3,]= as.numeric( colsumabs[temp]   )
  net[4,]= as.numeric( colsum[temp]   )
  tt=as.numeric(net[1,])
  #tt=c(1,200)
  indofleftright=tt>=166
  net[5,][indofleftright]="Right"
  net[5,][!indofleftright]="Left"
  net[6,]=sum(as.numeric(net[4,]))
  net[7,]=sum(as.numeric(net[3,]))
  #allpossibletemp=expand.grid(temp, temp)
  # for (j in 1:dim(allpossibletemp)[1]) {
  #   aaaa=allpossibletemp[j,]
  #   aaaa=unlist(aaaa, use.names=FALSE  )
  #   if (connectivitvals[ aaaa[1], aaaa[2] ] !=0 )  net[3,which(net[1,]==aaaa[1] )]= abs(connectivitvals[ aaaa[1], aaaa[2] ]) + abs(as.numeric(net[3,which(net[1,]==aaaa[1] )]))    
  # }
  subnetsresults[[i]]=net 
}



#for (i in 1:length(subnetsresults)) {
#  net=subnetsresults[i]
#  print(net[[1]][1:2,])
#}


for (i in 1:length(subnetsresults)) {
  net=subnetsresults[i]
  cat( i,'th sub-net: the summation of all edges in this sub-net is' ,sum(as.numeric(net[[1]][4,])), 'and summation of absolut values of all edges in this subnet is', sum(as.numeric(net[[1]][3,])),'\n')
  cat(  'the fsirst row is the Region #, second row is the name of Region, the third row is the sum of absulote values of the edges of each region, and the last row is the sum of edges of each region \n')
  print(net)
  cat( '\n \n \n')
}


capture.output(subnetsresults, file = "subnet.txt")


write.csv(subnetsresults, row.names = T)




#install.packages("xlsx")
library(xlsx)


for (i in 1:length(subnetsresults)){
  net=subnetsresults[[i]]
  write.xlsx2(net, "subnets.xlsx", sheetName =  paste0(i), append=TRUE )
}


