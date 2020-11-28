library(plyr)
library(reshape2)
library(ggplot2)

RemoveBackground_single <- function(background_profile_set, sig_profile,boundary){
  
  # Range of background
  centroid_background <- rowMeans(background_profile_set)
  sd_background <- apply(background_profile_set,1,sd)
  boundary_background <- centroid_background+boundary*sd_background
  
  
  diff_all_boundary <- sig_profile-boundary_background
  diff_all <- sig_profile-centroid_background
  
  diff_all[which(diff_all_boundary<0)] <- 0
  return(diff_all)
  
}

cos_similarity <- function(v1,v2){
  v1v2 <- sum(v1*v2)
  v1_length <- sqrt(sum(v1*v1))
  v2_length <- sqrt(sum(v2*v2))
  return(v1v2/v1_length/v2_length)
}

bootstrapGenomesfun <- function(genomes){
  
  return(apply(genomes, 2, function(x) rmultinom(1, sum(x), x)))
}

##############################################################################################
sub_catalogue=read.table("~/sub_catalogue_All.txt",sep="\t",header=TRUE)
muttype_template <- read.table("~/MutationType_template.txt", sep = "\t", header = T, as.is = T)

control_profile=sub_catalogue[,3:11]
row.names(control_profile)=sub_catalogue$MutationType
A3G_profile=sub_catalogue[,12:23]
row.names(A3G_profile)=sub_catalogue$MutationType

control_profile_percentage=control_profile #Percentage counts for control group#
control_profile_percentage[,1:dim(control_profile_percentage)[2]]=control_profile_percentage[,1:dim(control_profile_percentage)[2]]/colSums(control_profile_percentage[,1:dim(control_profile_percentage)[2]])[col(control_profile_percentage[,1:dim(control_profile_percentage)[2]])]

A3G_profile_percentage=A3G_profile #Percentage counts for A3G group#
A3G_profile_percentage[,1:dim(A3G_profile_percentage)[2]]=A3G_profile_percentage[,1:dim(A3G_profile_percentage)[2]]/colSums(A3G_profile_percentage[,1:dim(A3G_profile_percentage)[2]])[col(A3G_profile_percentage[,1:dim(A3G_profile_percentage)[2]])]

##Frobenuis distance##
control_profile_centroid=rowMeans(control_profile_percentage)
A3G_profile_centroid=rowMeans(A3G_profile_percentage)
clonecentroidDist=round(norm(as.matrix(control_profile_centroid-A3G_profile_centroid),"f"),digits=4)

control_percentage_all=NULL #for control subclone#
for (g in 1:10000){ #bootstrap resampling count of control samples, 10000 bootstrap#
  Control_FD=control_profile[,sample(1:dim(control_profile)[2],dim(control_profile)[2],replace=T)]
  Con_FD=bootstrapGenomesfun(Control_FD)
  Con_FD_percentage=Con_FD
  Con_FD_percentage[,1:dim(Con_FD_percentage)[2]]=Con_FD_percentage[,1:dim(Con_FD_percentage)[2]]/colSums(Con_FD_percentage[,1:dim(Con_FD_percentage)[2]])[col(Con_FD_percentage[,1:dim(Con_FD_percentage)[2]])]
  Con_FD_per=rowMeans(Con_FD_percentage)
  control_percentage_all=cbind(control_percentage_all,Con_FD_per)
}
control_percentage_all_diff <- control_percentage_all-control_profile_centroid
controlFrobeniusDist <- apply(control_percentage_all_diff,2, function(x) norm(as.matrix(x),"f"))
controlFrobeniusDist_threshold <- quantile(controlFrobeniusDist, probs = 0.95)

subclone_percentage_all=NULL #for A3G subclone#
for (j in 1:10000){ #bootstrap resampling count of A3G samples, 10000 bootstrap#
  RepCompound_FD=A3G_profile[,sample(1:dim(A3G_profile)[2],dim(A3G_profile)[2],replace=T)]
  FD=bootstrapGenomesfun(RepCompound_FD)
  FD_percentage=FD
  FD_percentage[,1:dim(FD_percentage)[2]]=FD_percentage[,1:dim(FD_percentage)[2]]/colSums(FD_percentage[,1:dim(FD_percentage)[2]])[col(FD_percentage[,1:dim(FD_percentage)[2]])]
  FD_per=rowMeans(FD_percentage)
  subclone_percentage_all=cbind(subclone_percentage_all,FD_per)
}
subclone_percentage_all_diff <- subclone_percentage_all-A3G_profile_centroid
SubcloneFrobeniusDist <- apply(subclone_percentage_all_diff,2, function(x) norm(as.matrix(x),"f"))
SubcloneFrobeniusDist_threshold <- quantile(SubcloneFrobeniusDist, probs = 0.95)

totalerror <- data.frame(cbind(controlFrobeniusDist,SubcloneFrobeniusDist),row.names = NULL)
names(totalerror) <- c("clone","subclone")
totalerror_melt <- melt(totalerror)
names(totalerror_melt) <- c("flag","FrobeniusDist")

binw <- 0.0005
mypalette <- c("lightpink2","seagreen2")
g1 <-ggplot(totalerror_melt, aes(x=FrobeniusDist,fill=flag)) + geom_histogram(binwidth=binw,alpha=1, position="identity")+scale_fill_manual(values=mypalette)

g1 <- g1+ geom_vline(aes(xintercept=controlFrobeniusDist_threshold), colour="red", linetype="dashed")
g1 <- g1+ geom_vline(aes(xintercept=SubcloneFrobeniusDist_threshold), colour="green", linetype="dashed")

g1 <- g1+ annotate("segment", x = clonecentroidDist, xend = clonecentroidDist, y = 200, yend = 0, colour="blue", size=1, arrow=arrow())
g1 <-g1 +theme(axis.text.x=element_text(size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))

#Signature generation#
result_profile=NULL
for(j in 1:dim(A3G_profile_percentage)[2]){
  a <- RemoveBackground_single(control_profile_percentage,A3G_profile_percentage[,j],1.65)
  result_profile <- cbind(result_profile,a)
}
min_b <- 1
max_b <- 0
for(j in 1:dim(A3G_profile_percentage)[2]){
  if((j+1)<=dim(A3G_profile_percentage)[2]){
    for(k in (j+1):dim(A3G_profile_percentage)[2]){
      c <- cos_similarity(result_profile[,j],result_profile[,k])
      min_b <- ifelse(c<min_b,c,min_b)
      max_b <- ifelse(c>max_b,c,max_b)
      print(paste0("min_b:",min_b,"; max_b:",max_b))
    }
  }
}
stability_sig <- c(min_b,max_b)

centroid_background_percentage=rowMeans(control_profile_percentage)
sd_background_percentage=apply(control_profile_percentage,1,sd)
boundary_background_percentage=centroid_background_percentage+1.65*sd_background_percentage #p=0.1#

diff_mean_all_percentage=NULL
bootstrapCompoundCount=NULL
bootstrapCompoindCountSu=NULL
g=c()
for (bt_num in 1:10000){
  RepCompound=A3G_profile[,sample(1:dim(A3G_profile)[2],dim(A3G_profile)[2],replace=T)]
  d=bootstrapGenomesfun(RepCompound)  #bootstrap resampling count of A3G samples, 10000 bootstrap#
  bootstrapCompoundCount=cbind(bootstrapCompoundCount,d)
  f=sum(d)/dim(d)[2]
  g=c(g,f)
  
  d_percentage=d
  d_percentage[,1:dim(d_percentage)[2]]=d_percentage[,1:dim(d_percentage)[2]]/colSums(d_percentage[,1:dim(d_percentage)[2]])[col(d_percentage[,1:dim(d_percentage)[2]])]
  
  diff_all_boundary_percentage <- rowMeans(d_percentage)-boundary_background_percentage
  diff_all_percentage <- rowMeans(d_percentage)-centroid_background_percentage
  
  diff_all_percentage[which(diff_all_boundary_percentage<0)] <- 0
  diff_mean_all_percentage <- cbind(diff_mean_all_percentage,diff_all_percentage)
  
}
bootstrapCompoundCountSum=matrix(g,nrow=1)
#average mutation burden of replicate is similar to orignal average mutation burden.
#here we used replicate's mutaiton burden * difference of percentge of substitutions to calculate the number of mutaiton of each substitution#

diff_mean_all=NULL
for (i in 1:10000){
  h=diff_mean_all_percentage[,i]*bootstrapCompoundCountSum[1,i]
  diff_mean_all=cbind(diff_mean_all,h)
}
diff_mean_all=as.data.frame(diff_mean_all)
diff_centroid_all=rowMeans(diff_mean_all)
diff_quantile_all_sd=apply(diff_mean_all,1,sd)

A3G_sig=data.frame(cbind(diff_centroid_all,diff_quantile_all_sd))
names(A3G_sig) <- c("centroid","sd")

A3G_sig$MutationType <- rownames(A3G_sig)
muts_basis_melt_summary[is.na(muts_basis_melt_summary)] <- 0
muts_basis_melt_summary <- muts_basis_melt_summary[order(muts_basis_melt_summary$Mutation),]

muts_basis_melt_summary$percentage <- muts_basis_melt_summary[,"centroid"]/sum(muts_basis_melt_summary[,"centroid"])
muts_basis_melt_summary$percentage_sd <- muts_basis_melt_summary[,"sd"]/sum(muts_basis_melt_summary[,"centroid"])

write.table(muts_basis_melt_summary,paste0("A3G", "_exposure.txt"),sep = "\t",row.names = F, col.names = T, quote = F)

mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")

p <- ggplot(data=muts_basis_melt_summary, aes(x=MutationType, y=centroid,fill=Mutation))+ geom_bar(stat="identity",position="dodge", width=.7)+xlab("Substitution Types")+ylab("Count")
p <- p+scale_x_discrete(limits = as.character(muts_basis_melt_summary$MutationType))+ggtitle(paste0("A3G_RB","_exposure","(stability=",round(max_b,2),")"))
p <- p+scale_fill_manual(values=mypalette)
p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=4.5,colour = "black"),
             axis.text.y=element_text(size=10,colour = "black"),
             axis.title.x = element_text(size=15),
             axis.title.y = element_text(size=15),
             plot.title = element_text(size=10),
             panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.border = element_rect(colour = "black", fill=NA))


