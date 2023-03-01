#This script is to merge the rehions in LD
#to call it:
#Rscript --vanilla ./merge_regions.r and-t10k.interchrom.hap.ld.gz 100000 | gzip -c > and-t10k.interchrom.hap.ld.merged.tab.gz
args<-commandArgs(trailing=T)
tab<-read.table(gzfile(args[1]),header=T,colClasses="numeric")
di<-as.numeric(args[2])
um1r<-paste(tab[,1],tab[,2],sep="_")
um1<-unique(um1r)
m1<-matrix(ncol=4,nrow=length(um1))
colnames(m1)<-c("CHR1","POSI","POSF","RID")
ne<-0
for (i in 1:length(um1)){
	#print(i)
	if (i==1){
		chrp<-as.numeric(unlist(strsplit(um1[i],split="_"))[1])
		posp<-as.numeric(unlist(strsplit(um1[i],split="_"))[2])
		id<-which(um1r==um1[i])
	}else{
		chr<-as.numeric(unlist(strsplit(um1[i],split="_"))[1])
		pos<-as.numeric(unlist(strsplit(um1[i],split="_"))[2])
		if (chrp==chr & (pos-posp)<=di){
			id<-c(id,which(um1r==um1[i]))
			chrp<-chr
			posp<-pos
		}else{
			ne<-ne+1
			m1[ne,]<-c(chrp,tab[id[1],2],tab[id[length(id)],2],paste(id,collapse=";"))
			chrp<-chr
			posp<-pos
			id<-which(um1r==um1[i])
		}			
	}
}
for (i in 1:ne){
	#print(i)
	id<-as.numeric(unlist(strsplit(m1[i,"RID"],split=";")))
	tab1<-tab[id,]
	tab2<-c()
	for (k in 1:11){
		if (sum(tab1[,3]==k)>0){
			tmpid<-sort(tab1[tab1[,3]==k,4],index.return=T)
			if (length(tab2)==0){
				tab2<-cbind(rep(k,sum(tab1[,3]==k)),tmpid$x,tab1[tab1[,3]==k,6][tmpid$ix])
			}else{
				tab2<-rbind(tab2,cbind(rep(k,sum(tab1[,3]==k)),tmpid$x,tab1[tab1[,3]==k,6][tmpid$ix]))
			}
		}
			
	}
	for (k in 1:nrow(tab2)){
		if (k==1){
			chrp<-tab2[k,1]
			posp<-tab2[k,2]
			posi<-posp
			posf<-posp
		}else{
			chr<-tab2[k,1]
			pos<-tab2[k,2]
			if (chrp==chr & (pos-posp)<=di){
				posf<-pos
				posp<-pos
			}else{
				write(paste(m1[i,1],m1[i,2],m1[i,3],chrp,posi,posf,sep="\t"),stdout())
				chrp<-chr
				posp<-pos
				posi<-pos
				posf<-pos
			}
			if (k==nrow(tab2)){
				write(paste(m1[i,1],m1[i,2],m1[i,3],chrp,posi,posf,sep="\t"),stdout())
			}
		}
	}
}
