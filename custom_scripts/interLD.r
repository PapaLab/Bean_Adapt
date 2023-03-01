#This script is for compute the Linkage Disequilibrium
#find overlaps between ribbon links comparing each row of the two files
find_overlap_fast<-function(x,y,slop=0,outstd=F){
	#out<-vector(mode = "list", length = 100000)
	out<-list()
	x<-as.matrix(x);x[,2]<-x[,2]-slop;x[,3]<-x[,3]+slop;x[,5]<-x[,5]-slop;x[,6]<-x[,6]+slop
	y<-as.matrix(y);y[,2]<-y[,2]-slop;y[,3]<-y[,3]+slop;y[,5]<-y[,5]-slop;y[,6]<-y[,6]+slop
	j<-0
	pos<-c()
	for (i in 1:nrow(x)){
		#print(i)
		a<-unlist(x[i,])	
		#check correspondace between chromosomes
		tmp<-y[(y[,1]==a[1] & y[,4]==a[4]),]
		if (length(tmp)==0){next}
		if(length(tmp)>6){
			tmp1<-tmp[((a[2]<=tmp[,2]&a[3]>=tmp[,2])|(a[2]<=tmp[,3]&a[3]>=tmp[,3])|(a[2]>=tmp[,2]&a[3]<=tmp[,3])|(a[2]<=tmp[,2]&a[3]>=tmp[,3])),]
		}else{
			tmp1<-tmp[((a[2]<=tmp[2]&a[3]>=tmp[2])|(a[2]<=tmp[3]&a[3]>=tmp[3])|(a[2]>=tmp[2]&a[3]<=tmp[3])|(a[2]<=tmp[2]&a[3]>=tmp[3]))]
		}
		if (length(tmp1)==0){next}
		if(length(tmp1)>6){
			tmp2<-tmp1[((a[5]<=tmp1[,5]&a[6]>=tmp1[,5])|(a[5]<=tmp1[,6]&a[6]>=tmp1[,6])|(a[5]>=tmp1[,5]&a[6]<=tmp1[,6])|(a[5]<=tmp1[,5]&a[6]>=tmp1[,6])),]
		}else{
			tmp2<-tmp1[((a[5]<=tmp1[5]&a[6]>=tmp1[5])|(a[5]<=tmp1[6]&a[6]>=tmp1[6])|(a[5]>=tmp1[5]&a[6]<=tmp1[6])|(a[5]<=tmp1[5]&a[6]>=tmp1[6]))]		
		}
		if (length(tmp2)==0){next}
		if(outstd){
			write(paste(c(a,tmp2),collapse="\t"),stdout())
		}else{
			#print(a)
			#print(tmp2)
			#print(length(tmp2))
			if(length(tmp2)<=6){
				j<-j+1
				out[[j]]<-c(a,tmp2)
				pos[j]<-i
			}else{
				for (w in 1:nrow(tmp2)){
					j<-j+1
					out[[j]]<-c(a,tmp2[w,])
					pos[j]<-i
				}
			}
		}
	}
	if (outstd==F){
		names(out)<-as.character(pos)
		return(out)
	}
}


fix_chr_names<-function(x){
	y<-x[,1]
	z<-x[,4]
	for (i in 1:nrow(x)){
		if (x[i,1]<10){y[i]<-paste("Chr0",as.character(x[i,1]),sep="")}else{y[i]<-paste("Chr",as.character(x[i,1]),sep="")}
		if (x[i,4]<10){z[i]<-paste("Chr0",as.character(x[i,4]),sep="")}else{z[i]<-paste("Chr",as.character(x[i,4]),sep="")}
	}
	x[,1]<-y
	x[,4]<-z
	return(x)	
}

#
#o<-find_overlap_fast(band,beuand,slop=100000)
#oo<-find_overlap_fast(bmes,beumes,slop=100000)
#am<-matrix(unlist(o),ncol=12,byrow=T)[,1:6]
#hist((as.vector(matrix(unlist(oo),ncol=12,byrow=T)[,c(1,4)])),breaks=10,freq=F,ylim=c(0,1))
#barplot(table(as.vector(matrix(unlist(oo),ncol=12,byrow=T)[,c(1,4)])))
#barplot(table(as.vector((am)[,c(1,4)]))/(nrow(am)*2))

require("RCircos")
#rw is the min width of the regions that are connected; the source and the sink regions should be at least rw bp to be considered.
rw<-500000
g<-read.table("pvulg-cent-pericent.genome",h=T)

eur<-read.table(gzfile("eur-t25k.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(eur)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
feur<-((eur[,3]-eur[,2])>=rw) & ((eur[,6]-eur[,5])>=rw)
beur<-as.data.frame(eur[feur,])

ame<-read.table(gzfile("ame-t25k.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(ame)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
fame<-((ame[,3]-ame[,2])>=rw) & ((ame[,6]-ame[,5])>=rw)
bame<-as.data.frame(ame[fame,])

and<-read.table(gzfile("and-t25k.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(and)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
fand<-((and[,3]-and[,2])>=rw) & ((and[,6]-and[,5])>=rw)
band<-as.data.frame(and[fand,])

euand<-read.table(gzfile("eu-and-t25k.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(euand)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
feuand<-((euand[,3]-euand[,2])>=rw) & ((euand[,6]-euand[,5])>=rw)
beuand<-as.data.frame(euand[feuand,])

mes<-read.table(gzfile("mes-t25k.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(mes)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
fmes<-((mes[,3]-mes[,2])>=rw) & ((mes[,6]-mes[,5])>=rw)
bmes<-as.data.frame(mes[fmes,])

eumes<-read.table(gzfile("eu-mes-t25k.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(eumes)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
feumes<-((eumes[,3]-eumes[,2])>=rw) & ((eumes[,6]-eumes[,5])>=rw)
beumes<-as.data.frame(eumes[feumes,])

#find overlaps between mesoamerica and european mesoamerican individuals (beumes:583, bmes:5286) ##(beumes:234, bmes:2280)
o<-find_overlap_fast(beumes,bmes)
#389 overlaps, but only 286 are unique (one link can overlap multiple links in the other group) ##34 overlaps, but only 32 are unique (one link can overlap multiple links in the other group)
#and now extract private links of beumes (297) ##and now extract private links of beumes (202)
popm<-beumes[-as.numeric(unique(names(o))),]
#and now extract shared links of beumes (286) ##and now extract shared links of beumes (32)
posm<-beumes[as.numeric(unique(names(o))),]
##find overlaps between andes and european andean individuals (beuand: 299, band:8447)#find overlaps between andes and european andean individuals (beuand: 275, band:4971)
o1<-find_overlap_fast(beuand,band)
##274 overlaps, but only 211 are unique (one link can overlap multiple links in the other group)#121 overlaps, but only 104 are unique (one link can overlap multiple links in the other group)
##and now extract private links of beuand (88) #and now extract private links of beuand (171) 
popa<-beuand[-as.numeric(unique(names(o1))),]
##and now extract shared links of beuand (211)#and now extract shared links of beuand (104)
posa<-beuand[as.numeric(unique(names(o1))),]



#and now the same groups using regions under selection
rw1<-50000
haand<-read.table(gzfile("and-t25k-hapflk-allk.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(haand)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
hafand<-((haand[,3]-haand[,2])>=rw1) & ((haand[,6]-haand[,5])>=rw1)
haband<-as.data.frame(haand[hafand,])

haeuand<-read.table(gzfile("eu-and-t25k-hapflk-allk.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(haeuand)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
hafeuand<-((haeuand[,3]-haeuand[,2])>=rw1) & ((haeuand[,6]-haeuand[,5])>=rw1)
habeuand<-as.data.frame(haeuand[hafeuand,])

hames<-read.table(gzfile("mes-t25k-hapflk-allk.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(hames)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
hafmes<-((hames[,3]-hames[,2])>=rw1) & ((hames[,6]-hames[,5])>=rw1)
habmes<-as.data.frame(hames[hafmes,])

haeumes<-read.table(gzfile("eu-mes-t25k-hapflk-allk.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(haeumes)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
hafeumes<-((haeumes[,3]-haeumes[,2])>=rw1) & ((haeumes[,6]-haeumes[,5])>=rw1)
habeumes<-as.data.frame(haeumes[hafeumes,])

hiand<-read.table(gzfile("and-t25k-hapflk-intersk.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(hiand)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
hifand<-((hiand[,3]-hiand[,2])>=rw1) & ((hiand[,6]-hiand[,5])>=rw1)
hiband<-as.data.frame(hiand[hifand,])

hieuand<-read.table(gzfile("eu-and-t25k-hapflk-intersk.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(hieuand)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
hifeuand<-((hieuand[,3]-hieuand[,2])>=rw1) & ((hieuand[,6]-hieuand[,5])>=rw1)
hibeuand<-as.data.frame(hieuand[hifeuand,])

himes<-read.table(gzfile("mes-t25k-hapflk-intersk.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(himes)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
hifmes<-((himes[,3]-himes[,2])>=rw1) & ((himes[,6]-himes[,5])>=rw1)
hibmes<-as.data.frame(himes[hifmes,])

hieumes<-read.table(gzfile("eu-mes-t25k-hapflk-intersk.interchrom.hap.ld.merged100kb.tab.gz"))
colnames(hieumes)<-c("Chromosome","chromStart","chromEnd","Chromosome.1","chromStart.1","chromEnd.1")
hifeumes<-((hieumes[,3]-hieumes[,2])>=rw1) & ((hieumes[,6]-hieumes[,5])>=rw1)
hibeumes<-as.data.frame(hieumes[hifeumes,])

#find overlaps between mesoamerica and european mesoamerican individuals (habeumes:137, habmes:162)##(habeumes:119, habmes:153)
hao<-find_overlap_fast(habeumes,habmes)
#66 overlaps, but only 65 are unique (one link can overlap multiple links in the other group)
#and now extract private links of beumes (72)
hapopm<-habeumes[-as.numeric(unique(names(hao))),]
#and now extract shared links of beumes (65)
haposm<-habeumes[as.numeric(unique(names(hao))),]
#find overlaps between andes and european andean individuals (habeuand: 18, haband:221) ##(beuand: 35, band:147)
hao1<-find_overlap_fast(habeuand,haband)
#17 overlaps, but only 17 are unique (one link can overlap multiple links in the other group)
#and now extract private links of beuand (1) 
hapopa<-habeuand[-as.numeric(unique(names(hao1))),]
#and now extract shared links of beuand (17)
haposa<-habeuand[as.numeric(unique(names(hao1))),]

#find overlaps between mesoamerica and european mesoamerican individuals (hibeumes:74, hibmes:106)##(habeumes:82, habmes:97)
hio<-find_overlap_fast(hibeumes,hibmes)
#44 overlaps, but only 44 are unique (one link can overlap multiple links in the other group)
#and now extract private links of beumes (30)
hipopm<-hibeumes[-as.numeric(unique(names(hio))),]
#and now extract shared links of beumes (44)
hiposm<-hibeumes[as.numeric(unique(names(hio))),]
#find overlaps between andes and european andean individuals (hibeuand: 17, hiband:121)##(beuand: 20, band:76)
hio1<-find_overlap_fast(hibeuand,hiband)
#17 overlaps, but only 17 are unique (one link can overlap multiple links in the other group)
#and now extract private links of beuand (0) 
hipopa<-hibeuand[-as.numeric(unique(names(hio1))),]
#and now extract shared links of beuand (17)
hiposa<-hibeuand[as.numeric(unique(names(hio1))),]


#do circos
chr.exclude <- NULL;
cyto.info <- g
tracks.inside <- 10
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside);

pdf("circos_plots_4groups_100kbmerging_500kbsize.m0.05_link.pdf",width=12)
par(mfrow=c(2,4))
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(bmes), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(bmes)))
text(x=0,y=2.5,labels="MES",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(beumes), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(beumes)))
text(x=0,y=2.5,labels="EU-MES",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(popm), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(popm)))
text(x=0,y=2.5,labels="Private (EU-MES)",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(posm), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(posm)))
text(x=0,y=2.5,labels="Shared (EU-MES)",cex=1.5)


RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(band), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(band)))
text(x=0,y=2.5,labels="AND",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(beuand), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(beuand)))
text(x=0,y=2.5,labels="EU-AND",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(popa), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(popa)))
text(x=0,y=2.5,labels="Private (EU-AND)",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(posa), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(posa)))
text(x=0,y=2.5,labels="Shared (EU-AND)",cex=1.5)

dev.off()


#plot regions under selection
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside);

pdf("circos_plots_4groups_100kbmerging_50kbsize_hapflk-allk.m0.05_link.pdf",width=12)
par(mfrow=c(2,4))
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(habmes), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(habmes)))
text(x=0,y=2.5,labels="MES",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(habeumes), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(habeumes)))
text(x=0,y=2.5,labels="EU-MES",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hapopm), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hapopm)))
text(x=0,y=2.5,labels="Private (EU-MES)",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(haposm), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(haposm)))
text(x=0,y=2.5,labels="Shared (EU-MES)",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(haband), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(haband)))
text(x=0,y=2.5,labels="AND",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(habeuand), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(habeuand)))
text(x=0,y=2.5,labels="EU-AND",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hapopa), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hapopa)))
text(x=0,y=2.5,labels="Private (EU-AND)",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(haposa), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(haposa)))
text(x=0,y=2.5,labels="Shared (EU-AND)",cex=1.5)

dev.off()



RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside);

pdf("circos_plots_4groups_100kbmerging_50kbsize_hapflk-intersk.m0.05_link.pdf",width=12)
par(mfrow=c(2,4))
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hibmes), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hibmes)))
text(x=0,y=2.5,labels="MES",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hibeumes), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hibeumes)))
text(x=0,y=2.5,labels="EU-MES",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hipopm), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hipopm)))
text(x=0,y=2.5,labels="Private (EU-MES)",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hiposm), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hiposm)))
text(x=0,y=2.5,labels="Shared (EU-MES)",cex=1.5)


RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hiband), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hiband)))
text(x=0,y=2.5,labels="AND",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hibeuand), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hibeuand)))
text(x=0,y=2.5,labels="EU-AND",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hipopa), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hipopa)))
text(x=0,y=2.5,labels="Private (EU-AND)",cex=1.5)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hiposa), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hiposa)))
text(x=0,y=2.5,labels="Shared (EU-AND)",cex=1.5)
dev.off()




o2<-find_overlap_fast(bame,beur)
#11648 overlaps
#and now extract private links of bame (7736) 
popame<-bame[-as.numeric(unique(names(o2))),]
#and now extract shared links of bame (6088)
posame<-bame[as.numeric(unique(names(o2))),]
#and now extract private links of beur (6088) 
popeur<-beur[-as.numeric(unique(names(o2))),]
#and now extract shared links of beur (6813) 
poseur<-beur[as.numeric(unique(names(o2))),]


#plot of private sites
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside);

pdf("circos_plots_6groups_onlyprivate_100kbmerging_50kbsize.m0.05_link_last.pdf",width=12)
par(mfcol=c(2,3))

#private america vs private eur
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(popame), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(popm)))
text(x=0,y=2.7,labels="AMERICA (n=7,736)",cex=1)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(popeur), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(popa)))
text(x=0,y=2.7,labels="EUROPE (n=6,813)",cex=1)

#global within europe
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(popm), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(popm)))
text(x=0,y=2.7,labels="EU_M (A,n=297)",cex=1)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(popa), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(popa)))
text(x=0,y=2.7,labels="EU_A (A,n=88)",cex=1)

#1-selection allk, this  set includes also the intersection dataset
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hapopm), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hapopm)))
text(x=0,y=2.7,labels="EU_M (S,n=72)",cex=1)

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
RCircos.Link.Plot(fix_chr_names(hapopa), track.num=1, by.chromosome=F,lineWidth=rep(2, nrow(hapopa)))
text(x=0,y=2.7,labels="EU_A (S,n=1)",cex=1)

dev.off()


##RCircos.Link.Plot(fix_chr_names(hiposa), track.num=1, by.chromosome=F,lineWidth=rep(5, nrow(hiposa)))

####excluding chromosomes
#RCircos.Set.Core.Components(cyto.info, chr.exclude=3:11,tracks.inside, tracks.outside)
#RCircos.Set.Plot.Area();
#RCircos.Chromosome.Ideogram.Plot()
#RCircos.Ribbon.Plot(b[b[,1]==1 & b[,4]==2,][1:10,], track.num=1, by.chromosome=F,twist=T)

write.table(popm,"Private_EU_M.m0.05.txt",quote=F,row.names=F)
write.table(popa,"Private_EU_A.m0.05.txt",quote=F,row.names=F)
write.table(hapopm,"Private_EU_M_extended.m0.05.txt",quote=F,row.names=F)
write.table(hapopa,"Private_EU_A_extended.m0.05.txt",quote=F,row.names=F)
write.table(hipopm,"Private_EU_M_restricted.m0.05.txt",quote=F,row.names=F)
write.table(hipopa,"Private_EU_A_restricted.m0.05.txt",quote=F,row.names=F)
