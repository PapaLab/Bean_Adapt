#GeneticLength(cM) Physicallength(Kb) Kb/cMPericentromere Kb/cMeuchromarticarms Kb/cMperchromosome
#Chr01 84.0 52183.5 5210 278 651.5
#Chr02 127.6 49033.7 3084 233 384.2
#Chr03 116.9 52218.6 3452 262 445.6
#Chr04 94.0 45793.2 4701 164 486.7
#Chr05 90.8 40237.5 2342 134 443.0
#Chr06 70.8 31973.2 6102 239 451.3
#Chr07 105.4 51698.4 9179 233 489.7
#Chr08 114.0 59634.6 6913 208 554.7
#Chr09 94.6 37399.6 3322 352 394.0
#Chr10 60.2 43213.2 5388 267 732.1
#Chr11 78.5 50203.6 5877 232 638.9
#Mean 94.3 46689.9 5052 237 515.6

#recombination rate kb/cM
schmutz_rec<-c(651.5,384.2,445.6,486.7,443.0,451.3,489.7,554.7,394.0,732.1,638.9)
#convert recrate in cM/Mb
schmutz_rec<-1/schmutz_rec*1000
#blair rec rate cM/Mb
blair_rec<-c(2.82,2.87,2.19,1.85,2.04,2.62,2.09,1.66,1.24,1.67,1.35)


get_input<-function(data=a,chr="all",n=-1,th=-1,avgr=1,posid=4){
if (chr=="all"){c1<-data}else{c1<-data[data[,1]==chr,]}
#first ind 10000 snp random
if (n>0){
	s1<-c1[sort(sample(nrow(c1),n)),c(1,3,posid)]
	}else{
	s1<-c1[,c(1,3,posid)]
}
s1[,3]<-as.character(s1[,3])
#perform thinnink in bp
if (th>0 & chr=="all"){print("Thinning is possible on separate chr only");return(-1)}
if (th>0){
	for (w in 1:nrow(s1)){
		if (w==1){
			start<-s1[w,2]
		}else{
			if ((s1[w,2]-start)<th){
				s1[w,3]<-"?"
			}else{
				start<-s1[w,2]			
			}
		}		
	}
}
s1[,2]<-s1[,2]/1000000/100*avgr
#s1[,2]<-s1[,2]/1000000/100*2.13#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0189597
#ancestry conversion
s1<-s1[s1[,3]!="?",]
s1[grep("A",s1[,3]),3]<-"0"
s1[grep("M",s1[,3]),3]<-"1"
#last fixing
s1[,3]<-as.numeric(s1[,3])
s1[,1]<-as.numeric(sub("Chr","",s1[,1]))
out<-list()
out[[1]]<-s1
out[[2]]<-(table(s1[,3])/nrow(s1))
return(out)
}



library(junctions)
a<-read.table(gzfile("/media/benazzo/My Passport/progetti/beans_modern/chromopainter/selection_scan/alternative/new_thresholds/individual_attr.txt.gz"),h=F)
#filter centromeric regions
aa<-a
aa[(aa[,1]=="Chr01")&(aa[,3]>=12200000)&(aa[,3]<=19900000),4:231]<-"?"
aa[(aa[,1]=="Chr02")&(aa[,3]>=5400000)&(aa[,3]<=10000000),4:231]<-"?"
aa[(aa[,1]=="Chr03")&(aa[,3]>=14800000)&(aa[,3]<=16900000),4:231]<-"?"
aa[(aa[,1]=="Chr04")&(aa[,3]>=15700000)&(aa[,3]<=22200000),4:231]<-"?"
aa[(aa[,1]=="Chr05")&(aa[,3]>=15300000)&(aa[,3]<=22700000),4:231]<-"?"
aa[(aa[,1]=="Chr06")&(aa[,3]>=2600000)&(aa[,3]<=2700000),4:231]<-"?"
aa[(aa[,1]=="Chr07")&(aa[,3]>=16700000)&(aa[,3]<=30300000),4:231]<-"?"
aa[(aa[,1]=="Chr08")&(aa[,3]>=24300000)&(aa[,3]<=38200000),4:231]<-"?"
aa[(aa[,1]=="Chr09")&(aa[,3]>=1500000)&(aa[,3]<=5800000),4:231]<-"?"
aa[(aa[,1]=="Chr10")&(aa[,3]>=30600000)&(aa[,3]<=31300000),4:231]<-"?"
aa[(aa[,1]=="Chr11")&(aa[,3]>=16000000)&(aa[,3]<=17100000),4:231]<-"?"

c1<-get_input(data=a,chr="Chr01",n=10000,posid=14)
c2<-get_input(data=a,chr="Chr02",n=10000,posid=14)
c3<-get_input(data=a,chr="Chr03",n=10000,posid=14)
c4<-get_input(data=a,chr="Chr04",n=10000,posid=14)
c5<-get_input(data=a,chr="Chr05",n=10000,posid=14)
c6<-get_input(data=a,chr="Chr06",n=10000,posid=14)
c7<-get_input(data=a,chr="Chr07",n=10000,posid=14)
c8<-get_input(data=a,chr="Chr08",n=10000,posid=14)
c9<-get_input(data=a,chr="Chr09",n=10000,posid=14)
c10<-get_input(data=a,chr="Chr10",n=10000,posid=14)
c11<-get_input(data=a,chr="Chr11",n=10000,posid=14)
ca<-rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)



#P1167 P1168 P1169 P1171 P1172 P1173 P1174 P1177 P1179 P1180 P1181 P1182 P1183 P1184 P1185 P1188 P1189 P1190 P1191 P1192 P1193 P1194 P1195 P1196 P1197 P1198 P1199 P1200 P1201 P1202 P1203 P1204 P1205 P1206 P1207 P1208 P1214 PECE179 WGS1 WGS10 WGS12 WGS13 WGS14 WGS15 WGS16 WGS17 WGS2 WGS3 WGS6 WGS7 WGS8 WGS9 P1100 P1101 P1102 P1103 P1104 P1107 P1108 P1109 P1111 P1112 P1113 P1114 P1115 P1116 P1117 P1118 P1119 P1120 P1121 P1122 P1123 P1124 P1125 P1126 P1127 P1128 P1129 P1130 P1131 P1132 P1133 P1134 P1135 P1136 P1137 P1138 P1139 P1140 P1141 P1142 P1143 P1144 P1145 P1147 P1148 P1149 P1150 P1151 P1152 P1153 P1154 P1156 P1157 P1158 P1159 P1160 P1161 P1162 P1163 P1164 P1165 P1166
#load the previous line
b<-unlist(read.table("clipboard",sep=" ",h=F,colClasses="character"))
b<-sub("P","",b)
#and<-read.table("/media/benazzo/My Passport/progetti/beans_modern/chromopainter/time_admixture/eu_and.ind.txt",colClasses="character")[,1]
#mes<-read.table("/media/benazzo/My Passport/progetti/beans_modern/chromopainter/time_admixture/eu_mes.ind.txt",colClasses="character")[,1]

#ind having p admix >=10%
#and<-read.table("/media/benazzo/My Passport/progetti/beans_modern/chromopainter/time_admixture/eu_and_10perc.ind.txt",colClasses="character")[,1]
#mes<-read.table("/media/benazzo/My Passport/progetti/beans_modern/chromopainter/time_admixture/eu_mes_10perc.ind.txt",colClasses="character")[,1]
idand<-seq(1,228,by=2)[which(b%in%and)]+3
idmes<-seq(1,228,by=2)[which(b%in%mes)]+3


#final analysis in andean inds
risand<-matrix(ncol=length(idand),nrow=3)
row.names(risand)<-c("250","500","1000")
colnames(risand)<-b[which(b%in%and)]
st<-list()
for (x in 1:length(idand)){
print(x)
ca<-list()
st[[1]]<-get_input(data=a,chr="Chr01",n=50000,th=1000,avgr=schmutz_rec[1],posid=idand[x])
st[[2]]<-get_input(data=a,chr="Chr02",n=50000,th=1000,avgr=schmutz_rec[2],posid=idand[x])
st[[3]]<-get_input(data=a,chr="Chr03",n=50000,th=1000,avgr=schmutz_rec[3],posid=idand[x])
st[[4]]<-get_input(data=a,chr="Chr04",n=50000,th=1000,avgr=schmutz_rec[4],posid=idand[x])
st[[5]]<-get_input(data=a,chr="Chr05",n=50000,th=1000,avgr=schmutz_rec[5],posid=idand[x])
st[[6]]<-get_input(data=a,chr="Chr06",n=50000,th=1000,avgr=schmutz_rec[6],posid=idand[x])
st[[7]]<-get_input(data=a,chr="Chr07",n=50000,th=1000,avgr=schmutz_rec[7],posid=idand[x])
st[[8]]<-get_input(data=a,chr="Chr08",n=50000,th=1000,avgr=schmutz_rec[8],posid=idand[x])
st[[9]]<-get_input(data=a,chr="Chr09",n=50000,th=1000,avgr=schmutz_rec[9],posid=idand[x])
st[[10]]<-get_input(data=a,chr="Chr10",n=50000,th=1000,avgr=schmutz_rec[10],posid=idand[x])
st[[11]]<-get_input(data=a,chr="Chr11",n=50000,th=1000,avgr=schmutz_rec[11],posid=idand[x])
for (i in 1:length(st)){
	if (st[[i]][[2]][1] >= 0.05 | st[[i]][[2]][1] <= 0.95){
		if (length(ca)==0){ca<-st[[i]][[1]]}else{ca<-rbind(ca,st[[i]][[1]])}	
	}
}
ph<-sum(ca[,3])/nrow(ca)
risand[1,x]<-estimate_time_haploid(ancestry_matrix = ca, N = 250,freq_ancestor_1 = ph)$time
risand[2,x]<-estimate_time_haploid(ancestry_matrix = ca, N = 500,freq_ancestor_1 = ph)$time
risand[3,x]<-estimate_time_haploid(ancestry_matrix = ca, N = 1000,freq_ancestor_1 = ph)$time
}

#final analysis in mesoamerican inds
rismes<-matrix(ncol=length(idmes),nrow=3)
row.names(rismes)<-c("250","500","1000")
colnames(rismes)<-b[which(b%in%mes)]
st<-list()
for (x in 1:length(idmes)){
print(x)
ca<-list()
st[[1]]<-get_input(data=a,chr="Chr01",n=50000,th=1000,avgr=schmutz_rec[1],posid=idmes[x])
st[[2]]<-get_input(data=a,chr="Chr02",n=50000,th=1000,avgr=schmutz_rec[2],posid=idmes[x])
st[[3]]<-get_input(data=a,chr="Chr03",n=50000,th=1000,avgr=schmutz_rec[3],posid=idmes[x])
st[[4]]<-get_input(data=a,chr="Chr04",n=50000,th=1000,avgr=schmutz_rec[4],posid=idmes[x])
st[[5]]<-get_input(data=a,chr="Chr05",n=50000,th=1000,avgr=schmutz_rec[5],posid=idmes[x])
st[[6]]<-get_input(data=a,chr="Chr06",n=50000,th=1000,avgr=schmutz_rec[6],posid=idmes[x])
st[[7]]<-get_input(data=a,chr="Chr07",n=50000,th=1000,avgr=schmutz_rec[7],posid=idmes[x])
st[[8]]<-get_input(data=a,chr="Chr08",n=50000,th=1000,avgr=schmutz_rec[8],posid=idmes[x])
st[[9]]<-get_input(data=a,chr="Chr09",n=50000,th=1000,avgr=schmutz_rec[9],posid=idmes[x])
st[[10]]<-get_input(data=a,chr="Chr10",n=50000,th=1000,avgr=schmutz_rec[10],posid=idmes[x])
st[[11]]<-get_input(data=a,chr="Chr11",n=50000,th=1000,avgr=schmutz_rec[11],posid=idmes[x])
for (i in 1:length(st)){
	if (st[[i]][[2]][1] >= 0.05 | st[[i]][[2]][1] <= 0.95){
		if (length(ca)==0){ca<-st[[i]][[1]]}else{ca<-rbind(ca,st[[i]][[1]])}	
	}
}
ph<-sum(ca[,3])/nrow(ca)
rismes[1,x]<-estimate_time_haploid(ancestry_matrix = ca, N = 250,freq_ancestor_1 = ph)$time
rismes[2,x]<-estimate_time_haploid(ancestry_matrix = ca, N = 500,freq_ancestor_1 = ph)$time
rismes[3,x]<-estimate_time_haploid(ancestry_matrix = ca, N = 1000,freq_ancestor_1 = ph)$time
}







