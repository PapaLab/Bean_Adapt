#This script runs chromopainter 
#get id codes for each ind for each file
for i in $(ls ./results/*.gz);do NAME=$(zcat $i|head -1|cut -f 4); echo "$i $NAME"; done > ind_ids.txt
grep -w -f /media/benazzo/dati1/beans_modern/finestructure/selected_genes/trees/1kb_windows_refined/eu_individuals.txt ind_ids.txt
# To pull out the samples and to compute the statistics
python zpaste.py $(grep -w -f /media/benazzo/dati1/beans_modern/finestructure/selected_genes/trees/1kb_windows_refined/eu_individuals.txt ind_ids.txt | cut -f 1 -d " " |awk '{printf "%s ",$0}END{printf "\n"}'| sed 's/.$//')|cut -f  $(cat <(echo "1") <(echo "2") <(echo "3") <(seq 6 7 1077) <(seq 7 7 1078)|sort -n|awk '{printf "%s,",$0}END{printf "\n"}'| sed 's/.$//')| awk 'BEGIN{printf "CHR\tSTART\tEND\tFAND\tFMES\tFMISS\n"}{tota=0;totm=0;totq=0;for (i=4;i<=NF;i++){if ($i=="A1"||$i=="A2"||$i=="A3"||$i=="AND"){tota=tota+1};if ($i=="M1"||$i=="M2"||$i=="MES"){totm=totm+1};if ($i=="?"){totq=totq+1}}printf "%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,tota/(NF-3),totm/(NF-3),totq/(NF-3)}'|awk '{if ($4>=0.811){printf "%s\tAND-->MES\n",$0}else if($5>=0.688){printf "%s\tMES-->AND\n",$0}}' > AND.MES.introgression.raw
cat AND.MES.introgression.raw|sed '1d'|/media/benazzo/dati/progs_genetica/bedtools2-master/bin/bedtools slop -b 2500 -g ../../../refgenome_annotations/Pvulgaris_442_v2.1.genome |/media/benazzo/dati/progs_genetica/bedtools2-master/bin/bedtools merge -d 10000 -i - |awk '{printf "%s\t%s\n",$0,$3-$2}'>put_reg_sl2.5k_mer10k.bed
/media/benazzo/dati/progs_genetica/bedtools2-master/bin/bedtools intersect -wa -wb -a put_reg_sl2.5k_mer10k.bed -b <(sed '1d' AND.MES.introgression.raw) |/media/benazzo/dati/progs_genetica/bedtools2-master/bin/bedtools groupby -g 1,2,3,4,11 -c 8,9,10,11 -o mean,mean,mean,count|awk '{if ($9>2)print $0}' > sel_reg_sl2.5k_mer10k_h2snps.bed




#extract individual components
#colnames
P1167 P1168 P1169 P1171 P1172 P1173 P1174 P1177 P1179 P1180 P1181 P1182 P1183 P1184 P1185 P1188 P1189 P1190 P1191 P1192 P1193 P1194 P1195 P1196 P1197 P1198 P1199 P1200 P1201 P1202 P1203 P1204 P1205 P1206 P1207 P1208 P1214 PECE179 WGS1 WGS10 WGS12 WGS13 WGS14 WGS15 WGS16 WGS17 WGS2 WGS3 WGS6 WGS7 WGS8 WGS9 P1100 P1101 P1102 P1103 P1104 P1107 P1108 P1109 P1111 P1112 P1113 P1114 P1115 P1116 P1117 P1118 P1119 P1120 P1121 P1122 P1123 P1124 P1125 P1126 P1127 P1128 P1129 P1130 P1131 P1132 P1133 P1134 P1135 P1136 P1137 P1138 P1139 P1140 P1141 P1142 P1143 P1144 P1145 P1147 P1148 P1149 P1150 P1151 P1152 P1153 P1154 P1156 P1157 P1158 P1159 P1160 P1161 P1162 P1163 P1164 P1165 P1166
#get european genotypes 
python zpaste.py $(grep -w -f /media/benazzo/dati1/beans_modern/finestructure/selected_genes/trees/1kb_windows_refined/eu_individuals.txt ind_ids.txt | cut -f 1 -d " " |awk '{printf "%s ",$0}END{printf "\n"}'| sed 's/.$//')|cut -f  $(cat <(echo "1") <(echo "2") <(echo "3") <(seq 6 7 1077) <(seq 7 7 1078)|sort -n|awk '{printf "%s,",$0}END{printf "\n"}'| sed 's/.$//')|gzip -c > individual_attr.txt.gz
#get genotypes for american admixed inds
python zpaste.py $(grep -v -w -f /media/benazzo/dati1/beans_modern/finestructure/selected_genes/trees/1kb_windows_refined/eu_individuals.txt ind_ids.txt | cut -f 1 -d " " |awk '{printf "%s ",$0}END{printf "\n"}'| sed 's/.$//')|cut -f  $(cat <(echo "1") <(echo "2") <(echo "3") <(seq 6 7 279) <(seq 7 7 280)|sort -n|awk '{printf "%s,",$0}END{printf "\n"}'| sed 's/.$//')|gzip -c > individual_attr_amer_not_pure.txt.gz
