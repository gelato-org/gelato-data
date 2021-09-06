
#### assemble update version GeLaTo
#Bioinformatic analysis
# Chiara Barbieri
# 30 March 2021

## run FST
# in the Server

## now with new script from Epifania


# prune

~/programs/plink --bfile GelatoHO_mergedSetMarchMaMi --indep-pairwise 200 25 0.4 --out x.tmp
~/programs/plink --bfile GelatoHO_mergedSetMarchMaMi --extract x.tmp.prune.in --recode12 --out GelatoHO_mergedSetMarchMaMi.pruned

#### plink --het for inbreeding on a pruned continental set

~/programs/plink --file GelatoHO_mergedSetMarchMaMi.pruned --het --allow-no-sex --maf 0.05 --out GelatoHO_mergedSetAugustBED.pruned




## R elaboration

perpopRED<-read.table("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perpopRED.txt",   sep = "\t", header=T, as.is=T)

library(ggplot2)
library(reshape)


# assign language name and language family from glottocode

languages<-read.csv("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/languages.csv", header=T, as.is = T)
colnames(languages)[1]<-"glottocodeBase"

perpopRED$glottolog.node1<-languages$Family_ID [match(perpopRED$glottocodeBase, languages$glottocodeBase)]
perpopRED$glottolog.NAME<-languages$Name [match(perpopRED$glottolog.node1, languages$glottocodeBase)]


### create the fst matrix symmetric
### create the list of Fst pairs LONG format and add information


infoID<-read.csv("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/infoGelatoMarch2020matchfam.csv",header=T, as.is=T)

FstList<-read.table("/Users/chiarabarbieri/switchdrive/GeLaTo/fst.csv",as.is=T,  sep=" ")
dim(FstList)
# [1] 91806     1
fstdouble<-FstList
fstdouble$V1<-FstList$V2
fstdouble$V2<-FstList$V1
# make it symmetric
FstList$case<-"single"
fstdouble$case<-"double" # to mark them for future analysis when i do not need repeated pairs
FstList<-rbind(FstList, fstdouble)
colnames(FstList)<-c("Pop1", "Pop2", "Fst", "case")
FstList$popslistemp<-paste0(FstList$Pop1,FstList$Pop2,sep="")

dim(FstList)
[1] 183612      5

fst2<-as.matrix(cast(FstList, Pop1~Pop2, value="Fst" ))
diag(fst2)<-0
write.table(fst2,"matrixfstALL.txt", sep="\t")


 dim(fst2)
#[1] 429 429
# with all the set of Gelato, including pops for which i do not have the glottocode

fstShapeMatrixMaMi<-fst2[which(rownames(fst2)%in%perpopRED$PopName),which(colnames(fst2)%in%perpopRED$PopName)] # only the set of  404 pops

isSymmetric.matrix(fstShapeMatrixMaMi) # check

write.table(fstShapeMatrixMaMi,"fstShapeMatrixMaMi.txt", sep = "\t")


#******************************************************************** 
#### FREEZE MAMI DATASET 4284 INDIVIDUALS
# HETEROZYGOSITY

het<-read.table("GelatoHO_mergedSetMarchMaMiAUTOSOMAL.pruned.het",header=T)

het$proportionHeteroz<-(het$N.NM.-het$O.HOM.)/ het$N.NM.
homozyg<-c()
variancehomozyg<-c()
proportionHeterozy<-c()

for (i in 1:nrow(perpopRED)){
  homozyg[i]<-mean(het[which(het$FID==perpopRED$PopName[i]),][,6])
  variancehomozyg[i]<-var(het[which(het$FID==perpopRED$PopName[i]),][,6])
  proportionHeterozy[i]<-mean(het[which(het$FID==perpopRED$PopName[i]),][,7])
  
}

perpopRED$homozyg<-homozyg
perpopRED$variancehomozyg<-variancehomozyg
perpopRED$proportionHeterozy<-proportionHeterozy


#********************************************************************
### create Melt file for pairs 
#******************************************************************** 


#******************************************************
# elaborate results FST

#FstList<-read.csv("/Volumes/MANNAIA/GeLaTo/HumOrigins/FST/Fst_mat.csv",as.is=T, header=T)

minimuminfoPOP<-perpopRED[,c("PopName","glottolog.NAME")]

colnames(minimuminfoPOP)<-c("Pop1","family1")   
FstList1<-merge(FstList,minimuminfoPOP)   # automatically excludes the populations in GeLaTo that are not represented in the MarchMami popset (no glottocode)
colnames(minimuminfoPOP)<-c("Pop2","family2")   
FstListREDinfo<-merge(FstList1,minimuminfoPOP)   # automatically excludes the populations in GeLaTo that are not represented in the MarchMami popset (no glottocode)


# FstListREDinfo<- FstList1[which(FstList1$Fst!=0),] # check for FST values ==0
FstListREDinfo$Fst[which(FstListREDinfo$Fst<0)]<-0   # replace the FST negative with a zero
dim(FstListREDinfo)
[1] 162812      7


minimuminfoPOP<-perpopRED[,c("PopName","glottocodeBase")]

colnames(minimuminfoPOP)<-c("Pop2","glottocodeBase2")
FstListREDinfo<-merge(FstListREDinfo,minimuminfoPOP)
colnames(minimuminfoPOP)<-c("Pop1","glottocodeBase1")
FstListREDinfo<-merge(FstListREDinfo,minimuminfoPOP)

minimuminfoPOP<-perpopRED[,c("PopName","glottolog.node1")]  # the code of the highest node (language family)
colnames(minimuminfoPOP)<-c("Pop2","glottocodeFamily2")
FstListREDinfo<-merge(FstListREDinfo,minimuminfoPOP)
colnames(minimuminfoPOP)<-c("Pop1","glottocodeFamily1")
FstListREDinfo<-merge(FstListREDinfo,minimuminfoPOP)


minimuminfoPOP<-perpopRED[,c("PopName","geographicRegion")]

colnames(minimuminfoPOP)<-c("Pop2","region2")
FstListREDinfo<-merge(FstListREDinfo,minimuminfoPOP)
colnames(minimuminfoPOP)<-c("Pop1","region1")
FstListREDinfo<-merge(FstListREDinfo,minimuminfoPOP)


minimuminfoPOP<-perpopRED[,c("PopName","lat", "lon")]

colnames(minimuminfoPOP)<-c("Pop1","lat1", "lon1")
FstListREDinfo<-  merge(FstListREDinfo,minimuminfoPOP)
colnames(minimuminfoPOP)<-c("Pop2","lat2", "lon2")
FstListREDinfo<-  merge(FstListREDinfo,minimuminfoPOP)



#******************************************************************** 
### geographic distance

library(fields)		#for geographic distances on coordinates
perpopREDgeo<-perpopRED[which(abs(as.numeric(perpopRED$lat))>0),]
lista<-cbind(as.numeric(perpopRED$lon),as.numeric(perpopRED$lat))
GEOdistances<-round(rdist.earth(lista, miles=FALSE))
rownames(GEOdistances)<-perpopRED$PopName
colnames(GEOdistances)<-perpopRED$PopName
library(reshape)
GEOmelt<-melt(GEOdistances)
colnames(GEOmelt)<- c("Pop1","Pop2","GEOdist")
GEOmelt$popslistemp<-paste0(GEOmelt$Pop1,GEOmelt$Pop2)
dim(GEOmelt)
[1] 163216      4

FstListREDinfo<-merge(FstListREDinfo,GEOmelt[,3:4])

#******************************************************************** 
# FST LINERARIZED

FstListREDinfo$FstLinear<-FstListREDinfo$Fst/(1-FstListREDinfo$Fst)


#******************************************************************** 
# Language Family Pairs

withinfam<-c()
for (i in 1:nrow(FstListREDinfo)){
  if (FstListREDinfo$family1[i]==FstListREDinfo$family2[i]){
    withinfam[i]<-FstListREDinfo$family1[i]
  }
  else {
    withinfam[i]<-"DIVERSE"
  }
}

FstListREDinfo$FAMILY<-withinfam

FstListREDinfo$SameFamily<-"YES"
FstListREDinfo$SameFamily[which(FstListREDinfo$FAMILY=="DIVERSE")]<-"NO"
FstListREDinfo$REGION<-"DIVERSE"
FstListREDinfo$REGION[which(FstListREDinfo$region2==FstListREDinfo$region1)]<-FstListREDinfo$region1[which(FstListREDinfo$region2==FstListREDinfo$region1)]
###


#++++++++++++++++++++++++++++++++++++++++++++
# Median FST global and within region, per each pop 

MedianFST<-c()
MedianFSTregion<-c()
for (i in 1:nrow(perpopRED)){
  tempblock<-FstListREDinfo[c(which(FstListREDinfo$Pop1==perpopRED$PopName[i]),which(FstListREDinfo$Pop2==perpopRED$PopName[i])),]
  MedianFST[i]<-median(tempblock$Fst)
  regiontarget<-perpopRED$geographicRegion[i]
  regionpop<-perpopRED$PopName[which(perpopRED$geographicRegion==regiontarget)]
  MedianFSTregion[i]<-median(tempblock$Fst[which(tempblock$Pop1%in%regionpop&tempblock$Pop2%in%regionpop)])
}
perpopRED$medianFST<-MedianFST
perpopRED$medianFSTregion<-MedianFSTregion




### ANALYSIS SESSION FST 
FstListGlotto_infowithinREgion<-FstListREDinfo[which(FstListREDinfo$region1==FstListREDinfo$region2),]

## exclude drifted pops or the Fst averages will be higher than normal
DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName
#DRIFTONI<-perpopRED[which(perpopRED$averageFSTregion>0.1&perpopRED$averageFST>0.1),]$PopName
DRIFTONI
[1] "Rennell_and_Bellona" "Ju_hoan_North"       "San"                 "Algonquin"           "Baining_Marabu"     
[6] "Chukchi"             "Eskimo_Sireniki"     "Itelmen"             "Koryak"              "Nivh"               
[11] "Onge"                "She"                 "Lahu"                "Nanai"               "Karitiana"          
[16] "Surui"               "Nganasan"            "Nganasan_UstAvam"    "Nganasan_Volochanka"

FstListGlotto_infowithinREgionNoDrif<-FstListGlotto_infowithinREgion[-c(which(FstListGlotto_infowithinREgion$Pop2%in%DRIFTONI), which(FstListGlotto_infowithinREgion$Pop1%in%DRIFTONI)),]
FstListGlottoIBD_infoNoDrift<-FstListREDinfo[-c(which(FstListREDinfo$Pop2%in%DRIFTONI), which(FstListREDinfo$Pop1%in%DRIFTONI)),]


#average FST within 1000 km exclude driftoni
# proportion FST adjusted for the median of the neighbors
#EXCLUDE DRIFTONI

radius<-1000

perpopRED$MedianFSTAdjustedNeighbors<-NA
perpopRED$numberofNeighbors<-NA

for (i in 1:nrow(perpopRED)){
  target<-perpopRED$PopName[i]
  tempblock<-FstListGlottoIBD_infoNoDrift[which(FstListGlottoIBD_infoNoDrift$Pop1==target),]
  tempneighbors<-tempblock[which(tempblock$GEOdist<radius),]
  neighborsnames<-unique(c(tempneighbors$Pop1,tempneighbors$Pop2))
  neighborsnames<-neighborsnames[-which(neighborsnames==target)]
  perpopRED$numberofNeighbors[i]<-length(neighborsnames)
  perpopRED$MedianFSTAdjustedNeighbors[i]<-  median(tempneighbors$Fst)
}
perpopRED$MedianFSTAdjustedNeighbors<-as.numeric(perpopRED$MedianFSTAdjustedNeighbors)
perpopRED$MedianFSTAdjustedNeighbors[which(perpopRED$PopName%in%DRIFTONI)]<-0.1  # mark the Drifted pops

# proportion FST, proportion heterozygosity, adjusted for the median of the neighbors
#EXCLUDE DRIFTONI

radius<-1000

perpopRED$proportionHeterozyAdjustedNeighbors<-c()

for (i in 1:nrow(perpopRED)){
  target<-perpopRED$PopName[i]
  tempblock<-FstListREDinfo[c(which(FstListREDinfo$Pop1==target),which(FstListREDinfo$Pop2==target)),]
  tempneighbors<-tempblock[which(tempblock$GEOdist<radius),]
  neighborsnames<-unique(c(tempneighbors$Pop1,tempneighbors$Pop2))
  neighborsnames<-neighborsnames[-which(neighborsnames==target)]
  perpopneighbors<-perpopRED[which(perpopRED$PopName%in%neighborsnames),]
  perpopRED$proportionHeterozyAdjustedNeighbors[i]<-  perpopRED$proportionHeterozy[i]/median(perpopneighbors$proportionHeterozy)
}


perpopRED$proportionFstAdjustedNeighbors<-c()

for (i in 1:nrow(perpopRED)){
  target<-perpopRED$PopName[i]
  tempblock<-FstListREDinfo[c(which(FstListREDinfo$Pop1==target),which(FstListREDinfo$Pop2==target)),]
  tempneighbors<-tempblock[which(tempblock$GEOdist<radius),]
  neighborsnames<-unique(c(tempneighbors$Pop1))
  neighborsnames<-neighborsnames[-which(neighborsnames==target)]
  neighborsnames<-neighborsnames[!neighborsnames%in%DRIFTONI]  # exclude Driftoni from the median value
  
  perpopneighbors<-perpopRED[which(perpopRED$PopName%in%neighborsnames),]
  perpopRED$proportionFstAdjustedNeighbors[i]<-  perpopRED$averageFST[i]/median(perpopneighbors$averageFST)
}

perpopRED$proportionFstAdjustedNeighbors<-as.numeric(perpopRED$proportionFstAdjustedNeighbors)

write.table(perpopRED, "/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perpopREDMaMi.txt", row.names = F,  sep = "\t", quote = F)


# ***************************************************
### ADD TMRCA from NE Calculated with IBD
# ***************************************************

# TMRCA = linearizedFST * 2Ne * generationsyears

dim(FstListREDinfo)
[1] 162812     20

generationsyears=29
perpopREDNe<-perpopRED[which(perpopRED$USEforNe_calculation=="YES"),]
possiblepopswithNE<-perpopREDNe$PopName


FstListREDinfo$TMRCA_doubleNe<-rep(NA,nrow(FstListREDinfo))
FstListREDinfo$TMRCA_doubleNe_5<-rep(NA,nrow(FstListREDinfo))
FstListREDinfo$TMRCA_doubleNe_95<-rep(NA,nrow(FstListREDinfo))

FstListREDinfoNOne<-FstListREDinfo[-which(FstListREDinfo$Pop1%in%possiblepopswithNE&FstListREDinfo$Pop2%in%possiblepopswithNE),]
FstListREDinfoYESne<-FstListREDinfo[which(FstListREDinfo$Pop1%in%possiblepopswithNE&FstListREDinfo$Pop2%in%possiblepopswithNE),]

for (i in  1:nrow(FstListREDinfoYESne)){
  pop1<-FstListREDinfoYESne$Pop1[i]
  pop2<-FstListREDinfoYESne$Pop2[i]
  Ne1<-perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm
  Ne2<-perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm
  FstListREDinfoYESne$TMRCA_doubleNe[i]<-FstListREDinfoYESne$FstLinear[i] * (Ne1+Ne2) * generationsyears
  FstListREDinfoYESne$TMRCA_doubleNe_5[i]<-FstListREDinfoYESne$FstLinear[i] * (  perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm_5perc+perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm_5perc) * generationsyears
  FstListREDinfoYESne$TMRCA_doubleNe_95[i]<-FstListREDinfoYESne$FstLinear[i] * (  perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm_95perc+perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm_95perc) * generationsyears
  }

FstListREDinfo<-rbind(FstListREDinfoYESne,FstListREDinfoNOne)
# write.table(FstListREDinfo, "/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/FstListREDinfo.txt", sep="\t", row.names = F, quote=F) 

