
#### assemble update version GeLaTo
#Bioinformatic analysis
# Chiara Barbieri
# May 2022

## run FST
# in the Server

## now with new script from Epifania




## R elaboration

#perpopRED<-read.table("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perpopRED.txt",   sep = "\t", header=T, as.is=T)
perpopRED<-read.table("/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/perpopREDMaMi2021.txt",   sep = "\t", header=T, as.is=T)

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
#******************************************

## FST pairwise file in MASTERFILE

perpopMASTER<-read_xlsx("/Users/chiarabarbieri/switchdrive/GeLaTo2022/perpopMASTER_653pops2022.xlsx")

FstList<-read.csv("fst_GelatoHO_mergedSetFeb2022_.txt",as.is=T, header=T, sep = "\t") 
## pairwise FST wide format, generated by Epifania Arango with the script https://github.com/epifaniarango/Fst_forLargeDatasets


dim(FstList)
[1] 212878      3
length(unique(FstList$Pop1))
[1] 652  ## ok because i have the asymmetrix matrix, i do not compare each population to itself.

fstdouble<-FstList
fstdouble$Pop1<-FstList$Pop2
fstdouble$Pop2<-FstList$Pop1
# make it symmetric
FstList$case<-"single"
fstdouble$case<-"double" # to mark them for future analysis when i do not need repeated pairs
FstList<-rbind(FstList, fstdouble)

FstList$popslistemp<-paste0(FstList$Pop1,FstList$Pop2,sep="")
dim(FstList)
[1] 425756      5


### now add information for each population in the pair

FstListinfo<-FstList

minimuminfoPOP<-perpopMASTER[,c("PopName","glottocodeBase")] # assign glottocode to each pop in the pair (if available)

colnames(minimuminfoPOP)<-c("Pop2","glottocodeBase2")
FstListinfo<-merge(FstListinfo,minimuminfoPOP)
colnames(minimuminfoPOP)<-c("Pop1","glottocodeBase1")
FstListinfo<-merge(FstListinfo,minimuminfoPOP)

minimuminfoPOP<-perpopMASTER[,c("PopName","glottolog.node1")]  # the code of the highest node (language family)
colnames(minimuminfoPOP)<-c("Pop2","glottocodeFamily2")
FstListinfo<-merge(FstListinfo,minimuminfoPOP)
colnames(minimuminfoPOP)<-c("Pop1","glottocodeFamily1")
FstListinfo<-merge(FstListinfo,minimuminfoPOP)

minimuminfoPOP<-perpopMASTER[,c("PopName","glottolog.NAME")]  

colnames(minimuminfoPOP)<-c("Pop1","family1")   
FstListinfo<-merge(FstListinfo,minimuminfoPOP)   # automatically excludes the populations in GeLaTo that are not represented in the MarchMami popset (no glottocode)
colnames(minimuminfoPOP)<-c("Pop2","family2")   
FstListinfo<-merge(FstListinfo,minimuminfoPOP)   # automatically excludes the populations in GeLaTo that are not represented in the MarchMami popset (no glottocode)

# replace the FST negative with a zero

FstListinfo$FST[which(FstListinfo$FST<0)]<-0   
#************************


minimuminfoPOP<-perpopMASTER[,c("PopName","geographicRegion")]

colnames(minimuminfoPOP)<-c("Pop2","region2")
FstListinfo<-merge(FstListinfo,minimuminfoPOP)
colnames(minimuminfoPOP)<-c("Pop1","region1")
FstListinfo<-merge(FstListinfo,minimuminfoPOP)


minimuminfoPOP<-perpopMASTER[,c("PopName","lat", "lon")]

colnames(minimuminfoPOP)<-c("Pop1","lat1", "lon1")
FstListinfo<-  merge(FstListinfo,minimuminfoPOP)
colnames(minimuminfoPOP)<-c("Pop2","lat2", "lon2")
FstListinfo<-  merge(FstListinfo,minimuminfoPOP)



#******************************************************************** 
### geographic distance

library(fields)		#for geographic distances on coordinates
perpopREDgeo<-perpopMASTER[which(abs(as.numeric(perpopMASTER$lat))>0),]
lista<-cbind(as.numeric(perpopMASTER$lon),as.numeric(perpopMASTER$lat))
GEOdistances<-round(rdist.earth(lista, miles=FALSE))
rownames(GEOdistances)<-perpopMASTER$PopName
colnames(GEOdistances)<-perpopMASTER$PopName
library(reshape)
GEOmelt<-melt(GEOdistances)
colnames(GEOmelt)<- c("Pop1","Pop2","GEOdist")
GEOmelt$popslistemp<-paste0(GEOmelt$Pop1,GEOmelt$Pop2)
dim(GEOmelt)
[1] 426409      4

FstListinfo$GEOdist<-GEOmelt$GEOdist[ match(FstListinfo$popslistemp,GEOmelt$popslistemp)]

#******************************************************************** 
# FST LINERARIZED

FstListinfo$FstLinear<-FstListinfo$FST/(1-FstListinfo$FST)


#******************************************************************** 
# Language Family Pairs

withinfam<-c()
for (i in 1:nrow(FstListinfo)){
  if (FstListinfo$family1[i]==FstListinfo$family2[i]){
    withinfam[i]<-FstListinfo$family1[i]
  }
  else {
    withinfam[i]<-"DIVERSE"
  }
}

FstListinfo$FAMILY<-withinfam

FstListinfo$SameFamily<-"YES"
FstListinfo$SameFamily[which(FstListinfo$FAMILY=="DIVERSE")]<-"NO"
FstListinfo$REGION<-"DIVERSE"
FstListinfo$REGION[which(FstListinfo$region2==FstListinfo$region1)]<-FstListinfo$region1[which(FstListinfo$region2==FstListinfo$region1)]
###

## freeze pairwise info file MASTER
write.table(FstListinfo, "MASTER_PairwiseFstListinfo.txt", sep="\t", row.names = F, quote=F) 


FstListinfo<-read.table("MASTER_PairwiseFstListinfo.txt", header=T, sep="\t")

#***********************************
# FstListinfo<- replace with the version that includes only the matches for GeLaTo MaMi
#***********************************
#*
### Fst list pair, LONG format. WORKING SUBSET of GeLaTo with reduced list of 397 populations

# Median FST global and within macro region, for each pop 

MedianFST<-c()
MedianFSTregion<-c()
for (i in 1:nrow(perpopRED)){
  tempblock<-FstListinfo[c(which(FstListinfo$Pop1==perpopRED$PopName[i])),]
  MedianFST[i]<-median(tempblock$FST)
  regiontarget<-perpopRED$geographicRegion[i]
  regionpop<-perpopRED$PopName[which(perpopRED$geographicRegion==regiontarget)]
  MedianFSTregion[i]<-median(tempblock$FST[which(tempblock$Pop2%in%regionpop)])
}
perpopRED$medianFST<-MedianFST
perpopRED$medianFSTregion<-MedianFSTregion




### ANALYSIS SESSION FST 
FstListGlotto_infowithinREgion<-FstListinfo[which(FstListinfo$region1==FstListinfo$region2),]

## exclude drifted pops or the Fst averages will be higher than normal
DRIFTONI<-perpopRED[which(perpopRED$medianFSTregion>0.1&perpopRED$medianFST>0.1),]$PopName
#DRIFTONI<-perpopRED[which(perpopRED$averageFSTregion>0.1&perpopRED$averageFST>0.1),]$PopName
DRIFTONI
[1] "Baining_Marabu"      "Chukchi"             "Eskimo_Sireniki"     "Itelmen"             "Ju_hoan_North"      
[6] "Karitiana"           "Koryak"              "Lahu"                "Nanai"               "Nganasan"           
[11] "Nganasan_UstAvam"    "Nganasan_Volochanka" "Nivh"                "Onge"                "Rennell_and_Bellona"
[16] "San"                 "She"                 "Surui"              

FstListGlotto_infowithinREgionNoDrif<-FstListGlotto_infowithinREgion[-c(which(FstListGlotto_infowithinREgion$Pop2%in%DRIFTONI), which(FstListGlotto_infowithinREgion$Pop1%in%DRIFTONI)),]
FstListGlottoIBD_infoNoDrift<-FstListinfo[-c(which(FstListinfo$Pop2%in%DRIFTONI), which(FstListinfo$Pop1%in%DRIFTONI)),]


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
  perpopRED$MedianFSTAdjustedNeighbors[i]<-  median(tempneighbors$FST)
}
perpopRED$MedianFSTAdjustedNeighbors<-as.numeric(perpopRED$MedianFSTAdjustedNeighbors)
#perpopRED$MedianFSTAdjustedNeighbors[which(perpopRED$PopName%in%DRIFTONI)]<-0.1  # mark the Drifted pops



# ***************************************************
### ADD TMRCA from NE Calculated with IBD
# ***************************************************



# import the list of IBD fragments merged for all the chromosome.
# the data was processed with refinedIBD and the available tool for gap merging.

ibd<-read.table("all.refinedIBD_invariant_sites.Merged", as.is=T)
colnames(ibd)<-c("firstID","firstHapIndex","secondID","secondHapIndex", "chromosome", "start","end","LOD","length")

# add the population name associated to each individual ID, for both individuals involved in the sharing
ibd$source1<- infoID$PopName [match(ibd$firstID, infoID$Sample_ID)] 
ibd$source2<- infoID$PopName [match(ibd$secondID, infoID$Sample_ID)]

ibdSamePops<-ibd[which(ibd$source1==ibd$source2),] # chose only the fragments shared within the same population
ibdSamePops<-ibdSamePops[which(ibdSamePops$length>2),] # i filter for fragments larger than 2 cM to reduce noise
# IMPORTANT:
ibdSamePops<-ibdSamePops[-which(ibdSamePops$length>100),] # i exclude exceptionally long fragments that are outliers (0.2% of the total fragments)



# save separate files with the list of IBD blocks shared within populations
pops<-table(ibdSamePops$source1)
popNames<-unlist(labels(pops))
npops<-length(popNames)

for (i in 1:npops){
  
  target<-popNames[i]
  temp<-ibdSamePops[which(ibdSamePops$source1==target),][,1:8]
  system(paste("mkdir", target))
  write.table(temp, paste0(target,"/",target,".ibd"),col.names = F, row.names = F, sep = "\t", quote = F)
}
write.table(popNames[1:npops], "poplist.txt", row.names = F, col.names = F, quote = F)


#*****************************
# RUN IBDNe
#*****************************

### in the terminal, in the same folder where you downloaded the ibdne file (here called "ibdne.07May18.6a4.jar", but check for latest version)
# check also the appropriate map file

mypops=$( cat poplist.txt )
for target in  $mypops; do
cat ${target}/${target}.ibd | java -jar ibdne.07May18.6a4.jar map=allchromosomeGRCh37.map out=${target}/${target}IBDNE nthreads=8
done
####


# check to exclude the populations for which the IBDNe could not run (maybe not enough IBD segments overall) - VALID FOR LARGE LIST OF POPULATIONS
popNames2<-NA
for (i in 1:length(popNames)){
  target<-popNames[i]
  if(file.exists(paste0(target,"/",target,"IBDNE.ne"))){
    popNames2 [i]<- target
  }
}

popNames2<-popNames2[!is.na(popNames2)]
npops2<-length(popNames2)
tableNE<-matrix(NA, npops2, 5)
tableNE[,1]<-popNames2


## keep the first 50 generation for each population, skip the last two generations
# in the same loop, create a table for all the values of variation in size through generation for each population, and plot the BSP like plot for each population

npops2<-length(popNames2)
tableNE<-c()
for (i in 1:npops2){
  target<-popNames2[i]
  temp<-read.table(paste0(target,"/",target,"IBDNE.ne"), header=T, as.is = T)[5:51,] #from  3 to 50 generations ago
  temp$PopNames<-target
  tableNE<-rbind (tableNE, temp)
  pdf(paste0(target,"NePLOT.pdf"),useDingbats = FALSE)
  plot(temp[,1],temp[,2],type="l",log="y",xlim=c(0,50),ylim=c(100,1e10),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",lwd=2.5,lty=32,col="blue",main=target)
  abline(h=c(1e2,1e3,1e4,1e5,1e6,1e7),v=seq(20,80,20),col="lightgray",lty="dotted")
  lines(temp[,1],temp[,3],lwd=1)
  lines(temp[,1],temp[,4],lwd=1)
  dev.off()
  
}


colnames(tableNE)<-c("generation","mean","5perc","95perc", "PopName")
### harmonic mean of the Ne in the last 50 generations
library("psych")

for (i in 1:nrow(perpopRED2)){
  target<-perpopRED2$PopName[i]
  perpopRED2$harmonic2cm[i]<-harmonic.mean(tableNE[which(tableNE$PopName==target),2])
  perpopRED2$harmonic2cm_5perc[i]<-harmonic.mean(tableNE[which(tableNE$PopName==target),3])
  perpopRED2$harmonic2cm_95perc[i]<-harmonic.mean(tableNE[which(tableNE$PopName==target),4])
  
}

# difference between min and max Ne as an indication of variation of size in time
for (i in 1:nrow(perpopRED2)){
  target<-perpopRED2$PopName[i]
  perpopRED2$differenceNe2cM[i]<-max(tableNE[which(tableNE$PopName==target),2])-min(tableNE[which(tableNE$PopName==target),2])
}


# correlation coefficient slope between 20 and 3 gen ago as an indication of who expanded or collapsed in recent time
# use Spearmann for test

for (i in 1:nrow(perpopRED2)){
  target<-perpopRED2$PopName[i]
  testtemp<-cor.test(c(1:47),tableNE[which(tableNE$PopName==target),2])
  perpopRED2$correlationTestEstimate[i]<-testtemp$estimate
  perpopRED2$correlationTestPvalue[i]<-round(testtemp$p.value,digits = 4)
  
}

# now bin these results in increase or decrease in size

perpopRED2$variationSize<-"QuiteConstant"
perpopRED2$variationSize[which(perpopRED2$correlationTestEstimate>0.3)]<-"Decline"
perpopRED2$variationSize[which(perpopRED2$correlationTestEstimate<(-0.3))]<-"Expansion"
perpopRED2$variationSize[which(perpopRED2$correlationTestPvalue>0.05)]<-"NotSignificant"


write.table(perpopRED2, "perpopRED2_testNe.txt", row.names = F, quote = F, sep="\t")

# Manual screen populations for relatively constant Ne through time, not too broad confidence intervals and not too high Ne
# use the indicators and compare with the pop size variation Plots

#*************************
# ******************************************************
# TMRCA = linearizedFST * 2Ne * generationsyears

dim(FstListinfo)
[1] 157212     22

generationsyears=29
perpopREDNe<-perpopRED[which(perpopRED$USEforNe_calculation=="YES"),]
possiblepopswithNE<-perpopREDNe$PopName



FstListinfo$TMRCA_doubleNe<-rep(NA,nrow(FstListinfo))
FstListinfo$TMRCA_doubleNe_5<-rep(NA,nrow(FstListinfo))
FstListinfo$TMRCA_doubleNe_95<-rep(NA,nrow(FstListinfo))

FstListREDinfoNOne<-FstListinfo[-which(FstListinfo$Pop1%in%possiblepopswithNE&FstListinfo$Pop2%in%possiblepopswithNE),]
FstListREDinfoYESne<-FstListinfo[which(FstListinfo$Pop1%in%possiblepopswithNE&FstListinfo$Pop2%in%possiblepopswithNE),]

for (i in  1:nrow(FstListREDinfoYESne)){
  pop1<-FstListREDinfoYESne$Pop1[i]
  pop2<-FstListREDinfoYESne$Pop2[i]
  Ne1<-perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm
  Ne2<-perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm
  FstListREDinfoYESne$TMRCA_doubleNe[i]<-FstListREDinfoYESne$FstLinear[i] * (Ne1+Ne2) * generationsyears  # instead of (Ne1+Ne2)/2 * 2
  FstListREDinfoYESne$TMRCA_doubleNe_5[i]<-FstListREDinfoYESne$FstLinear[i] * (  perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm_5perc+perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm_5perc) * generationsyears
  FstListREDinfoYESne$TMRCA_doubleNe_95[i]<-FstListREDinfoYESne$FstLinear[i] * ( perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm_95perc+perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm_95perc) * generationsyears
}

FstListinfo<-rbind(FstListREDinfoYESne,FstListREDinfoNOne)
# write.table(FstListinfo, "/Users/chiarabarbieri/Library/Mobile Documents/com~apple~CloudDocs/GeLaTo/FstListREDinfo.txt", sep="\t", row.names = F, quote=F) 
# write.table(FstListREDinfo,"FstListREDinfo_MaMi2022.txt", row.names = F,  sep = "\t", quote = F)

#### another way to do TMRCA   https://genome.cshlp.org/content/23/9/1514

# formula for T generations ago: log(1-fst)/log(1-(1/(2*Ne)))

FstListinfo$TMRCA_Bhatia<-rep(NA,nrow(FstListinfo))


FstListREDinfoNOne<-FstListinfo[-which(FstListinfo$Pop1%in%possiblepopswithNE&FstListinfo$Pop2%in%possiblepopswithNE),]
FstListREDinfoYESne<-FstListinfo[which(FstListinfo$Pop1%in%possiblepopswithNE&FstListinfo$Pop2%in%possiblepopswithNE),]

for (i in  1:nrow(FstListREDinfoYESne)){
  pop1<-FstListREDinfoYESne$Pop1[i]
  pop2<-FstListREDinfoYESne$Pop2[i]
  Ne1<-perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm
  Ne2<-perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm
  Ne<-((Ne1+Ne2)/2)
  FstListREDinfoYESne$TMRCA_Bhatia[i]<-log(1-FstListREDinfoYESne$FstLinear[i])/log(1-(1/(2*Ne))) * generationsyears
  }

FstListinfo<-rbind(FstListREDinfoYESne,FstListREDinfoNOne)



#### HARMONIC MEAN TMRCA ####
## test harmonic mean of Ne following reviewer's suggestion
library("psych")

FstListinfo$TMRCA_harmonicNe<-rep(NA,nrow(FstListinfo))
FstListinfo$TMRCA_harmonicNe_5<-rep(NA,nrow(FstListinfo))
FstListinfo$TMRCA_harmonicNe_95<-rep(NA,nrow(FstListinfo))

FstListREDinfoNOne<-FstListinfo[-which(FstListinfo$Pop1%in%possiblepopswithNE&FstListinfo$Pop2%in%possiblepopswithNE),]
FstListREDinfoYESne<-FstListinfo[which(FstListinfo$Pop1%in%possiblepopswithNE&FstListinfo$Pop2%in%possiblepopswithNE),]

for (i in  1:nrow(FstListREDinfoYESne)){
  pop1<-FstListREDinfoYESne$Pop1[i]
  pop2<-FstListREDinfoYESne$Pop2[i]
  Ne1<-perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm
  Ne2<-perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm
  FstListREDinfoYESne$TMRCA_harmonicNe[i]<-FstListREDinfoYESne$FstLinear[i] * (harmonic.mean(c(Ne1,Ne2))*2) * generationsyears
   FstListREDinfoYESne$TMRCA_harmonicNe_5[i]<-FstListREDinfoYESne$FstLinear[i] * (harmonic.mean(c(perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm_5perc,perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm_5perc))*2) * generationsyears
   FstListREDinfoYESne$TMRCA_harmonicNe_95[i]<-FstListREDinfoYESne$FstLinear[i] * (harmonic.mean( c( perpopREDNe[which(perpopREDNe$PopName==pop1),]$harmonic2cm_95perc,perpopREDNe[which(perpopREDNe$PopName==pop2),]$harmonic2cm_95perc))*2) * generationsyears
}
FstListinfo<-rbind(FstListREDinfoYESne,FstListREDinfoNOne)


write.table(FstListinfo,"FstListREDinfo_MaMi2022_July.txt", row.names = F,  sep = "\t", quote = F)
write.table(perpopRED,"PerpopRED_MaMi2022_July.txt", row.names = F,  sep = "\t", quote = F)



