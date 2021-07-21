setwd("~/Box Sync/MPL/Lesser DOB sponges/coevolution_analyses/pipeline4/02_input_spongetree/")
library(dplyr)
tab18<-read.table("18SAccessionReport.tsv",header = F);colnames(tab18)<-c("r18s","Sequence_ID","rel")
tab18<-tab18 %>% select(Sequence_ID,r18s)
meta18<-read.csv("source_modifiers_18s.csv",header = T) %>% select(Sequence_ID,Organism)
tabc<-read.table("COIAccessionReport.tsv",header=F);colnames(tabc)<-c("cox1","Sequence_ID","rel")
tabc<-tabc%>% select(Sequence_ID,cox1)
metaC<-read.csv("source.coi.out2.csv",header=T)  %>% select(Sequence_ID,Organism)

metaC<-full_join(tabc,metaC) %>% select(Sequence_ID, Organism, cox1)
meta18<-full_join(tab18,meta18) %>% select(Sequence_ID, Organism, r18s)
newgb<-full_join(metaC,meta18)
nsamp<-nrow(newgb)
levels(newgb$Organism)
nspp<-length(unique(newgb$Organism))

unique(newgb$Organism[grep("\\ sp",newgb$Organism)])

`%notin%` <- Negate(`%in%`)
#read metadat
meta16<-read.csv("../01_inputs_DADA/metaDOB.csv",header = T,stringsAsFactors = F) %>% 
filter(Species_consensus %notin% c("seawater","sediment"))   %>% 
#filter(project=="DOB")   %>% #remove nonDOB samples
filter(depth=="50fsw")   %>% 
filter(discard!='y')   %>% #remove previously examined, outliers etc #53 DOB 50fsw water samples
filter(DNAid %notin% c("KY265","KY284","KY285","KY286","KY287"))   %>%
select(DNAid,Species_consensus)
nrow(meta16)

metatr<-read.csv("treemeta.csv",header=T,stringsAsFactors = F) %>% select(tip,Species_consensus) %>% filter(Species_consensus %notin% c("seawater","sediment"))
metas<-full_join(meta16,metatr,by=c("DNAid"="tip")); metas[is.na(metas)]<-"missing"; colnames(metas)<-c("DNAid","meta16sp","tree_sp")
dim(metas)
metas[metas$meta16sp!=metas$tree_sp,]

#compare species names used in metadatas and submitted to in GB
metas2<-left_join(metas,newgb, by=c("DNAid"="Sequence_ID"))
c2<-metas2 %>%  filter(Organism!=tree_sp)
c2[300:400,]
metas2 %>% filter(meta16sp!="missing") %>%  filter(tree_sp!="missing") %>% filter(meta16sp!=tree_sp)
