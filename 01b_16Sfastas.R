mainDir<-"~/Box Sync/MPL/Lesser DOB sponges/coevolution_analyses/pipeline4/"
setwd(mainDir)
outDir<-"01_results/"
ifelse(!dir.exists(file.path(mainDir, outDir)), dir.create(file.path(mainDir, outDir)), FALSE)
system("ls")
## read in same filtered ps object created for 16S analyses
load(paste0(outDir,"ps4.rds"))
load(paste0(outDir,"ps4.seqs.rds"))
ps4
ASVseqs.ps4

tax_table(ps4)

classes<-character();num<-numeric(); 
for (i in 1:length(unique(tax_table(ps4)[,"Class"])))
  
{
  cla<-as.character(unique(tax_table(ps4)[,"Class"])[i])
  if(is.na(cla)==F){
    clalab<-as.character(tax_table(ps4)[which(tax_table(ps4)[,"Class"]==cla),"KPC"])[1]
    classes<-c(classes,clalab)
    keep<-taxa_names(ps4)[which(tax_table(ps4)[,"Class"]==cla)]
    psi<-prune_taxa(keep,ps4)
    num<- c(num,ntaxa(psi))
  }
}
nas<-taxa_names(tax_table(ps4)[which(is.na(tax_table(ps4)[,"Class"])),])
totalcounts.orderNA<-sum(taxa_sums(prune_taxa(nas,ps4)))
totalcounts.claknown<-sum(taxa_sums(prune_taxa(!(taxa_names(ps4) %in% nas),ps4)))
totalcounts.orderNA/(totalcounts.orderNA+totalcounts.claknown)
totalcounts.orderNA/sum(taxa_sums(ps4))
# "Order": about 25% of bacterial abundance not assigned to Order level
# "Class": about 4% of bact abundance not assigned to Class level

keep<-taxa_names(ps4)[which(tax_table(ps4)[,"Phylum"]=="Poribacteria")]
psPori<-prune_taxa(keep,ps4)
psPori<-prune_samples(names(which(sample_sums(psPori)!=0)),psPori)
psPori.ra <- transform_sample_counts(psPori, function(OTU) OTU/sum(OTU)) #convert to relative 
ps.Pori.m<-merge_samples(psPori.ra,"Species_consensus",fun=mean)
# repair factors after merge:
sample_data(ps.Pori.m)$Species_consensus <-sample_names(ps.Pori.m)
#sample_data(ps.Pori.m)
#taxa_sums(ps.Pori.m)
#sample_sums(ps.Pori.m)
#otu_table(ps.Pori.m)
plot_heatmap(ps.Pori.m) #order sample.order by sponge phylo and taxa.order by bact tree


Ctab<-data.frame(classes,num)
nrow(Ctab)
nrow(Ctab[which(Ctab$num>10),])
hist(Ctab$num[which(Ctab$num>10)],breaks=50)
head(Ctab)

Ctab2<-Ctab[which(Ctab$num>10),]
Ctab2[order(Ctab2$num),]
dim(Ctab2)
#retrieve seqs for each Order with sufficient members
#ASVseqs.ps4

library(ShortRead)
for (i in 1:nrow(Ctab2)){
  myC<- as.character(Ctab2$classes[i])
  print(myC)
  seqsIDsKeep<-taxa_names(ps4)[which(tax_table(ps4)[,"KPC"] == myC)]
  print(length(seqsIDsKeep))
  #create fasta for each...
  seqsKeep<-ASVseqs.ps4[seqsIDsKeep]
  db_out <- data.frame(ids=names(seqsKeep),seqs=seqsKeep)
  fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
  # write spongy fasta for filtering
  writeFasta(fasta, file = paste0(myC,"_seqs.fasta"))
  
}
system("mkdir bact_classes_fasta")
system("mv *fasta bact_classes_fasta")
#todo: go back and capture the Class-less bacteria Phylum, eg. Poribacteria

for (i in as.character(unique(tax_table(ps4)[nas,"KP"])))
{ 
  ps.i<-prune_taxa(taxa_names(ps4)[tax_table(ps4)[,"KP"]==i ],ps4)
   print(i)
  seqsIDsKeep<-taxa_names(ps.i)
  print(length(seqsIDsKeep))
  #create fasta for each...
  if (length(seqsIDsKeep)>10) {
  seqsKeep<-ASVseqs.ps4[seqsIDsKeep]
  db_out <- data.frame(ids=names(seqsKeep),seqs=seqsKeep)
  fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
  # write spongy fasta for filtering
  writeFasta(fasta, file = paste0(i,"_seqs.fasta"))
  }

}
system("mkdir bact_smallphyla_fasta")
system("mv *.fasta bact_smallphyla_fasta")


phys<-character();num<-numeric(); 
for (i in 1:length(unique(tax_table(ps4)[,"KP"])))
  
{
  phyi<-as.character(unique(tax_table(ps4)[,"KP"])[i])
  if(is.na(phyi)==F){
    phyilab<-as.character(tax_table(ps4)[which(tax_table(ps4)[,"KP"]==phyi),"KP"])[1]
    phys<-c(phys,phyilab)
    keep<-taxa_names(ps4)[which(tax_table(ps4)[,"KP"]==phyi)]
    psi<-prune_taxa(keep,ps4)
    num<- c(num,ntaxa(psi))
  }
}
nas<-taxa_names(tax_table(ps4)[which(is.na(tax_table(ps4)[,"Phylum"])),])
totalcounts.orderNA<-sum(taxa_sums(prune_taxa(nas,ps4)))
totalcounts.orderknown<-sum(taxa_sums(prune_taxa(!(taxa_names(ps4) %in% nas),ps4)))
totalcounts.orderNA/(totalcounts.orderNA+totalcounts.orderknown)
totalcounts.orderNA/sum(taxa_sums(ps4))
# "Order": about 25% of bacterial abundance not assigned to Order level
# "Class": about 4% of bact abundance not assigned to Class level
# "Phylum" : about 1% of bac abundance not assigned to a Phylum
Ptab<-data.frame(phys,num)
nrow(Ptab)
nrow(Ptab[which(Ptab$num>10),])
hist(Ptab$num[which(Ptab$num>10)],breaks=50)
head(Ptab)

Ptab2<-Ptab[which(Ptab$num>10),]
dim(Ptab2)
Ptab2[order(Ptab2$num),]

#retrieve seqs for each Order with sufficient members
#ASVseqs.ps4

for (i in 1:nrow(Ptab2)){
  myP<- as.character(Ptab2$phys[i])
  print(myP)
  seqsIDsKeep<-taxa_names(ps4)[which(tax_table(ps4)[,"KP"] == myP)]
  print(length(seqsIDsKeep))
  #create fasta for each...
  seqsKeep<-ASVseqs.ps4[seqsIDsKeep]
  db_out <- data.frame(ids=names(seqsKeep),seqs=seqsKeep)
  fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
  # write spongy fasta for filtering
  writeFasta(fasta, file = paste0(myP,"_seqs.fasta"))
  
}
system("mkdir bact_phyla_fasta")
system("mv *.fasta bact_phyla_fasta")

#run 16S trees with RAxML and calc PD distances for all pairs
# run class trees for: Proteobacteria, Planctomycetes, Chloroflexi (≥ 2000 members in phylum)
# discard taxonomic-NA fastas


  