library(phytools)
library(methods)
library(phyloseq)
mainDir<-"~/Box Sync/MPL/Lesser DOB sponges/coevolution_analyses/pipeline4/"
setwd(mainDir)
outDir<-"08_pglmm/"
ifelse(!dir.exists(file.path(mainDir, outDir)), dir.create(file.path(mainDir, outDir)))
start0<-Sys.time()
load("06_swadj_results/ps.adj.hostspecies.ra.RDS") #loads 'ps.adj.hostspecies.ra'
ps.hostspecies.ra<-ps.adj.hostspecies.ra

load("02_tree_results/spongespeciesTimeTree.OGless.Rds") #loads cladecols cladelist colsF colsO timetr treemeta.u
ctr<-timetr
ctr$tip.label[!(ctr$tip.label %in% sample_names(ps.hostspecies.ra))] 
sample_names(ps.hostspecies.ra)[!(sample_names(ps.hostspecies.ra) %in% ctr$tip.label )] 

saved_nodelabels<-ctr$node.label
ctr$node.label<-NULL
tr2<-force.ultrametric(ctr)
tr2$edge.length[tr2$edge.length<=0]<-min(tr2$edge.length[tr2$edge.length>0])

#save spongespecies ps object and matching sponge tree
keep<-sample_names(ps.hostspecies.ra)[(sample_names(ps.hostspecies.ra) %in% tr2$tip.label )]
pstmp<-prune_samples(keep,ps.hostspecies.ra) #excludes enviro and tree-missing sponge samples
ps.spongespecies<-prune_taxa(taxa_sums(pstmp)>0,pstmp) #discard rare taxa
hosttr<-keep.tip(tr2,tip=which(tr2$tip.label %in% sample_names(ps.spongespecies)))
save(ps.spongespecies,hosttr,file=paste0(outDir,'hosttr4mcmc.Rds'))

### H2: ASVs diversify or diverge more within a host than between hosts.

#ways to test if  microbes diversifying more within host (family, genus, species) than between? 
# divergence and diversification of specialists...
# divergence == average of PDs among ASVs unique to a host (null: all PDs between hosts)
# diversification == number of nodes with ASR 0.95MLE for a single host (null: num/prop of all nodes recieving 0.95 for any host)
## phylo GLMM on each family

