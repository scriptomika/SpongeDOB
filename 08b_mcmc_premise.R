start0<-Sys.time()
library(phytools)
library(ape)
library(phyloseq)
library(methods)
load("hosttr4mcmc.Rds")# ps.spongespecies hosttr

args = commandArgs(trailingOnly=TRUE)
#call script: Rscript script.R ../rooted_trees/RAxML_result.Bacteria_Acidobacteria_seqs
#Phylogenetic covariance matrices based on the bacterial and host phylogenies were then generated (Hadfield et al 2014 Am Nat; Pollock etal 2018):
# --Phylogenetic covariance matrices based on the bacterial and host phylogenies were generated using the function inverseA on each host tree. 
# --The Kronecker product of the resulting matrices was then computed for use as the ‘coevolutionary’ covariance matrix. 
# --The Kronecker product of each phylogenetic covariance matrix and an identity matrix was computed for use as 
#     -microbial identity x host phylogeny interaction effect
#     -microbial phylogeny × host identity interaction effect

#trfiles <- list.files(path="bact_trees/MLtrees", pattern="RAxML_result*", full.names=TRUE, recursive=FALSE)

#have bash wrapper loop through EACH bacterial family, to do:
#fam<-"Acidobacteria"
#trb<-read.tree("bact_trees/rootedtrees/RAxML_result.Bacteria_Acidobacteria_seqs")
print(args[1])
trb<-read.tree(args[1])
trb
tmp<- unlist(strsplit(args[1],"\\."))[6]
fam<-unlist(strsplit(tmp,"_seqs"))[1]
print(fam)
	og<-grep("Outgroup",trb$tip.label)
if (length(og)!=0){
	print(paste0("Outgroup: ",trb$tip.label[og]))
	rtrb<-root(trb,outgroup=og)
	##root then drop outgroups
	##plot(rtrb)
	rtrb<-drop.tip(rtrb,tip=og)
	##trb<-drop.tip(trb,tip=c("OutgroupA_sp51056","OutgroupA_sp50834","OutgroupT_sp40784","OutgroupV_sp51192"))
	}else{rtrb<-trb} 	#subtrees are rooted and outgroups pruned already

#only proceed with analysis for files not yet completed:
if(!file.exists(paste0(fam,'_mcmc_solutions.txt'))){

#drop enviro-only tips from microbial tree
pruned.bacttree<-keep.tip(rtrb,tip=which(rtrb$tip.label %in% taxa_names(ps.spongespecies)))

taxon_data <- merge_phyloseq(ps.spongespecies,pruned.bacttree)

sample_data(taxon_data)$sample_sum <- sample_sums(taxon_data)

#kill run if count table empty
if(all(sample_sums(taxon_data) < 10)==TRUE) {
	stop(paste0(fam," not present in sponge samples after adjustments."))
	}else{
n.pruned <- prune_samples(sample_sums(taxon_data) >= 10, taxon_data) # 
pruned.hosttree <- drop.tip(hosttr,hosttr$tip.label[!hosttr$tip.label %in% sample_data(n.pruned)$Species_consensus])

sample_data(n.pruned)$Species_consensus[!sample_data(n.pruned)$Species_consensus %in% pruned.hosttree$tip.label] <- NA
sample_data(n.pruned)$
Species_consensus <- droplevels(sample_data(n.pruned)$Species_consensus)

c.pruned <- prune_samples(!is.na(sample_data(n.pruned)$Species_consensus), n.pruned)
pruned <- filter_taxa(c.pruned, function(x) any(x>0),TRUE)

otutable <- as.matrix(as.data.frame(otu_table(pruned)))

library(MCMCglmm)
library(reshape2)
assocs <- melt(otutable,as.is=T)
assocs <- data.frame(count=assocs$value,otu=assocs$Var2,sample=assocs$Var1)
assocs$count<-round(assocs$count)
# todo add addititonal data that could be in mixed model..would be better if phyloseq object was samples not mergedinto species
assocs <- merge(assocs,sample_data(n.pruned)[,c('locale','Species_consensus','sample_sum')],by.x='sample',by.y=0,all=F)
assocs$sample_sum<-round(assocs$sample_sum)

inv.bact.full <- inverseA(pruned.bacttree,scale=F)

inv.bact <- inv.bact.full$Ainv
inv.host.full <- inverseA(hosttr)

inv.host <- inv.host.full$Ainv


host.otuA<-as(kronecker(inv.host, inv.bact), "dgCMatrix")                   # coevolutionary effect
host.otuAS<-as(kronecker(inv.host, Diagonal(nrow(inv.bact))), "dgCMatrix")  # host evolutionary effect
host.otuSA<-as(kronecker(Diagonal(nrow(inv.host)), inv.bact), "dgCMatrix")  # parasite evolutionary effect

rownames(host.otuA)<-apply(expand.grid(rownames(inv.bact), rownames(inv.host)), 1, function(x){paste(x[2],x[1], sep=".")})
rownames(host.otuAS)<-rownames(host.otuSA)<-rownames(host.otuA)


##assocs$otu																 # non-phylogenetic main effect for bacteria
##assocs$Species_consensus												 # non-phylogenetic main effect for hosts
assocs$otu.phy<-assocs$otu                                 				     # phylogenetic main effect for bacteria
assocs$Species_consensus<-assocs$Species_consensus                   # phylogenetic main effect for hosts
assocs$Host.otu<-paste(assocs$Species_consensus, assocs$otu, sep=".")      # non-phylogentic interaction effect
assocs$Host.otu[is.na(assocs$Species_consensus)] <- NA
assocs$Host.otu.cophy<-paste(assocs$Species_consensus, assocs$otu, sep=".")  # phylogentic coevolution effect
assocs$Host.otu.cophy[is.na(assocs$Species_consensus)] <- NA
assocs$Host.otu.hostphy<-paste(assocs$Species_consensus, assocs$otu, sep=".") # phylogentic host evolutionary effect (specifies whether abundance is determined by an interaction between non-phylogenetic otu and the phylogenetic position of the host)
assocs$Host.otu.hostphy[is.na(assocs$Species_consensus)] <- NA
assocs$Host.otu.otuphy<-paste(assocs$Species_consensus, assocs$otu, sep=".") # phylogentic parasite evolutionary effect (specifies whether abundance is determined by an interaction between non-phylogenetic host species and the phylogenetic position of the otu)
assocs$Host.otu.otuphy[is.na(assocs$Species_consensus)] <- NA
assocs$geo.otu <- paste(assocs$locale, assocs$otu, sep=".")


randfacts <- c('otu.phy','otu','geo.otu','Host.otu.hostphy','Host.otu.otuphy','Host.otu','Host.otu.cophy')


rand <- as.formula(paste0('~ ',paste(randfacts, collapse=' + ')))



priorC <- list(B=list(mu=c(0,1), V=diag(c(1e+8,1e-6))), R=list(V=1, nu=0))

## priors for the random evolutionary effects (from Hadfield):
phypri<-lapply(1:length(randfacts), function(x){list(V=1, nu=1, alpha.mu=0, alpha.V=1000)})

## combine priors:
priorC$G<-phypri
names(priorC$G)<-paste("G", 1:length(randfacts), sep="")

#run mcmc 
#Zaneveld ran nitt=125000, thin=5,burnin=25000
startmc<-Sys.time() 
mc <- MCMCglmm(count ~ log(sample_sum),
               random = rand,
               family="poisson",
               data=assocs,
               ginverse=list(otu.phy=inv.bact, Host.otu.hostphy=host.otuAS, Host.otu.otuphy=host.otuSA, Host.otu.cophy=host.otuA),
               prior=priorC,
               nitt=12500,
               thin=5,
               burnin=2500,
               pr=T)
endmc<-Sys.time()
print('saving results')

save(mc,file=paste0(fam,'_mcmc_res.RData'), compress=T, compression_level=9)

pdf(file=paste0(fam,'_VCV_%03d.pdf'), onefile=F)
plot(mc$VCV)
dev.off()

print('summarizing results')

sm <- summary(mc, random=T)

print('saving summary')

save(sm, file=paste0(fam,'_mcmc_summary.RData'), compress=T, compression_level=9)

write.table(sm$Gcovariance, paste0(fam,'_mcmc_Gcovariance.txt'), sep='\t', quote=F)
write.table(sm$solutions, paste0(fam,'_mcmc_solutions.txt'), sep='\t', quote=F)
endfinal<-Sys.time()
print('Total time')
endfinal - start0

print('MCMC time')
endmc - startmc
}
}
# Binary models were fit with MCMCglmm using a single fixed effect of the
# log of the sequencing depth, ‘global’ random effects of host phylogeny,
# host identity, microbial phylogeny, and microbial identity, all
# combinations of host-by-microbe phylogenetic and identity random
# interaction effects, and a geographic area-by- microbial identity random
# interaction effect. Altogether this approach is similar to the models
# described in reference62. Our models were fit with a chain length of
# 1,250,000, thinning interval of 50, and burn-in of 250,000. After the
# model was fit, convergence was assessed by verifying that the Effective
# Sample Sizes (ESS) of all covariance terms were greater than 200.
# Intraclass correlation coefficients (ICCs) were calculated for each
# iteration, with 95% credible intervals calculated with HPDinterval.
# Factors with ICC lower credible bounds greater than 0.01 were considered
# significant.