library(tidyverse)
library(ggpubr)
library(purrr)
library(viridis)
library(reshape2)
library(wesanderson)
library(scales)
library(data.table)
library(ggExtra)
library(cowplot)
library(ggridges)
library(qqman)
library(ggrepel)
library(gt)
library(BiocManager)
BiocManager::install("GenomicRanges")
library(GenomicRanges)
BiocManager::install("genomation")
library(genomation)
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
library(VennDiagram)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
BiocManager::install("rrvgo")
library(rrvgo)
library(purrr)
library(jsonlite)
library(clusterProfiler)



###########################################################################
################### PCA ###########################################

library("tidyverse")
library('viridis')
library('ggpubr')


###### PCA of the goodeids
pca <- read_table("Goodeids.eigenvec", col_names = FALSE)
eigenval <- scan('Goodeids.eigenval')


#Get rid of duplicated column 
pca <- pca[,-1]
names(pca)[1] <- "species"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pca$species<-gsub('AS', 'Ameca splendens', pca$species)
pca$species<-gsub('AT', 'Ataeniobius toweri', pca$species)
pca$species<-gsub('CB', 'Crenicthys baileyi', pca$species)
pca$species<-gsub('GA', 'Goodea atripinnis', pca$species)
pca$species<-gsub('XR', 'Xenotaenia resolanae', pca$species)
pca$species<-gsub('XC', 'Xenoophorus captivus', pca$species)
pca$species<-gsub('IF', 'Ilyodon furcidens', pca$species)
pca$species<-gsub('CL', 'Characodon lateralis', pca$species)
pca$species<-gsub('GM', 'Girardinichthys multiradiatus', pca$species)

pve <- data.frame(PC = 1:9, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a <- a + ylab("Percentage variance explained") + theme_light()

b <- ggplot(pca, aes(PC1, PC2, col = species, label=species)) + geom_point(size = 3, alpha=0.6,
                                                                           position=position_jitter(width=0.06,height=0.06))
b <- b + scale_colour_brewer(palette = "Set1")
b <- b + coord_equal() + theme_light()
b <-b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
b
### PC2 vs PC3 - PC1 seems to explain technical variation.
PC23 <- ggplot(pca, aes(PC2, PC3, col = species, label=species)) + geom_point(size = 3, alpha=0.6,
                                                                           position=position_jitter(width=0.06,height=0.06))
PC23 <- PC23 + scale_colour_brewer(palette = "Set1")
PC23 <- PC23 + coord_equal() + theme_light()
PC23 <- PC23 + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))
PC23

### PC1 vs PC3
c <- ggplot(pca, aes(PC1, PC3, col = species, label=species)) + geom_point(size = 3, alpha=0.6,
                                                                           position=position_jitter(width=0.06,height=0.06))
c <- c + scale_colour_brewer(palette = "Set1")
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

###PC1 vs PC4
d <- ggplot(pca, aes(PC1, PC4, col = species, label=species)) + geom_point(size = 3, alpha=0.6,
                                                                           position=position_jitter(width=0.06,height=0.06))
d <- d + scale_colour_brewer(palette = "Set1")
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)"))

tiff("Biogeography_PCA_subsetted_data.png", units="in", width=10, height=7, res=300)
ggarrange(b, PC23, nrow = 1, ncol =2)
dev.off()





##################### seperate chromosomes
eigenvec_files <- list.files(path="chroms", pattern="*.eigenvec", full.names=TRUE, recursive=FALSE)
eigenval_files <- list.files(path="chroms", pattern="*.eigenval", full.names=TRUE, recursive=FALSE)
PCA_elements<-list(eigenvec_files, eigenval_files)



pca <- read_table("chroms/Goodeid.NC_024353.1_RagTag.PCA.eigenvec", col_names = FALSE)
eigenval <- scan("chroms/Goodeid.NC_024353.1_RagTag.PCA.eigenval")


#Get rid of duplicated column 
pca <- pca[,-1]
names(pca)[1] <- "species"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pca$species<-gsub('AS', 'Ameca splendens', pca$species)
pca$species<-gsub('AT', 'Ataeniobius toweri', pca$species)
pca$species<-gsub('CB', 'Crenicthys baileyi', pca$species)
pca$species<-gsub('GA', 'Goodea atripinnis', pca$species)
pca$species<-gsub('XR', 'Xenotaenia resolanae', pca$species)
pca$species<-gsub('XC', 'Xenoophorus captivus', pca$species)
pca$species<-gsub('IF', 'Ilyodon furcidens', pca$species)
pca$species<-gsub('CL', 'Characodon lateralis', pca$species)
pca$species<-gsub('GM', 'Girardinichthys multiradiatus', pca$species)

pve <- data.frame(PC = 1:9, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a <- a + ylab("Percentage variance explained") + theme_light()

b <- ggplot(pca, aes(PC1, PC2, col = species, label=species)) + geom_point(size = 3, alpha=0.6,
                                                                           position=position_jitter(width=0.06,height=0.06))
b <- b + scale_colour_brewer(palette = "Set1")
b <- b + coord_equal() + theme_light()
b <-b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))


### PC1 vs PC3
c <- ggplot(pca, aes(PC1, PC3, col = species, label=species)) + geom_point(size = 3, alpha=0.6,
                                                                           position=position_jitter(width=0.06,height=0.06))
c <- c + scale_colour_brewer(palette = "Set1")
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

###PC1 vs PC4
d <- ggplot(pca, aes(PC1, PC4, col = species, label=species)) + geom_point(size = 3, alpha=0.6,
                                                                           position=position_jitter(width=0.06,height=0.06))
d <- d + scale_colour_brewer(palette = "Set1")
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)"))

png("PCA_LG23.png", units="in", width=10, height=7, res=300)
ggarrange(a, b, c, d, nrow = 2, ncol =2 )
dev.off()

###########################################################################
#############       Estimating divergence times   ###################
library(bppr)
library(coda)
library(ape)
if (!any(rownames(installed.packages()) == "MCMCtreeR")) install.packages("MCMCtreeR")
library(MCMCtreeR, quietly = TRUE, warn.conflicts = FALSE)
library(data.table)
library(ggtree)
if (!any(rownames(installed.packages()) == "paleotree")) install.packages("paleotree")
if (!any(rownames(installed.packages()) == "strap")) install.packages("strap")
library(strap)
library(phangorn)

#mcmc.goodeids <- read_table("mcmc_bpp-goodieds.txt")
#Goodeid_tree <- ape::read.tree("SVDquartets/Goodeid_copy.tre")
mcmc.goodeids <- read_table("Mcmc_Goodeids_SecondRun.txt")
Goodeid_tree <- ape::read.tree("Phylogenomics_ChapteR3/SVDquartets/ASTRAL_goodeids.tre")
plot(Goodeid_tree)
Goodeid_tree
# Calculate posterior means
apply(mcmc.goodeids, 2, mean)
# Calculate 95% CIs:
t(apply(mcmc.goodeids, 2, quantile, probs=c(.025,.975)))
#coda package to calculate HPD 95% 
coda::HPDinterval(coda::as.mcmc(mcmc.goodeids))

bppr::hominids$mcmc

#Xiphophorus divergence time (3.5e-9) from Powell et al
#Probably better to use given Poecilid mutation rate 
goodeid_divergence<-msc2time.r(mcmc.goodeids, u.mean = 3.5e-9,
                               u.sd=.1e-9, g.mean = 1, g.sd = 0.5)
#Apply the above parameters to estimate mean divergence times.
apply(goodeid_divergence, 2, mean)

options(scipen=999)
#Put divergence times estimates in dataset.
goodeid_ages_ne<-coda::HPDinterval(coda::as.mcmc(goodeid_divergence))
goodeid_ages_ne<-as.data.frame(goodeid_ages_ne)
setDT(goodeid_ages_ne, keep.rownames = TRUE)

#Filter to just get divergence times, not population sizes.
goodeid_ages <- goodeid_ages_ne %>%
  dplyr::filter(str_detect(rn, "^t_"))

##  Plot the density of gene trees used for astral.
#Read in the multi tree
MultiTreeProteinCodingGeneTrees<-read.tree("Phylogenomics_ChapteR3/Goodeid_gene_trees.contact.treefile")
#Plot the densitree for the multi tree
densityTree(MultiTreeProteinCodingGeneTrees,use.edge.length=FALSE,type="phylogram",nodes="inner")

#densiTree(MultiTreeProteinCodingGeneTrees, type = "cladogram",
#         direction = "rightwards")



timetree <- mcmc2multiphylo(Goodeid_tree, goodeid_divergence, "t_", thin=0.00001)
gcon <- Goodeid_tree$tip.label
timetree
phangorn::densiTree(MultiTreeProteinCodingGeneTrees, col="black", alpha=0.04, consensus = gcon, label.offset=.01,
                    direction="leftwards", scaleX = TRUE)
write.nexus(timetree, file = "Goodeid_timetree_bppr", translate = TRUE)

#MCMC.tree.plot(timetree_goodeid)
#timetree_goodeid <- ape::read.nexus("Goodeid_timetree_bppr")
#timetree_goodeid<-timetree_goodeid[[1]]
#timetree_goodeid<-makeNodeLabel(timetree_goodeid, method="number", prefix="Node")

# Plotting geological times.
Geological_times <- read.csv('Goodeid_Ages_GeologicalTime2.csv') # Import geological ages
Geological_times<-Geological_times %>% remove_rownames %>% column_to_rownames(var="Species")

tree_l <- DatePhylo(Goodeid_tree, Geological_times, method="equal", rlen=1)


geoscalePhylo(tree=ladderize(tree_l,right=FALSE),ages=Geological_times,ranges=TRUE,units=c("Period","Epoch"),boxes = "Epoch",
              cex.tip=0.8,cex.ts=0.05,cex.age=0.9, erotate = 0, arotate=0, x.lim=c(25,0),width=3)
plot(timetree[[1]])


# Plot densitree for the SVDQuartets

btrees<-read.nexus(file="Phylogenomics_ChapteR3/SVDquartets/Goodeid_phylogeny_filtered_dataset.tre",
                   tree.names=NULL)

ggdensitree(btrees, alpha=.3, colour='black') + 
  geom_tiplab(size=3) + hexpand(.35)

###########################################################################
##########################     G-phocs    ##############################

#Estimating migration rate parameters
mcmc.gphocs <- read_table("G-phocs.mcmc.txt")
#Estimate confidence intervals 
mcmc.gphocs_HPD<-coda::HPDinterval(coda::as.mcmc(mcmc.gphocs))
mcmc.gphocs_HPD<-as.data.frame(mcmc.gphocs_HPD)
g_phocs_scales_estimates<-msc2time.r(mcmc.gphocs, u.mean = 3.5e-9,
                                     u.sd=.1e-9, g.mean = 1, g.sd = 0.5)

#Get means for scales estimates
apply(g_phocs_scales_estimates, 2, mean)
apply(mcmc.goodeids, 2, mean)

################# PLOT densitree for bootstrap trees.
library(phangorn)
BiocManager::install("ggtree")
library(treeio)
library(ggtree)
SVDquarts_nexus<-ape::read.nexus(file='Phylogenomics_ChapteR3/SVDquartets/Goodeid_phylogeny_filtered_dataset.tre')
IQTREE_genetrees<-ape::read.tree(file='Goodeid_gene_trees.contact.treefile')


#as.phylo(SVDquarts_nexus)
#consensus_tree <- read.nexus('Phylogenomics_ChapteR3/SVDquartets/Dated_Goodeids_filteredData')

ggdensitree(SVDquarts_nexus, alpha=.3, colour='steelblue') + 
  geom_tiplab(size=3) + xlim(0, 45)



png('Cloudogram.png', width=7, height=5, units='in', res=300)
densiTree(btrees, type="cladogram")
dev.off()

#densityTree(IQTREE_genetrees,type="cladogram",nodes="intermediate")
densityTree(SVDquarts_nexus,type="cladogram",nodes="intermediate")
#ASTRAL and gene tree concordance analysis.


##########################################################################
############### Dstatistic ######################################


Dmin <- read.table('Phylogenomics_ChapteR3/Goodeids_Dtrios_Dmin.txt', header=T)
Dmin<-within(Dmin, Interaction <- paste(P2, P3, sep='_'))
fdr_pvalues<-p.adjust(Dmin$p.value, method = 'BH', n = 56)
Dmin <- cbind(Dmin, fdr_pvalues)
Dmin_ordered <- Dmin[order(Dmin$Interaction, -abs(Dmin$Dstatistic) ), ]
#Dmin_rmdups<-Dmin_ordered[ !duplicated(Dmin_ordered$Interaction), ]  

mean(Dmin$Dstatistic)
Dmin %>%
  filter(fdr_pvalues < 0.05) %>%
  summarise(mean(Dstatistic), mean(f4.ratio), n = n())


Dmin <- filter(Dmin, fdr_pvalues < 0.05)
#mean(Dmin$Dstatistic)
#Dmin <- filter(Dmin, f4.ratio > 0.01)
#Dmin <- filter(Dmin_rmdups, Dstatistic > 0.1)
Dmin<-within(Dmin, Interaction <- paste(P1,P2, P3, sep='_'))

ggplot(Dmin, aes(x=Interaction, y=Dstatistic, colour=p.value))+
  geom_point(size=3)+theme_minimal()+
  scale_colour_gradientn(colours = pal)+
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size=6))


mean(Dmin$Dstatistic)
mean(Dmin$f4.ratio)


#################################################
### August 2021 
Dmin$P1<-gsub("_"," ",as.character(Dmin$P1))
Dmin$P2<-gsub("_"," ",as.character(Dmin$P2))
Dmin$P3<-gsub("_"," ",as.character(Dmin$P3))

Dmin_plot<-ggplot(data = Dmin, aes(x=P2, y=P3, fill=Dstatistic)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "pink", 
                       midpoint = 0.03, space = "Lab", 
                       name="D-statistic (minimum)\n") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(size=12))+
  coord_fixed()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=8),
    legend.direction = "vertical")

Dmin_plot

f4_ratio<-ggplot(data = Dmin, aes(x=P2, y=P3, fill=f4.ratio)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "pink", 
                       midpoint = 0.015, space = "Lab", 
                       name="f4-ratio\n") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(size=12))+
  coord_fixed()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=8),
    legend.direction = "vertical")

f4_ratio


Dmin_plot+theme(legend.direction = "vertical")

png('Phylogenomics_ChapteR3/f4ratio_Goodeids.png', width=7, height=5, units='in', res=300)
f4_ratio
dev.off()

png('Phylogenomics_ChapteR3/Dmin_Goodeids.png', width=7, height=5, units='in', res=300)
Dmin_plot
dev.off()

### Fbranch statistic data analysis -- where is gene flow actually prevalent.
Fbranch_pvalues<-read_table('Phylogenomics_ChapteR3/FbranchASTRALGoodeids_pvals.txt')

Fbranch_pvals_long<-pivot_longer(Fbranch_pvalues, 
             cols = Outgroup:Ameca_splendens,
             names_to = "Nodes",
             values_to = "Pvalues")

Fbranch_pvals_long$Pvalues_adjust<-p.adjust(Fbranch_pvals_long$Pvalues, method = 'BH', n = 126)
Confident_fbranch<-filter(Fbranch_pvals_long, Pvalues_adjust < 0.05)


############### Gene flow results - f_d across genome ###################
#Things to do, quantify levels of gene flow in the hypothesised direction.
#Questions to ask:
###### Is there more gene flow from AT into monophyly or between GM and monophyly (how much has been preserved?)
###### Do some species show much less gene flow on average than others (purifying selection hinting?)
##### What are the genes being exchanged and are they different between GM and AT?

#First thing to do is load all the files up.
Gir_Ata_XR_Introgression<-read.table('Phylogenomics_ChapteR3/Dinvestigate_conservative_attempt2/Girardinichthys_Ataenobius_Xenotaenia_localFstats__250_100.txt', header=T)
Gir_Ata_Ame_Introgression<-read.table('Phylogenomics_ChapteR3/Dinvestigate_conservative_attempt2/Girardinichthys_Ataenobius_Ameca_localFstats__250_100.txt', header=T)
Gir_Ata_Char_Introgression<-read.table('Phylogenomics_ChapteR3/Dinvestigate_conservative_attempt2/Girardinichthys_Ataenobius_Characodon_localFstats__250_100.txt', header=T)
Gir_Ata_Good_Introgression<-read.table('Phylogenomics_ChapteR3/Dinvestigate_conservative_attempt2/Girardinichthys_Ataenobius_Goodea_localFstats__250_100.txt', header=T)
Gir_Ata_Ilyo_Introgression<-read.table('Phylogenomics_ChapteR3/Dinvestigate_conservative_attempt2/Girardinichthys_Ataenobius_Ilyodon_localFstats__250_100.txt', header=T)
Gir_Ata_XC_Introgression<-read.table('Phylogenomics_ChapteR3/Dinvestigate_conservative_attempt2/Girardinichthys_Ataenobius_Xenoophorus_localFstats__250_100.txt', header=T)


# First we want to figure out the direction of gene flow using the fdm statistic. 

ggplot(Gir_Ata_XC_Introgression, aes(x=f_dM))+
  geom_histogram()+
  geom_vline(xintercept = 0)

ggplot(Gir_Ata_Ilyo_Introgression, aes(x=f_dM))+
  geom_histogram()+
  geom_vline(xintercept = 0)


### Next thing to do is to get the gene names
# files<-c(Gir_Ata_XR_Introgression, Xen_Ame_Gir_Introgression, Xen_Good_Ata_Introgression, 
#   Xen_Good_Gir_Introgression, Xen_Xeno_Ata_Introgression, Xen_Xeno_Gir_Introgression)
# 
# # for(file in files){
# #   file$Gene_Identifier<-paste(file$chr,"-", file$windowEnd)
# #   file$transcriptID<-GM_gtf$V9[match(file$Gene_Identifier, GM_gtf$Gene_Identifier)]
# #   file$Annotation<-annot$Preferred_name[match(file$transcriptID, annot$query_name)]
# # }

Gir_Ata_XR_Introgression %>% filter(f_d > 0) %>% summarise(mean(f_d))
Gir_Ata_Ame_Introgression %>% filter(f_d > 0) %>% summarise(mean(f_d))
Gir_Ata_Char_Introgression %>% filter(f_d > 0) %>% summarise(mean(f_d))
Gir_Ata_Good_Introgression%>% filter(f_d > 0) %>% summarise(mean(f_d))
Gir_Ata_Ilyo_Introgression %>% filter(f_d > 0) %>% summarise(mean(f_d))
Gir_Ata_XC_Introgression %>% filter(f_d > 0) %>% summarise(mean(f_d))


Gir_Ata_XR_Introgression <- Gir_Ata_XR_Introgression %>% filter(f_d > 0) %>%filter(f_d > quantile(f_d, .90))
Gir_Ata_Ame_Introgression <- Gir_Ata_Ame_Introgression %>% filter(f_d > 0) %>% filter(f_d > quantile(f_d, .90))
Gir_Ata_Char_Introgression <- Gir_Ata_Char_Introgression %>% filter(f_d > 0) %>% filter(f_d > quantile(f_d, .90))
Gir_Ata_Good_Introgression <- Gir_Ata_Good_Introgression%>% filter(f_d > 0) %>% filter(f_d > quantile(f_d, .90))
Gir_Ata_Ilyo_Introgression <- Gir_Ata_Ilyo_Introgression %>% filter(f_d > 0) %>% filter(f_d > quantile(f_d, .90))
Gir_Ata_XC_Introgression <- Gir_Ata_XC_Introgression %>% filter(f_d > 0) %>% filter(f_d > quantile(f_d, .90))

mean(Gir_Ata_XR_Introgression$f_d)
mean(Gir_Ata_Ame_Introgression$f_d)
mean(Gir_Ata_Char_Introgression$f_d)
mean(Gir_Ata_Good_Introgression$f_d)
mean(Gir_Ata_Ilyo_Introgression$f_d)
mean(Gir_Ata_XC_Introgression$f_d)

Gir_Ata_XR_Introgression_GR<-makeGRangesFromDataFrame(Gir_Ata_XR_Introgression, keep.extra.columns=FALSE,ignore.strand=FALSE)
Gir_Ata_Ame_Introgressio_GR<-makeGRangesFromDataFrame(Gir_Ata_Ame_Introgression, keep.extra.columns=FALSE,ignore.strand=FALSE)
Gir_Ata_Char_Introgression_GR<-makeGRangesFromDataFrame(Gir_Ata_Char_Introgression, keep.extra.columns=FALSE,ignore.strand=FALSE)
Gir_Ata_Good_Introgression_GR<-makeGRangesFromDataFrame(Gir_Ata_Good_Introgression, keep.extra.columns=FALSE,ignore.strand=FALSE)
Gir_Ata_Ilyo_Introgression_GR<-makeGRangesFromDataFrame(Gir_Ata_Ilyo_Introgression, keep.extra.columns=FALSE,ignore.strand=FALSE)
Gir_Ata_XC_Introgression_GR<-makeGRangesFromDataFrame(Gir_Ata_XC_Introgression, keep.extra.columns=FALSE,ignore.strand=FALSE)

#gff_file <- system.file("~/Desktop","Goodeid_Cyprinodon_analysis", "GM.filtered.Omega.2.gtf", package="GenomicFeatures")
gmxdb <- makeTxDbFromGFF("augustus.filter.gff3", format="gff")
gmxdb
#GetOverlaps
subsetByOverlaps(Gir_Ata_XR_Introgression_GR,gmxdb)
#Ameca_splendens and Ataenobious toweri
Gir_Ata_XR_Introgression_GR_OverlapGenes<-transcriptsByOverlaps(gmxdb, Gir_Ata_XR_Introgression_GR)
Gir_Ata_XR_Introgression_GR_OverlapGenes<-as.data.frame(Gir_Ata_XR_Introgression_GR_OverlapGenes)
Gir_Ata_XR_Introgression_GR_OverlapGenes$Annotation<-annot$Preferred_name[match(Gir_Ata_XR_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]
Gir_Ata_XR_Introgression_GR_OverlapGenes$Description<-annot$X.4[match(Gir_Ata_XR_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]

#Goodea_atripinnis and Ataenobious toweri
Gir_Ata_Ame_Introgressio_GR_OverlapGenes<-transcriptsByOverlaps(gmxdb, Gir_Ata_Ame_Introgressio_GR)
Gir_Ata_Ame_Introgressio_GR_OverlapGenes<-as.data.frame(Gir_Ata_Ame_Introgressio_GR_OverlapGenes)
Gir_Ata_Ame_Introgressio_GR_OverlapGenes$Annotation<-annot$Preferred_name[match(Gir_Ata_Ame_Introgressio_GR_OverlapGenes$tx_name, annot$query_name)]
Gir_Ata_Ame_Introgressio_GR_OverlapGenes$Description<-annot$X.4[match(Gir_Ata_Ame_Introgressio_GR_OverlapGenes$tx_name, annot$query_name)]


#Xenoophorus_captivus and Ataenobious toweri
Gir_Ata_Char_Introgression_GR_OverlapGenes<-transcriptsByOverlaps(gmxdb, Gir_Ata_Char_Introgression_GR)
Gir_Ata_Char_Introgression_GR_OverlapGenes<-as.data.frame(Gir_Ata_Char_Introgression_GR_OverlapGenes)
Gir_Ata_Char_Introgression_GR_OverlapGenes$Annotation<-annot$Preferred_name[match(Gir_Ata_Char_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]
Gir_Ata_Char_Introgression_GR_OverlapGenes$Description<-annot$X.4[match(Gir_Ata_Char_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]


###Now for the Girardinicthys comparisons
#Ameca_splendens and GM
Gir_Ata_Good_Introgression_GR_OverlapGenes<-transcriptsByOverlaps(gmxdb, Gir_Ata_Good_Introgression_GR)
Gir_Ata_Good_Introgression_GR_OverlapGenes<-as.data.frame(Gir_Ata_Good_Introgression_GR_OverlapGenes)
Gir_Ata_Good_Introgression_GR_OverlapGenes$Annotation<-annot$Preferred_name[match(Gir_Ata_Good_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]
Gir_Ata_Good_Introgression_GR_OverlapGenes$Description<-annot$X.4[match(Gir_Ata_Good_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]


#Goodea_atripinnis and GM
Gir_Ata_Ilyo_Introgression_GR_OverlapGenes<-transcriptsByOverlaps(gmxdb, Gir_Ata_Ilyo_Introgression_GR)
Gir_Ata_Ilyo_Introgression_GR_OverlapGenes<-as.data.frame(Gir_Ata_Ilyo_Introgression_GR_OverlapGenes)
Gir_Ata_Ilyo_Introgression_GR_OverlapGenes$Annotation<-annot$Preferred_name[match(Gir_Ata_Ilyo_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]
Gir_Ata_Ilyo_Introgression_GR_OverlapGenes$Description<-annot$X.4[match(Gir_Ata_Ilyo_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]

#Xenoophorus and GM
Gir_Ata_XC_Introgression_GR_OverlapGenes<-transcriptsByOverlaps(gmxdb, Gir_Ata_XC_Introgression_GR)
Gir_Ata_XC_Introgression_GR_OverlapGenes<-as.data.frame(Gir_Ata_XC_Introgression_GR_OverlapGenes)
Gir_Ata_XC_Introgression_GR_OverlapGenes$Annotation<-annot$Preferred_name[match(Gir_Ata_XC_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]
Gir_Ata_XC_Introgression_GR_OverlapGenes$Description<-annot$X.4[match(Gir_Ata_XC_Introgression_GR_OverlapGenes$tx_name, annot$query_name)]


genes_Gir_Ata_XC <- Gir_Ata_XC_Introgression_GR_OverlapGenes$tx_name
genes_Gir_Ata_Ilyo <- Gir_Ata_Ilyo_Introgression_GR_OverlapGenes$tx_name
genes_Gir_Ata_Good <- Gir_Ata_Good_Introgression_GR_OverlapGenes$tx_name
genes_Gir_Ata_Char <- Gir_Ata_Char_Introgression_GR_OverlapGenes$tx_name
genes_Gir_Ata_Ame <- Gir_Ata_Ame_Introgressio_GR_OverlapGenes$tx_name
genes_Gir_Ata_XR <- Gir_Ata_XR_Introgression_GR_OverlapGenes$tx_name



list_introgressed_genes<-list(XC_AT=genes_Gir_Ata_XC,GA_AT=genes_Gir_Ata_Good,AS_AT=genes_Gir_Ata_Ame,
                              AT_CH=genes_Gir_Ata_Char, AT_IT=genes_Gir_Ata_Ilyo,
                              AT_XR=genes_Gir_Ata_XR)
library(ggvenn)

png('Introgressed_genes_VENNDIAGRAM_GoodeidsGIR.png', width=7, height=5, units='in', res=300)
ggvenn(
  list_introgressed_genes_Gir,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4, text_size = 3)
dev.off()

png('Introgressed_genes_VENNDIAGRAM_GoodeidsATA.png', width=7, height=5, units='in', res=300)
ggvenn(
  list_introgressed_genes_Ata,
  fill_color = c("#CD534CFF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4, text_size = 3)
dev.off()


#Now get the middle genes where there is complete overlap.
#list_overlapacrossall_GM<-calculate.overlap()
#list_overlapacrossall_GM$a5
#list_overlapacrossall_AT<-calculate.overlap(list_introgressed_genes_Ata)

#get dataframe of all the genes so that you can plot.
Dataframe_overlap_introgressed_genes<-data.frame(do.call(cbind, list_introgressed_genes))

#Get overlap of all the genes across the different comparisons.
Overlapping_introgressed_genes<-data.frame(Reduce(intersect, list(genes_Gir_Ata_XC,genes_Gir_Ata_Ilyo, genes_Gir_Ata_Ame,
                       genes_Gir_Ata_Char, genes_Gir_Ata_Good, genes_Gir_Ata_XR)))

# Annotate the overlap of all 6 comparisons.
names(Overlapping_introgressed_genes)[1] <- "gene_id"
Overlapping_introgressed_genes$Annotation<-annot$Preferred_name[match(Overlapping_introgressed_genes$gene_id, annot$query_name)]
Overlapping_introgressed_genes$Description<-annot$X.4[match(Overlapping_introgressed_genes$gene_id, annot$query_name)]


#Get the complete overall of genes.
Genes_Overlapping_All_GIR<-subset(Gir_Ata_XC_Introgression_GR_OverlapGenes, tx_name %in% list_overlapacrossall_GM$a5)
Genes_Overlapping_All_ATA<-subset(Gir_Ata_Char_Introgression_GR_OverlapGenes, tx_name %in% list_overlapacrossall_AT$a5)







#Very conservative gene set
gene.df <- bitr(Overlapping_introgressed_genes$Annotation, fromType = "SYMBOL",
                 toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                 OrgDb = organism)
 
ego <- enrichGO(gene          = gene.df$ENTREZID,
                 OrgDb         = organism,
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.5)
Genes_Overlapping_Introgression_geneontology<-as.data.frame(ego)


png(file="Phylogenomics_ChapteR3/dotplotFDGeneFlowHighConfidenceGO.png",width=7,height=8)
dotplot(ego, showCategory=20, color="pvalue", orderBy = "x")
dev.off()




# Very lenient gene set -- can we recover the initial cohort of genes that were introgressed.
lenient_introgressed.df <- bitr(Gir_Ata_XR_Introgression_GR_OverlapGenes$Annotation, fromType = "SYMBOL",
                toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                OrgDb = organism)

ego_lenient_intro <- enrichGO(gene          = lenient_introgressed.df$ENTREZID,
                OrgDb         = organism,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.5)
Genes_XRLenient_Introgression_geneontology<-as.data.frame(ego_lenient_intro)


png(file="Phylogenomics_ChapteR3/dotplotFDGeneFlowHighConfidenceGO.png",width=7,height=8)
dotplot(ego, showCategory=20, color="pvalue", orderBy = "x")
dev.off()



##  What about the chromosomal position of introgression -- is it localised to particular chromosomes.
GmAtXR_Chr<-read.table('Phylogenomics_ChapteR3/Dinvestigate_attempt3_newcords/Girardinichthys_Ataenobius_Xenotaenia_localFstats__250_100.txt', header=T)
GmAtAme_Chr<-read.table('Phylogenomics_ChapteR3/Dinvestigate_attempt3_newcords/Girardinichthys_Ataenobius_Ameca_localFstats__250_100.txt', header=T)
GmAtChar_Chr<-read.table('Phylogenomics_ChapteR3/Dinvestigate_attempt3_newcords/Girardinichthys_Ataenobius_Characodon_localFstats__250_100.txt', header=T)
GmAtGood_Chr<-read.table('Phylogenomics_ChapteR3/Dinvestigate_attempt3_newcords/Girardinichthys_Ataenobius_Goodea_localFstats__250_100.txt', header=T)
GmAtIlyo_Chr<-read.table('Phylogenomics_ChapteR3/Dinvestigate_attempt3_newcords/Girardinichthys_Ataenobius_Ilyodon_localFstats__250_100.txt', header=T)
GmAtXC_Chr<-read.table('Phylogenomics_ChapteR3/Dinvestigate_attempt3_newcords/Girardinichthys_Ataenobius_Xenoophorus_localFstats__250_100.txt', header=T)


#Convert chromosome names to something more readable.
GmAtXR_Chr$Chromosome_Names<-Gm_scaffold_names$V1[match(GmAtXR_Chr$chr, Gm_scaffold_names$V2)]
GmAtAme_Chr$Chromosome_Names<-Gm_scaffold_names$V1[match(GmAtAme_Chr$chr, Gm_scaffold_names$V2)]
GmAtChar_Chr$Chromosome_Names<-Gm_scaffold_names$V1[match(GmAtChar_Chr$chr, Gm_scaffold_names$V2)]
GmAtGood_Chr$Chromosome_Names<-Gm_scaffold_names$V1[match(GmAtGood_Chr$chr, Gm_scaffold_names$V2)]
GmAtIlyo_Chr$Chromosome_Names<-Gm_scaffold_names$V1[match(GmAtIlyo_Chr$chr, Gm_scaffold_names$V2)]
GmAtXC_Chr$Chromosome_Names<-Gm_scaffold_names$V1[match(GmAtXC_Chr$chr, Gm_scaffold_names$V2)]


#Plot genome-wide f_d for all.


##########################################################################
#################### Estimating Divergence (dxy) genome-wide ######################################
######
library(data.table)

Dxy_goodeids<-fread('Phylogenomics_ChapteR3/pixy_dxy.txt.gz', sep='\t')
# Filter dxy with less than 500 sites to only get robust estimates.
# Filter Crenicthys baileyi 
Dxy_goodeids <- Dxy_goodeids %>%
  na.omit() %>%
  filter(no_sites > 1000)


# rename start and end windows to be more easily understood.

names(Dxy_goodeids)[names(Dxy_goodeids) == "window_pos_1"] <- "start"
names(Dxy_goodeids)[names(Dxy_goodeids) == "window_pos_2"] <- "end"

# Gm scaffolds to real chromosome names

Gm_scaffold_names <- read.csv("Phylogenomics_ChapteR3/OldGm_chromandscaffolds.txt",
                              header=FALSE)

#Change the column names in dxy datasets to match

Dxy_goodeids$chromosome <- sub("_RagTag", "", Dxy_goodeids$chromosome)

Dxy_goodeids$Chromosome_Names<-Gm_scaffold_names$V1[match(Dxy_goodeids$chromosome, Gm_scaffold_names$V2)]


# Plot dxy across chromosomes for different comparisons

Dxy_goodeids %>%
  filter(pop1 == 'Crenicthys_baileyi') %>%
  mutate(Chromosome_Names = factor(Chromosome_Names, levels = c(1:24, "X", "Y"))) %>%
  ggplot(aes(x=Chromosome_Names, y=avg_dxy, color=pop2))+
  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom="pointrange", position=position_dodge(width=.75))+
  scale_colour_viridis(option='D', discrete=T)+
  theme_classic()+
  labs(y=expression(D[XY]),
       x="Chromosomes",
       colour="Species (against C.baileyi)")


# Do the same thing but look at the split of GM from all others
Dxy_goodeids %>%
  filter(pop1 == 'Girardinicthys_multiradiatus') %>%
  mutate(Chromosome_Names = factor(Chromosome_Names, levels = c(1:24, "X", "Y"))) %>%
  ggplot(aes(x=Chromosome_Names, y=avg_dxy, color=pop2))+
  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom="pointrange", position=position_dodge(width=.75))+
  scale_colour_viridis(option='D', discrete=T)+
  theme_classic()+
  labs(y=expression(D[XY]),
       x="Chromosomes",
       colour="Species (against G.multiradiatus)")




# Shows reduced divergence on the X and Y across all species when compared to C.baileyi...
# Also see a small increase in chromosome 14?
# Reduced dxy on sex chromosomes likely due to population contractions -- affect X chromosomes more than Y chromosome.


Dxy_goodeids %>%
  filter(pop1 == 'Girardinicthys_multiradiatus' & pop2=="Characodon_lateralis" & Chromosome_Names == "14") %>%
  ggplot(aes(x = (start + end)/2, y = avg_dxy, colour=pop2))+
  geom_jitter()+
  xlab("Genomic Position (MB)")+
  ylab("Statistic Value (dxy)")+
  theme_classic()+
  scale_colour_viridis(option='D', discrete=T)


# Get genes by creating a variable to match with gtf file. 
#Dxy_goodeids_againstCrenicthys$Gene_Identifier<-paste(Dxy_goodeids_againstCrenicthys$chromosome,"-", Dxy_goodeids_againstCrenicthys$end)
#Dxy_goodeids_againstCrenicthys$transcriptID<-GM_gtf$V9[match(Dxy_goodeids_againstCrenicthys$Gene_Identifier, GM_gtf$Gene_Identifier)]

#Categorise species by biology to extract potential candidates
# Dxy_goodeids %>%Dxy_goodeids
#   mutate(Trophotaenia_character = case_when(pop2 %in% c("Ameca_splendens","Goodea_atripinnis","Girardinicthys_multiradiatus",
#                                                       "Xenotaenia_resolanae", "Ilyodon_furcidens","Characodon_lateralis",
#                                                       "Xenoophorus_captivus", "Ataenobious_toweri") ~ "Trophotaenia",
#                                           pop2 %in% c("Ataenobious_toweri") ~ "No_trophotaenia"))

Dxy_goodeids <- Dxy_goodeids %>% 
  mutate(Chromosome_categories = case_when(Chromosome_Names %in% c(1:24) ~ "Autosomes",
                                            Chromosome_Names %in% c("X") ~ "X chromosome",
                                            Chromosome_Names %in% c("Y") ~ "Y chromosome"))
  
# What is the distribution of dxy for aut v.s. sex chromosomes amongst species pairs.

Dxy_Line_AM_XC <- Dxy_goodeids %>%
  filter(pop1 == 'Ameca_splendens' & pop2=="Xenoophorus_captivus" & Chromosome_categories == c("X chromosome", "Autosomes")) %>%
  group_by(Chromosome_categories)

t.test(avg_dxy ~ Chromosome_categories, Dxy_Line_AM_XC)

Dxy_Line_AT_GM <- Dxy_goodeids %>%
  filter(pop1 == 'Ataenobious_toweri' & pop2=="Girardinicthys_multiradiatus" & Chromosome_categories == c("X chromosome", "Autosomes")) %>%
  group_by(Chromosome_categories) %>%
  summarize(mean_x = mean(avg_dxy))


Dxy_Line_AT_GM <- Dxy_goodeids %>%
  filter(pop1 == 'Ataenobious_toweri' & pop2=="Girardinicthys_multiradiatus" & Chromosome_categories == c("X chromosome", "Autosomes")) %>%
  group_by(Chromosome_categories) %>%
  summarize(mean_x = mean(avg_dxy))



t.test(avg_dxy ~ Chromosome_categories, Dxy_Line_AT_GM)

Dxy_Line_IF_XR <- Dxy_goodeids %>%
  filter(pop1 == 'Xenotaenia_resolanae' & pop2=="Ilyodon_furcidens" & Chromosome_categories == c("X chromosome", "Autosomes")) %>%
  group_by(Chromosome_categories) %>%
  summarize(mean_x = mean(avg_dxy))

t.test(avg_dxy ~ Chromosome_categories, Dxy_Line_IF_XR)

#Plot the difference between Autosomes and X chromosome
# XC v.s AS
png(file="Phylogenomics_ChapteR3/DensityPlot_AutX_Dxy_ASXC.png",width=4,height=3, units='in', res=300)
Dxy_goodeids %>%
  filter(pop1 == 'Ameca_splendens' & pop2=="Xenoophorus_captivus" & 
           Chromosome_categories == c("X chromosome", "Autosomes")) %>%
    ggplot(aes(x=avg_dxy, fill=Chromosome_categories,
               colour = Chromosome_categories))+
    geom_density(aes(y =..density..), 
                 alpha=0.5)+
  theme_pubr()+
    #facet_wrap(~Chromosome_categories, scales="free")+
    geom_vline(data = Dxy_Line_AM_XC, aes(xintercept = mean_x, colour = Chromosome_categories), 
               size=0.75, linetype='dashed')+
    scale_fill_viridis(option='D', discrete=T)+
    scale_colour_viridis(option='D', discrete=T)+
    labs(x=expression(D[XY]),
         fill="Chromosome",
         colour="Chromosome")+
  xlim(0,0.05)
dev.off()    

#AT v.s. GM
png(file="Phylogenomics_ChapteR3/DensityPlot_AutX_Dxy_ATGM.png",width=4,height=3, units='in', res=300)
Dxy_goodeids %>%
  filter(pop1 == 'Ataenobious_toweri' & pop2=="Girardinicthys_multiradiatus" & 
           Chromosome_categories == c("X chromosome", "Autosomes")) %>%
  ggplot(aes(x=avg_dxy, fill=Chromosome_categories,
             colour = Chromosome_categories))+
  geom_density(aes(y =..density..), 
               alpha=0.5)+
  theme_pubr()+
  #facet_wrap(~Chromosome_categories, scales="free")+
  geom_vline(data = Dxy_Line_AT_GM, aes(xintercept = mean_x, colour = Chromosome_categories), 
             size=0.75, linetype='dashed')+
  scale_fill_viridis(option='D', discrete=T)+
  scale_colour_viridis(option='D', discrete=T)+
  labs(x=expression(D[XY]),
       fill="Chromosome",
       colour="Chromosome")+
  xlim(0,0.05)
dev.off()

#IF v.s. XR
png(file="Phylogenomics_ChapteR3/DensityPlot_AutX_Dxy_IFXR.png",width=4,height=3, units='in', res=300)
Dxy_goodeids %>%
  filter(pop1 == 'Xenotaenia_resolanae' & pop2=="Ilyodon_furcidens" & 
           Chromosome_categories == c("X chromosome", "Autosomes")) %>%
  ggplot(aes(x=avg_dxy, fill=Chromosome_categories,
             colour = Chromosome_categories))+
  geom_density(aes(y =..density..), 
               alpha=0.5)+
  theme_pubr()+
  #facet_wrap(~Chromosome_categories, scales="free")+
  geom_vline(data = Dxy_Line_IF_XR, aes(xintercept = mean_x, colour = Chromosome_categories), 
             size=0.75, linetype='dashed')+
  scale_fill_viridis(option='D', discrete=T)+
  scale_colour_viridis(option='D', discrete=T)+
  labs(x=expression(D[XY]),
       fill="Chromosome",
       colour="Chromosome")+
  xlim(0,0.05)
dev.off()

### Plot chromosomal dxy for all species to find 
Dxy_goodeids %>%
  filter(pop1 == 'Crenicthys_baileyi') %>%
  ggplot(aes(x=Chromosome_Names, y=avg_dxy, colour=pop2))+
  geom_jitter(size=1.75,alpha=0.7)+
  scale_colour_viridis(option='D', discrete=TRUE)+labs(y=expression('d'[XY]))+
  theme_minimal_hgrid()+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=13,face="bold"))+
  facet_grid(cols = vars(Chromosome_Names),
             space = "free_x",
             scales = "free_x",
             switch = "x") +
  labs(x = "Chromosome") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  theme(axis.line.y=element_line(colour = 'grey'))





# Use this grouping characteristic to collect average  
Dxy_Viviparity_character<-Dxy_goodeids_againstCrenicthys_groupsadded %>%
  dplyr::select(transcriptID, Trophotaenia_character, avg_dxy) %>%
  pivot_wider(names_from = Trophotaenia_character, values_from = avg_dxy, values_fn=mean)


Difference_DXY_againstCrenichtys$Annotation<-annot$Preferred_name[match(Difference_DXY_againstCrenichtys$transcriptID, annot$query_name)]
Difference_DXY_againstCrenichtys$Description<-annot$X.4[match(Difference_DXY_againstCrenichtys$transcriptID, annot$query_name)]

####How many show positive difference between trophotaenia and non-trophotaenia species?
Tropho_genes_potentially_dxy<-Difference_DXY_againstCrenichtys %>%
  top_frac(.05, Diff_dxy)





gene_dxy.df <- bitr(Tropho_genes_potentially_dxy$Annotation, fromType = "SYMBOL",
                toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                OrgDb = organism)

ego_dxy <- enrichGO(gene          = gene_dxy.df$ENTREZID,
                OrgDb         = organism,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1)

TROPHO_DXY_ALL<-as.data.frame(ego_dxy)
dotplot(ego, showCategory=30)


#For Yoli:

Dxy_goodeids %>%
  filter(pop1 == "Goodea_atripinnis" & pop2=="Xenoophorus_captivus" & Chromosome_categories == c("X chromosome", "Autosomes")) %>%
  group_by(Chromosome_categories) %>%
  summarize(mean_x = mean(avg_dxy))


#########################################################################
#####################    Evolutionary rate  #########################
##########################################################################
species_list_goodeids <- c('transcriptID', 'transcriptID2', 'Ameca_splendens', 'Ataenobius_toweri', 'Crenicthys_baileyi',
                           'Characodon_lateralis', 'Goodea_atripinnis','Girardinicthys_multiradiatus',
                           'Ilyodon_furcidens','Xenoophorus_captivus','Xenotaenia_resolanae')

#Omega - READ file in.
Goodeid_TerminalBranch_Omega<- read.delim(gzfile("FitMG94.Phylogenetics.LineageSpecificOmega.results.gz"))
#Then transpose the table so that all species have their own column rather gene.
Goodeid_TerminalBranch_Omega_tranpose <- data.table::transpose(Goodeid_TerminalBranch_Omega)
#Make sure genes are not rownames but their own columns.
rownames(Goodeid_TerminalBranch_Omega_tranpose) <- colnames(Goodeid_TerminalBranch_Omega_tranpose$V1)
#Name columns based on species_names.
colnames(Goodeid_TerminalBranch_Omega_tranpose)<-species_list_goodeids
#Check to see if it looks right.
head(Goodeid_TerminalBranch_Omega_tranpose)
#Get rid of additional genename column. 
Goodeid_TerminalBranch_Omega_tranpose = subset(Goodeid_TerminalBranch_Omega_tranpose, select = -c(transcriptID2) )


#Synonymous sub rate 
Goodeid_TerminalBranch_dS<- read.delim(gzfile("FitMG94.Phylogenetics.LineageSpecificOmega.dS.results.gz"))
List_col_names_dS<-colnames(Goodeid_TerminalBranch_dS)
Goodeid_TerminalBranch_dS_tranpose <- data.table::transpose(Goodeid_TerminalBranch_dS, keep.names = "rn")
rownames(Goodeid_TerminalBranch_dS_tranpose) <- colnames(Goodeid_TerminalBranch_dS_tranpose$V1)
colnames(Goodeid_TerminalBranch_dS_tranpose)<-species_list_goodeids
head(Goodeid_TerminalBranch_dS_tranpose)
#Goodeid_TerminalBranch_dS_tranpose<-t(Goodeid_TerminalBranch_dS)
#Goodeid_TerminalBranch_dS_tranpose<-as.data.frame(Goodeid_TerminalBranch_dS_tranpose)

#Nonsynonymous sub rate 
Goodeid_TerminalBranch_dN<- read.delim(gzfile("FitMG94.Phylogenetics.LineageSpecificOmega.dN.results.gz"))
List_col_names_dN<-colnames(Goodeid_TerminalBranch_dN)
Goodeid_TerminalBranch_dN_tranpose <- data.table::transpose(Goodeid_TerminalBranch_dN, keep.names = "rn")
rownames(Goodeid_TerminalBranch_dN_tranpose) <- colnames(Goodeid_TerminalBranch_dN_tranpose$V1)
colnames(Goodeid_TerminalBranch_dN_tranpose)<-species_list_goodeids
#Check the top of the file.
head(Goodeid_TerminalBranch_Omega_tranpose)

#### Get the positions of the guppy genome. 
# 
# Guppy_coordinates<-fread('Gmultiradiatus_to_guppy.genes.bed', header=F)
# #Add column names
# colnames(Guppy_coordinates) <- c('Chromosomes', 'start', 'end','transcriptID')
# #Add suffix to elements in 'gene'
# Guppy_coordinates$transcriptID<-lapply(Guppy_coordinates$transcriptID, paste0, ".t1")
# #Join together guppy coordinates with omega table.
# Omega_GuppyCoordinates_Goodeid<-merge(Goodeid_TerminalBranch_Omega_tranpose, Guppy_coordinates, by='transcriptID')
# 

#Pivot dataset for omega from wide to long.
Goodeid_omega_terminal<-pivot_longer(Goodeid_TerminalBranch_Omega_tranpose,
             cols=Ameca_splendens:Xenotaenia_resolanae,
             names_to="Species",
             values_to='Omega')

#Filter anything larger than dn/ds of 10.
Goodeid_omega_terminal<-Goodeid_omega_terminal %>%
  dplyr::filter(Omega < 10)

#Add the mid point for the plot. 
# Goodeid_omega_terminal<-Goodeid_omega_terminal %>%
#   mutate(mid = (start+end)/2)

##### Put in the Guppy coordinates since they have linkage groups
# Guppy_old_newname<-read_csv('Guppy_oldname_newname.csv')
# #Match and add to main dataframe.
# Goodeid_omega_terminal$Linkage_group<-Guppy_old_newname$New_col_names[match(Goodeid_omega_terminal$Chromosomes, Guppy_old_newname$Old_col_names)]
# Guppy_coordinates$Linkage_group <- Guppy_old_newname$New_col_names[match(Guppy_coordinates$Chromosomes, Guppy_old_newname$Old_col_names)]
# #Make sure omega is a double.
# Goodeid_omega_terminal$Omega<-as.double(Goodeid_omega_terminal$Omega)


##Average omega 
averaged_omega_over_genes<-Goodeid_omega_terminal %>%
  dplyr::select(transcriptID, Omega, mid, Linkage_group) %>%
  group_by(transcriptID) %>%
  summarise(Ave_omega=mean(Omega)) %>%
  ungroup()


#Goodeid_omega_terminal$<-averaged_omega_over_genes$Ave_omega[match(Goodeid_omega_terminal$transcriptID, averaged_omega_over_genes$transcriptID)]


### Get distributions of Omega for different species. This gives an indication of what type of selection may be dominating.
##  We may have an expectation that species with lower effective population sizes have higher dN/dS, and vice versa.

ggplot(Goodeid_omega_terminal, aes(x=Omega, fill=Species))+
  geom_density(aes(y = ..count..), alpha=0.1)+
  scale_fill_viridis(option='D', discrete=T)+
  geom_vline(xintercept = 1.0, linetype='dashed')+
  theme_bw()

ggplot(Goodeid_omega_terminal, aes(y=Omega, x=Species, fill=Species))+
  geom_boxplot(width=0.1)+
  scale_fill_viridis(option='D', discrete=T)+
  geom_vline(xintercept = 1.0, linetype='dashed')+
  theme_bw()+
  coord_flip()


# Do sexually dimorphic species have greater divergence in trophotaenia genes and in genes under positive selection?

# Classify species into dimoprhic, intermediate dimorphism and monomorphic based on summary of previous behavioural work.
Goodeid_omega_terminal<-Goodeid_omega_terminal %>%
  mutate(Sexual_Dimorphism = case_when(Species %in% c("Ameca_splendens","Girardinicthys_multiradiatus", "Characodon_lateralis") ~ "Dimorphic",
                              Species %in% c("Crenicthys_baileyi", "Goodea_atripinnis", "Ilyodon_furcidens") ~ "Monomorphic",
                              Species %in% c('Ataenobius_toweri', "Xenoophorus_captivus", "Xenotaenia_resolanae") ~ "Intermediate"))

# Pivot to get classification
Goodeid_omega_terminal_SD<-Goodeid_omega_terminal %>%
  group_by(Sexual_Dimorphism, transcriptID) %>%
  mutate(SD_mean_omega = mean(Omega))


# Now we want to get only the trophotaenia specific genes.
# Add in trophotaenia genes from Manfred paper -- is there a significant enrichment of positively selected genes found in trophotaenia (more than expected by chance?)
Trophotaenia_data_Manfred <- read.csv('Phylogenomics_ChapteR3/GMManfredPaper_SupplementTab8_TrophotaeniaRNAseq.csv')
# filter them by their adjusted p-value
Trophotaenia_genes_Manfred <- Trophotaenia_data_Manfred %>%
  filter(p.adj.tr.rest < 0.05)

#Now we want to match using the name of annotations.
Goodeid_omega_terminal$Annotation<-annot$Preferred_name[match(Goodeid_omega_terminal$transcriptID, annot$query_name)]
Goodeid_omega_terminal$Description<-annot$X.4[match(Goodeid_omega_terminal$transcriptID, annot$query_name)]
# Now match with Trophotaenia genes from Manfred
# Make sure case sensitivity is adjusted.
Goodeid_omega_terminal[[5]]<-tolower(Goodeid_omega_terminal[[5]])
Trophotaenia_data_Manfred[[37]]<-tolower(Trophotaenia_data_Manfred[[37]])
Goodeid_omega_terminal_troSD=Goodeid_omega_terminal[Goodeid_omega_terminal$Annotation %in% Trophotaenia_data_Manfred$symbol, ]
#Finally, we want gene-wise means for each category of sexual dimorphism.
Goodeid_omega_terminal_troSD_summary <- Goodeid_omega_terminal_troSD %>%
  group_by(transcriptID, Sexual_Dimorphism)
  summarise(SexDmean_omega = mean(as.double(Omega)))

# Test for significant differences and follow-up with post-hoc if significant.  
kruskal.test(SexDmean_omega ~ Sexual_Dimorphism, data = Goodeid_omega_terminal_troSD_summary)
# There is a significant result:
#Kruskal-Wallis chi-squared = 22.669, df = 2, p-value = 1.196e-05
# Between who: posthoc test is wilcox
pairwise.wilcox.test(Goodeid_omega_terminal_troSD_summary$SexDmean_omega, 
                     Goodeid_omega_terminal_troSD_summary$Sexual_Dimorphism,
                     p.adjust.method = "BH")

# So we observe a very significant difference in evolution of trophotaenia genes between monomorphic/dimorphic and intermediate, but not between dimorphic and monomorphic.
# 

ggplot(Goodeid_omega_terminal_troSD_summary, aes(x=Sexual_Dimorphism, y=SexDmean_omega, fill=Sexual_Dimorphism))+
  geom_boxplot(width=0.3)+
  theme_bw()+
  labs(x="Sexual Dimorphism",
       y=expression(d[N]/d[S]),
       fill="Sexual Dimorphism")

# Do we find the same thing for the genomic background though?
Goodeid_omega_terminal_control <- Goodeid_omega_terminal %>%
  group_by(transcriptID, Sexual_Dimorphism) %>%
  summarise(SexDmean_omega = mean(as.double(Omega)))

kruskal.test(SexDmean_omega ~ Sexual_Dimorphism, data = Goodeid_omega_terminal_control)

ggplot(Goodeid_omega_terminal_control, aes(x=Sexual_Dimorphism, y=SexDmean_omega, fill=Sexual_Dimorphism))+
  geom_boxplot(width=0.3)+
  theme_bw()


gene_outlierdnds.df <- bitr(Gene_candidates_Trophotaenia$Annotation, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                    OrgDb = organism)

ego_dnds <- enrichGO(gene          = gene_outlierdnds.df$ENTREZID,
                    OrgDb         = organism,
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1)

TROPHO_DNDS_ALL<-as.data.frame(ego_dnds)
dotplot(ego_dnds, showCategory=30)


#########################################################################
#####################    Branch-site test Absrel   #########################
##########################################################################
library(purrr)
library(tidyverse)
library(jsonlite)
BiocManager::install("org.Hs.eg.db")

#path <- "~/Desktop/MOUNTPOINT/"
#files <- dir(path, pattern = "*.txt")


#Read in results from ABSREL
ABSREL_results<- read_table('ABSREL_ALL_GENES_TESTED.2.txt', col_names = F)

#Once you've read in results, name the columns something intuitive.
# Total genes tested: 20,138
ABSREL_results<-ABSREL_results %>% 
  dplyr::rename(
    gene = X1,
    node = X2,
    baseline_omega = X3,
    num_rate_class = X4,
    tested = X5,
    prop_sites_selected = X6,
    LRT = X7,
    uncorrected_p = X8,
    corrected_p = X9
    
  )

#Now only retrieve the significant results (corrected pvalue)
ABSREL_results<-ABSREL_results %>%
  dplyr::filter(corrected_p <= 0.05) 

#left with 1735 records, but we want to know how many genes
Number_of_genes_under_pos_sel_internalbranch<-unique(ABSREL_results$gene)
length(Number_of_genes_under_pos_sel_internalbranch)
### 1,482 genes show evidence of positive selection on an internal branch.
# As a percentage:
(1482/20138)*100
(1293/20138)*100
(168/20138)*100
(21/20138)*100

#On branch as percentage
#Leading to GA_AS
(108/20138)*100
(73/20138)*100
(97/20138)*100
(431/20138)*100
(404/20138)*100
(579/20138)*100

#Now we want to get their annotations.
#Need to add suffix to match values. 
ABSREL_results$gene<-lapply(ABSREL_results$gene, paste0, ".t1")
ABSREL_results$Annotation<-annot$Preferred_name[match(ABSREL_results$gene, annot$query_name)]
ABSREL_results$Description<-annot$X.4[match(ABSREL_results$gene, annot$query_name)]

organism = "org.Dr.eg.db"
#organism = "org.Hs.eg.db"

#How many genes under selection in more than one branch?
Repeated_positive_sel_genes<-ABSREL_results %>% 
  group_by(gene, Annotation, Description) %>% 
  #filter(n()>1) %>%
  mutate(n = n(),
            mean_prop_sites_under_selection=mean(prop_sites_selected))

#Repeated_positive_sel_genes<-distinct_at(Repeated_positive_sel_genes, vars(gene))
Repeated_positive_sel_genes %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(n == 2)


absrel_gene.df <- bitr(ABSREL_results$Annotation, fromType = "SYMBOL",
                toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                OrgDb = organism)

absrel_ego <- enrichGO(gene          = absrel_gene.df$ENTREZID,
                OrgDb         = organism,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1)

ABSREL_GO<-as.data.frame(absrel_ego)
ABSREL_GO <- filter(ABSREL_GO, p.adjust < 0.05)
#Plot and save dotplot with annotations
png('dotplotABSRELPhylogenomics.png', width=7, height=5, units='in', res=300)
dotplot(absrel_ego, showCategory=20)
dev.off()

#Very quickly we find genes evolving that are important for DNA damage, reproductive processes,
# cillium movement etc.


# How many of these genes are involved in mammalian pregnancy or canonical immune related pathways.
# How many are known imprinted genes?
mouse_pregnant_genes <- read_tsv("Genes_mammalian_pregnancy(1).txt")
adaptive_immune_genes <- read_tsv("adaptive_immune_genes(1).txt")
imprinted_genes <- read_tsv("Imprinted_genes.txt")

mouse_pregnant_genes[[2]]<-tolower(mouse_pregnant_genes[[2]])
adaptive_immune_genes[[2]]<-tolower(adaptive_immune_genes[[2]])
imprinted_genes[[2]]<-tolower(imprinted_genes[[2]])
ABSREL_results[[10]]<-tolower(ABSREL_results[[10]])

intersect(mouse_pregnant_genes$Symbol, ABSREL_results$Annotation)
intersect(adaptive_immune_genes$Symbol, ABSREL_results$Annotation)
intersect(imprinted_genes$Symbol, ABSREL_results$Annotation)

mkk <- enrichKEGG(gene = gene.df$ENTREZID,
                   organism = 'dre',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
KEGG_absrel<-as.data.frame(mkk)

#absrel_WP<-enrichWP(gene.df$ENTREZID, organism = "Danio rerio") 
#WP_absrel <- as.data.frame(absrel_WP)


#How many genes show evidence of positive selection on each branch (is it uniform?)
Node_3_results<- Repeated_positive_sel_genes %>%
  filter(node == "Node3") # 108 genes
Node_6_results<-Repeated_positive_sel_genes %>%
  filter(node == "Node6") #73 genes
Node_8_results<-Repeated_positive_sel_genes %>%
  filter(node == "Node8") #97 genes
Node_9_results<-Repeated_positive_sel_genes %>%
  filter(node == "Node9") # 431 genes
Node_11_results<-Repeated_positive_sel_genes %>%
  filter(node == "Node11") #404 genes
Node_14_results<-Repeated_positive_sel_genes %>%
  filter(node == "Node14") # 579 genes

unique(ABSREL_results$node)
#intersect(Gir_Ata_XR_Introgression_GR_OverlapGenes$Annotation, Node_11_results$Annotation)

# What genes are repeatedly evolving under positive selection?
######
# DNA damage response under positive selection. 
# Hypothesis: DNA damage induced by volcanic activity?


ABSREL_proteins_name_STRING<-ABSREL_results$Annotation
write(ABSREL_proteins_name_STRING,'ABSREL_proteins_selection.txt')


# Are positive selected proteins distributed non-randomly genome-wide.
#ABSREL_results$Linkage_group<-Guppy_coordinates$Linkage_group[match(ABSREL_results$gene, Guppy_coordinates$transcriptID)]




### Some stats -- what is the proportion of positively selected 



##################################################################################################
#############      MEME - estimating codons under selection in internal branches  ###################
###################################################################################################

#Read in MEME file 
MEME_results<- read_table('MEME_Sites_Under_Selection_Internal_Branch_Goodeids.txt', col_names = F)

#Once you've read in results, name the columns something intuitive.
MEME_results<-MEME_results %>% 
  dplyr::rename(
    gene = X1,
    site = X2,
    alpha = X3,
    beta_neg = X4,
    prop_beta_neg = X5,
    beta_pos = X6,
    prop_beta_pos = X7,
    LRT = X8,
    p_value = X9,
    num_branches_pos_sel = X10,
    Total_branch_len = X11, 
    MEME_logl = X12, 
    FEL_logl = X13)


#left with 1735 records, but we want to know how many genes
Number_of_genes_under_pos_sel_internalbranch_MEME<-unique(MEME_results$gene)
length(Number_of_genes_under_pos_sel_internalbranch_MEME)


#Add annotations:
MEME_results$Annotation<-annot$Preferred_name[match(MEME_results$gene, annot$query_name)]
MEME_results$Description<-annot$X.4[match(MEME_results$gene, annot$query_name)]


gene_MEME.df <- bitr(MEME_results$Annotation, fromType = "SYMBOL",
                toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                OrgDb = organism)

ego_m <- enrichGO(gene          = gene_MEME.df$ENTREZID,
                OrgDb         = organism,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)

MEME_ALL<-as.data.frame(ego_m)
dotplot(ego_m, showCategory=30)


No_sites_under_selection_meme<-MEME_results %>% group_by(gene) %>% count(gene)
No_sites_under_selection_meme$Annotation<-annot$Preferred_name[match(No_sites_under_selection_meme$gene, annot$query_name)]
No_sites_under_selection_meme$Description<-annot$X.4[match(No_sites_under_selection_meme$gene, annot$query_name)]

################################## ABSREL -- Detecting ATOWERI genes trophotaenia loss ##############################
atoweri_genes<- read_table('Phylogenomics_ChapteR3/ALLSITES.ATOWERISELECTION.tsv.sites_under_selection.txt', col_names = F)

#Once you've read in results, name the columns something intuitive.
# Total genes tested: 20,138
atoweri_genes<-atoweri_genes %>% 
  dplyr::rename(
    gene = X1,
    node = X2,
    baseline_omega = X3,
    num_rate_class = X4,
    tested = X5,
    prop_sites_selected = X6,
    LRT = X7,
    uncorrected_p = X8,
    corrected_p = X9
    
  )

atoweri_genes<-atoweri_genes %>%
  dplyr::filter(corrected_p <= 0.05) 

atoweri_genes$gene<-lapply(atoweri_genes$gene, paste0, ".t1")
atoweri_genes$Annotation<-annot$Preferred_name[match(atoweri_genes$gene, annot$query_name)]
atoweri_genes$Description<-annot$X.4[match(atoweri_genes$gene, annot$query_name)]

#Length convert into percentage
(1177/20138)*100


# Add in trophotaenia genes from Manfred paper -- is there a significant enrichment of positively selected genes found in trophotaenia (more than expected by chance?)
Trophotaenia_data_Manfred <- read.csv('Phylogenomics_ChapteR3/GMManfredPaper_SupplementTab8_TrophotaeniaRNAseq.csv')
# filter them by their adjusted p-value
Trophotaenia_genes_Manfred <- Trophotaenia_data_Manfred %>%
  filter(p.adj.tr.rest < 0.05)
  
# get symbol consisted
atoweri_genes[[10]]<-tolower(atoweri_genes[[10]])
Trophotaenia_data_Manfred[[37]]<-tolower(Trophotaenia_data_Manfred[[37]])

# Genes that are in trophotaenia and also under positive selection in AT.
list_postive_selected_genesAT_expressedinTR<-intersect(Trophotaenia_data_Manfred$symbol, atoweri_genes$Annotation)
length(list_postive_selected_genesAT_expressedinTR)

# What is the fold change of these genes?
# First we build on creating a dataset with these genes, and add fold change data for each gene.
atoweri_genes_trophotaenia_genes <- atoweri_genes %>%
  filter(Annotation %in% list_postive_selected_genesAT_expressedinTR)


intersect(atoweri_genes_trophotaenia_genes$Annotation, ABSREL_results$Annotation)
# Now add in the other component.
#Got to remove slc37a1 since there are so many duplicates.
atoweri_genes_trophotaenia_genes$logFC.tr.rest<-Trophotaenia_data_Manfred$logFC.tr.rest [match(atoweri_genes_trophotaenia_genes$Annotation, Trophotaenia_data_Manfred$symbol)]



AT_candidate_genes.df <- bitr(list_postive_selected_genesAT_expressedinTR, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                     OrgDb = organism)

at_candidate_ego <- enrichGO(gene          = AT_candidate_genes.df$ENTREZID,
                  OrgDb         = organism,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1)

AT_CANDIDATE_GO<-as.data.frame(at_candidate_ego)
# Does this overlap with any high impact mutations that may deactivate any of these genes?
Potentially_lost_genes_AT_list<-intersect(Variants_HOMOALT_CB_Viviparity$Annotation, list_postive_selected_genesAT_expressedinTR)
AT_CANDIDATE_GO

atoweri_genes_trophotaenia_genes <- atoweri_genes %>%
  filter(Annotation %in% list_postive_selected_genesAT_expressedinTR)

##################################  Variant analysis -- gene loss ##############################
#Input table
High_impact_variants<-read_table('Only_High_Impact_simple_vars_nomissing.tsv')

##Annotate table
High_impact_variants$Gene<-lapply(High_impact_variants$Gene, paste0, ".t1")
High_impact_variants$Annotation<-annot$Preferred_name[match(High_impact_variants$Gene, annot$query_name)]
High_impact_variants$Description<-annot$X.4[match(High_impact_variants$Gene, annot$query_name)]

#Now get all that are homozygous alt for CB that are homozygous ref for all other species
#These are genes that 'activated' in Goodeids but were previously non-coding loci.


Variants_HOMOALT_CB_Viviparity<-High_impact_variants %>%
  filter(CB == '0/0', GM == '0/0',XR == '0/0', XC == '0/0', AS == '0/0',
         AT == '1/1', IF == '0/0', CL == '0/0', GA == '0/0')

length(unique(Variants_HOMOALT_CB_Viviparity$Gene))
intersect(Variants_HOMOALT_CB_Viviparity$Gene, ABSREL_results$gene)

gene_HIGHIMPACTVARs.df <- bitr(Variants_HOMOALT_CB_Viviparity$Annotation, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                     OrgDb = organism)

ego_HI <- enrichGO(gene          = gene_HIGHIMPACTVARs.df$ENTREZID,
                  OrgDb         = organism,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.50)

Gene_ontology_highimpactvarsCB<-data.frame(ego_HI)


# Find genes in A.toweri overlap analysis
Variants_Loss_genes_AT <- Variants_HOMOALT_CB_Viviparity %>%
  filter(Annotation %in% Potentially_lost_genes_AT_list)




#Tranpose data for everything against CB
Variants_HOMOALT_CB_Viviparity_long<-pivot_longer(Variants_HOMOALT_CB_Viviparity,
                                        cols=GM:GA,
                                        names_to="Species",
                                        values_to='Genotypes')


ggplot(Variants_HOMOALT_CB_Viviparity_long, aes(x=Type, fill=Type))+
  geom_bar(width=0.5)+theme_bw()+
  labs(x='Genotypes',
       y='Count (High-impact mutations specific to Goodeinae)',
       fill='Types of mutation')+
  coord_flip()




#How important are what changes?

#Split Type into separate types so that it can be easily plotted.
High_impact_variants<-High_impact_variants %>% 
  mutate(Type = strsplit(as.character(Type), "&")) %>% 
  unnest(Type)

#Tranpose the dataset
High_impact_variants_long<-pivot_longer(High_impact_variants,
               cols=GM:GA,
               names_to="Species",
               values_to='Genotypes')

# Now count the number of 'types' of mutations. 
High_impact_variants_summary<-High_impact_variants_long %>%
  group_by(Gene) %>%
  summarise(N=n())
  
#High_impact_variants_summary$Annotation<-annot$Preferred_name[match(High_impact_variants_summary$Gene, annot$query_name)]


ggplot(High_impact_variants_long, aes(x=Genotypes, fill=Genotypes))+
  geom_bar(width=0.5)+theme_bw()+
  labs(x='Genotypes',
       y='Count',
       fill='Types of mutation')+
  facet_wrap(~Species)+
  coord_flip()


dotplot(ego_HI, showCategory=30)


Variants_HOMOALT_CB_Viviparity_Selection<-subset(Variants_HOMOALT_CB_Viviparity, Gene %in% ABSREL_results$gene)
Variants_HOMOALT_CB_Viviparity_Selection_Intogression_AmeAta<-subset(Variants_HOMOALT_CB_Viviparity_Selection, Gene %in% Gir_Ata_XR_Introgression_GR_OverlapGenes$tx_name)
Variants_HOMOALT_CB_Viviparity_Selection_Intogression_AmeGir<-subset(Variants_HOMOALT_CB_Viviparity_Selection, Gene %in% Gir_Ata_Good_Introgression_GR_OverlapGenes$tx_name)
Variants_HOMOALT_CB_Viviparity_Selection_Intogression_GoodGir<-subset(Variants_HOMOALT_CB_Viviparity_Selection, Gene %in% Gir_Ata_Ilyo_Introgression_GR_OverlapGenes$tx_name)
Variants_HOMOALT_CB_Viviparity_Selection_Intogression_GoodAta<-subset(Variants_HOMOALT_CB_Viviparity_Selection, Gene %in% Gir_Ata_Ame_Introgressio_GR_OverlapGenes$tx_name)
Variants_HOMOALT_CB_Viviparity_Selection_Intogression_XenoGir<-subset(Variants_HOMOALT_CB_Viviparity_Selection, Gene %in% Gir_Ata_Char_Introgression_GR_OverlapGenes$tx_name)
Variants_HOMOALT_CB_Viviparity_Selection_Intogression_XenoAta<-subset(Variants_HOMOALT_CB_Viviparity_Selection, Gene %in% Gir_Ata_Char_Introgression_GR_OverlapGenes$tx_name)








# 
# ############### Estimates of relaxed selection ######
# RELAX_dataset<-fread('RELAX_survey_output_allgenes.txt.gz', sep=',', header = T)
# 
# RELAX_dataset <- RELAX_dataset %>%
#   filter(corrected_pval < 0.05) %>%
#   mutate(Selection_intensity = case_when(K > 1 ~ "Intensified", 
#                                          K < 1 ~ "Relaxed"))
# 
# means_RELAX_dataset<-RELAX_dataset %>%
#   group_by(Branch) %>%
#   mutate('Branch_K_mean'= mean(K),
#          'Branch_K_median'= median(K))
# 
# means_RELAX_dataset$K[means_RELAX_dataset$K == 0] <- 0.01
# 
# ggplot(RELAX_dataset, aes(x=Selection_intensity))+
#   geom_bar(aes(fill = Branch))+
#   facet_wrap(~Branch,  scales = "free")
# 
# ggplot(means_RELAX_dataset, aes(x=K))+
#   geom_freqpoly(bins=50, aes(colour = Branch))+
#   facet_wrap(~Branch,  scales = "free_y")+theme_bw()+
#   geom_vline(xintercept = 1, linetype='dashed')+
#   geom_vline(data = means_RELAX_dataset, mapping = aes(xintercept = Branch_K_mean), colour='red', linetype='dashed')+
#   #geom_vline(data = means_RELAX_dataset, mapping = aes(xintercept = Branch_K_median), colour='blue')+ #median
#   scale_x_continuous(trans='log',labels = label_number(accuracy = 1))+
#   labs(x="K (selection intensity parameter)")+
#   theme(legend.position = "none")
# 
# 
# RELAX_dataset$Gene<-lapply(RELAX_dataset$Gene, paste0, ".t1")
# RELAX_dataset$Annotation<-annot$Preferred_name[match(RELAX_dataset$Gene, annot$query_name)]
# RELAX_dataset$Description<-annot$X.4[match(RELAX_dataset$Gene, annot$query_name)]
# 

