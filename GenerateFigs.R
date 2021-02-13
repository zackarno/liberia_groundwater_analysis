library(tidyverse)
library(phyloseq)
library(ggplot2)
library(broom)
library(dplyr)
library(magrittr)
library(Hmisc)
library(vegan)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(directlabels)
load('../uclust_otus/z_seqs/base_files/base_file_factory/phylo_liberia.phyloseq')#data1
source("functions/utils.R")


sample_sums(data1)
set.seed(1)
# tidy up data ------------------------------------------------------------

data_pruned<-prune_taxa(taxa_sums(data1) > 0, data1)


original_tax_table<-as(tax_table(data_pruned),"matrix") %>% data.frame()
cleaned_tax_table<-original_tax_table %>%
  rownames_to_column() %>%
  filter(Kingdom %in% c("Bacteria","Archaea"),
         Family !="mitochondria"|is.na(Family),
         Class!="Chloroplast"|is.na(Class)) %>%
  column_to_rownames("rowname") %>% as.matrix()

tax_table(data_pruned)<-cleaned_tax_table


# rarefy ------------------------------------------------------------------

data_rarefined<-rarefy_even_depth(data_pruned,rngseed = 711)
sample_sums(data_rarefined)


# Rename well column ------------------------------------------------------

sample_data(data_rarefined)$well_no <- as.character(gsub("..$","",get_variable(data_rarefined, "Z_ID")))

# mege duplicates ---------------------------------------------------------

data_samples_merged<-merge_samples(data_rarefined,"well_no")
#add refseq data to merge
data_samples_merged<-merge_phyloseq(data_samples_merged,refseq(data_rarefined))


# Pruning & Normalizing------------------------------------------------

# normalize to 100 %
data_samples_normalized = transform_sample_counts(data_samples_merged, function(x) 1E2 * x/sum(x))

abundance_threshold<-0

df_taxa_pruned_to_threshold<-prune_taxa(taxa_sums(data_samples_normalized) > 1, data_samples_transformed)

#normalize to 100 %
df_taxa_pruned_renormalized= transform_sample_counts(df_taxa_pruned_to_threshold, function(x) 1E2*x/sum(x))
sample_sums(df_taxa_pruned_renormalized)

sample_data(df_taxa_pruned_renormalized)$well_no<-as.factor(sample_data(df_taxa_pruned_renormalized)$well_no)


df_normalized_phy_glom<-tax_glom(df_taxa_pruned_renormalized,taxrank="Phylum")

(major_phylum_bp1<-plot_bar(df_normalized_phy_glom, "well_no", "Abundance","Phylum", title="Phylum (OTUS>1% Abundance)")+
    ylab("Relative Abundance")+xlab("Sample") +theme_bw()+
    theme(axis.title.x = element_text(size=16),
          axis.text.x = element_text(angle = 70, hjust = 1,size=8, colour = "black"),
          axis.title.y = element_text(size=16),
          axis.text.y = element_text(angle = 90, hjust = 0.5,size=8,colour="black"),
          plot.title = element_text(size=16, hjust=0.5,colour="black",face="bold"))+
    scale_fill_manual(values=one_percent_colors,na.value="#555555ff")+
    ggtitle("Microbial Community Structure (OTUs > 1%)"))






# filter to top 15 --------------------------------------------------------

# it becomes an interesting question on whether we are interesting in looking at the top 15 most frequent
# phyla or if we want to look at the top 15 most frequent biota by Phyla

# it seems there is only really 17 phyla when you glom after filtering by 1 % does it really make sense to choose top 15 from a list of only 17? should I make choose top 10? top 5?
#choose 1 for phyla
graph_top_15_by<-c("phylum_glom","total")[1]
if(graph_top_15_by=="total"){
top_15_names<-names(sort(taxa_sums(df_taxa_pruned_renormalized), decreasing=TRUE)[1:15])}else{
top_15_names<-names(sort(taxa_sums(df_normalized_phy_glom), decreasing=TRUE)[1:15])}

top_15<-prune_taxa(top_15_names,df_taxa_pruned_renormalized)
top_15_renormalized= transform_sample_counts(top_15, function(x) 1E2*x/sum(x))
####paper figure-  HERES THE TICKET####
major_phylum_bp1<-plot_bar(top_15_renormalized, "well_no", "Abundance","Phylum", title="Phylum (OTUS>1% Abundance)")+
  ylab("Relative Abundance")+xlab("Sample") +theme_bw()+
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_text(angle = 70, hjust = 1,size=8, colour = "black"),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(angle = 90, hjust = 0.5,size=8,colour="black"),
        legend.text=element_text(size=8),
        legend.title = element_text(size=9),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(colour = "black", fill=alpha("white",0.7)),
        plot.margin=margin(0,0,0,0,"cm"),
        plot.title = element_text(size=16, hjust=0.5,colour="black",face="bold"),
        legend.position=c(0.3,0.4))+
  scale_fill_manual(values=colv3,na.value="#555555ff")+
  ggtitle("Microbial Community Structure")
major_phylum_bp1


# setwd("C:/Users/Zack/Dropbox/Grad School/Liberia/thesis_writing/Thesis_Drafts/figures/pngs/microbes")
# png(file="09_01_StackedBar.png",units="in",width=6, height=5,res=300)
major_phylum_bp1
dev.off()



windows();major_phylum_bp1

########write abundance table
Phylumfac=factor(tax_table(top_15_renormalized)[,"Phylum"])
phytab = apply(otu_table(top_15_renormalized), MARGIN = 1, function(x) {
  tapply(x, INDEX = Phylumfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
phytab
phytab2<-phytab[, order(as.integer(colnames(phytab)))]


land_use_metadata<-read.csv('../uclust_otus/z_seqs/random_forest/all_wells_new_landuse_meta.csv', stringsAsFactors = F) %>%
  rename(well_no="Z_ID") %>%
  mutate(well_no=as.factor(well_no))
lu %>% dim()


all_metadata<-sample_data(df_taxa_pruned_renormalized) %>% data.frame() %>% left_join(land_use_metadata, by=c("well_no"))

#convert to mmol for CoDa transformations
library(robCompositions)
major_ion_mmol_matrix<-convert_major_ions_to_mmol(all_metadata)

library(robCompositions)
#do roobust ilr cluster analysis on major ion data-- a few types
clst5<-clustCoDa(major_ion_mmol_matrix,k=5,method="ward.D2",vals=TRUE,scale="robust",transformation="pivotCoord")
clst4<-clustCoDa(major_ion_mmol_matrix,k=4,method="ward.D2",vals=TRUE,scale="robust",transformation="pivotCoord")
clst3<-clustCoDa(major_ion_mmol_matrix,k=3,method="ward.D2",vals=TRUE,scale="robust",transformation="pivotCoord")
clst2<-clustCoDa(major_ion_mmol_matrix,k=2,method="ward.D2",vals=TRUE,scale="robust",transformation="pivotCoord")


dist_matrix_trace<-dist(scale(log10(all_metadata[,trace_elements])),method="euclidean")
trace_res<-hclust(dist_matrix_trace, method="ward.D")
trec_res_clust4<-cutree(trace_res,4)
trec_res_clust3<-cutree(trace_res,3)

all_metadata<- all_metadata %>%
  mutate(trace_clust4=trec_res_clust4,
         trace_clust3=trec_res_clust3,
         coda_clust5=clst5$cluster,
         coda_clust4=clst4$cluster,
         coda_clust3=clst3$cluster,
         cod_clust2=clst2$cluster,
         geo_class1= ifelse(Geo_New_3==1, "Crystalline","Sedimentary"),
         geo_class2= case_when(
           Geo_New_2==1~"Igneous",
           Geo_New_2==2~"Metamorphic",
           TRUE~"Sedimentary"
         )

         )



metadata_env<-sample_data(all_metadata)
sample_names(metadata_env)<-sample_data(metadata_env)$well_no %>% as.character()
#extract phyloseq parts for new object
OTU<-otu_table(df_taxa_pruned_renormalized)
TAX<-tax_table(df_taxa_pruned_renormalized)
PHY<-phy_tree(df_taxa_pruned_renormalized)
REF<-refseq(df_taxa_pruned_renormalized)


phylo_with_all<-phyloseq(OTU,TAX,PHY,REF,metadata_env)


set.seed(1)
NMDS1<-metaMDS(data.frame(otu_table(phylo_with_all)),trace=FALSE,try=1000,trymax=1000)

windows();stressplot(MDS1)

####simple NMDS plot- extract data from NMDS####
#extract species scores
names(NMDS1)
species_scores<-as.data.frame(NMDS1$species)
species_scores$species <- rownames(species.scores)
names(species_scores)[c(1, 2)] <- c("x", "y")
species_scores$z <- NA
sample_scores<-as.data.frame(NMDS1$points)
sample_scores$sample_no<-rownames(sample_scores)
names(sample_scores)[c(1, 2)] <- c("x", "y")
sample_scores$z <- NA
#merge species scores and tax table
#by denvovo to get phylyum to symbolize in ggplot


tax_table_df=as.data.frame(phylo_with_all@tax_table@.Data)
tax_table_df$species<-rownames(tax_table_df)
species_scores2<-left_join(species_scores$species,tax_table_df, by="species")
#merge sample scores and envdata by Z_ID to get env data to symbolize points
sample_scores2<-metadata_env %>% left_join(sample_scores, by=c("well_no"="sample_no"))


###plot by land use and Kingdom ####
#### lets plot fig 39/fig 40 partitioned at Domain Level

####plot domain by landuse####
colv3<-c( "#386CB0", "#BEAED4", "#FDC086", "#d9f4c7ff","#E7298A","#6dcc29ff", "#fc913eff", "#666666", "#41841cff", "#3bc8ecff", "#7570B3",   "#fcfe42ff","#BF5B17","#333333ff", "#e70b1eff")

NMDS.domain.LC<-ggplot(data=species_scores2,aes(x=x,y=y))+geom_point(aes(color=Kingdom),size=3)+scale_color_manual(values=c("#E7298A","#BEAED4"))+
  #geom_point(data=sample.scores,aes(x=x,y=y),colour="black",shape=8,size=4)+theme_bw()+
  #geom_segment(data=vec.sig,aes(x=0,xend=x,y=0,yend=y),
  #arrow = arrow(length = unit(0.5, "cm")),colour="black",inherit_aes=FALSE) +
  #geom_text(data=vec.sig,aes(x=x,y=y-.05,label=parameter),size=3)+theme_bw()+
  labs(x="NMDS1",y="NMDS2")+geom_label_repel(data=sample_scores2,aes(x=x,y=y,label=well_no,fill=land_use.y),alpha=0.6,size=3, fontface = 'bold',color='black', point.padding = unit(0,"lines"),box.padding=unit(0,"lines"))+scale_fill_manual(name="Land Use",values=c("blue","red","green"),labels=c("Agriculture","City","Scattered Forest\n& Ag."))+theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=8, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5,size=8,colour="black"),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(colour = "grey", fill=alpha("white",0.6)),
        plot.margin=margin(0,0,0,0,"cm"),
        plot.title = element_text(size=12,hjust=0.5))+
  #legend.position=c(-.18,0.6))+
  #guides(fill=guide_legend(ncol=2))+
  ggtitle("NMDS Plot: Taxa and Samples")
# windows();NMDS.domain.LC
# setwd("C:/Users/Zack/Dropbox/Grad School/Liberia/thesis_writing/Thesis_Drafts/figures/pngs/microbes")
# png(file="09_01_NMDS_domain_LC.png",units="in",width=6.9, height=4,res=300)

# dev.off()
# d9.MDS
###double check lables#####
NMDS_landcover<-ggplot(data=species_scores2,aes(x=x,y=y))+geom_point(aes(color=Kingdom),size=3)+scale_color_manual(values=c("#E7298A","#BEAED4"))+
  labs(x="NMDS1",y="NMDS2")+
  geom_label_repel(data=sample_scores2,aes(x=x,y=y,label=well_no,fill=land_use.y),
                   alpha=0.6,size=3, fontface = 'bold',color='black', point.padding = unit(0,"lines"),box.padding=unit(0,"lines"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=8, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5,size=8,colour="black"),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(colour = "grey", fill=alpha("white",0.6)),
        plot.margin=margin(0,0,0,0,"cm"),
        plot.title = element_text(size=12,hjust=0.5))+
  #legend.position=c(-.18,0.6))+
  #guides(fill=guide_legend(ncol=2))+
  ggtitle("NMDS Plot: Taxa and Samples")
NMDS_landcover



####plot by simple geo1#####
NMDS_geo_class1<-ggplot(data=species_scores2,aes(x=x,y=y))+geom_point(aes(color=Kingdom),size=3)+scale_color_manual(values=c("#E7298A","#BEAED4"))+
  labs(x="NMDS1",y="NMDS2")+
  geom_label_repel(data=sample_scores2,aes(x=x,y=y,label=well_no,fill=geo_class1),
                   alpha=0.6,size=3, fontface = 'bold',color='black', point.padding = unit(0,"lines"),box.padding=unit(0,"lines"))+
  scale_fill_manual(name="Land Use",values=c("blue","red","green"),labels=c("Crystalline","Sedimentary"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=8, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5,size=8,colour="black"),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(colour = "grey", fill=alpha("white",0.6)),
        plot.margin=margin(0,0,0,0,"cm"),
        plot.title = element_text(size=12,hjust=0.5))+
  #legend.position=c(-.18,0.6))+
  #guides(fill=guide_legend(ncol=2))+
  ggtitle("NMDS Plot: Taxa and Samples")
# windows();ggNMDS.geo1
# setwd("C:/Users/Zack/Dropbox/Grad School/Liberia/thesis_writing/Thesis_Drafts/figures/pngs/microbes")
# png(file="09_01_NMDS_domain_Geo1.png",units="in",width=5.9, height=4,res=300)
# ggNMDS.geo1
# dev.off()

NMDS_geo_class1
###SIMPLE GEO 2 - with dikes#####
NMDS_geo_class2<-ggplot(data=species_scores2,aes(x=x,y=y))+geom_point(aes(color=Kingdom),size=3)+scale_color_manual(values=c("#E7298A","#BEAED4"))+
  #geom_point(data=sample.scores,aes(x=x,y=y),colour="black",shape=8,size=4)+theme_bw()+
  #geom_segment(data=vec.sig,aes(x=0,xend=x,y=0,yend=y),
  #arrow = arrow(length = unit(0.5, "cm")),colour="black",inherit_aes=FALSE) +
  #geom_text(data=vec.sig,aes(x=x,y=y-.05,label=parameter),size=3)+theme_bw()+
  labs(x="NMDS1",y="NMDS2")+geom_label_repel(data=sample_scores2,aes(x=x,y=y,label=well_no,fill=geo_class2),alpha=0.6,size=3, fontface = 'bold',color='black', point.padding = unit(0,"lines"),box.padding=unit(0,"lines"))+scale_fill_manual(name="Bedrock Geology",values=c("blue","red","green"),labels=c("Metamorphic","Igneous","Sedimentary"))+theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=8, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5,size=8,colour="black"),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(colour = "grey", fill=alpha("white",0.6)),
        plot.margin=margin(0,0,0,0,"cm"),
        plot.title = element_text(size=12,hjust=0.5))+
  ggtitle("NMDS Plot: Taxa and Samples")

# vector fitting ----------------------------------------------------------
envr_data_cols<-c("Ca", "pop_2km", "elev_m", "dist_river",
  "pH_f", "cond_l", "WLC", "Mg", "Na", "K", "HCO3", "SO4", "Cl",
  "NO3", "SiO2", "Al", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
  "Zn", "Sr", "Mo", "Cd", "Sb", "Ba", "Pb", "U", "As", "Se",
  "cond_f",
  "DEM_interp" ,
  "dist_coast")

metadata_merged<-sample_data(df_taxa_pruned_renormalized) %>% data.frame()
env_vars<-c("DEM_interp",
  "elev_m", "dist_coast", "dist_river",
  "temp_c", "pH_f", "cond_l", "pH_l", "cond_f",
  "SWL", "TD", "WLC",  "Pop_1km", "chloride_rain",
  "sodium_rain", "Chirps_5m_precip", "NDVI_1_yr",
   "Ca", "Mg", "Na", "K",
  "HCO3", "SO4", "Cl", "NO3", "SiO2", "Al", "V", "Cr", "Mn", "Fe",
  "Co", "Ni", "Cu", "Zn", "Sr", "Mo", "Cd", "Sb", "Ba", "Pb", "U",
  "As", "Se")

metadata_merged %>% colnames() %>% dput()
asdf<-metadata_merged %>% select(env_vars) %>% butteR::get_na_response_rates()
# metadata_merged %>% select(env_vars) %>% View()
asdf


windows();metadata_merged %>%
  select(envr_data_cols) %>%
  rownames_to_column() %>%
  pivot_longer(-rowname, names_to= "param", values_to = "vals") %>%
  ggplot()+geom_histogram(aes(x=vals)) + facet_wrap(~param)


metadata_merged %>%
  select(env_vars) %>% filter_all(any_vars(.==0))

env_data_to_vecfit<-metadata_merged %>%
  select(env_vars) %>%
  mutate_all(.funs = function(x){scale(log10(x))})


envfit_1<-envfit(NMDS1,env_data_to_vecfit,permu=1000)
envfit_1_df<-data.frame(envfit_1$vectors$arrows,r2=envfit_1$vectors$r,p=envfit_1$vectors$pvals)

envfit_sig<-envfit_1_df %>%
  rownames_to_column() %>%
  filter(p<=0.05)

# extract arrows for biplot over NMDS
envfit_sig<-envfit_sig %>%
  rename(parameter="rowname") %>%
  mutate(x=NMDS1 * sqrt(r2),
         y=NMDS2 *sqrt(r2))




# plot nmds with taxa and vectors fit -------------------------------------


ggenvfit<-ggplot(data=species_scores2,aes(x=x,y=y))+geom_point(aes(color=Phylum),size=3)+
  scale_color_manual(values=vcols,guide=guide_legend(ncol = 1))+
  theme_bw()+
  geom_segment(data=envfit_sig,aes(x=0,xend=x,y=0,yend=y),
               arrow = arrow(length = unit(0.5, "cm")),colour="black",inherit_aes=FALSE) +
  geom_text(data=envfit_sig,aes(x=x,y=y-.05,label=parameter),size=3)+
  theme_bw()+
  labs(x="NMDS1",y="NMDS2")+
  geom_label_repel(data=sample_scores2,aes(x=x,y=y,
                                           label=well_no,
                                           fill=factor(e_coli_PA)),
                   alpha=0.7,fontface="bold",color="black",size=3,point.padding = unit(0,"lines"),
                   box.padding=unit(0,"lines"))+
  scale_fill_manual(name="E. coli field\ntest result",values=c("lightblue","red"),labels=c("Absent","Present"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=8, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5,size=8,colour="black"),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(colour = "grey", fill=alpha("white",0.7)),
        plot.margin=margin(0,0,0,0,"cm"),
        plot.title = element_text(size=12,margin=margin(b=0),hjust=0.5))+
  ggtitle("NMDS Plot With Taxa, Samples, and Significant Vector Fits")

ggenvfit

NMDS1
envfit_sig
td_oridsurf<-ordisurf(NMDS1,metadata_merged$TD,plot=TRUE)#WLC

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}


td_oridsurf_xyz<-extract.xyz(td_oridsurf)

# species_scores2$z<-NA
# sample_scores$z<-NA

envfit_sig$z<-NA



envfit_sig$TD_fill_tag<-ifelse(envfit_sig$parameter== "TD","strong","weak")


depth_spline<-ggplot(data=td_oridsurf_xyz,aes(x=x,y=y,z=z))+stat_contour(aes(colour=..level..),binwidth=2)+
  geom_point(data=species_scores2,aes(x=x,y=y))+
  labs(x = "NMDS1", y = "NMDS2")+
  geom_label_repel(data=sample_scores2,
                   aes(x=x,y=y,label=well_no),
                   fontface="bold",
                   color="black",
                   size=3,
                   point.padding = unit(0,"lines"),
                   box.padding=unit(0,"lines"))+
  geom_segment(data=envfit_sig,aes(x=0,xend=x,y=0,yend=y),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE) +
  geom_label_repel(data=envfit_sig,aes(x=x*1.1,y=y*1.1,label=parameter,fill=TD_fill_tag, alpha=TD_fill_tag))+scale_alpha_discrete(range=c(0.8,0.3),guide=FALSE)+scale_fill_manual(values=c("red","white"),guide=FALSE)+
  #geom_text(data=vec.sig,aes(x=x,y=y,label=parameter),size=4)+
  scale_colour_gradient(high = "red", low = "black")+
  theme_bw()
depth_spline
SiO2_spline.z<-direct.label(depth_spline,
                            list("top.points", rot=75, cex=1))
envfit_sig$z<-NA

spline_graph<-ggplot(data=td_oridsurf_xyz,aes(x=x,y=y,z= z))+stat_contour(aes(colour=..level..),binwidth=2)+
  geom_point(data=species_scores2,aes(x=x,y=y))+labs(x = "NMDS1", y = "NMDS2")+
  geom_label_repel(data=sample_scores2$z,aes(x=x,y=y,label=well_no),fontface="bold",color="black",size=3,point.padding = unit(0,"lines"),box.padding=unit(0,"lines"))+
  geom_segment(data=envfit_sig,aes(x=0,xend=x,y=0,yend=y),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE) +
  geom_label_repel(data=envfit_sig,aes(x=x*1.1,y=y*1.1,label=parameter,fill=TD_fill_tag, alpha=TD_fill_tag))+scale_alpha_discrete(range=c(0.8,0.3),guide=FALSE)+scale_fill_manual(values=c("red","white"),guide=FALSE)+
  #geom_text(data=vec.sig,aes(x=x,y=y,label=parameter),size=4)+
  scale_colour_gradient(high = "red", low = "black")+
  theme_bw()

