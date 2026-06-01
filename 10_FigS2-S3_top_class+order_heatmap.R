library(pheatmap)
library(tidyverse)
library(magrittr)
library(forcats)
library(vegan)
# Load data ---------------------------------------------------------------
asv <- read.delim("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/arc_unrarefild/ASV_Arc_60cm.txt",
                  header = TRUE)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/env_manage_220812.csv", 
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/tax.txt",
                  header = TRUE)
# clean raw data ----------------------------------------------------------
env1 <- env %>%
  filter(env$sample_id %in% colnames(asv))

tax1 <- tax %>% 
  filter(tax$asv_id %in% asv$asv_id)
# 01 Top Class  ---------------------------------------------
top_class <- asv %>% 
  pivot_longer(-asv_id, 
               names_to = 'sample_id', 
               values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Class)%>%
  summarise(ClassAbun = sum(Abun)) %>%
  mutate(RelAbun = 100*ClassAbun/sum(ClassAbun)) %>%
  ungroup()%>%
  arrange(desc(RelAbun))

topallclass <- top_class$Class[-3]
topallclass

plot_data1 <- asv %>% 
  pivot_longer(-asv_id, 
               names_to = 'sample_id', 
               values_to = 'Abun')%>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(site,
           dna_type,
           Phylum,
           Class) %>% 
  summarise(ClassAbun = sum(Abun)) %>%
  ungroup %>% 
  group_by(site,
           dna_type)%>%
  mutate(RelAbun = 100*ClassAbun/sum(ClassAbun)) %>%
  ungroup()

any_na <- any(is.na(plot_data1))

# merge site&dna_type and transfer wider data
plot_data2 <- plot_data1 %>%
  filter(Class %in% topallclass) %>%
  mutate(tmp = paste0(.$site, .$dna_type)) %>%
  pivot_wider(id_cols = Class, 
              names_from = tmp, 
              values_from = RelAbun) %>% 
  mutate_at(vars(-Class), ~replace_na(., 0))

plot_data3 <- as.matrix(plot_data2[2:9]) %>% `rownames<-`(plot_data2$Class)

row_anno <- tax %>% 
  select(Phylum, 
         Class) %>%
  distinct() %>%
  filter(Class %in% topallclass)

#row_anno[row_anno$Class=='c__NA','Phylum'] <- 'Unassigned'
row_anno %<>% distinct()
tmp2 <- data.frame(row_anno$Phylum) 
rownames(tmp2) <- row_anno$Class
row_anno <- tmp2
names(row_anno) <- 'Phylum'
unique(row_anno$Phylum)
my_colour = list(
  Phylum = RColorBrewer::brewer.pal(8, 'Set1')%>% `names<-`(unique(row_anno$Phylum))
)

drows = vegdist(plot_data3, 
                method="bray", 
                binary=FALSE, 
                diag=FALSE, 
                upper=FALSE, 
                na.rm = F)

dcols = vegdist(t(plot_data3), 
                method="bray", 
                binary=FALSE, 
                diag=FALSE, 
                upper=FALSE, 
                na.rm = F)

pdf(file = '/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/02 🌼Chile Archaea/heatmap/Top_class.pdf',
    height = 5,
    width = 10)
pheatmap(plot_data3, 
         scale='row', 
         cutree_cols = 4,
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         #clustering_distance_rows = drows,
         clustering_distance_cols = dcols,
         clustering_method = "ward.D2",
         annotation_row = row_anno, 
         annotation_colors  = my_colour
)
dev.off()

# 02 Top order  ---------------------------------------------
top_order <- asv %>% 
  pivot_longer(-asv_id, 
               names_to = 'sample_id', 
               values_to = 'Abun') %>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(Order)%>%
  summarise(OrderAbun = sum(Abun)) %>%
  mutate(RelAbun = 100*OrderAbun/sum(OrderAbun)) %>%
  ungroup()%>%
  arrange(desc(RelAbun))

topallorder <- top_order$Order[2:16]
topallorder

plot_data1 <- asv %>% 
  pivot_longer(-asv_id, 
               names_to = 'sample_id', 
               values_to = 'Abun')%>% 
  filter(Abun>0) %>%
  left_join(tax1) %>%
  left_join(env1) %>%
  group_by(site,
           dna_type,
           Phylum,
           Class,
           Order) %>% 
  summarise(OrderAbun = sum(Abun)) %>%
  ungroup %>% 
  group_by(site,
           dna_type)%>%
  mutate(RelAbun = 100*OrderAbun/sum(OrderAbun)) %>%
  ungroup()

any_na <- any(is.na(plot_data1))

# merge site&dna_type and transfer wider data
plot_data2 <- plot_data1 %>%
  filter(Order %in% topallorder) %>%
  mutate(tmp = paste0(.$site, .$dna_type)) %>%
  pivot_wider(id_cols = Order, 
              names_from = tmp, 
              values_from = RelAbun) %>% 
  mutate_at(vars(-Order), ~replace_na(., 0))
# Method1
plot_data3 <- as.matrix(plot_data2[2:9]) %>% `rownames<-`(plot_data2$Order)
# # Method2
# tmp1 <- plot_data2$Genus
# plot_data3 <- as.matrix(plot_data2[2:9])
# rownames(plot_data3) <- tmp1

row_anno <- tax %>% 
  select(Phylum, 
         Order) %>%
  distinct() %>%
  filter(Order %in% topallorder)

#row_anno[row_anno$Class=='c__NA','Phylum'] <- 'Unassigned'
row_anno %<>% distinct()

tmp2 <- data.frame(row_anno$Phylum) 

rownames(tmp2) <- row_anno$Order
row_anno <- tmp2
names(row_anno) <- 'Phylum'

my_colour = list(
  Phylum = RColorBrewer::brewer.pal(7, 'Set1')%>% `names<-`(unique(row_anno$Phylum))
)

drows = vegdist(plot_data3, 
                method="bray", 
                binary=FALSE, 
                diag=FALSE, 
                upper=FALSE, 
                na.rm = F)

dcols = vegdist(t(plot_data3), 
                method="bray", 
                binary=FALSE, 
                diag=FALSE, 
                upper=FALSE, 
                na.rm = F)

pdf(file = '/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/02 🌼Chile Archaea/heatmap/Top_order.pdf',
    height = 5,
    width = 10)
pheatmap(plot_data3, 
                  scale='row', 
                  cutree_cols = 2,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE, 
                  #clustering_distance_rows = drows,
                  clustering_distance_cols = dcols,
                  clustering_method = "ward.D2",
                  annotation_row = row_anno, 
                  annotation_colors  = my_colour
                  )
dev.off()
