library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(rdacca.hp)
library(RColorBrewer)
# load data ---------------------------------------------------------------
asv <- read.delim("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/file_A.txt",
                  header = TRUE,
                  row.names = 1)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/env_manage_220812.csv", 
                row.names = 1,
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/tax.txt",
                  row.names = 1,
                  header = TRUE)
# 01. top ASV for all sample-------------------------------------------------------------
asv <- asv %>% 
  filter(Kingdom == 'Archaea')

env1 <- env %>%
  filter(rownames(.) %in% colnames(asv),
         depth < 60) %>% 
  data.frame()

filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]

filt2 <- rowSums(asv1) > 0#删除全部是0的asv
asv2 <- asv1[filt2,]

filt3 <- colSums(asv2) > 200#delete samples < 200 reads
asv3 <- asv2[,filt3]

filt4 <- rownames(tax) %in% rownames(asv3)
tax1 <-tax[filt4,]# 

env2 <- env1 %>% 
  filter(row.names(.) %in% colnames(asv3)) %>% 
  mutate(group_col = paste(site_dna,
                           depth_label,
                           sep = "_"))# 可以发现有25个样本被删除，因为reads 达不到200


  group_info <- data.frame(col_name = rownames(env2),
                         group = env2$group_col)

asv4 <- asv3 %>%
  mutate(id = rownames(.)) %>%  # Add row identifier
  pivot_longer(-id, names_to = "col_name", values_to = "value") %>%
  left_join(group_info, 
            by = "col_name") %>% 
  group_by(id, group) %>%  # Group by row number and group
  summarise(value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = "group", values_from = "value") %>%
  as.data.frame()

rownames(asv4) <- asv4[,1]
asv5 <- asv4[,-1]

tasv <- t(asv5)
re_t <- 100*tasv/rowSums(tasv)
asv6 <- t(re_t)

summary(asv6)
colSums(asv6)#check if all sum = 100

asv7 <- asv6 [,unique(env2$group_col)]
asv8 <- asv7[match(rownames(tax1), rownames(asv7)),] %>% 
  as.data.frame()

rankRA <- asv8 %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance)) 
  
top_50_asv <- rankRA %>%
  slice_head(n = 50) %>% 
  select(-mean_abundance)

filt5 <- rownames(tax1) %in% rownames(top_50_asv)
tax2 <- tax1[filt5,]

# 为 asv-table 和 tax-table 添加行名
asv_table <- top_50_asv %>%
  rownames_to_column(var = "asv_id")
tax_table <- tax2 %>%
  rownames_to_column(var = "asv_id")

# 合并 ASV table 和 tax table
merged_data <- left_join(asv_table, tax_table, by = "asv_id")
# 将数据转换为长格式
long_data <- merged_data %>%
  pivot_longer(
    cols = -c(asv_id, Kingdom, Phylum, Class, Order, Family, Genus),
    names_to = "sample_id",
    values_to = "relative_abundance"
  ) %>%
  mutate(asv_anno = paste(Class,Order,Family,Genus,asv_id,sep = "_")
  )
# 按照门（Phylum）、目（Order）和 asv-id 排序数据
sorted_data <- long_data %>%
  arrange(Phylum,Class,Order,Family,Genus,asv_id)
# 生成包含 asv-id 和目（Order）信息的新纵坐标
sorted_data <- sorted_data %>%
  mutate(y_label = asv_anno)
unique(long_data$sample_id)

tmp_s <- c("AZ_iDNA_0_5cm","AZ_iDNA_5_10cm","AZ_iDNA_10_20cm","AZ_iDNA_20_40cm","AZ_iDNA_40_60cm",
           "AZ_eDNA_0_5cm","AZ_eDNA_5_10cm","AZ_eDNA_10_20cm","AZ_eDNA_20_40cm","AZ_eDNA_40_60cm",
           "SG_iDNA_0_5cm","SG_iDNA_5_10cm","SG_iDNA_10_20cm","SG_iDNA_20_40cm","SG_iDNA_40_60cm",
           "SG_eDNA_0_5cm","SG_eDNA_5_10cm","SG_eDNA_10_20cm","SG_eDNA_20_40cm","SG_eDNA_40_60cm",
           "LC_iDNA_0_5cm","LC_iDNA_5_10cm","LC_iDNA_10_20cm","LC_iDNA_20_40cm","LC_iDNA_40_60cm",
           "LC_eDNA_0_5cm","LC_eDNA_5_10cm","LC_eDNA_10_20cm","LC_eDNA_20_40cm","LC_eDNA_40_60cm",
           "NB_iDNA_0_5cm","NB_iDNA_5_10cm","NB_iDNA_10_20cm","NB_iDNA_20_40cm","NB_iDNA_40_60cm",
           "NB_eDNA_0_5cm","NB_eDNA_5_10cm","NB_eDNA_10_20cm","NB_eDNA_20_40cm","NB_eDNA_40_60cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,
                            levels = tmp_s,
                            labels = c("AZ_iDNA_0_5cm","AZ_iDNA_5_10cm","AZ_iDNA_10_20cm","AZ_iDNA_20_40cm","AZ_iDNA_40_60cm",
                                       "AZ_eDNA_0_5cm","AZ_eDNA_5_10cm","AZ_eDNA_10_20cm","AZ_eDNA_20_40cm","AZ_eDNA_40_60cm",
                                       "SG_iDNA_0_5cm","SG_iDNA_5_10cm","SG_iDNA_10_20cm","SG_iDNA_20_40cm","SG_iDNA_40_60cm",
                                       "SG_eDNA_0_5cm","SG_eDNA_5_10cm","SG_eDNA_10_20cm","SG_eDNA_20_40cm","SG_eDNA_40_60cm",
                                       "LC_iDNA_0_5cm","LC_iDNA_5_10cm","LC_iDNA_10_20cm","LC_iDNA_20_40cm","LC_iDNA_40_60cm",
                                       "LC_eDNA_0_5cm","LC_eDNA_5_10cm","LC_eDNA_10_20cm","LC_eDNA_20_40cm","LC_eDNA_40_60cm",
                                       "NB_iDNA_0_5cm","NB_iDNA_5_10cm","NB_iDNA_10_20cm","NB_iDNA_20_40cm","NB_iDNA_40_60cm",
                                       "NB_eDNA_0_5cm","NB_eDNA_5_10cm","NB_eDNA_10_20cm","NB_eDNA_20_40cm","NB_eDNA_40_60cm"),
                            ordered = TRUE))
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
summary(filtered_data)# relative_abundance Max.   :57.16146
# 创建一个包含要阴影的样本范围的数据框（没用，使用Illustrator吧😁）
# samples-to-highlight <- c("AZiDNA0-5cm","AZiDNA5-10cm","AZiDNA10-20cm","AZiDNA20-40cm",
#                             "LCiDNA0-5cm","LCiDNA10-20cm","LCiDNA20-40cm","LCiDNA5-10cm")
#color_palette <- brewer.pal(n = length(unique(long_data$Phylum)), name = "Set2")
# 创建命名颜色向量
unique(filtered_data$Phylum)
custom_color_palette <- c("p__Thermoplasmatota"="#c82423",
                          "p__Crenarchaeota"="#2878B5")

pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/02 🌼Chile Archaea/bubble plot/top50asv.pdf",
    width = 18,
    height = 10)

ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", 
                                  size = 0.5),
        panel.grid.major = element_line(size = 0.1, 
                                        linetype = "solid", 
                                        color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), 
                        breaks = c(1,10,20,40,60),
                        limits = c(0,60)) +
  labs(x = "Sample", y = "", title = "top 50 ASV") +
  guides(size = guide_legend(title = "Relative Abundance"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
