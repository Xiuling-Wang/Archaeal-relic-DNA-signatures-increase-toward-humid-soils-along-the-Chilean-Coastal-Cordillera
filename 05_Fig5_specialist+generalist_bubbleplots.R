library(tidyverse)
library(vegan)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(rdacca.hp)
library(RColorBrewer)
library(labdsv)#Indicator Value
asv <- read.delim("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/file_A.txt",
                  header = TRUE,
                  row.names = 1)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/env_2024.csv", 
                row.names = 1,
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/tax.txt",
                  row.names = 1,
                  header = TRUE)
asv <- asv %>% 
  filter(Kingdom == 'Archaea')
# 01 iDNA Specialist ---------------------------------------------------------
env1 <- env %>%
  filter(rownames(.) %in% colnames(asv),
         depth < 60,
         dna_type == "iDNA") %>% 
  data.frame()

filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]

filt2 <- rowSums(asv1) > 0#删除全部是0的asv
asv2 <- asv1[filt2,]

filt3 <- colSums(asv2) > 200#delete samples < 200 reads
asv3 <- asv2[,filt3]

filt4 <- rownames(tax) %in% rownames(asv3)
tax1 <-tax[filt4,]# 

# new asv (grouped)
group_info <- data.frame(col_name = rownames(env1),
                         group = env1$site_dna_depth)

asv2.5 <- asv3 %>%
  mutate(id = rownames(.)) %>%  # Add row identifier
  pivot_longer(-id, names_to = "col_name", values_to = "value") %>%
  left_join(group_info, by = "col_name") %>%  # Join with group-info
  group_by(id, group) %>%  # Group by row number and group
  summarise(value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = "group", values_from = "value") %>% 
  as.data.frame()

rownames(asv2.5) <- asv2.5[[1]]
asv2.6 <- asv2.5[,-1]
asv2.7 <- asv2.6 [,unique(env1$site_dna_depth)]
asv2.8 <- asv2.7[match(rownames(tax1), rownames(asv2.7)),]
# trans and relative abundance🤩
asv2.8_t <- t(asv2.8)
total_reads <- rowSums(asv2.8_t)
relative_abundance <- 100*asv2.8_t / total_reads
relative_abundance_i <- t(relative_abundance) %>% 
  as.data.frame()
relative_abundance_i$mean_ra <- rowMeans(relative_abundance_i) 
#set threshold
filt_asv_i <- relative_abundance_i %>%
  filter(mean_ra > 0.1) %>% #keep over 0.001 asv
  select(-"mean_ra")
filt_asv_i_t <- t(filt_asv_i)
#create group for analyse
group_info2 <- env1 %>%
  select(site_dna_depth, site) %>%
  distinct(site_dna_depth, .keep_all = TRUE) %>% 
  as.data.frame()
rownames(group_info2) <- group_info2[,1]
group_info2 <- group_info2[,-1,drop = F]

group_info2$site <- as.factor(group_info2$site)
levels_info <- levels(group_info2$site)
#indval analyse
iva <- indval(filt_asv_i_t,group_info2$site)
# relfrq = relative frequency of species in classes
# relabu = relative abundance of species in classes
# indval =  the indicator value for each species
# maxcls = the class each species has maximum indicator value for
# indcls = the indicator value for each species to its maximum class
# pval = the probability of obtaining as high an indicator values as observed over the
# specified iterations
gr <- levels_info[iva$maxcls[iva$pval <= 0.05 & iva$indcls > 0.8]]
#gr <- iva$maxcls[iva$pval <= 0.05 & iva$indcls > 0.8]
iv <- iva$indcls[iva$pval <= 0.05 & iva$indcls > 0.8]
pv <- iva$pval[iva$pval <= 0.05 & iva$indcls > 0.8]
fr <- apply(filt_asv_i_t > 0, 2, sum)[iva$pval <= 0.05 & iva$indcls > 0.8]

indvalsummary <- data.frame(group = gr, 
                            indval = iv, 
                            pvalue = pv, 
                            freq = fr)
indvalsummary <- indvalsummary[order(indvalsummary$fr, -indvalsummary$indval, decreasing = TRUE),]
indvalsummary
tab_s1 <- indvalsummary[order(indvalsummary$group), ]

tab_s1

tax_tmp <- tax1 %>% 
  filter(rownames(.) %in% rownames(tab_s1)) %>% 
  as.data.frame()

TableS1 <- merge(tab_s1,
                 tax_tmp,
                 by = "row.names", 
                 all.x = TRUE) %>% 
  arrange(Phylum)

head(TableS1)
write.table(TableS1, 
            file = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/02 🌼Chile Archaea/😁Figure in Paper/Table/Table_S_indvalsummary.txt", 
            sep = "\t", 
            row.names = T)

filt4 <- rownames(filt_asv_i) %in% rownames(indvalsummary)
specialist_i <- filt_asv_i[filt4,]
filt5 <- rownames(tax1) %in% rownames(specialist_i)
tax2 <- tax1[filt5,]
# 为 asv-table 和 tax-table 添加行名
asv_table <- specialist_i %>%
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
unique(sorted_data$sample_id)
tmp123 <- c("AZiDNA0_5cm","AZiDNA5_10cm","AZiDNA10_20cm","AZiDNA20_40cm","AZiDNA40_60cm",
            "SGiDNA0_5cm","SGiDNA5_10cm","SGiDNA10_20cm","SGiDNA20_40cm","SGiDNA40_60cm",
            "LCiDNA0_5cm","LCiDNA5_10cm","LCiDNA10_20cm","LCiDNA20_40cm","LCiDNA40_60cm",
            "NBiDNA0_5cm","NBiDNA5_10cm","NBiDNA10_20cm","NBiDNA20_40cm","NBiDNA40_60cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZiDNA0_5cm","AZiDNA5_10cm","AZiDNA10_20cm","AZiDNA20_40cm","AZiDNA40_60cm",
                                       "SGiDNA0_5cm","SGiDNA5_10cm","SGiDNA10_20cm","SGiDNA20_40cm","SGiDNA40_60cm",
                                       "LCiDNA0_5cm","LCiDNA5_10cm","LCiDNA10_20cm","LCiDNA20_40cm","LCiDNA40_60cm",
                                       "NBiDNA0_5cm","NBiDNA5_10cm","NBiDNA10_20cm","NBiDNA20_40cm","NBiDNA40_60cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
unique(filtered_data$Phylum)
custom_color_palette <- c("p__Crenarchaeota"="#1f78b4",
                          "p__Thermoplasmatota"="#f1393b")
summary(filtered_data)
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/02 🌼Chile Archaea/IDA_Indicator species/Specialist (iDNA pool).pdf",
    width = 13,
    height = 6)
ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0, 10), 
                        breaks = c(1,5,10,15,25),
                        limits = c(0.01,26)) +
  labs(x = "Sample", y = "", title = "Specialist (iDNA pool)") +
  guides(size = guide_legend(title = "Relative Abundance (%)"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
# 02 iDNA Generalist ---------------------------------------------------------
asv_table <- asv2.8
# 计算 ASV 在样本中出现的频率
asv_freq <- asv_table %>%
  apply(1, function(x) sum(x > 0)) %>%
  as.data.frame() %>%
  rename(frequency = 1) %>%
  rownames_to_column(var = "ASV")
# 计算 ASV 在所有站点中的相对丰度
asv_relab <- asv_table %>%
  apply(1, function(x) sum(x)) %>%
  as.data.frame() %>%
  rename(relative_abundance = 1) %>%
  rownames_to_column(var = "ASV") %>%
  mutate(relative_abundance = relative_abundance / sum(relative_abundance) * 100)
# 合并两个数据框
asv_freq_relab <- inner_join(asv_freq, asv_relab, by = "ASV")
# 筛选出高度广泛分布的 ASV
generalist_asv <- asv_freq_relab %>%
  filter(frequency / ncol(asv_table) > 0.65 & relative_abundance > 0.05)

filt6 <- rownames(filt_asv_i) %in% generalist_asv$ASV
generalist_tab <- filt_asv_i[filt6,]

filt7 <- rownames(tax1) %in% rownames(generalist_tab)
tax3 <- tax1[filt7,]
# 为 asv-table 和 tax-table 添加行名
asv_table <- generalist_tab %>%
  rownames_to_column(var = "asv_id")
tax_table <- tax3 %>%
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
unique(sorted_data$sample_id)
tmp123 <- c("AZiDNA0_5cm","AZiDNA5_10cm","AZiDNA10_20cm","AZiDNA20_40cm","AZiDNA40_60cm",
            "SGiDNA0_5cm","SGiDNA5_10cm","SGiDNA10_20cm","SGiDNA20_40cm","SGiDNA40_60cm",
            "LCiDNA0_5cm","LCiDNA5_10cm","LCiDNA10_20cm","LCiDNA20_40cm","LCiDNA40_60cm",
            "NBiDNA0_5cm","NBiDNA5_10cm","NBiDNA10_20cm","NBiDNA20_40cm","NBiDNA40_60cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZiDNA0_5cm","AZiDNA5_10cm","AZiDNA10_20cm","AZiDNA20_40cm","AZiDNA40_60cm",
                                       "SGiDNA0_5cm","SGiDNA5_10cm","SGiDNA10_20cm","SGiDNA20_40cm","SGiDNA40_60cm",
                                       "LCiDNA0_5cm","LCiDNA5_10cm","LCiDNA10_20cm","LCiDNA20_40cm","LCiDNA40_60cm",
                                       "NBiDNA0_5cm","NBiDNA5_10cm","NBiDNA10_20cm","NBiDNA20_40cm","NBiDNA40_60cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
filtered_data <- sorted_data
unique(filtered_data$Phylum)

summary(filtered_data)
# 创建一个包含要阴影的样本范围的数据框（没用，使用AI吧😁）
custom_color_palette <- c("p__Crenarchaeota"="#1f78b4",
                          "p__Thermoplasmatota"="#f1393b")
# 使用自定义颜色调色板
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/02 🌼Chile Archaea/IDA_Indicator species/Generalist (iDNA pool).pdf",
    width = 13,
    height = 2)

ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), 
                        breaks = c(1,5,10,20,40),
                        limits = c(0,40)) +
  labs(x = "Sample", y = "", title = "Generalist (iDNA pool)") +
  guides(size = guide_legend(title = "Relative Abundance (%)"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()

rm(list = ls())
# 03 eDNA Specialist ---------------------------------------------------------
asv <- read.delim("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/file_A.txt",
                  header = TRUE,
                  row.names = 1)
env <- read.csv("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/env_2024.csv", 
                row.names = 1,
                header = TRUE)
tax <- read.table("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/data/02 🌼Chile Archaea/tax.txt",
                  row.names = 1,
                  header = TRUE)
asv <- asv %>% 
  filter(Kingdom == 'Archaea')
#
env1 <- env %>%
  filter(rownames(.) %in% colnames(asv),
         depth < 60,
         dna_type == "eDNA") %>% 
  data.frame()

filt1 <- colnames(asv) %in% rownames(env1)
asv1 <- asv[,filt1]

filt2 <- rowSums(asv1) > 0# clean asv
asv2 <- asv1[filt2,]

filt3 <- colSums(asv2) > 200# clean sample
asv3 <- asv2[,filt3]

filt4 <- rownames(tax) %in% rownames(asv3)
tax1 <-tax[filt4,]
# new asv (grouped)
group_info <- data.frame(col_name = rownames(env1),
                         group = env1$site_dna_depth)

asv2.5 <- asv3 %>%
  mutate(id = rownames(.)) %>%  # Add row identifier
  pivot_longer(-id, names_to = "col_name", values_to = "value") %>%
  left_join(group_info, by = "col_name") %>%  # Join with group-info
  group_by(id, group) %>%  # Group by row number and group
  summarise(value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = "group", values_from = "value") %>% 
  as.data.frame()

rownames(asv2.5) <- asv2.5[,1]
asv2.6 <- asv2.5[,-1]
asv2.7 <- asv2.6 [,unique(env1$site_dna_depth)]
asv2.8 <- asv2.7[match(rownames(tax1), rownames(asv2.7)),]
# trans and relative abundance important😅
asv2.8_t <- t(asv2.8)
total_reads <- rowSums(asv2.8_t)
relative_abundance <- 100*asv2.8_t / total_reads
relative_abundance_i <- t(relative_abundance) %>% 
  as.data.frame()
relative_abundance_i$mean_ra <- rowMeans(relative_abundance_i) 
#set threshold
filt_asv_i <- relative_abundance_i %>%
  filter(mean_ra > 0.1) %>% #keep over 0.001 asv
  select(-"mean_ra")
filt_asv_i_t <- t(filt_asv_i)
#create group for analyse
group_info2 <- env1 %>%
  select(site_dna_depth, site) %>%
  distinct(site_dna_depth, .keep_all = TRUE) %>% 
  as.data.frame()
rownames(group_info2) <- group_info2[,1]
group_info2 <- group_info2[,-1,drop = F]

group_info2$site <- as.factor(group_info2$site)
levels_info <- levels(group_info2$site)
#indval analyse
iva <- indval(filt_asv_i_t,group_info2$site)
gr <- levels_info[iva$maxcls[iva$pval <= 0.05 & iva$indcls > 0.8]]
#gr <- iva$maxcls[iva$pval <= 0.05 & iva$indcls > 0.8]
iv <- iva$indcls[iva$pval <= 0.05 & iva$indcls > 0.8]
pv <- iva$pval[iva$pval <= 0.05 & iva$indcls > 0.8]
fr <- apply(filt_asv_i_t > 0, 2, sum)[iva$pval <= 0.05 & iva$indcls > 0.8]

indvalsummary <- data.frame(group = gr, 
                            indval = iv, 
                            pvalue = pv, 
                            freq = fr)
indvalsummary <- indvalsummary[order(indvalsummary$fr, -indvalsummary$indval, decreasing = TRUE),]
indvalsummary
tab_s1 <- indvalsummary[order(indvalsummary$group), ]
tab_s1

tax_tmp <- tax1 %>% 
  filter(rownames(.) %in% rownames(tab_s1)) %>% 
  as.data.frame()

TableS1 <- merge(tab_s1,
                 tax_tmp,
                 by = "row.names", 
                 all.x = TRUE) %>% 
  arrange(Phylum)

head(TableS1)
write.table(TableS1, 
            file = "/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/02 🌼Chile Archaea/😁Figure in Paper/Table/Table_S_indvalsummary_e.txt", 
            sep = "\t", 
            row.names = T)

filt4 <- rownames(filt_asv_i) %in% rownames(indvalsummary)
specialist_i <- filt_asv_i[filt4,]
filt5 <- rownames(tax1) %in% rownames(specialist_i)
tax2 <- tax1[filt5,]
# 为 asv-table 和 tax-table 添加行名
asv_table <- specialist_i %>%
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
unique(sorted_data$sample_id)
tmp123 <- c("AZeDNA0_5cm","AZeDNA5_10cm","AZeDNA10_20cm","AZeDNA20_40cm","AZeDNA40_60cm",
            "SGeDNA0_5cm","SGeDNA5_10cm","SGeDNA10_20cm","SGeDNA20_40cm","SGeDNA40_60cm",
            "LCeDNA0_5cm","LCeDNA5_10cm","LCeDNA10_20cm","LCeDNA20_40cm","LCeDNA40_60cm",
            "NBeDNA0_5cm","NBeDNA5_10cm","NBeDNA10_20cm","NBeDNA20_40cm","NBeDNA40_60cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,levels = tmp123,
                            labels = c("AZeDNA0_5cm","AZeDNA5_10cm","AZeDNA10_20cm","AZeDNA20_40cm","AZeDNA40_60cm",
                                       "SGeDNA0_5cm","SGeDNA5_10cm","SGeDNA10_20cm","SGeDNA20_40cm","SGeDNA40_60cm",
                                       "LCeDNA0_5cm","LCeDNA5_10cm","LCeDNA10_20cm","LCeDNA20_40cm","LCeDNA40_60cm",
                                       "NBeDNA0_5cm","NBeDNA5_10cm","NBeDNA10_20cm","NBeDNA20_40cm","NBeDNA40_60cm"),
                            ordered = TRUE)
  )
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
filtered_data <- sorted_data
unique(filtered_data$Phylum)
custom_color_palette <- c("p__Crenarchaeota"="#1f78b4",
                          "p__Thermoplasmatota"="#f1393b")
summary(filtered_data)
str(filtered_data)
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/02 🌼Chile Archaea/IDA_Indicator species/Specialist (eDNA pool).pdf",
    width = 13,
    height = 11)
ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0.1, 10), 
                        breaks = c(1,5,10,25,50),
                        limits = c(0.1,50)) +
  labs(x = "Sample", y = "", title = "Specialist (eDNA pool)") +
  guides(size = guide_legend(title = "Relative Abundance (%)"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
# 04 eDNA Generalist ---------------------------------------------------------
asv_table <- asv2.8
# 计算 ASV 在样本中出现的频率
asv_freq <- asv_table %>%
  apply(1, function(x) sum(x > 0)) %>%
  as.data.frame() %>%
  rename(frequency = 1) %>%
  rownames_to_column(var = "ASV")
# 计算 ASV 在所有站点中的相对丰度
asv_relab <- asv_table %>%
  apply(1, function(x) sum(x)) %>%
  as.data.frame() %>%
  rename(relative_abundance = 1) %>%
  rownames_to_column(var = "ASV") %>%
  mutate(relative_abundance = relative_abundance / sum(relative_abundance) * 100)
# 合并两个数据框
asv_freq_relab <- inner_join(asv_freq, asv_relab, by = "ASV")
# 筛选出高度广泛分布的 ASV
generalist_asv <- asv_freq_relab %>%
  filter(frequency / ncol(asv_table) > 0.65 & relative_abundance > 0.05)

filt6 <- rownames(filt_asv_i) %in% generalist_asv$ASV
generalist_tab <- filt_asv_i[filt6,]

filt7 <- rownames(tax1) %in% rownames(generalist_tab)
tax3 <- tax1[filt7,]
# 为 asv-table 和 tax-table 添加行名
asv_table <- generalist_tab %>%
  rownames_to_column(var = "asv_id")
tax_table <- tax3 %>%
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
unique(sorted_data$sample_id)

tmp123 <- c("AZeDNA0_5cm","AZeDNA5_10cm","AZeDNA10_20cm","AZeDNA20_40cm","AZeDNA40_60cm",
            "SGeDNA0_5cm","SGeDNA5_10cm","SGeDNA10_20cm","SGeDNA20_40cm","SGeDNA40_60cm",
            "LCeDNA0_5cm","LCeDNA5_10cm","LCeDNA10_20cm","LCeDNA20_40cm","LCeDNA40_60cm",
            "NBeDNA0_5cm","NBeDNA5_10cm","NBeDNA10_20cm","NBeDNA20_40cm","NBeDNA40_60cm")
# 添加一个用于排序的变量
sorted_data <- sorted_data %>%
  mutate(order_var = match(y_label, unique(y_label)),
         sample_id = factor(sample_id,
                            levels = tmp123,
                            labels = c("AZeDNA0_5cm","AZeDNA5_10cm","AZeDNA10_20cm","AZeDNA20_40cm","AZeDNA40_60cm",
                                       "SGeDNA0_5cm","SGeDNA5_10cm","SGeDNA10_20cm","SGeDNA20_40cm","SGeDNA40_60cm",
                                       "LCeDNA0_5cm","LCeDNA5_10cm","LCeDNA10_20cm","LCeDNA20_40cm","LCeDNA40_60cm",
                                       "NBeDNA0_5cm","NBeDNA5_10cm","NBeDNA10_20cm","NBeDNA20_40cm","NBeDNA40_60cm"),
                            ordered = TRUE)
  )
#view(sorted_data)
filtered_data <- sorted_data[sorted_data$relative_abundance != 0,]
filtered_data <- sorted_data
unique(filtered_data$Phylum)

summary(filtered_data)
#view(filtered_data)
# 创建一个包含要阴影的样本范围的数据框（没用，使用AI吧😁）
custom_color_palette <- c("p__Crenarchaeota"="#1f78b4",
                          "p__Thermoplasmatota"="#f1393b")
# 使用自定义颜色调色板
pdf("/Users/xwang/Library/Mobile Documents/com~apple~CloudDocs/001 👩🏻‍🎓PhD Project/phd_es16s/reports/02 🌼Chile Archaea/IDA_Indicator species/Generalist (eDNA pool).pdf",
    width = 13,
    height = 3)

ggplot(filtered_data,
       aes(x = sample_id,
           y = order_var, 
           size = relative_abundance,
           color = Phylum)) +
  geom_point(alpha = 0.9) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # 修改vjust参数以调整x轴标签与坐标轴之间的距离
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_line(size = 0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  scale_y_continuous(labels = unique(filtered_data$y_label), 
                     breaks = unique(filtered_data$order_var), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = custom_color_palette) +
  scale_size_continuous(range = c(0, 10), 
                        breaks = c(1,5,10,25,50),
                        limits = c(0.01,50)) +
  labs(x = "Sample", y = "", title = "Generalist (eDNA pool)") +
  guides(size = guide_legend(title = "Relative Abundance (%)"),
         color = guide_legend(title = "Phylum", 
                              override.aes = list(size = 5)))
dev.off()
rm(list = ls())
