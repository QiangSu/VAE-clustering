library(Rtsne)

count <- read.csv("C:/Users/Qiang/Desktop/shenzhen Uui/project/GaussF/N7_spike_splitting_output.csv", header = TRUE, fileEncoding = "UTF-8")
nrow(count)
#count_uni = unique(count) # 去除重复数据
#count_uni = na.omit(count_uni) #去除缺失值
head(count)
data = count[,2:26] # 获取性状信息
head(data)

sample = count[,1] # 获取品种信息
head(sample)
#write.table(sample,'C:/Users/Qiang/Desktop/new_pipeline/drug_Injury/tSNA_subtoxicity_sample.txt', sep = '\t', col.names = NA, quote = FALSE)

set.seed(321) # 设置随机数种子
tsne_out = Rtsne(data, check_duplicates=FALSE, dims = 2, pca = T, max_iter = 1000, theta = 0.4, perplexity = 20, verbose = F) # 进行t-SNE降维分析
head(tsne_out)
plot(tsne_out$Y)
write.table(tsne_out$Y,'C:/Users/Qiang/Desktop/shenzhen Uui/project/GaussF/tSNE_USF2_201_kmer_splitting_output.txt', sep = '\t', col.names = NA, quote = FALSE)




count <- read.csv("C:/Users/Qiang/Desktop/shenzhen Uui/project/GaussF/N7_spike_splitting_output.csv", fileEncoding = "UTF-8")
count
count_uni = unique(count) # 去除重复数据
count_uni = na.omit(count_uni) #去除缺失值
#head(count_uni)
data = count_uni[,2:26] # 获取性状信息
#head(data)
sample = count_uni[,1] # 获取品种信息
#head(sample)
#write.table(sample,'C:/Users/Qiang/Desktop/new_pipeline/drug_Injury/tSNA_subtoxicity_sample.txt', sep = '\t', col.names = NA, quote = FALSE)

library(umap)
drug.umap = umap::umap(data)
drug.umap
#head(drug.umap$layout)
plot(drug.umap$layout)
write.table(drug.umap$layout,'C:/Users/Qiang/Desktop/shenzhen Uui/project/GaussF/umap_N7_spike_splitting_output.txt', sep = '\t', col.names = NA, quote = FALSE)

