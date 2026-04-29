####Importing Microarray Data ####
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install('GEOquery')
library(GEOquery)

dir.create('GSE79768 data' , showWarnings = FALSE)

dataset <- getGEO('GSE79768' , GSEMatrix = TRUE, getGPL = FALSE, destdir = 'GSE79768 data') 

exprSet <- dataset[[1]]

geneexpr_matrix <- exprs(exprSet)

sample_metadata <- pData(exprSet)

dim(geneexpr_matrix)

head(sample_metadata)

####Preparation of Data ####
colnames(geneexpr_matrix)

colnames(sample_metadata)

unique(sample_metadata$characteristics_ch1.3)

?ifelse
?grepl

groups <-ifelse(grepl('Sinus Rhythm' , sample_metadata$characteristics_ch1.3, ignore.case = TRUE),'SR', 'AF')

geneexpr_matrix <- log2(geneexpr_matrix + 1)

table(groups)

boxplot(geneexpr_matrix, outline=FALSE, las=2, main='samples boxplot')

#### PCA ####
?prcomp

GSEPCA <- prcomp(t(geneexpr_matrix), scale. = TRUE)

plot(GSEPCA$x [,1], GSEPCA$x [,2], col = as.factor(groups),
pch = 19, xlab = 'PC1' , ylab = 'PC2' , asp = 1) 

legend('topright', legend = unique(groups), col = 1:length(unique(groups)), pch = 19)

#### CASE and CONTROL ####
groups <- factor(groups , levels = c('SR', 'AF'))
table(groups)

#### mean calculation ####
?rowMeans

case_mean <- rowMeans(geneexpr_matrix[, groups == 'AF'])
control_mean <- rowMeans(geneexpr_matrix[, groups == 'SR'])

#### log2FC ####
log2FC <- case_mean - control_mean

##### p-value calculation ####
pvalue <- apply(geneexpr_matrix, 1, function(x) 
  {t.test(x[groups == 'AF'] , x[groups == 'SR']) $p.value }) 

####fdr####
FDR <- p.adjust(pvalue, method = 'fdr')

####Stats Table ####
?data.frame
statstable <- data.frame(gene = rownames(geneexpr_matrix), mean.case = case_mean,
  mean.control = control_mean , log2FC = log2FC , p.value = pvalue , FDR = FDR)

#### significance ####
statstable <- statstable[order(statstable$FDR), ]
head (statstable)

#### probe ids > gene names ####
annonfile <- getGEO('GPL570' , AnnotGPL = TRUE) 

annontable <- Table(annonfile)
head(annontable) 
colnames(annontable)

statstable$probe.id <- rownames(statstable)

statstable.final <- merge(statstable , annontable [, c('ID', 'Gene symbol', 'Gene title')] , 
                          by.x = 'probe.id' , by.y= 'ID' , all.x = TRUE)

#### sort by FDR ####
statstable.finalfdr <- statstable.final [order(statstable.final$FDR),]
head(statstable.finalfdr)

significant.genes <- statstable.finalfdr[statstable.finalfdr$FDR <0.05 ,]

head(significant.genes, 20)

dim(significant.genes)

write.csv(significant.genes, 'C:/Users/user/Downloads/12107951 (1)/Project/DEA_results/DEA_resultssignificantgenes.csv', row.names = FALSE)
write.csv(statstable.finalfdr,'C:/Users/user/Downloads/12107951 (1)/Project/DEA_results/Statstablefinalfdr.csv', row.names = FALSE )

####save top20 significant ####

t20 <- head(statstable.finalfdr , 20)

write.csv(t20, 'C:/Users/user/Downloads/12107951 (1)/Project/DEA_Results/Top20sigDEA', row.names = FALSE)

#### Senescence dataset ####

install.packages("readxl")
library('readxl')

senescence_genes <- read_excel(file.choose('All Senescent Proteins from IMR90 with Alias Names.xlsx'))
senescence_genes <- senescence_genes[[1]]
head(senescence_genes)

senescence_genes <- unlist(strsplit(senescence_genes, ','))



#### cleaning gene lists ####
?toupper
?trimws
senescence_genes <- toupper(trimws(senescence_genes))
senescence_genes <- senescence_genes[!is.na(senescence_genes) & senescence_genes != ""]
senescence_genes <- unique(senescence_genes)

####Significant DE gene selection ####
dea.genes <- statstable.finalfdr$'Gene symbol'[statstable.finalfdr$FDR < 0.05]
dea.genes <- toupper(trimws(dea.genes))
dea.genes <- dea.genes[!is.na(dea.genes) & dea.genes != ""]
dea.genes <- unique(dea.genes)

#### Overlap analysis ####

Overlap_genes <- intersect(dea.genes , senescence_genes)

length(Overlap_genes)
head(Overlap_genes, 20)
summary(Overlap_genes)

statstable.finalfdr$gene_symbol_clean <- toupper(trimws(statstable.finalfdr$`Gene symbol`))
statstable.finalfdr$gene_symbol_clean[is.na(statstable.finalfdr$gene_symbol_clean)] <- ""

dea.overlap <- statstable.finalfdr[
  statstable.finalfdr$gene_symbol_clean %in% Overlap_genes,
]

length(Overlap_genes)
nrow(dea.overlap)
length(unique(dea.overlap$gene_symbol_clean))
#### Overlap Analysis Results ####
write.csv(Overlap_genes, 'C:/Users/user/Downloads/12107951 (1)/Project/resultsSen.DEA.overlap.csv' , row.names = FALSE)

top20_DEA.SEN_genes <- head(dea.overlap , 20)

write.csv(top20_DEA.SEN_genes, 'C:/Users/user/Downloads/12107951 (1)/Project/results/top20_DEA.SEN_genes' , 
          row.names = FALSE)

#### top 20 table new #### 
top20_clean <- statstable.finalfdr[order(statstable.finalfdr$FDR),]
top20_clean <- top20_clean[!duplicated(toupper(trimws(top20_clean$'Gene symbol'))),]

top20_clean <- top20_clean[!is.na(top20_clean$'Gene symbol') & top20_clean$'Gene symbol' != '',] 

top20_clean <- head(top20_clean , 20)

write.csv(top20_clean,'Top20_DEAres_CLean',  row.names = FALSE)


####  setting up volcano plot ####
install.packages("ggplot2")
install.packages("dplyr")
install.packages("ggrepel")

library(ggplot2) 
library(ggrepel)
library(dplyr)

DEAresults <- statstable.finalfdr

p_cutoff <- 0.05
log2FC_cutoff <- 0.1


DEAresults$genesymbols <- DEAresults$'Gene symbol'

DEAresults$negLogFDR <- -log10(DEAresults$FDR)


senescence_genes <- toupper(trimws(senescence_genes))
DEAresults$gene_symbol <- toupper(trimws(DEAresults$'Gene symbol' ))

#### volcano plot gene categorization ####
DEAresults$category <- 'Not significant'

DEAresults$category[DEAresults$FDR < p_cutoff & DEAresults$log2FC > log2FC_cutoff
] <- 'UPREGULATED'

DEAresults$category[DEAresults$FDR < p_cutoff & DEAresults$log2FC < -log2FC_cutoff
] <- 'DOWNREGULATED'

DEAresults$category[DEAresults$FDR <p_cutoff & abs(DEAresults$log2FC) > log2FC_cutoff
& DEAresults$gene_symbol %in% senescence_genes
] <- 'SENESCENT'

DEAresults$category <- factor(DEAresults$category , 
levels = c('Not significant', 'UPREGULATED' , 'DOWNREGULATED' , 'SENESCENT') )

T20label <- DEAresults %>% 
filter(FDR < p_cutoff & abs(log2FC) > log2FC_cutoff) %>%
arrange(FDR) %>% 
slice_head(n=20)

#### volcano plot ####
Volcano_plot <- ggplot(DEAresults, aes(x = log2FC, y= negLogFDR, colour = category))+
geom_point(alpha= 0.75, size = 2)+
scale_color_manual(values = c('Not significant' = 'grey',
'UPREGULATED' = 'red', 'DOWNREGULATED' = 'blue' , 'SENESCENT' = 'darkgreen'))+
geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = 'dashed')+
geom_hline(yintercept = -log10(p_cutoff), linetype = 'dashed')+
geom_text_repel(data = T20label, aes(label = gene_symbol),
size = 3, max.overlaps = 25, box.padding =0.3, point.padding = 0.2 )+
theme_minimal(base_size = 12)+
labs(title = 'DEA_Sen Volcano Plot', x='Log2FC', y= '-log10(FDR)', color= 'Category')
print(Volcano_plot)


##### heatmap #####

install.packages("pheatmap")
library(pheatmap)
library(dplyr)

top20_DEA <- head(top20_clean, 20)

top20_probes <- top20_DEA$probe.id


heatmap.data <- geneexpr_matrix[top20_probes,]

heatmap_scaled_data <- t(scale(t(heatmap.data)))

annotation_columns <- data.frame(Group = groups)
rownames(annotation_columns) <- colnames(heatmap_scaled_data)

rownames(heatmap_scaled_data) <-top20_clean$'Gene symbol'


pheatmap(heatmap_scaled_data , annotation_col = annotation_columns, show_rownames = TRUE, fontsize_row = 4,
cluster_rows = TRUE, cluster_cols = TRUE,  show_colnames = FALSE, main = 'Top 20 DEA genes ')


####Violin plots ####
install.packages("ggstatsplot")
library(ggstatsplot)
 
genes <- geneexpr_matrix["1565483_at", ]
 
 vbplot <- data.frame(
group = factor(groups), 
express_value = as.numeric(genes)
)
 
 ggbetweenstats(
data = vbplot,
x = group,
y = express_value,
type = "parametric",
pairwise.comparisons = FALSE,
title = "EGFR expression, AF vs SR", xlab = "group",  
  ylab = "Expression Level"
)
 
genes2 <- geneexpr_matrix['"1555786_s_at"']

vbplot <- data.frame(
group = factor(groups),
express_value = as.numeric(genes) 
)

ggbetweenstats(data = vbplot,
x = group,
y = express_value,
type = "parametric", 
pairwise.comparisons = FALSE,
title = "CRLS1 expression, AF vs SR",  
xlab = "group",  
ylab = "Expression Level" )


####pathway analysis ####
install.packages("gprofiler2")
library(gprofiler2)
library(dplyr)
library(ggplot2)

dea_genes <- DEAresults %>%
filter(!is.na(`Gene symbol`), `Gene symbol`!= '')

sig_genes_list <- DEAresults$'Gene symbol'[DEAresults$FDR <0.05]
sig_genes_list <- toupper(trimws(sig_genes_list))
sig_genes_list <- sig_genes_list[!is.na(sig_genes_list) & sig_genes_list != '']
sig_genes_list <- unique(sig_genes_list)
sig_genes_list <-head(sig_genes_list,500)


length(sig_genes_list)
head(sig_genes_list, 20)

?gost


Path_res <- gost(query = sig_genes_list , organism = 'hsapiens' , sources = c('GO:BP')
, significant = TRUE, correction_method = 'fdr')    


write.csv(Path_res$result, 'gprofiler_top500_GOenrichment.csv',
row.names = FALSE)
head(Path_res$result)
Path_res_clean <- Path_res$result[, c(
'term_name' , 'p_value', 'intersection_size', 'source')]

write.csv(Path_res_clean, 'gprofiler_top500_Goenrichment.csv',
row.names = FALSE)


#### Go result bar plot Sorting ####
library(ggplot2)

Go_results <- Path_res_clean

top_terms <- Go_results[order(Go_results$p_value), ]

top_terms <- top_terms[!grepl(
  "regulation|metabolic process|biological process",
  top_terms$term_name
), ]

top_terms <- head(top_terms, 10)

top_terms$negLogP <- -log10(top_terms$p_value)

#### GO result bar plot plotting ####

?reorder

ggplot(top_terms,
aes(x = reorder(term_name, negLogP),
y = negLogP)) +
geom_bar(stat = "identity", fill = "darkblue") +
coord_flip() +
theme_minimal() +
labs(title = "Top GO Biological Pathways",
x = "GO Term",
y = "-log10(p-value)"
)

#### corelation ####

egfr_probeid <- c('1565483_at' , '1565484_x_at')

egfr_exprdata <- geneexpr_matrix[egfr_probeid,]
head(egfr_probeid)
head(egfr_exprdata)

egfr_cor <- cor(t(egfr_exprdata) , method = 'pearson',)

egfr_corvalue <- cor(t(egfr_exprdata)) [1,2] 
head(egfr_corvalue)

#### top 20 table new #### 
top20_clean <- statstable.finalfdr[order(statstable.finalfdr$FDR),]
top20_clean <- top20_clean[!duplicated(toupper(trimws(top20_clean$'Gene symbol'))),]

top20_clean <-  head(top20_clean , 20)

top20_clean <- top20_clean[!is.na(top20_clean$'Gene symbol') & top20_clean$'Gene symbol' != '',] 

top20_clean <- head(top20_clean , 20)

write.csv(top20_clean,'Top20_DEAres_CLean',  row.names = FALSE)