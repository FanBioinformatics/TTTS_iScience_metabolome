##################### #####################
########## ##########
##### TTTS twins 
#####
#####
#####
library(readxl);
options(stringsAsFactors = F);

## data after clean up
meta<-read_xlsx('2022 03 metabolomics_all data.xlsx', sheet = 2);
meta<-as.data.frame(meta);  ## 416 rows
meta.rnames<-paste(meta$`Compound Name`, meta$`RT <=1.2, ?, deID'ed by Q, replaced by RP data`, sep='_');

# meta<-meta[which(meta$`Compound Name`!='.'), ]; ## 233
meta.mat<-meta[, 15:46]; ## remove T2_52D and T7_121DR
# pos<-which(colnames(meta.mat)%in%c('T2_52D', 'T7_121DR'));
# meta.mat<-meta.mat[, -pos];

samples<-read.table('sampleinfor.txt', as.is = T, sep='\t', header=T);
cnames<-as.character(unlist(lapply(sapply(colnames(meta.mat), strsplit, split='_'), function(x)x[1])));
colnames(meta.mat)<-cnames;
rownames(samples)<-samples$ID;
colnames(meta.mat)<-cnames;

meta.mat<-meta.mat[, samples$ID];
rownames(meta.mat)<-meta.rnames;

### check for missing values
meta.mat<-as.matrix(meta.mat);
null.pos<-apply(meta.mat, 1, function(x){
  length(which(is.na(x)))/ncol(meta.mat);
})

## remove items with more than 50% NAs first
# pos<-which(null.pos>=.3);
# meta.mat<-meta.mat[-pos, ];

library(EnhancedVolcano); library(limma); library(Biobase);
eSet <- ExpressionSet(assayData = as.matrix(apply(meta.mat+1, c(1, 2), log2)), 
                      phenoData = as(samples, "AnnotatedDataFrame"));
design.trt=model.matrix(~status, data = samples);
fit <- lmFit(eSet, design.trt);
fit <- eBayes(fit)
df <- topTable(fit, coef=2, n=2000);

meta$logFC<-log2(meta$FC);
hh1<-EnhancedVolcano(meta, pCutoff=0.001,
                    lab = meta$`Compound Name`, FCcutoff=2,
                    x = 'logFC', title ='', subtitle = '',
                    y = 'P all T vs Ctr');

hh2<-EnhancedVolcano(df, pCutoff=0.001,
                     lab = rownames(df), FCcutoff=2,
                     x = 'logFC', title ='', subtitle = '',
                     y = 'adj.P.Val')+theme(aspect.ratio=1)
library(cowplot)
plot_grid(hh2, '', '', '', ncol=2);
dev.copy2pdf(file='0_differential vocanoplot_all case vs. ctrl.pdf', family='ArialMT');

########## top differential metabolites
library(ComplexHeatmap); library(gplots); library(viridis);

diffs<-meta.mat[rownames(df)[which(abs(df$logFC)>=2&df$adj.P.Val<0.001)], ]; ## 28
write.table(df[rownames(diffs), ], file='0_allT vs. normal.txt', 
                                   col.names = T, row.names = T, sep='\t', quote = F);

diffs.z<-t(apply(log2(diffs+1), 1, function(x){
  (x-mean(x[1:10], na.rm=T))/sd(x, na.rm = T)
}))

### not showing names with @
tmp.name<-diffs.z; p1<-grep('@', rownames(diffs.z)); tmp.name<-tmp.name[-p1, ]; rownames(tmp.name)<-gsub('^._', '', rownames(tmp.name));

Heatmap(tmp.name, col = plasma(1024), 
row_split = meta[which(rownames(meta.mat)%in%rownames(diffs.z)[-p1]), 'SubPathway'],
  row_title_rot = 0, column_split = samples[colnames(diffs.z), 'status'],
top_annotation=HeatmapAnnotation(Subtype=samples[colnames(diffs.z), 'status'], 
  col = list(Subtype=c('1Control'='black', '2Case'='blue'))));

dev.copy2pdf(file='0_heatmap of top differential metabolites.pdf', family='ArialMT');

#### pca using case-control differential metabolites
library(factoextra);
diffs.z[is.na(diffs.z)]<-0;
res.pca <- prcomp(t(diffs.z), scale. = F);
t1<-fviz_eig(res.pca);

t2<-fviz_pca_ind(res.pca, repel = TRUE,
                 col.ind=rep(c('Control', 'Case'), times=c(10, 22)))
par(mfrow=c(2, 2), pty='s')
# plot_grid(t1, t2, '', '', ncol=2);
plot_grid(t1, t2, '', '', ncol=2);

library(car);
print('')

print(t3<-dataEllipse(res.pca$x[, 1], res.pca$x[, 2], 
            groups = as.factor(samples[rownames(res.pca$x), 'status']),
            group.labels = c('Contrl', 'Case'),
            level=.95, fill=TRUE, fill.alpha=0.1))

dev.copy2pdf(file='0_PCA_differential metabolites from all T to C.pdf', family='ArialMT');

#### #### #### #### #### #### #### 
#### pca using the most variable metabolites to see sub-clusters within cases
meta.sd<-apply(meta.mat[, -c(1:10)], 1, sd, na.rm=T);
meta50<-names(sort(meta.sd, decreasing = T)[1:50]);
library(ellipse);
meta50.mat<-meta.mat[meta50, ];
diffs.z<-t(apply(log2(meta50.mat+1), 1, function(x){
  (x-mean(x, na.rm=T))/sd(x, na.rm = T)
}))
diffs.z[is.na(diffs.z)]<-0;
res.pca2 <- prcomp(t(diffs.z), scale. = F);
tt1<-fviz_eig(res.pca2);

tt2<-fviz_pca_ind(res.pca2, repel = TRUE,
                 col.ind=rep(c('Control', 'Case'), times=c(10, 22)))

# plot_grid(t1, t2, '', '', ncol=2);
par(mfrow=c(2, 2), pty='s');
plot_grid(tt1, tt2, '', '', ncol=2);
print('')

################## cluster 1 vs cluster 2 in cases
cl1<-paste('T', c(3, 4, 6, 13, 14, 15, 16), sep='');
cl2<-setdiff(colnames(meta.mat)[grep('T', colnames(meta.mat))], c(cl1, 'T22'));

samples$cluster<-samples$status;
samples[cl1, 'cluster']<-'2Cluster1';
samples[cl2, 'cluster']<-'2Cluster2';
samples$tmp<-samples$cluster; samples$tmp[32]<-'1Control';

library(cowplot)
plot_grid(t1, t2, tt1, tt2, ncol=2);
dev.copy2pdf(file='0_PCA using differential genes.pdf', family='ArialMT');

print(tt3<-dataEllipse(res.pca2$x[, 1], res.pca2$x[, 2], 
                 groups = as.factor(samples[rownames(res.pca2$x), 'tmp']),
                 group.labels = c('Contrl', 'Cluster1', 'Cluster2'),
                 level=.95, fill=TRUE, fill.alpha=0.1))
dev.copy2pdf(file='0_ellipse pca overlay.pdf', family='ArialMT');

### need to get rid of some NA to remove compounds 
eSet <- ExpressionSet(assayData = as.matrix(apply(meta.mat[, c(cl1, cl2)]+1, c(1, 2), log2)))
design.trt=model.matrix(~rep(c('cl1', 'cl2'), times=c(7, 14)));
fit <- lmFit(eSet, design.trt);
fit <- eBayes(fit)
df <- topTable(fit, coef=2, n=2000);

hh1<-EnhancedVolcano(df, pCutoff=0.05,
                     lab = rownames(df), FCcutoff=2,
                     x = 'logFC', title ='', subtitle = '',
                     y = 'adj.P.Val')+theme(aspect.ratio = 1)+ylim(c(0, 5))

plot_grid(hh1,'', '', '', ncol=2);
dev.copy2pdf(file='0_differential vocanoplot cl1 vs. cl2.pdf', family='ArialMT');

## ggplot volcano with labels
volcano.label<-read_xlsx('20220803/2022 08 03 edits_TTTS metab fig.xlsx');
volcano.label<-as.data.frame(volcano.label);

df$label2<-''; df[volcano.label$`Supp Table 4 file label`, 'label2']<-volcano.label$`Volcano Label`;
label.df<-df[which(df$label2!=''), ];

library(ggrepel);

EnhancedVolcano(df, x='logFC', y='adj.P.Val', 
                lab = '', pCutoff=0.05, ylim = c(0, 6), xlim = c(-5, 10),
                title = '', subtitle = '', FCcutoff = 2,colAlpha=1,
                col = c("grey50", "grey50", "grey50", "brown")) +
  geom_text_repel(data = label.df, mapping = aes(x=label.df$logFC, 
                  y=-log10(label.df$adj.P.Val)), xlim = c(6, 9),
                  segment.colour = "grey",
                  label = label.df$label2, direction = "y")+
                  theme(aspect.ratio = 1)
dev.copy2pdf(file='0_vocalnoplot_labels.pdf', family = 'ArialMT', width=7, height=7)


diff.1<-df[which(abs(df$logFC)>=2 & df$adj.P.Val<0.001), ]; ## 47
diffs.1<-meta.mat[rownames(df)[which(abs(df$logFC)>=2&df$adj.P.Val<0.001)], ]; ## 
diffs.z<-t(apply(log2(diffs.1+1), 1, function(x){
  (x-mean(x[1:10], na.rm=T))/sd(x, na.rm = T)
}))

#### T22 is too similar to control, is not included for calculation
col.cl<-rep('Ctrl', 32);
names(col.cl)<-colnames(meta.mat);
col.cl[cl1]<-'Cl1'; col.cl[cl2]<-'Cl2'; 
col.split<-col.cl;
col.split[which(col.split%in%c('Ctrl', 'Cl1'))]<-'Others';

### not showing names with @
tmp.name<-diffs.z; p1<-grep('@', rownames(diffs.z)); tmp.name<-tmp.name[-p1, ]; rownames(tmp.name)<-gsub('^._', '', rownames(tmp.name));

Heatmap(tmp.name, col = plasma(1024), 
        row_split = meta[which(meta.rnames%in%rownames(diffs.z)[-p1]), 'SubPathway'],
        row_title_rot = 0, column_split = col.split,
        top_annotation=HeatmapAnnotation(Subtype=col.cl, 
        col = list(Subtype=c(Ctrl='black', Cl1='steelblue', Cl2='blue', Case='light blue'))));
dev.copy2pdf(file='0_heatmap of differential metabolites cl1 vs. cl2.pdf', family='ArialMT', width=16);

diff.names<-rownames(diffs.z);
diff.df<-df[diff.names, ];
write.table(diff.df, file='0_differential metabolites cl1 vs. cl2.txt', 
            col.names = T, row.names = T, sep='\t', quote = F);

### visualize with picked labels
supp.meta<-read_xlsx('figure panels/20220731/2022 07 28 edits_TTTS metab fig.xlsx',
                     sheet = 1);
supp.meta<-as.data.frame(supp.meta);
xx<-setdiff(supp.meta$logFC, meta$`Compound Name`);
pos.orig<-c(217, 337, 383, 72, 340, 59, 129, 115)
old2new<-cbind(diffs=xx, new=meta.rnames[pos.orig]);
rownames(old2new)<-xx;

xx2<-setdiff(supp.meta$logFC, xx);
pos<-which(meta$`Compound Name`%in%xx2)
xx2.meta<-meta[pos, ]
xx2.meta$rnames<-meta.rnames[pos]
rownames(xx2.meta)<-xx2.meta$`Compound Name`;

supp.meta.df<-rbind(old2new, cbind(diffs=rownames(xx2.meta), new=xx2.meta$rnames))
rownames(supp.meta.df)<-supp.meta.df[, 1]
supp.meta.df<-supp.meta.df[supp.meta$logFC, ]

rownames(supp.meta)<-supp.meta.df[, 2];
cmm.pos<-intersect(rownames(supp.meta), rownames(tmp.name));
supp.meta<-supp.meta[cmm.pos, ];

Heatmap(tmp.name[cmm.pos, ], col = plasma(1024), 
        row_split = supp.meta$`Fig Pathway`,
        row_labels = supp.meta$`Heatmap Label for Fig`,
        row_title_rot = 0, column_split = col.split,
        top_annotation=HeatmapAnnotation(Subtype=col.cl, 
        col = list(Subtype=c(Ctrl='black', Cl1='steelblue', Cl2='blue', Case='light blue'))));
dev.copy2pdf(file='00_heatmap of differential metabolites cl1 vs. cl2.pdf', family='ArialMT', width=16);


########### ########### ########### ########### ########### 
########### ########### 
########### 
###########  expression RNA
########### ########### ########### ########### 
library(readxl);
options(stringsAsFactors = F)
cf.exp<-as.data.frame(read_xlsx('filtered.count.matrix.xlsx'));
rownames(cf.exp)<-cf.exp[, 1];
cf.exp<-cf.exp[, -1];

exo.exp<-as.data.frame(read_xlsx('filtered.VSD.GEP.mtx.xlsx'));
rownames(exo.exp)<-exo.exp[, 1];
exo.exp<-exo.exp[, -1];

nor.exp<-exo.exp[rownames(cf.exp), ]; ## 13369    31
nor.cf.exp<-nor.exp[, 1:15];
cf.vt<-rep(c('1Control', '2Case'), times=c(6, 9))
nor.exo.exp<-nor.exp[, 16:31];
exo.vt<-rep(c('1Control', '2Case'), times=c(6, 10))

####### pca + ellipse
library(factoextra);
## use sd 5000
sd<-sort(apply(nor.cf.exp, 1, sd, na.rm=T), decreasing = T);
sd5000<-names(sd)[1:1000];
tmp<-nor.cf.exp[sd5000, ]
res.pca <- prcomp(t(tmp), scale. = F);
t1<-fviz_eig(res.pca)

t2<-fviz_pca_ind(res.pca, repel = TRUE,
                 col.ind=rep(c('Control', 'Case'), times=c(6, 9)))
# plot_grid(t1, t2, '', '', ncol=2);
library(car);
par(mfrow=c(2, 2));
dataEllipse(res.pca$x[, 1], res.pca$x[, 2], 
                groups = as.factor(cf.vt),
                group.labels = c('Contrl', 'Case'),
                level=.95, fill=TRUE, fill.alpha=0.1)

sd<-sort(apply(nor.exo.exp, 1, sd, na.rm=T), decreasing = T);
sd5000<-names(sd)[1:1000];
tmp<-nor.exo.exp[sd5000, ]
res.pca <- prcomp(t(tmp), scale. = F);
tt1<-fviz_eig(res.pca)

tt2<-fviz_pca_ind(res.pca, repel = TRUE,
                 col.ind=rep(c('Control', 'Case'), times=c(6, 10)))
# plot_grid(t1, t2, '', '', ncol=2);
library(car);
par(mfrow=c(2, 2));
dataEllipse(res.pca$x[, 1], res.pca$x[, 2], 
            groups = as.factor(exo.vt),
            group.labels = c('Contrl', 'Case'),
            level=.95, fill=TRUE, fill.alpha=0.1);

################ ################ ################ ################ 
################ deg vocano using this result
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

exo.deg<-read_xlsx('exosomalRNA_TTTS.vs.exosomalRNA_BCM.ttest.result.xlsx', col_names = T)
exo.deg<-as.data.frame(exo.deg); rownames(exo.deg)<-exo.deg[, 1];
h1<-EnhancedVolcano(exo.deg, x='log2FoldChange', y='Qvalue', lab=exo.deg$Geneid, pCutoff=0.1, FCcutoff=log2(1.5), ylim=c(0, 3))+theme(aspect.ratio = 1)

cf.deg<-read_xlsx('cfRNA_TTTS.vs.cfRNA_BCM.ttest.result.xlsx', col_names = T)
cf.deg<-as.data.frame(cf.deg); rownames(cf.deg)<-cf.deg[, 1];
h2<-EnhancedVolcano(cf.deg, x='log2FoldChange', y='Qvalue', lab=cf.deg$Geneid, pCutoff=0.1, FCcutoff=log2(1.5), ylim=c(0, 3))+theme(aspect.ratio = 1)

library(cowplot);
plot_grid(h1, h2, '', '', ncol=2)
ggsave2(file='0_exo_cf_volcano.pdf', family='ArialMT');

####
exo.sig<-rownames(exo.deg)[which(exo.deg$Qvalue<0.1&exo.deg$absFoldChange>=1.5)]
write.csv(exo.deg[exo.sig, ], file='0_allT vs. normal_significant_exoRNA.csv');

cf.sig<-rownames(cf.deg)[which(cf.deg$Qvalue<0.1&cf.deg$absFoldChange>=1.5)]
write.csv(cf.deg[cf.sig, ], file='0_allT vs. normal_significant_cfRNA.csv');

##### split directions and plot piechart for diff genes
deg.list<-list(exoRNA.up=exo.sig[which(exo.deg[exo.sig, 'Direction']==1)],
               exoRNA.down=exo.sig[which(exo.deg[exo.sig, 'Direction']=='-1')],
               cfRNA.up=cf.sig[which(cf.deg[cf.sig, 'Direction']==1)],
               cfRNA.down=cf.sig[which(cf.deg[cf.sig, 'Direction']=='-1')])

library(eulerr);
par(mfrow=c(2, 2), pty='s')
plot(euler(deg.list, shape = "ellipse"), quantities = TRUE)
sig.labels<-intersect(deg.list$exoRNA.up, deg.list$cfRNA.up);
### [1] "TRAK2"    "SH3PXD2A" "TAF15"    "RAB13"    "SMARCC1"  "NCL"      "MAP4K4"  
dev.copy2pdf(file='0_venn_exo_cf_DEGs.pdf', family='ArialMT')

################ ################ ################ ################ ################ ################ ################ ################ 
################ ################ ################ ################ 
################ ################ GSEA analysis
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
library(clusterProfiler);
library(org.Hs.eg.db);
library(ggnewscale); library(enrichplot); library(wordcloud); library(xlsx)

all.gsea<-read.gmt('msigdb.v7.5.1.entrez.gmt.txt');
#### remove gene sets with out clear information in their names
## start with module_
# pos<-grep('^MODULE_', all.gsea[, 1]);
# all.gsea<-all.gsea[-pos, ];

## exo-RNA
original_gene_list <- exo.deg$log2FoldChange;
names(original_gene_list) <- exo.deg$Geneid;
gene_list<-na.omit(original_gene_list);
gene_list = sort(gene_list, decreasing = TRUE); 

name2id<-bitr(names(gene_list), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db');
name2id<-name2id[which(!duplicated(name2id[, 1])), ];
name2id<-name2id[which(!duplicated(name2id[, 2])), ]; ## 12596
rownames(name2id)<-name2id$ALIAS;

gene.fi<-gene_list;
gene.fi<-gene_list[name2id[, 1]]; names(gene.fi)<-as.character(name2id[, 2]);
gene.fi = sort(gene.fi, decreasing = TRUE); 

exo.gsea <- GSEA(gene.fi, TERM2GENE = all.gsea);
write.csv(exo.gsea@result, file='0_exoRNA_GSEA_all.sig.csv')

## cf-RNA
original_gene_list <- cf.deg$log2FoldChange;
names(original_gene_list) <- cf.deg$Geneid;
gene_list<-na.omit(original_gene_list);
gene_list = sort(gene_list, decreasing = TRUE); 

name2id<-bitr(names(gene_list), fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db');
name2id<-name2id[which(!duplicated(name2id[, 1])), ];
name2id<-name2id[which(!duplicated(name2id[, 2])), ]; ## 12596
rownames(name2id)<-name2id$ALIAS;

gene.fi<-gene_list;
gene.fi<-gene_list[name2id[, 1]]; names(gene.fi)<-as.character(name2id[, 2]);
gene.fi = sort(gene.fi, decreasing = TRUE); 

cf.gsea <- GSEA(gene.fi, TERM2GENE = all.gsea);
write.csv(cf.gsea@result, file='0_cfRNA_GSEA_all.sig.csv')

library(readxl);
all.expr<-read_xlsx('filtered.VSD.GEP.mtx.xlsx');
all.expr<-as.data.frame(all.expr); rownames(all.expr)<-all.expr[, 1]; all.expr<-all.expr[, -1];

cf.expr<-all.expr[, grep('_cfRNA_', colnames(all.expr))]; ## control 6, case 9
exo.expr<-all.expr[, grep('_exosomalRNA_', colnames(all.expr))]; ## control 6, case 10

library(RColorBrewer);
coul <- rev(colorRampPalette(brewer.pal(8, "PiYG"))(1024));

library(ComplexHeatmap);
pos<-grep('BCM131', colnames(cf.expr));
h1<-Heatmap(t(scale(t(cf.expr[cf.sig, -pos]))), col=coul, 
        row_split = factor(cf.deg[cf.sig, 'Direction'], levels=c("-1", '1')), show_column_names = F,
        column_split = rep(c('Ctrl', 'TTTS'), times=c(6-1, 9)));
dev.copy2pdf(file='00_cfRNA_deg heatmap.pdf', family='ArialMT');

pos<-c(10, 11, 12);
h2<-Heatmap(t(scale(t(exo.expr[exo.sig, -pos]))), col=coul, show_column_names = F,
        row_split = factor(exo.deg[exo.sig, 'Direction'], levels=c("1", '-1')), 
        column_split = rep(c('Ctrl', 'TTTS'), times=c(6, 10-3)));
dev.copy2pdf(file='00_exoRNA_deg heatmap.pdf', family='ArialMT');


######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## ######## ######## ######## ######## ######## 
######## ######## exoRNA
sort.exo.gsea<-exo.gsea@result[order(exo.gsea@result$NES), ]
tmp.sort<-rbind(head(sort.exo.gsea, 15), tail(sort.exo.gsea, 15))

par(mfrow=c(2, 2), pty='s')
barplot(tmp.sort$NES, names.arg = tmp.sort$ID, horiz = T, las=2, xlim = c(-4, 3), border = F)
barplot(sort.exo.gsea[grep('HALLMARK_', sort.exo.gsea$ID), ]$NES, names.arg = sort.exo.gsea[grep('HALLMARK_', sort.exo.gsea$ID), ]$ID, horiz = T, las=2, xlim = c(-4, 3), border = F)
barplot(sort.exo.gsea[grep('KEGG_', sort.exo.gsea$ID), ]$NES, names.arg = sort.exo.gsea[grep('KEGG_', sort.exo.gsea$ID), ]$ID, horiz = T, las=2, xlim = c(-4, 3), border = F)

dev.copy2pdf(file='0_barchart_top_tail_15gene sets_KEGG_hallmark_exoRNA.pdf', family='ArialMT');

######## ######## cfRNA
sort.cf.gsea<-cf.gsea@result[order(cf.gsea@result$NES), ]
tmp.sort<-rbind(head(sort.cf.gsea, 15), tail(sort.cf.gsea, 15))

par(mfrow=c(2, 2), pty='s')
barplot(tmp.sort$NES, names.arg = tmp.sort$ID, horiz = T, las=2, xlim = c(-4, 3), border = F)
barplot(sort.cf.gsea[grep('HALLMARK_', sort.cf.gsea$ID), ]$NES, names.arg = sort.cf.gsea[grep('HALLMARK_', sort.cf.gsea$ID), ]$ID, horiz = T, las=2, xlim = c(-4, 3), border = F)
barplot(sort.cf.gsea[grep('KEGG_', sort.cf.gsea$ID), ]$NES, names.arg = sort.cf.gsea[grep('KEGG_', sort.cf.gsea$ID), ]$ID, horiz = T, las=2, xlim = c(-4, 3), border = F)

dev.copy2pdf(file='0_barchart_top_tail_15gene sets_KEGG_hallmark_cfRNA.pdf', family='ArialMT');


######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## ######## ######## ######## ######## ######## visualize a few terms using gsea

p1 <- gseaplot2(exo.gsea, geneSetID = 'KEGG_DILATED_CARDIOMYOPATHY', title = 'KEGG_DILATED_CARDIOMYOPATHY')
p2 <- gseaplot2(exo.gsea, geneSetID = 'HALLMARK_OXIDATIVE_PHOSPHORYLATION', title = 'HALLMARK_OXIDATIVE_PHOSPHORYLATION')
p5<-gseaplot2(exo.gsea, geneSetID ='REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION', title='Translation elongation')

p3 <- gseaplot(cf.gsea, geneSetID ='HALLMARK_MYOGENESIS', title = 'HALLMARK_MYOGENESIS');
p4 <- gseaplot(cf.gsea, geneSetID ='HALLMARK_HYPOXIA', title = 'HALLMARK_HYPOXIA');

cowplot::plot_grid(plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4]), '', '', '', ncol=2)
dev.copy2pdf(file='0_GSEA examples.pdf', family='ArialMT')

cowplot::plot_grid(plot_grid(p1, p5, p2, ncol=3, 
                  labels=LETTERS[1:3]), '', ncol=1)
dev.copy2pdf(file='0_supp GSEA exo examples.pdf', family='ArialMT', width=11, height=7)

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## ######## ######## ######## ######## ######## visualize barcharts for GSEAs
exo.rec<-read_xlsx('figure panels/20220731/2022 07 28 edits_TTTS RNAseq fig.xlsx', sheet = 2)
exo.rec<-as.data.frame(exo.rec)[c(1:10, 30:40), ];

col <- exo.rec$`Figure Label`
val <- exo.rec$NES
df <- data.frame(col, val)
df<-df[order(df$val), ];

# create graph
h1<-ggplot(data = df, aes(x = reorder(col, val), y = val)) +
  geom_bar(stat = 'identity', aes(fill = exo.rec$qvalues)) +
  scale_fill_gradientn(colours = c('blue', 'grey'))+
  coord_flip()+NoGrid()+theme(legend.position = 'bottom')

exo.keg<-read_xlsx('figure panels/20220731/2022 07 28 edits_TTTS RNAseq fig.xlsx', 
                   sheet = 3)
exo.keg<-as.data.frame(exo.keg)[c(1:10, 17:26), ];

col <- exo.keg$`Figure Label`
val <- exo.keg$NES
df <- data.frame(col, val)

# create graph
h2<-ggplot(data = df, aes(x = reorder(col, val), y = val)) +
  geom_bar(stat = 'identity', aes(fill = exo.keg$qvalues)) +
  scale_fill_gradientn(colours = c('blue', 'grey'))+
  coord_flip()+NoGrid()+theme(legend.position = 'bottom')


cf.rec<-read_xlsx('figure panels/20220731/2022 07 28 edits_TTTS RNAseq fig.xlsx', 
                  sheet = 4)
cf.rec<-as.data.frame(cf.rec)[c(1:10, 32:39), ];

col <- cf.rec$`Figure Label`
val <- cf.rec$NES
df <- data.frame(col, val)
df<-df[order(df$val), ];

# create graph
h3<-ggplot(data = df, aes(x = reorder(col, val), y = val)) +
  geom_bar(stat = 'identity', aes(fill = cf.rec$qvalues)) +
  scale_fill_gradientn(colours = c('blue', 'grey'))+
  coord_flip()+NoGrid()+theme(legend.position = 'bottom')


cf.keg<-read_xlsx('figure panels/20220731/2022 07 28 edits_TTTS RNAseq fig.xlsx', 
                   sheet = 5)
cf.keg<-as.data.frame(cf.keg);

col <- cf.keg$`Fig Label`
val <- cf.keg$NES
df <- data.frame(col, val)

# create graph
h4<-ggplot(data = df, aes(x = reorder(col, val), y = val)) +
  geom_bar(stat = 'identity', aes(fill = cf.keg$qvalues)) +
  scale_fill_gradientn(colours = c('blue', 'grey'))+
  coord_flip()+NoGrid()+theme(legend.position = 'bottom')

plot_grid(h1, h2, h3, h4, ncol=2, axis = 'l', align = 'hv', byrow = F)
dev.copy2pdf(file='figure panels//0_reactome+kegg picked.pdf', family='ArialMT',
             width=9, height=7);

######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## ######## ######## 
############### word cloud for met pathways 
library(wordcloud);
word.c<-read_xlsx('20220803/2022 08 03 edits_TTTS metab fig.xlsx', sheet = 2);
word.c<-as.data.frame(word.c);

words<-word.c$`pathway term used`;
words<-words[which(!is.na(words))];

wds<-sapply(words, strsplit, split=' ')
wds<-as.character(unlist(wds))
wds.tly<-table(wds)
wds.df<-data.frame(word=names(wds.tly), freq=as.numeric(wds.tly))[-3, ]

wordcloud(words = wds.df$word, freq = wds.df$freq, min.freq = 1, 
          max.words=200, random.order=FALSE, rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))

################################# ######## ######## ######## ######## ######## 
######## ######## ######## ######## ######## 
######## ######## data clean date 20221101

###### all T vs. controls
## ggplot volcano with labels
volcano.label<-read_xlsx('20221028/2022 10 26 supp table_diff metab all cases vs ctrl clean.xlsx', sheet=2);
volcano.label<-as.data.frame(volcano.label);

huh<-sapply(volcano.label$Metabolite, function(x)cbind(x, agrep(x, rownames(df)), rownames(df)[agrep(x, rownames(df))]));
huh<-do.call('rbind', huh);
huh.df<-cbind(huh, df[as.numeric(huh[, 2]), ]);
## pos 2,107, 123,88,55,117, 72,67,41, 142, 136, 59, 134, 106,57
pos<-c(2,107, 123,88,55,117, 72,67,41, 142, 136, 59, 134, 106,57)
volcano.label<-huh.df[which(as.numeric(huh.df$V2)%in%pos), ];

df$label2<-''; df[volcano.label$V3, 'label2']<-volcano.label$x;
label.df<-df[which(df$label2!=''), ];

library(ggrepel);

EnhancedVolcano(df, x='logFC', y='adj.P.Val', 
                lab = '', pCutoff=0.05, xlim=c(-15, 30),
                title = '', subtitle = '', FCcutoff = 2,colAlpha=1,
                col = c("grey50", "grey50", "grey50", "brown")) +
  geom_text_repel(data = label.df, mapping = aes(x=label.df$logFC, 
                  y=-log10(label.df$adj.P.Val)), xlim = c(25, 29),
                  segment.colour = "grey",
                  label = label.df$label2, direction = "y")+
  theme(aspect.ratio = 1)
dev.copy2pdf(file='20221028/0_vocalnoplot_labels_original.pdf', family = 'ArialMT', width=7, height=7);


###### cl1 vs. cl2
## ggplot volcano with labels
volcano.label<-read_xlsx('20221028/2022 10 26 supp table_diff metab c1 vs c2 clean.xlsx', sheet=2);
volcano.label<-as.data.frame(volcano.label);
rownames(volcano.label)<-volcano.label$Metabolite;
rnames<-volcano.label$`Volcano Label`;

huh<-sapply(volcano.label$Metabolite, function(x)cbind(x, volcano.label[x, 2], 
            agrep(x, rownames(df)), rownames(df)[agrep(x, rownames(df))]));
huh<-do.call('rbind', huh);
huh.df<-cbind(huh, df[as.numeric(huh[, 3]), ]);
## pos 
pos<-c(21,4,22,1,38,11,26,5,35,6, 8,46,18,51,23,32,52);
rpos<-paste(rnames, pos, sep='');
volcano.label<-huh.df[which(paste(huh.df$V2, huh.df$V3, sep='')%in%rpos), ];

df$label2<-''; df[volcano.label$V4, 'label2']<-volcano.label$x;
label.df<-df[which(df$label2!=''), ];

library(ggrepel);
EnhancedVolcano(df, x='logFC', y='adj.P.Val', 
                lab = '', pCutoff=0.05, ylim=c(0, 6),
                title = '', subtitle = '', FCcutoff = 2,colAlpha=1,
                col = c("grey50", "grey50", "grey50", "brown")) +
  geom_text_repel(data = label.df, mapping = aes(x=label.df$logFC, 
                                                 y=-log10(label.df$adj.P.Val)), xlim = c(6, 9),
                  segment.colour = "grey",
                  label = label.df$label2, direction = "y")+
  theme(aspect.ratio = 1)
dev.copy2pdf(file='20221028/0_vocalnoplot_labels_cl1 vs. cl2_original.pdf', family = 'ArialMT', width=7, height=7);

############### ############### ############### ############### ############### 
############### rank c and t average signals
ctrls<-apply(meta.mat[, 1:10], 1, median, na.rm=T);
cases<-apply(meta.mat[, -c(1:10)], 1, median, na.rm=T);

ctrls<-sort(ctrls, decreasing = T); 
cases<-sort(cases, decreasing = T);

df.ctrl<-data.frame(name=names(ctrls), signal=as.numeric(ctrls));
df.case<-data.frame(name=names(cases), signal=as.numeric(cases));

write.table(df.ctrl, file='20221028/control average signal rank decreasing.txt', sep='\t', quote = F, 
            row.names = F);
write.table(df.case, file='20221028/cases average signal rank decreasing.txt', sep='\t', quote = F, 
            row.names = F);








