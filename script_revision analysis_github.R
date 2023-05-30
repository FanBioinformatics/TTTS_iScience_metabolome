##################### #####################
########## ##########
##### check for sex difference
#####
#####
#####

library(readxl);
options(stringsAsFactors = F);
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

df.case.sex<-as.data.frame(read_xlsx('sex difference/TTTS metab cohort_fetal sex .xlsx', sheet = 1));
df.case.sex$`TTTS ID`<-paste('T', df.case.sex$`TTTS ID`, sep='');
rownames(df.case.sex)<-df.case.sex$`TTTS ID`;
df.case.sex$sex<-NA;
df.case.sex$sex[which(df.case.sex$`Female fetus`==0)]<-'male'
df.case.sex$sex[which(df.case.sex$`Female fetus`==1)]<-'female'

df.contrl.sex<-as.data.frame(read_xlsx('sex difference/TTTS metab cohort_fetal sex .xlsx', sheet = 2));
df.contrl.sex$`TTTS ID`<-paste('C', df.contrl.sex$`BCM ID`, sep='');
df.contrl.sex$sex<-NA;
df.contrl.sex$sex[which(df.contrl.sex$`Female fetus`==0)]<-'male';
df.contrl.sex$sex[which(df.contrl.sex$`Female fetus`==1)]<-'female';

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### difference in T3S
library(EnhancedVolcano); library(limma); library(Biobase);
eSet <- ExpressionSet(assayData = as.matrix(apply(meta.mat[, rownames(df.case.sex)]+1, 
        c(1, 2), log2)), phenoData = as(df.case.sex, "AnnotatedDataFrame"));
design.trt=model.matrix(~sex, data = df.case.sex);

fit <- lmFit(eSet, design.trt);
fit <- eBayes(fit)
df.case <- topTable(fit, coef=2, n=2000);
write.csv(df.case, file='sex difference/Sex difference in TTTS.csv')

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### difference in contrl
rownames(df.contrl.sex)<-df.contrl.sex$`BCM ID`;
df.contrl.sex<-df.contrl.sex[which(!is.na(df.contrl.sex$sex)), ];
eSet <- ExpressionSet(assayData = as.matrix(apply(meta.mat[, df.contrl.sex$`BCM ID`]+1, 
        c(1, 2), log2)), phenoData = as(df.contrl.sex, "AnnotatedDataFrame"));
design.trt=model.matrix(~sex, data = df.contrl.sex);

fit <- lmFit(eSet, design.trt);
fit <- eBayes(fit)
df.contrl <- topTable(fit, coef=2, n=2000);
write.csv(df.contrl, file='sex difference/Sex difference in controls.csv')

##### plot results
EnhancedVolcano(df.case, pCutoff=0.05, y='P.Value', ylim=c(0, 3),
                     lab = rownames(df.case), FCcutoff=0,
                     x = 'logFC', title ='Sex difference in TTTS', subtitle = '')+
                     theme(aspect.ratio=1)
dev.copy2pdf(file='sex difference/volcano_Sex difference in TTTS.pdf', family='ArialMT');

EnhancedVolcano(df.contrl, pCutoff=0.05, ylim=c(0, 2),
                     lab = rownames(df.contrl), FCcutoff=0, y='P.Value',
                     x = 'logFC', title ='Sex difference in Controls', subtitle = '')+
                     theme(aspect.ratio=1);
dev.copy2pdf(file='sex difference/volcano_Sex difference in controls.pdf', family='ArialMT');

######### overlapping none for all!

intersect(rownames(df.case)[which(df.case$logFC>0 & df.case$P.Value<0.05)],
          rownames(df.contrl)[which(df.contrl$logFC>0 & df.contrl$P.Value<0.05)])

intersect(rownames(df.case)[which(df.case$logFC<0 & df.case$P.Value<0.05)],
          rownames(df.contrl)[which(df.contrl$logFC<0 & df.contrl$P.Value<0.05)])

intersect(rownames(df.case)[which(df.case$logFC>0 & df.case$P.Value<0.05)],
          rownames(df.contrl)[which(df.contrl$logFC<0 & df.contrl$P.Value<0.05)])

intersect(rownames(df.case)[which(df.case$logFC<0 & df.case$P.Value<0.05)],
          rownames(df.contrl)[which(df.contrl$logFC>0 & df.contrl$P.Value<0.05)])


############## ## ## ## ## ## ## revision 
## covariates 
library(readxl);
covars<-read_xlsx('TTTS metab ms covariates.xlsx', sheet = 1);
covars<-as.data.frame(covars);
rownames(covars)<-covars$`Study ID`;
tmp.mat<-meta.mat[, rownames(covars)];

cc1<-apply(tmp.mat, 1, function(x){
  h1<-cor.test(x, covars$`Mat age`);
  c(cor=h1$estimate, p=h1$p.value);
})

cc1<-t(cc1);
cc1<-as.data.frame(cc1);
cc1$FDR<-p.adjust(cc1[, 2]);

cc2<-apply(tmp.mat, 1, function(x){
  h2<-cor.test(x, covars$GA);
  c(cor=h2$estimate, p=h2$p.value)
  
})
cc2<-t(cc2);
cc2<-as.data.frame(cc2);
cc2$FDR<-p.adjust(cc2[, 2]);

par(mfrow=c(2, 2), pty='s');
plot(density(cc1$FDR), main='Correlations between Mat. age \nand metabolites(FDR distribution)')
plot(density(cc2$FDR), main='Correlations between GA \nand metabolites(FDR distribution)')

dev.copy2pdf(file='0_revision_covariants.pdf', family='ArialMT');


