# FROM NATHAN
library(dplyr)
## ck

unique_id_var = "individual"


sce <- sce_subset
sce$log10_total_counts <- log10(sce$total_counts)
mat <- as(SingleCellExperiment::counts(sce), "dgCMatrix")
#norm_mat <- sctransform::vst(mat)
vst_out <- sctransform::vst(
  umi = mat,
  cell_attr = as.data.frame(SummarizedExperiment::colData(sce)),
  latent_var = c("log10_total_counts", "pc_mito"),
  batch_var = NULL,
  latent_var_nonreg = NULL,
  n_genes = 2000,
  n_cells = NULL,
  return_gene_attr = TRUE,
  return_cell_attr = TRUE,
  method = "poisson",
  do_regularize = TRUE,
  residual_type = "pearson",
  show_progress = TRUE,
  return.only.var.genes = FALSE
  )
norm_mat_qn <- limma::normalizeQuantiles(norm_mat$y)
ridx <- rownames(sce) %in% rownames(norm_mat_qn)
sce <- sce[ridx, ]
SingleCellExperiment::normcounts(sce) <- norm_mat_qn



#fit <- limma::lmFit(SingleCellExperiment::normcounts(sce), model.matrix(~sce$sex))
#fit <- limma::lmFit(SingleCellExperiment::counts(sce), model.matrix(~sce$sex))
fit <- limma::lmFit(SingleCellExperiment::normcounts(pb), model.matrix(~pb$sex))
#fit <- limma::lmFit(SingleCellExperiment::counts(sce), model.matrix(~sce$diagnosis))
fit <- limma::eBayes(fit)

## or
head(coef(fit))
#contr <- limma::makeContrasts("sce$diagnosiscontrol", levels = colnames(coef(fit)))
#tmp <- contrasts.fit(fit, contrasts = contr)
#tt  <- topTable(fit, coef = "annot$dxcontrol", adjust = "fdr",number=10000000)
tt  <- limma::topTable(fit, coef = "pb$diagnosiscontrol", adjust = "fdr",number=10000000)
tt  <- limma::topTable(fit, coef = "pb$sexmale", adjust = "fdr",number=10000000)
#tt  <- limma::topTable(fit, adjust = "fdr",number=10000000)
tt$ensembl_gene_id <- rownames(tt)
chr_mappings <- map_ensembl_gene_id(tt[1:50,]$ensembl_gene_id, mappings = c("chromosome_name", "external_gene_name"))
dt <- dplyr::left_join(tt, chr_mappings)
table(chr_mappings$chromosome_name)
plot(table(chr_mappings$chromosome_name))


#table(chr_mappings$chromosome_name)
dt$padj <- dt$adj.P.Val
dt$gene <- dt$external_gene_name

de_genes <- dt %>%
  dplyr::filter(padj <= 0.05 & 2^abs(logFC) >= 2^0.25)

scFlow:::.volcano_plot(dt, fc_threshold = 2^0.25, pval_cutoff = 0.05, n_label = 12)

write.table(de_genes, "de_genes.tsv", row.names = F, sep = "\t")
write.table(dt, "de_genes.tsv", row.names = F, sep = "\t")

table(dt[1:100,]$chromosome_name)


#install_github('MacoskoLab/liger')
#detach("package:liger", unload=TRUE)
library(scFlow)
mysce = scFlow::read_sce("/var/final_sce")
aaa = mysce[,mysce$cluster_celltype=="Astro"]
#save(aaa,file="aaa.rda")
library(sctransform)
#devtools::install_github(repo = 'ChristophH/sctransform')
exp_dgC = as(aaa@assays@data[[1]], "CsparseMatrix")
normalized_data <- sctransform::vst(exp_dgC)$y
#save(normalized_data,file="normalized_data.rda")
rowS = Matrix::rowSums(exp_dgC)

library(Matrix)
normDat3=Matrix(normDat2)
allAnnot = data.frame(indv=aaa$individual,region=aaa$brain_region,mito=aaa$pc_mito,ribo=aaa$pc_ribo,PM_delay=aaa$PM_delay,diagnosis=aaa$diagnosis,braak=aaa$Braak_tangle_stage,gender=aaa$sex,RIN=aaa$RIN,age=aaa$age)
#load("/Users/natske/Downloads/allAnnot-6.rda")
allAnnot$braak = as.character(allAnnot$braak)
allAnnot$braak[allAnnot$braak=="'0'"] = 0
allAnnot$braak[allAnnot$braak=="I"] = 1
allAnnot$braak[allAnnot$braak=="I_Steve"] = 1
allAnnot$braak[allAnnot$braak=="I_or_II_Steve"] = 2
allAnnot$braak[allAnnot$braak=="III"] = 3
allAnnot$braak[allAnnot$braak=="IV"] = 4
allAnnot$braak[allAnnot$braak=="V"] = 5
allAnnot$braak[allAnnot$braak=="VI"] = 6
allAnnot$braak = as.numeric(allAnnot$braak)
#load("/Users/natske/Downloads/allDat.rda")
meanDat = matrix(0,nrow=dim(normDat3)[1],ncol=length(unique(allAnnot$indv)))
colnames(meanDat) = as.character(unique(allAnnot$indv))
rownames(meanDat) = rownames(normDat3)
#region="EC"
for(region in c("EC","SSC")){
  library(preprocessCore)
  library(limma)
  for(indv in as.character(unique(allAnnot$indv))){
    whichCells = allAnnot$indv==indv & allAnnot$region==region
    theData = normDat3[,whichCells]
    # NOTE: scTransform doesn't make each cell have equal amounts
    # hist(Matrix::colSums(theData))
    subAnnot = allAnnot[whichCells,]
    cpmData = normalize.quantiles(as.matrix(theData)) #t(t(theData)/colSums(theData))

    design <- model.matrix(~subAnnot$mito)#+subAnnot$ribo)
    fit <- lmFit(cpmData, design)
    fit <- eBayes(fit)
    vals = residuals(fit,cpmData)
    cpmData2 = normalize.quantiles(as.matrix(vals)) #t(t(theData)/colSums(theData))
    meanDat[,indv] = Matrix::rowSums(cpmData2)

    #meanDat[,indv] = Matrix::rowSums(cpmData)
  }
  normedDat = normalize.quantiles(meanDat)
  colnames(normedDat) = colnames(meanDat)
  rownames(normedDat) = rownames(meanDat)

  annot = unique(data.frame(indv=as.character(allAnnot$indv),
                            region=as.character(allAnnot$region),
                            gender=as.character(allAnnot$gender),
                            dx=as.character(allAnnot$diagnosis),
                            PM_delay=as.character(allAnnot$PM_delay),
                            age=as.character(allAnnot$age),
                            RIN=as.character(allAnnot$RIN),
                            braak=as.numeric(allAnnot$braak)))
  annot$PM_delay = as.numeric(as.character(annot$PM_delay))
  annot$RIN = as.numeric(as.character(annot$RIN))
  annot$PM_delay[is.na(annot$PM_delay)] = mean(annot$PM_delay,na.rm=TRUE)
  annot$RIN[is.na(annot$RIN)] = mean(annot$RIN,na.rm=TRUE)
  annot = annot[annot$region==region,]
  rownames(annot) = annot$indv
  annot = annot[colnames(normedDat),]
  #annot$indv==colnames(normedDat)

  library(limma)
  #Group <- factor(annot$dx+annot$gender)
  #Group <- factor(annot$gender)
  design <- model.matrix(~annot$braak+annot$gender+as.numeric(annot$RIN))
  #design <- model.matrix(~annot$gender)
  #design <- model.matrix(~annot$dx+annot$gender+as.numeric(annot$PM_delay))
  fit <- lmFit(normedDat, design)
  fit <- eBayes(fit)
  #tt  <- topTable(fit, coef = "annot$dxcontrol", adjust = "fdr",number=10000000)
  tt  <- topTable(fit, coef = "annot$braak", adjust = "fdr",number=10000000)
  tt[1:10,]
  write.csv(tt,file=sprintf("tt_%s_dx.csv",region))
  tt  <- topTable(fit, coef = "annot$gendermale", adjust = "fdr",number=10000000)
  tt[1:10,]
  write.csv(tt,file=sprintf("tt_%s_gender.csv",region))
}

gene="ENSG00000229807" # Xist
geneDat = data.frame(exp=meanDat[gene,],what=annot$gender)
library(ggplot2)
ggplot(geneDat) + geom_boxplot(aes(x=what,y=exp))
design <- model.matrix(~annot$dx+annot$gender+as.numeric(annot$PM_delay))
fit <- lmFit(normedDat, design)
fit <- eBayes(fit)
tt  <- topTable(fit, coef = "annot$dxcontrol", adjust = "fdr",number=10000000)
#tt  <- topTable(fit, coef = "annot$braak", adjust = "fdr",number=10000000)
tt[1:10,]
write.csv(tt,file=sprintf("tt_%s_dx.csv",region))
tt  <- topTable(fit, coef = "annot$gendermale", adjust = "fdr",number=10000000)
tt[1:10,]
write.csv(tt,file=sprintf("tt_%s_gender.csv",region))
gene="ENSG00000198719" # top hit from combiz
geneDat = data.frame(exp=meanDat[gene,],what=annot$dx)
normedDat["ENSG00000198719",]
