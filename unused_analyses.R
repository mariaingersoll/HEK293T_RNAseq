######this file contains random bits of code that didn't go into the final analysis

######after shared_DEGs05_correlation_new.pdf:
#Now doing GO analysis on the N only, NLS only, and shared DEGs (via Fisher Exact test)
  #First format your input files
iso2go_name <- read.table("./Data_Files/GRCh38_iso2go.tab", sep="\t")
colnames(iso2go_name)=paste(c("gene_ID", "GO"))

  #starting with just the shared ones because those were easiest to pull out from above; doesn't matter whether we do N or NLS for this because they have the same
head(shared_N)
shared_fisher <- shared_N
shared_fisher <- rownames_to_column(shared_fisher, var="gene_ID")
head(shared_fisher) 
shared_fisher <- shared_fisher %>%
  mutate(Fishers = 1) %>%
  dplyr::select(gene_ID, Fishers) %>%
  full_join(iso2go_name) %>%
  mutate(Fishers = if_else(is.na(Fishers), 0, 1)) %>%
  select(gene_ID, Fishers)
head(shared_fisher)
write.csv(shared_fisher, file="./Data_Files/shared_fisher_GO.csv", quote=F, row.names=FALSE)


  #now trying to isolate the just N and just NLS ones
N_DEGs05_DF <- as.data.frame(N_DEGs05) %>%
  mutate(gene_ID = N_DEGs05) %>%
  mutate(Fishers_N = 1)%>%
  dplyr::select(gene_ID, Fishers_N)
NLS_DEGs05_DF <- as.data.frame(NLS_DEGs05) %>%
  mutate(gene_ID = NLS_DEGs05) %>%
  mutate(Fishers_NLS = 1)%>%
  dplyr::select(gene_ID, Fishers_NLS)
fj <- full_join(N_DEGs05_DF, NLS_DEGs05_DF) %>%
  full_join(iso2go_name)
tail(fj)

fishers_N_only <- fj %>%
  mutate(Fishers_N = if_else(is.na(Fishers_N), 0, 1)) %>%
  select(gene_ID, Fishers_N)
write.csv(fishers_N_only, file="./Data_Files/fishers_N_only_GO.csv", quote=F, row.names=FALSE)

fishers_NLS_only <- fj %>%
  mutate(Fishers_NLS = if_else(is.na(Fishers_NLS), 0, 1)) %>%
  select(gene_ID, Fishers_NLS)
write.csv(fishers_NLS_only, file="./Data_Files/fishers_NLS_only_GO.csv", quote=F, row.names=FALSE)
  #didn't do this for the new Refseq, because the gene IDs aren't set up right in the iso2go

######### after Venns
#Create a term go_input_NFlag, a term go_input_NLSFlag, and a term go_input_NNLS that has the -logpval as a neg number if the l2fc is < 0 and a pos number if the l2fc is > 0
  #- These terms will go into the GO analysis
  #- These terms are including na.omit so anything without a padj value will be removed

  ##N to Flag
head(res_N_Flag)
res_N_Flag = res_N_Flag %>%
  tibble::rownames_to_column(var = "gene_ID")
head(res_N_Flag)
nrow(res_N_Flag)
#59048

res_N_Flag_name <- res_N_Flag 
colnames(res_N_Flag_name)[1] = "ISO"
res_N_Flag_name <- res_N_Flag_name%>%
  left_join(iso2gene_good)

go_input_NFlag = res_N_Flag %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(gene_ID, mutated_p_updown)
nrow(go_input_NFlag)
#10966
head(go_input_NFlag)
colnames(go_input_NFlag) <- c("gene", "pval")
head(go_input_NFlag)
write.csv(go_input_NFlag, file="./Data_Files/NFlag_GO_24.csv", quote=F, row.names=FALSE)

##NLS to Flag
head(res_NLS_Flag)
res_NLS_Flag = res_NLS_Flag %>%
  tibble::rownames_to_column(var = "gene_ID")
head(res_NLS_Flag)
nrow(res_NLS_Flag)
#59048

res_NLS_Flag_name <- res_NLS_Flag 
colnames(res_NLS_Flag_name)[1] = "ISO"
res_NLS_Flag_name <- res_NLS_Flag_name%>%
  left_join(iso2gene_good)

go_input_NLSFlag = res_NLS_Flag %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(gene_ID, mutated_p_updown)
nrow(go_input_NLSFlag)
#9822
head(go_input_NLSFlag)
colnames(go_input_NLSFlag) <- c("gene", "pval")
head(go_input_NLSFlag)
write.csv(go_input_NLSFlag, file="./Data_Files/NLSFlag_GO_24.csv", quote=F, row.names=FALSE)

##N to NLS
head(res_N_NLS)
res_N_NLS = res_N_NLS %>%
  tibble::rownames_to_column(var = "gene_ID")
head(res_N_NLS)
nrow(res_N_NLS)
#59048

res_N_NLS_name <- res_N_NLS 
colnames(res_N_NLS_name)[1] = "ISO"
res_N_NLS_name <- res_N_NLS_name%>%
  left_join(iso2gene_good)

go_input_NNLS = res_N_NLS %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(gene_ID, mutated_p_updown)
nrow(go_input_NNLS)
#59048
head(go_input_NNLS)
colnames(go_input_NNLS) <- c("gene", "pval")
head(go_input_NNLS)
write.csv(go_input_NNLS, file="./Data_Files/NNLS_GO_24.csv", quote=F, row.names=FALSE)
###not doing this for the new Refseq yet

#################
#Prepping files for GSEA input
#need the rld expression value for each gene for each sample you're comparing
head(rldpvals_NFlag)
head(rld_NFlag)
head(rldpvals_NLSFlag)
head(rld_NLSFlag)

#make the iso a new column with the label "description"
rld_NFlag_new <- rld_NFlag
rld_NFlag_new <- rownames_to_column(rld_NFlag_new, var="DESCRIPTION")
head(rld_NFlag_new)

rld_NLSFlag_new <- rownames_to_column(rld_NLSFlag, var="DESCRIPTION")
head(rld_NLSFlag_new)

##combine your rld files in one of the following ways
#rld_N_NLS <- rld_NFlag_new %>%
#  left_join(rld_NLSFlag_new)

#rld_N_NLS <- left_join(x=rld_NFlag_new, y=rld_NLSFlag_new)

rld_N_NLS <- full_join(x=rld_NFlag_new, y=rld_NLSFlag_new)
head(rld_N_NLS)

#Get shared DEGS from Venns with 1.5 and 2 lfc filters
resNFlag<- read.delim("./Data_Files/N_Flag_DE_24.txt")
head(resNFlag)
resNFlag$lfc_NFlag <- resNFlag$log2FoldChange
head(resNFlag)

resNLSFlag<- read.delim("./Data_Files/NLS_Flag_DE_24.txt")
head(resNLSFlag)
resNLSFlag$lfc_NLSFlag <- resNLSFlag$log2FoldChange
head(resNLSFlag)

#1.5
intersect05_1.5 <- intersect(N_DEGs05_lfc1.5, NLS_DEGs05_lfc1.5)
shared_N_1.5 <- as.data.frame(resNFlag)[rownames(resNFlag) %in% intersect05_1.5, ] %>% select(log2FoldChange)
N_name <- c("N_lfc")
colnames(shared_N_1.5) <- N_name
shared_NLS_1.5 <- as.data.frame(resNLSFlag)[rownames(resNLSFlag) %in% intersect05_1.5, ] %>% select(log2FoldChange)
NLS_name <- c("NLS_lfc")
colnames(shared_NLS_1.5) <- NLS_name

shared_both_lfc_1.5 <- cbind(shared_N_1.5, shared_NLS_1.5)
head(shared_both_lfc_1.5)
scatter_1.5 <- ggscatter(shared_both_lfc_1.5, x = "N_lfc", y = "NLS_lfc", 
                         add = "reg.line", conf.int = TRUE, 
                         cor.coef = TRUE, cor.method = "pearson",
                         xlab = "DEG LogFoldChange N to Flag", ylab = "DEG LogFoldChange NLS to Flag")+ theme(panel.grid.major = element_line(colour = "gray80", linewidth = 0.5, linetype="dotted")) 

#2
intersect05_2 <- intersect(N_DEGs05_lfc2, NLS_DEGs05_lfc2)
shared_N_2 <- as.data.frame(resNFlag)[rownames(resNFlag) %in% intersect05_2, ] %>% select(log2FoldChange)
N_name <- c("N_lfc")
colnames(shared_N_2) <- N_name
shared_NLS_2 <- as.data.frame(resNLSFlag)[rownames(resNLSFlag) %in% intersect05_2, ] %>% select(log2FoldChange)
NLS_name <- c("NLS_lfc")
colnames(shared_NLS_2) <- NLS_name

shared_both_lfc_2 <- cbind(shared_N_2, shared_NLS_2)
head(shared_both_lfc_2)
scatter_2 <- ggscatter(shared_both_lfc_2, x = "N_lfc", y = "NLS_lfc", 
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "pearson",
                       xlab = "DEG LogFoldChange N to Flag", ylab = "DEG LogFoldChange NLS to Flag")+ theme(panel.grid.major = element_line(colour = "gray80", linewidth = 0.5, linetype="dotted")) 

  #could do some type of scatter plot with the GSEA terms for the unique degs, like all the gsea terms from each set of unique degs 

#############
##after GSEA
#Fishers of the shared and unique DEGs
### Start with the shared DEGs (only doing BP)
# BP
input="shared_fisher_GO.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
# 30 GO terms at 10% FDR

resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as shared_fisher_GO_BP as 8x8 portrait


#now do the genes unique to N
input="fishers_N_only_GO.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
# 25 GO terms at 10% FDR

resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)

#saved as N_only_fisher_GO_BP as 8x7 landscape

#now do the genes unique to NLS
input="fishers_NLS_only_GO.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
# 23 GO terms at 10% FDR

resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)

#saved as NLS_only_fisher_GO_BP as 8x7 landscape




