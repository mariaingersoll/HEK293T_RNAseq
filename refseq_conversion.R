#get the list of the shared genes (no lfc cutoff) in REFSEQ_new format
head(shared_N_gsea_new)
nrow(shared_N_gsea_new)
#left_join it to have the numbers after the decimal point
sN_new1 <- shared_N_gsea_new %>%
  left_join(cds_seq2seq)
head(sN_new1)
#left_join it to have the "ISO" names
sN_new1 <- sN_new1 %>%
  left_join(cds_seq2iso_1)
head(sN_new1)
nrow(sN_new1)
#save that as a file
write.csv(sN_new1, file="../sN_new1.csv", quote=F)
length(unique(sN_new1$REFSEQ_new))
#1192
length(unique(sN_new1$REFSEQ))
#1249
length(unique(sN_new1$ISO))
#1249

#N-only genes (no lfc cutoff)
head(N_only_gsea_new)
nrow(N_only_gsea_new)
#161
oN_new1 <- N_only_gsea_new %>%
  left_join(cds_seq2seq)
head(oN_new1)
oN_new1 <- oN_new1 %>%
  left_join(cds_seq2iso_1)
head(oN_new1)
nrow(oN_new1)
#169
#save that as a file
write.csv(oN_new1, file="../oN_new1.csv", quote=F)
length(unique(oN_new1$REFSEQ_new))
#161
length(unique(oN_new1$REFSEQ))
#169
length(unique(oN_new1$ISO))
#169

#NLS-only genes (no lfc cutoff)
head(NLS_only_gsea_new)
nrow(NLS_only_gsea_new)
#509
oNLS_new1 <- NLS_only_gsea_new %>%
  left_join(cds_seq2seq)
head(oNLS_new1)
oNLS_new1 <- oNLS_new1 %>%
  left_join(cds_seq2iso_1)
head(oNLS_new1)
nrow(oNLS_new1)
#524
#save that as a file
write.csv(oNLS_new1, file="../oNLS_new1.csv", quote=F)
length(unique(oNLS_new1$REFSEQ_new))
#509
length(unique(oNLS_new1$REFSEQ))
#524
length(unique(oNLS_new1$ISO))
#524


######## 1.5
head(shared_N_gsea_1.5_new)
nrow(shared_N_gsea_1.5_new)
#217
sN_1.5_new1 <- shared_N_gsea_1.5_new %>%
  left_join(cds_seq2seq)
head(sN_1.5_new1)
#left_join it to have the "ISO" names
sN_1.5_new1 <- sN_1.5_new1 %>%
  left_join(cds_seq2iso_1)
head(sN_1.5_new1)
nrow(sN_1.5_new1)
#226
#save that as a file
write.csv(sN_1.5_new1, file="../sN_1.5_new1.csv", quote=F)
length(unique(sN_1.5_new1$REFSEQ_new))
#217
length(unique(sN_1.5_new1$REFSEQ))
#226
length(unique(sN_1.5_new1$ISO))
#226

#make a variable for the genes from the shared terms that overlap with the Enrichr rela kd
relakd <- c("SRP19", "POP7", "RTCB", "ACTB", "TOMM22", "PSMD8", "CDC20", "CHCHD2", "PSMD7", "GABARAPL2", "SUPT16H", "ATP1B3", "COX6B1", "RRAGA", "MTHFD1", "UQCRC1", "SRSF2", "SUCLG1", "UQCRC2", "SNRPF", "SELENOT", "SNRPB", "ANP32A", "PSMD14", "DDX1", "TXN", "PRDX3", "PRDX2", "PRDX4", "PCBP1", "CYC1", "HSPA9", "XRCC5", "DHCR24", "TBCA", "PRDX6", "LSM3", "AIMP2", "CNIH1", "UQCRQ", "BAMBI", "CYCS", "RNPS1", "COX7B", "LAMC1", "HMGB1", "YBX1", "RBM3", "MRPL3", "YWHAQ", "C1QBP", "DHX15", "UQCRFS1", "ELOB", "RAB8A", "COX8A", "PARP1", "HSBP1", "ZFR", "TUBG1", "MRPS7", "DDOST", "APRT", "PSMA2", "DAD1", "NDUFS6", "CKS2", "PABPC1", "SLC25A5", "ARF5", "VAMP3", "NDUFB8", "SF3B3", "ECHS1", "NDUFB6", "UQCR10", "COX5A", "SRM", "PSMB7", "ATP5F1B", "PSMB3", "POLR2C", "PSMB1", "EMC8", "NDUFA9", "TIMM8B", "GINS2", "SF3A1", "MARCKSL1", "TUBB4B", "HSPE1", "EIF3K", "SCD", "PSMC2", "DRG1", "ATP5PO", "RBMX")
#read in csv where you've added the Entrez symbols
sN_new2 <- read.csv(file="../sN_new1_symbol.csv", row.names = 1)
head(sN_new2)
#extract the relakd terms
sN_rel <- sN_new2 %>%
  subset(SYMBOL %in% relakd)
nrow(sN_rel)
#101

#now get the lfc values for the shared NLS
head(shared_NLS_gsea_new)
nrow(shared_NLS_gsea_new)
#1192
sNLS_new1 <- shared_NLS_gsea_new %>%
  left_join(cds_seq2seq)
head(sNLS_new1)
sNLS_new1 <- sNLS_new1 %>%
  left_join(cds_seq2iso_1)
head(sNLS_new1)
nrow(sNLS_new1)
#1249
#save that as a file
write.csv(sNLS_new1, file="../sNLS_new1.csv", quote=F)
length(unique(sNLS_new1$REFSEQ_new))
#1192
length(unique(sNLS_new1$REFSEQ))
#1249
length(unique(sNLS_new1$ISO))
#1249

#now add the sNLS lfc values sN_rel
sN_NLS_rel <- sN_rel %>%
  left_join(sNLS_new1)
nrow(sN_NLS_rel)
#101

#in order to make a heatmap, I will want the expression value across all the samples
#not just the average lfc per treatment
#and I might as well include the flag samples as well
rld_new_named <- read.csv(file="../mingers_rldall_new_named_24.csv", row.names = 1)

#make each symbol unique
sN_NLS_rel_un <- sN_NLS_rel %>%
  mutate(SYMBOL = make.names(SYMBOL, unique = TRUE))

sN_NLS_rel_rld <- sN_NLS_rel_un %>%
  left_join(rld_new_named) %>%
  select(-lfc_NFlag, -lfc_NLSFlag, -REFSEQ_new, -REFSEQ, -ISO) %>%
  column_to_rownames(var = "SYMBOL")

#make a heatmap
#graypurp <- c("gray50", "purple", "gray50", "purple", "gray50", "purple", "gray50", "purple", "gray50", "purple")

sN_rel_means = apply(sN_NLS_rel_rld, 1, mean)
explc_sN_rel = sN_NLS_rel_rld-sN_rel_means

sN_rel_hm <- pheatmap(as.matrix(explc_sN_rel),
         color=rev(colorRampPalette(brewer.pal(n=11, name="RdBu"))(20)),
         cluster_rows = T, show_rownames = T, fontsize_row=8, fontsize_col=10, treeheight_row = 10, treeheight_col = 10,
         cluster_cols = T, show_colnames = T, angle_col= "315", cellheight = 8, cellwidth = 14,
         breaks=c(seq(from=-1,to=1,by=0.1)),
         scale="none",
         legend=T)

##load in the entire TF_Perturbations_Followed_by_Expression dataset
TF_pert <- read.delim("../TF_Perturbations_Followed_by_Expression_table.txt")
TF_pert_top10 <- TF_pert[1:10,]
##the whole TF_pert_top10 could be a sup table

TF_pert_top10 <- TF_pert_top10 %>%
  mutate(Term = gsub("RELA KD HUMAN GSE1676 CREEDSID GENE 1787 UP","RELA KD HUMAN GSE1676", Term)) %>%
  mutate(Term = gsub("HNF1B OE HUMAN GSE3308 CREEDSID GENE 1750 UP","HNF1B OE HUMAN GSE3308", Term)) %>%
  mutate(Term = gsub("TP53 KD HUMAN GSE1676 CREEDSID GENE 1784 UP","TP53 KD HUMAN GSE1676", Term)) %>%
  mutate(Term = gsub("HNF1A OE HUMAN GSE3308 CREEDSID GENE 279 DOWN","HNF1A OE HUMAN GSE3308", Term)) %>%
  mutate(Term = gsub("DOT1L INHIBITION HUMAN GSE29828 CREEDSID GENE 1239 DOWN","DOT1L INHIBITION HUMAN GSE29828", Term)) %>%
  mutate(Term = gsub("NFKB1 INACTIVATION HUMAN GSE20667 CREEDSID GENE 2577 DOWN","NFKB1 INACTIVATION HUMAN GSE20667.2577", Term)) %>%
  mutate(Term = gsub("SOX4 KD HUMAN GSE4225 CREEDSID GENE 64 DOWN","SOX4 KD HUMAN GSE4225", Term)) %>%
  mutate(Term = gsub("MYB OE HUMAN GSE2815 CREEDSID GENE 1740 UP","MYB OE HUMAN GSE2815", Term)) %>%
  mutate(Term = gsub("NFKB1 INACTIVATION HUMAN GSE20667 CREEDSID GENE 2576 DOWN","NFKB1 INACTIVATION HUMAN GSE20667.2576", Term)) %>%
  mutate(Term = gsub("MYC KNOCKIN MOUSE GSE39443 CREEDSID GENE 1847 UP","MYC KNOCKIN MOUSE GSE39443", Term))

TF_pert_top10$Term <- factor(TF_pert_top10$Term, levels = c("RELA KD HUMAN GSE1676", "HNF1B OE HUMAN GSE3308", "TP53 KD HUMAN GSE1676", "HNF1A OE HUMAN GSE3308", "DOT1L INHIBITION HUMAN GSE29828", "NFKB1 INACTIVATION HUMAN GSE20667.2577", "SOX4 KD HUMAN GSE4225", "MYB OE HUMAN GSE2815", "NFKB1 INACTIVATION HUMAN GSE20667.2576", "MYC KNOCKIN MOUSE GSE39443"))

#make a bargraph
tf_pert_bar <- ggplot(data=TF_pert_top10, aes(x=Term, y=Combined.Score))+
  geom_bar(position="stack", stat="identity", color="black", alpha=0.6)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(color="black", fill="white"),
        axis.text = element_text(colour="black", size=10), 
        axis.title.x = ggtext::element_markdown(face="bold", size=12),
        axis.title.y = ggtext::element_markdown(face="bold",size=12),
        panel.grid.major = element_line(colour = "gray80", linewidth = 0.5, linetype="dotted")) 
#saved as tf_pert_bar as portrait 5x4
