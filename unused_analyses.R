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
