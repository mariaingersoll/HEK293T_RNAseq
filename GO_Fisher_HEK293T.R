#GO analysis with Fishers exact (1s and 0s)

setwd("/Users/mariaingersoll/Desktop/BU_Research/Davies-Gilmore/Mengrui_HEK293T/HEK293T_RNAseq/Data_Files")

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
#
resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as shared_fisher_GO_BP.pdf as portrait with 13x10 dimensions in Figures

#N only
# BP
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
#
resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)


