####GO analysis for MENGRUI HEK293T cells (Flag/control, N, NLS)
#GO input files were created from HEK293T_Mengrui.Rmd (ending around line 600)
#all _GO_24.csv files need to be in unix format to load into the GO commands, so open and save them in BBedit

#also, make sure all the GO MWU files required for the script are in your working directory

setwd("/Users/mariaingersoll/Desktop/BU_Research/Davies-Gilmore/Mengrui_HEK293T/HEK293T_RNAseq/Data_Files")

#### N to Flag

# CC
input="NFlag_GO_24.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#3 terms at 10% FDR
resultsCC=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as CC_NFlag_022224.pdf as landscape with 10x3 dimensions in Figures
resultsCC

#BP
input="NFlag_GO_24.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#3 terms at 10% FDR
resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as BP_NFlag_022224.pdf as landscape with 10x3 dimensions in Figures
resultsBP

#MF
input="NFlag_GO_24.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#3 terms at 10% FDR
resultsMF=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as MF_NFlag_022224.pdf as landscape with 10x3 dimensions in Figures
resultsMF

##########
#### NLS to Flag

# CC
input="NLSFlag_GO_24.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#8 terms at 10% FDR
resultsCC=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as CC_NLSFlag_022224.pdf as landscape with 10x3 dimensions in Figures
resultsCC

#BP
input="NLSFlag_GO_24.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#3 terms at 10% FDR
resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as BP_NLSFlag_022224.pdf as landscape with 10x3 dimensions in Figures
resultsBP

#MF
input="NLSFlag_GO_24.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#3 terms at 10% FDR
resultsMF=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as MF_NLSFlag_022224.pdf as landscape with 10x3 dimensions in Figures
resultsMF

#############
#### N to NLS

# CC
input="NNLS_GO_24.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#0 terms at 10% FDR


#BP
input="NNLS_GO_24.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#0 terms at 10% FDR

#MF
input="NNLS_GO_24.csv"
goAnnotations="GRCh38_iso2go.tab"
goDatabase="go.obo"
goDivision="MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#3 terms at 10% FDR
resultsMF=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as MF_NNLS_022224.pdf as landscape with 10x3 dimensions in Figures but not coloring the terms
resultsMF
