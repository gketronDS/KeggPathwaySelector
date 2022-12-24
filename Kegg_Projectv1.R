## ----setup, echo=FALSE--------------------------------------------------------
library(knitr)
options(width=80)
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(plotly)
library(KEGGREST)
library(png)
print("est runtime: 8min 30s")
## ----Accessing List of Genes--------------------------------------------------------
## Genes from : https://reader.elsevier.com/reader/sd/pii/S0896627322006006?token=0E7C634165845D34F9DB21F406C7E7BE7E48451562824D62061DFBF3B95EFA910A9550326C080AD66F23F890CB9E46DD&originRegion=us-east-1&originCreation=20221203005723
##Read from Excel
REST_genes<-read_excel("/Users/gabrielketron/Downloads/1-s2.0-S0896627322006006-mmc3.xlsx",sheet="REST (shared)")
excel_sheets("/Users/gabrielketron/Downloads/1-s2.0-S0896627322006006-mmc3.xlsx")
##Only take genes the study identified that impact Alzheimer's Disease.
Deg_REST<-REST_genes[(REST_genes$type=="DEG"),]
Deg_REST<- Deg_REST %>% 
  select(gene_names)
Deg_LIST<- list(Deg_REST)

## ----Get HSA values for Genes--------------------------------------------------------
##Turn list into characters for API to search
Deg_GENE <- Deg_LIST[[1]][[1]]
out_GENE<- list()
genehsa = list()
##Search for everything about the gene in KEGG
#gene_test <- keggFind("genes", Deg_GENE[[140]])
##Format result to get hsa values
#gene_test <- as.data.frame(gene_test) # Duplicate example data
#gene_test$row_names <- row.names(gene_test) # Apply row.names function
#genehsa <- gene_test$row_names[[1]]
for (x in 1:length(Deg_GENE)) {
  print(x / 293*100)
  print("% Finished - Finding HSA Values")
  gene_test <- keggFind("genes", Deg_GENE[[x]])
  ##Format result to get hsa values
  gene_test <- as.data.frame(gene_test) # Duplicate example data
  if(dim(gene_test)[1] != 0){
    gene_test$row_names <- row.names(gene_test) # Apply row.names function
    for (i in 1:length(gene_test[[2]])){
      strwrk <- gene_test[[2]][i]
      strwrk<-strsplit(strwrk, split = ":")
      if (strwrk[[1]][1]=="hsa"){
        out_GENE[[length(out_GENE)+1]]<-gene_test[[2]][i]
      }
    } 
    #for (y in 1:length(genehsa)){
     # out_GENE[[length(out_GENE)+1]] <- genehsa[[y]]
    #}
  }
  else{
    genehsa <- "hsa:0"
    out_GENE[[length(out_GENE)+1]] <- genehsa
  }
}
## add out-GENE as a column in Deg_REST and remove hsa: 0
out_GENE <- unique(out_GENE)
in_GENE <- out_GENE[out_GENE != "hsa:0"]
print("HSA Values Recieved")
## ----Get Pathways from Genes--------------------------------------------------------
## This section allows us to assign genes to pathways for later comparison. 
in_GENE3 <- t(in_GENE) 
pathquery2<-as.data.frame(matrix(in_GENE3, ncol=664, byrow = TRUE))
track_path<- list()
gene_track_id = list()
for (x in 1:length(pathquery2)){
  print(x / 664*100)
  print("% Finished - Finding Path Values")
  hold<-keggLink("pathway",in_GENE[x])
  for (y in 1:length(hold)){
    if (length(hold[y])!=0){
      if (!is.na(hold[y])){
        track_path[[length(track_path)+1]]<-hold[[y]]
        gene_track_id[[length(gene_track_id)+1]]<-in_GENE[x]
      }
    }
  }
}
track_path_1 <- as.data.frame(matrix(track_path))
gene_track_set <- as.data.frame(matrix(gene_track_id))
## Path_occurance counts how often each pathway comes up, as column "n".
path_occurance<-track_path_1 %>%
  count(V1)
#pathwayIN <- list() #Old attempt at track_path
pathwayID <- list() #Name of the pathway as a KEGG data structure.
pathwayNAME <- list() #Name in text
pathwayDRUGS <- list() #List of drugs that effect the path
pathwayDRUGL <- list() # # of drugs which effect the path
pathwayCOMPS <- list() # Compounds that effect the path
pathwayCOMPL <- list() # number of compounds that effect the path
pathwayREL <- list() # List of related paths
pathwayRELL <- list() # Number of related paths
pathwayREFL <- list() # Number of references in KEGG
pathwayPROT <- list() # List of proteins the pathway makes
pathwayGENEL <- list() # number of unique proteins the pathway makes
for (x in 1:length(path_occurance[[1]])){
  print(x / 273*100)
  print("% Finished - Finding Path Data")
  pathway_test <- keggGet(path_occurance[[1]][x])
  pathwayID[[length(pathwayID)+1]] <- pathway_test[[1]][["PATHWAY_MAP"]]
  pathwayNAME[[length(pathwayNAME)+1]] <-pathway_test[[1]][["NAME"]]
  #pathwayDRUGS[[length(pathwayDRUGS)+1]] <- pathway_test[[1]][["DRUG"]]
  pathwayDRUGL[[length(pathwayDRUGL)+1]] <- length(pathway_test[[1]][["DRUG"]])
  ##pathwayCOMPS[[length(pathwayCOMPS)+1]] <- pathway_test[[1]][["COMPOUND"]]
  pathwayCOMPL[[length(pathwayCOMPL)+1]] <- length(pathway_test[[1]][["COMPOUND"]])
  #pathwayREL[[length(pathwayREL)+1]] <- pathway_test[[1]][["REL_PATHWAY"]]
  pathwayRELL[[length(pathwayRELL)+1]] <- length(pathway_test[[1]][["REL_PATHWAY"]])
  pathwayREFL[[length(pathwayREFL)+1]] <- length(pathway_test[[1]][["REFERENCE"]])
  pathwayPROT[[length(pathwayPROT)+1]]<-pathway_test[[1]][["GENE"]]
  pathwayGENEL[[length(pathwayGENEL)+1]]<-length(pathway_test[[1]][["GENE"]])
}
pathWAY <- mutate(path_occurance, ID = pathwayID, NAME = pathwayNAME, 
                  DRUGL = pathwayDRUGL, COMPL = pathwayCOMPL, 
                  RELL = pathwayRELL, REF = pathwayREFL, 
                  PROT = pathwayGENEL)
pathWAY <- as.data.frame(lapply(pathWAY, unlist))
pathWAY <- pathWAY[order(-pathWAY$n),]
## ----Get Genes from Genes--------------------------------------------------------
in_GENE2 <- t(in_GENE) 
gene_query<-as.data.frame(matrix(in_GENE2, ncol=664, byrow = TRUE))
gene_NAME <- list() # Name of gene
gene_CDS <- list() #Hsa value of gene as a number
gene_path <- list() # list of pathways the gene is a part of
gene_PATHS <- list()# # of pathways the gene is a part of
gene_POS <- list() # gene position in the genome
gene_AA <- list() # amino acid sequence the gene encodes
gene_AL <- list() # Length of amino acid sequence
gene_DA <- list() # DNA sequence of the gene
gene_DL <- list() # Length of dna sequence
gene_ORTH <- list() #orthology reference number for later use
for (x in 1:length(gene_query)){
  print(x/664*100)
  print("% Finished - Finding Gene Data")
  gene_test2 <- keggGet(gene_query[[x]])
  gene_NAME[[length(gene_NAME)+1]] <- gene_test2[[1]][["SYMBOL"]]
  gene_CDS[[length(gene_CDS)+1]] <- gene_test2[[1]][["ENTRY"]]
  #gene_path[[length(gene_path)+1]] <- gene_test2[[1]][["PATHWAY"]]
  gene_PATHS[[length(gene_PATHS)+1]] <- length(gene_test2[[1]][["PATHWAY"]])
  strwrk<-gene_test2[[1]][["POSITION"]]
  strwrk<-strsplit(strwrk, split = ":")
  if (strwrk[[1]][1]=="X"){
    strwrk<- as.numeric(0)
  }
  strwrk<-as.numeric(strwrk[[1]][1])
  gene_POS[[length(gene_POS)+1]] <- strwrk
  gene_AA[[length(gene_AA)+1]] <- gene_test2[[1]][["AASEQ"]]
  #gene_AL[[length(gene_AL)+1]] <- gene_test2[[1]][["AASEQ"]]@ranges@width
  gene_DA[[length(gene_DA)+1]] <- gene_test2[[1]][["NTSEQ"]]
  #gene_DL[[length(gene_DL)+1]] <- gene_test2[[1]][["NTSEQ"]]@ranges@width
  gene_ORTH[[length(gene_ORTH)+1]] <-gene_test2[[1]][["ORTHOLOGY"]]
}
for (x in 1:length(gene_AA)){
  gene_AL[[length(gene_AL)+1]] <- gene_AA[[x]]@ranges@width
} 
for (x in 1:length(gene_DA)){
  gene_DL[[length(gene_DL)+1]] <- gene_DA[[x]]@ranges@width
}
gene_query <- gene_query<-as.data.frame(matrix(gene_query, nrow=664))
geneWAY <- mutate(gene_query, NAME = gene_NAME, CDS = gene_CDS, PATHS = gene_PATHS,
                  POS = gene_POS, DL = gene_DL)
geneWAY <- as.data.frame(lapply(geneWAY, unlist))
geneWAY <- geneWAY[order(-geneWAY$PATHS),]

## ----Data Visualizations and Trends--------------------------------------------------------
#geneWAY
## Relationship Between gene CDS and gene position
ggplot(geneWAY, aes(x=CDS,y=POS,color=POS)) + 
  geom_smooth()+
  geom_point(size=0.5)+
  ggtitle("Relationship Between gene \n CDS and gene position")
##Relationship Between CDS and Gene Paths
ggplot(geneWAY, aes(x=CDS,y=PATHS,color=POS)) + 
  geom_smooth()+
  geom_point(size=0.5)+
  ggtitle("Relationship Between CDS and Gene Paths")
##Relationship Between CDS and DNA Length
ggplot(geneWAY, aes(x=CDS,y=DL,color=POS)) +
  geom_smooth()+
  geom_point(size=0.5) +
  ggtitle("Relationship Between CDS and DNA Length")
##Relationship between Gene postion and gene paths
ggplot(geneWAY, aes(x=POS,y=PATHS,color=POS)) + 
  geom_smooth()+
  geom_point(size=0.5) +
  ggtitle("Relationship between Gene postion and gene paths")
##relationship between Gene position and DNA length
ggplot(geneWAY, aes(x=POS,y=DL,color=POS)) + 
  geom_smooth()+
  geom_point(size=0.5) +
  ggtitle("relationship between Gene position and DNA length")
##relationship between Gene paths and DNA length
ggplot(geneWAY, aes(x=PATHS,y=DL,color=POS)) + 
  geom_smooth()+
  geom_point(size=0.5) +
  ggtitle("relationship between Gene paths and DNA length")

#pathWAY
##Relationship between N and Related paths
ggplot(pathWAY, aes(x=n,y=RELL)) + 
  geom_smooth()+
  geom_point(size=0.5)+
  ggtitle("Relationship between contributing genes and related paths")
#relationship between N and number of proteins
lm_eqn <- function(df){
  m <- lm(PROT ~ n, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
ggplot(pathWAY, aes(x=n,y=PROT)) + 
  geom_smooth()+
  geom_smooth(method = "lm",se=FALSE, color = "Red" )+
  geom_text(x = 30, y = 1000, label = lm_eqn(pathWAY), parse = TRUE)+
  geom_point(size=0.5)+
  xlab("Number of Genes encoding the Pathway") +
  ylab("Number of Transcripted Proteins")
  ggtitle("Relationship Between Encoding Genes and Proteins") +
  theme_bw()
#Relationship between N and references
ggplot(pathWAY, aes(x=n,y=REF)) + 
  geom_smooth()+
  geom_point(size=0.5)+
  ggtitle("Relationship between contributing genes and references")
#Relationship Drugs and compounds
ggplot(pathWAY, aes(x=DRUGL,y=COMPL,color=n)) + 
  geom_smooth()+
  geom_point(size=0.5)+
  ggtitle("Relationship Drugs and compounds")
#relationship between Drugs and proteins
ggplot(pathWAY, aes(x=DRUGL,y=PROT,color=n)) + 
  geom_smooth()+
  geom_point(size=0.5)+
  ggtitle("relationship between Drugs and proteins")
#relationship between compounds and proteins
ggplot(pathWAY, aes(x=COMPL,y=PROT,color=n)) + 
  geom_smooth()+
  geom_point(size=0.5)+
  ggtitle("relationship between compounds and proteins")
#Final Results
#Picture of pathway with most genes w # of genes
strwrk <- pathWAY[[1]][[9]]
strwrk<-strsplit(strwrk, split = ":")
temp <- strwrk[[1]][2]
png<-keggGet(temp,"image")
t <- tempfile()
writePNG(png,t)
if (interactive()) browseURL(t)
print(pathWAY$NAME[[1]])
print(pathWAY$n[[1]])
print("Associated Genes")
#gene with most pathways w number of pathways
print(geneWAY$NAME[[1]])
print(geneWAY$PATHS[[1]])
print("Associated Pathways")
#Total number of unique pathways
print("289 Unique Genes ")
print(length(path_occurance[[1]]))
print("Unique Pathways")