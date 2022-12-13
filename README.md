# KeggPathwaySelector
Software for getting information on pathways related to a list of genes. 
Alzheimer’s Disease is a neurodegenerative disease that causes most dementia cases around the world. Recently, a team at the Stanford University School of Medicine identified 1063 genes associated with the transcription of proteins and RNA relevant in the formation of Neurofibrillary tangles (Otero-Garcia et al. 2022). My project took this list of relevant genes and compared it against the Kyoto Encyclopedia of Genes and Genomes’ R API to determine which metabolic pathways appeared most often, which gene played a role in the most amount of record pathways, and how the number of genes encoding a pathway related to protein and RNA production. 

1.	The "Ribosome - Homo sapiens (human)" pathway occurred most often, with 59 genes related to its functions. 

2.	The gene PRKACB, also known as CAFD2, PKA_C-beta, or PKACB, participated in 70 metabolic pathways.

3.	Participating genes in a pathway show a U-shaped relationship with the number of produced proteins. At a lower number of genes, the relationship is proportional, but as the number of encoding genes grow, the transcribed proteins decrease. The linear regression only has an Rsquared value of 0.265, supporting the U-shaped relationship hypothesis. 

References: 
[1] - M. Otero-Garcia, S. U. Mahajani, D. Wakhloo, W. Tang, Y.-Q. Xue, S. Morabito, J. Pan, J. Oberhauser, A. E. Madira, T. Shakouri, Y. Deng, T. Allison, Z. He, W. E. Lowry, R. Kawaguchi, V. Swarup, and I. Cobos, “Molecular signatures underlying neurofibrillary tangle susceptibility in alzheimer’s disease,” Neuron, vol. 110, no. 18, 2022. 
![image](https://user-images.githubusercontent.com/118032019/207422774-d242fda8-74ae-47a9-a0ee-ca5757718a00.png)

