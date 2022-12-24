# KeggPathwaySelector
Software for getting information on pathways related to a list of genes. Full paper can be found as Identifying_unexpected... in the file list.

BME 295 Data Science  12/06/22
Gabriel Ketron  55928130

Final Project
[KEGG Genes and Pathways]

	Alzheimer’s disease is a neurodegenerative disease that causes most dementia cases around the world. Ongoing investigation into the mechanisms 
	of Alzheimer's disease cover a spectrum of genomic, metabolic, and transcriptomic data. The aim of this paper is to demonstrate R software 
	developed to congregate relevant multiomic data to identify unexpected metabolic pathways which may impact Alzheimer's Disease.
This software queries a list of genes in the Kyoto Encyclopedia of Genes and Genomes’ R API to determine which metabolic pathways appeared most often, which gene played a role in the most amount of record pathways, and how the number of genes encoding a pathway related to protein and RNA production. By applying this software to Alzheimer's disease-related genes, this paper identifies the metabolic pathway that these genes contribute to the most, and by extension, which pathways may play a role in metabolic pathology.

Key Findings

1.	The "Ribosome - Homo sapiens (human)" pathway occurred most often, with 59 genes related to its functions. Covid-19 and two synapse metabolic pathways were also identified as possible correlations.

2.	The gene PRKACB, also known as CAFD2, PKA_C-beta, or PKACB, participated in 70 metabolic pathways.

3.	Participating genes in a pathway show a U-shaped relationship with the number of produced proteins. At a lower number of genes, the relationship is proportional, but as the number of encoding genes grow, the transcribed proteins decrease. The linear regression only has an Rsquared value of 0.265, supporting the U-shaped relationship hypothesis. 

References: 
[1] - M. Otero-Garcia, S. U. Mahajani, D. Wakhloo, W. Tang, Y.-Q. Xue, S. Morabito, J. Pan, J. Oberhauser, A. E. Madira, T. Shakouri, Y. Deng, T. Allison, Z. He, W. E. Lowry, R. Kawaguchi, V. Swarup, and I. Cobos, “Molecular signatures underlying neurofibrillary tangle susceptibility in alzheimer’s disease,” Neuron, vol. 110, no. 18, 2022. 

![image](https://user-images.githubusercontent.com/118032019/207422916-a432d1a9-a32c-45a9-b481-4c0d8e7a81a1.png)
