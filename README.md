# Epitranscriptome-Analysis

## Abstract:
Epitranscriptome describes the vast number of nucleotide modifications to the RNA molecules in cells. With the discoveries of novel antibody detection of specific modifications, less abundant RNA modifications can be purified and deciphered. These modifications have been shown to play a completely new level of gene expression regulation during translation [1]. The most studied modification is the N6-methyladnesoine (m6A) since it is the most abundant and has been shown to be dynamic. With already peak called bed files from m6A seq of various experiments with human kidney and brain cells, machine learning analysis were performed such as PCA and kmeans clustering. Results from these algorithms suggests that the m6A methylation is unique to each cell type and the landscape of these m6A methylation sites can be changed at various rates depending on the experimental conditions.

## Introduction:
With the advancements in next generation sequencing technologies and innovations, researchers are able to use it to probe and analyze gene expression, epigenetic marks, Epigenomics regulation, spatial temporal arrangement of the genome in nuclease and many more. Recently, one such advancement has been made to selectively map certain mRNA nucleotide modifications in the transcriptome of the cell. Almost more than 100 mRNA modifications have been known to researches for almost half a century, but the only tool available, before the boom of next generation sequencing analysis, was mass spectroscopy. The tool allowed researchers to detect different kinds of modifications on the transcriptome but was too ambiguous to map the modification and this was almost impossible for low abundant modifications, since almost all of tRNA and rRNA molecules (which make up 80% of transcriptome) are heavily methylated and even with advanced isolation techniques for mRNA, a small amount of tRNA or rRNA traces can lead to drastic bias in data (1). This problem was finally over come with two main advancements: Next generation sequencing and modification affinity specific immunoprecipitation. Even though this new method does not allow to precisely map modification to a particular nucleotide, it shows the relative abundance with respect to its genome via peak calling algorithms.

The sudden interest in the study of epitranscriptome modification worldwide is due to the discovery of its reversible and dynamic states throughout cell cycle and one such modifications, which has been shown to be reversible and most abundant, is the N6-methyl-adenosine (m6A). Certain molecular complexes are found to have the function to write (Methyltransferase), erase (demethylase) and read this modification and direct the mRNA to certain functions such as nuclear RNA export, control of protein translation, mRNA decay and splicing (2). Recently it has been specifically shown to regulate stem cell differentiation though heavy methylation of mRNA transcripts, responsible for pluripotent maintenance Transcription Factors. Knockout and depletion of METTL3 and METTL14 (methyltransferases) in cells leads to decreased variability and premature apoptosis of arrested stem cell. Furthermore, the modifications have been shown to be enriched around the stop and start codon of the mRNA molecule (2). The distribution of this m6A modification is extensively studied in Darminine et al (2) and they provide a complete analysis of how these exon specific distributions is conserved in all exons in all cell type since they are found to be heavily methylating the start and the stop codons of each exons. However they are no further analysis on how the overall transcriptome level m6A distribution differs with respect to various experimental conditions.

## Hypothesis:
According to the literature, the m6A peak distribution differs among cell types but does not offer any conclusion on how they might differ among various experimental conditions in the same cell line. Since the m6A methylation regulation is highly dynamic, I hypothesis that the methylation landscape shows considerable difference between conditions even though they are of the same cell line.

## Methods:
Two excel files were used for this experiment, where the main file contained the dataset for all identified m6A peaks in brain cells and hepG2 (Human liver) cell line under various conditions such as Carcinoma cancer, Heat Shocked, UV radiated, hepatocyte Growth Factor and Interferon expressed cells. This file contained approximately 26,000 rows and 49 columns, therefore unused fields were removed by querying in MS Access to reduce computational load. Then the files were exported to csv file formats due to complications with converting to tab delimited txt file formats. These csv files represented bed files, where each gene contained fields such as chromosome, transcription start and end positions on the chromosome, coding or noncoding field, and for each conditions, the number of peaks in each gene and the number of reads for each exons were populated. This file was loaded into the python script using load_m6A_seq.py and the second file contained genes that were only differentially expressed between two conditions at a time.

Part one of the analysis pipeline, which considers all genes, performs zero filtering, TPM normalization, Principle component analysis, clustering and plots a 2D plot of PC1 vs PC2 with clustering and a 3D plot of PC1, PC2 and PC3. For the second part of the pipeline, just the differentially expressed genes between experiments were considered. The differentially expressed genes were extracted from the all gene matrix by first removing any repeated genes that were differentially expressed in more than one comparison. Then these genes were used to extract the row information for each of these genes for the first matrix. Then it performs the above machine learning analysis with the TPM normalized matrix. The third part of the analysis considers just the genes that were not differentially expressed and performs the above-mentioned analysis. To obtain the genes that were not differentially expressed, the differentially expressed file from the second part of the analysis is used to remove these rows form the main all gene matrix.

## Results: 
The outputs from Part one of the m6A analysis pipeline contain the 2D and 3D PCA plots for all six experimental conditions. Starting with Human liver carcinoma cancer cells (HEPG2 cell line), the next five conditions included HEPG2 cells that were heat shocked by incubating the cells for 30 minutes at 43 degree Celsius, UV radiated at 254nm, incubating the cells in IFN Interferon solutions overnight to heighten their anti-viral defense mechanisms, incubated in Hepatocyte Growth factor solution overnight, which plays a major role in embryonic organ development and then finally the last conditions is with Human brain cells. The PCA plots in general without the clustering shows considerable relevance to the biology of each sample conditions. The m6A peaks variance in brain cells is projected apart from the rest of the conditions, which suggest that the m6A peak distribution might be completely different from the rest of the experimental conditions being analyzed. The same trend can also be observed in the 3D plot for all genes where the brain cell’s variance is plotted opposite to all the other conditions. This might also suggest that the m6A peak has a completely different landscape of distribution among different cell types and tissue types. This correlates with the literature that concludes that even though these distributions follow the same pattern of methylation on each exon such as near the start and end codons of the exons, the genes they occupy change considerably between cell types.

We can also observe that since all the other conditions are based on the HEPG2 cell line but undergo different conditioning during experiments, they all appear in the vicinity of each other. But the Interferon experimental condition appears to be clustered separately form the other HEPG2 experimental cells. This might be due to the fact that after incubating the IFN, the cell’s anti-viral response is heightened which change the m6A expression profile and projects it apart from others. In PCA clustering with 4 cluster groups, we can see that the large clustering group form 3 cluster analyzes has divided into two, which suggest that UV radiation DNA damage does not change the m6A methylation expression of the affected cell and this layer of regulation is unaffected by this type of damage. But the other two conditions, growth factor incubation and heat-shocked cells are grouped together. These two conditions seem to be closing towards interferon, which might suggest that prolonged exposure to either growth factor or higher temperature can slowly change the m6A methylation distribution.

In the second part of the analysis pipeline, the results suggest a more refined observation of the differential landscape of the m6A methylation variance between experimental conditions since we are focusing the machine learning algorithm to only work with differentially expressed exons rows. Compared to the all gene analysis observation, the variance between the brain cell’s m6a methylation variance and the HEPG2 cell lines are more apparent since they appear more distant from each other. In the four clustering PCA plot, we can observe that Human Liver carcinoma cells, UV radiated and Heat shocked cells appear to have the same variance in their m6A transcriptome expression. But the interferon and growth factor incubated cells are grouped separately on their own like the brain cell, which suggest that their differentially expressed m6A variance are very different from the HEPG2 cell line, but the variance in other non-differentially expressed exons are similar enough between growth factor induced cells and HEPG2 cell lines, that the overall landscape of m6A methylation is comparable to the original cell line. The 2D PCA plot and the 3D plot of the non-differentially m6A peak expressed exons are very similar to the normal expression patters, this might be due to the fact that the number of differentially expressed exons with m6A peak is approximately 600 exons but this is overcome by the similarly expressed exons which count for approximately 25,000 exons.

## Conclusion:
From the results of the m6A methylation analysis pipeline, it seems that the difference in m6A distribution is highly differential among cell types and the distribution can also be changed rapidly or slowly on the same cell line depending on the experimental conditions. This conclusion defeats the hypotheses where it I suggested that the distribution is highly sensitive to changes in the cell’s environment but the results show that these sites are most conserved during physical stress but can transform rapidly to chemical changes.

## References: 
1) Geula, S., S. Moshitch-Moshkovitz, D. Dominissini, A. A. Mansour, N. Kol, M. Salmon-Divon, V. Hershkovitz, E. Peer, N. Mor, Y. S. Manor, M. S. Ben-Haim, E. Eyal, S. Yunger, Y. Pinto, D. A. Jaitin, S. Viukov, Y. Rais, V. Krupalnik, E. Chomsky, M. Zerbib, I. Maza, Y. Rechavi, R. Massarwa, S. Hanna, I. Amit, E. Y. Levanon, N. Amariglio, N. Stern-Ginossar, N. Novershtern, G. Rechavi, and J. H. Hanna. "M6A mRNA methylation facilitates resolution of naive pluripotency toward differentiation." Science 347.6225 (2015): 1002-006. Web.
2) Dominissini, Dan, Sharon Moshitch-Moshkovitz, Schraga Schwartz, Mali Salmon-Divon, Lior Ungar, Sivan Osenberg, Karen Cesarkas, Jasmine Jacob-Hirsch, Ninette Amariglio, Martin Kupiec, Rotem Sorek, and Gideon Rechavi. "Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq." Nature 485.7397 (2012): 201-06. Web.


![alt text](https://github.com/Prdeeepg/Epitranscriptome-Analysis/blob/master/PCA2DAllgene.png)
![alt text](https://github.com/Prdeeepg/Epitranscriptome-Analysis/blob/master/PCA3DAllgene.png)


