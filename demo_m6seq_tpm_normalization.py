from load_m6Aseq_counts import load_m6Aseq_counts
import numpy as np
from prettytable import PrettyTable
from Print_Table import Table
from m6APeaks_and_ReadCount import PeaksAndReads
from zero_filter import zero_filter
from TPM_normalization import TPM_norm
from Norm_matrix_load import normtable
from PCA2Dand3D import PCA
import matplotlib.pyplot as pyplot
from kmeancluster_orignal import kmeans
x = PrettyTable()

def main():
    line = "----------------------------------------------------------------------------------------------------------"
    print line
    print "\nLoading the bed file:",

    matrix,gene_order,lib_order = load_m6Aseq_counts('Supplementary4.txt',delimiter =',') # load the txt file

    print "\n\nDimensions of Matrix = ", np.shape(matrix)
    print "Dimesnions of gene_order = ", np.shape(gene_order)
    print "Dimesnions of lib_order = ", np.shape(lib_order)

    print "                                 \nCOMPLETE\n"

    print line
    print "First 4 rows of the bed file (transposed): \n"

    Table(matrix, lib_order, gene_order)  # print 4 rows of the matrix

    print "\n\n", line
    print "Transcripts and Peaks Summary\n "

    PeaksAndReads(matrix) # prints read statistics

    print "\n", line
    print "Removing all zero filled row entries by looking at the no. of peaks identified"

    filtered_matrix, filtered_genes = zero_filter(matrix, gene_order) # filters all zero rows

    print "\n", line
    print "First 4 rows of the zero filtered matrix :\n"

    Table(filtered_matrix, lib_order, filtered_genes) # prints first 4 rows of the filtered matrix

    print "\n", line
    print "Exons and read statistics of filtered matrix :\n "

    PeaksAndReads(filtered_matrix) # prints filtered read statistics

    print line
    print "TPM Normalizing"

    NormMatrix, rangeofexon, rangeofgene = TPM_norm(filtered_matrix, Norm_Type="transcript")      # Calculates the TPM normalized matrix
    normtable(NormMatrix, filtered_genes)       # loads the normalized matrix to csv file

    pyplot.figure()
    pyplot.hist(rangeofgene[:,0],10, normed=True)
    pyplot.xlabel('Transcription length in bases')
    pyplot.ylabel('Relative Frequency')
    pyplot.title('Range of Transcripton Length for all genes')
    pyplot.show()

    pyplot.figure()
    pyplot.hist(rangeofexon,20, normed=True,label=["Human Liver Carcinoma", "Heatshock", "UV","Hepatocyte GF", "Interferons", "Brain"])
    pyplot.xlabel('Number of peaks identified on each gene')
    pyplot.ylabel('Relative Frequency')
    pyplot.title('Range of number of peaks in for all genes')
    pyplot.xlim(0,20000)
    pyplot.show()

    print line
    print "PCA Analysis "

    PCA(NormMatrix,'PCA2DAllgene.png','PCA3DAllgene.png')

    print line

    kmeans("PCA_Clustering_AllGene", 3, "PCA_Allgene_Clustering.png")
    kmeans("PCA_Clustering_AllGene", 4, "PCA_Allgene_Clustering_4.png")

if __name__ == '__main__':
        main()

