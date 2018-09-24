from load_m6Aseq_counts import load_m6Aseq_counts
import numpy as np
from prettytable import PrettyTable
from Print_Table import Table
from m6APeaks_and_ReadCount import PeaksAndReads
from zero_filter import zero_filter
from TPM_normalization import TPM_norm
from load_m6Aseq_Differential_counts import load_m6Aseq_differential_counts
from PCA2Dand3D import PCA
import matplotlib.pyplot as pyplot
from kmeancluster_orignal import kmeans
x = PrettyTable()


def main():

    print "----------------------------------------------------------------------------------------------------------"
    print "\nLoading the bed file:",

    diffmatrix,diffgene_order,difflib_order = load_m6Aseq_differential_counts('differential.txt',delimiter =',')
    matrix,gene_order,lib_order = load_m6Aseq_counts('Supplementary4.txt',delimiter =',')
    print np.shape(diffmatrix)
    print "     COMPLETE\n"
    print "----------------------------------------------------------------------------------------------------------"
    print "importing only differential expressed gene from matrix : \n"

    # Correcting for repeated transcripts in differential gene order.
    print np.shape(diffgene_order)
    myset = set(diffgene_order)
    diffgene_order = list(myset)
    print np.shape(diffgene_order)

    diffgene_order = set(gene_order).symmetric_difference(myset)
    diffgene_order = list(diffgene_order)
    print np.shape(diffgene_order)

    reads = np.arange(np.shape(diffgene_order)[0]*49).reshape(np.shape(diffgene_order)[0],49).astype('object')
    gene = np.arange(np.shape(diffgene_order)[0]).reshape(np.shape(diffgene_order)[0],1).astype('object')
    q = 0
    for i in range(0, np.shape(diffgene_order)[0]):
        for j in range(0, np.shape(gene_order)[0]):
            if diffgene_order[i] == gene_order[j]:
                gene[i] = gene_order[j]
                reads[i] = matrix[j,]

    print "reads = ", np.shape(reads)
    print reads[0:10,]
    print "gene order = ", np.shape(gene)

    print "----------------------------------------------------------------------------------------------------------"
    print "First 4 rows of the bed file (transposed): \n"

    #Table(matrix, lib_order, gene_order)
    Table (reads, lib_order, gene)

    print "\n\n--------------------------------------------------------------------------------------------------------"
    print "Exons and read statistics\n "

    PeaksAndReads(reads)

    print "\n---------------------------------------------------------------------------------------------------------"
    print "Removing all zero filled row entries by looking at the no. of peaks identified"

    #filtered_matrix, filtered_genes = zero_filter(matrix, gene_order)
    filtered_matrix, filtered_genes = zero_filter(reads, gene)

    print "\n----------------------------------------------------------------------------------------------------------"
    print "First 4 rows of the zero filtered matrix :\n"

    #Table(filtered_matrix, lib_order, filtered_genes)

    print "\n----------------------------------------------------------------------------------------------------------"
    print "Exons and read statistics of filtered matrix :\n "

    #PeaksAndReads(filtered_matrix)

    print "----------------------------------------------------------------------------------------------------------"
    print "TPM Normalizing"

    NormMatrix,rangeofexon, rangeofgene  = TPM_norm(filtered_matrix, Norm_Type="transcript")      # Calculates the TPM normalized matrix

    #normtable(NormMatrix, filtered_genes)       # loads the normalized matrix to csv file

    pyplot.figure()
    pyplot.hist(rangeofgene[:,0],10, normed=True)
    pyplot.xlabel('Counts')
    pyplot.ylabel('Relative Frequency')
    pyplot.title('Probability distribution for counts in ')
    pyplot.show()

    pyplot.figure()
    pyplot.hist(rangeofexon,20, normed=True)
    pyplot.xlabel('Counts')
    pyplot.ylabel('Relative Frequency')
    pyplot.title('Probability distribution for counts in ')
    pyplot.xlim(0,20000)
    pyplot.show()

    print "----------------------------------------------------------------------------------------------------------"
    print "PCA Analysis "

    PCA(NormMatrix,'PCA2DDiffgene.png','PCA3DDiffgene.png')



    print "----------------------------------------------------------------------------------------------------------"
    print " "

    kmeans("PCA_Clustering_NONDiffGene", 3, "PCA_NONDiffgene_Clustering.png")
    kmeans("PCA_Clustering_NONDiffGene", 4, "PCA_NONDiffgene_Clustering_4.png")

if __name__ == '__main__':
        main()

