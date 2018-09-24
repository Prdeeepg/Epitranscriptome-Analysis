import numpy as np
from prettytable import PrettyTable
x = PrettyTable()


def zero_filter(matrix, gene_order):
    c = [8, 15, 22, 29, 36, 43]
    gene_count = np.arange(np.shape(matrix)[0]).reshape(np.shape(matrix)[0], 1).astype('object')
    deleted_rows = 0
    for i in range(0, np.shape(matrix)[0]):
        try:
            gene_count[i] = float(matrix[i, 8]) + float(matrix[i, 15]) + float(matrix[i, 22]) + float(matrix[i, 29]) + float(matrix[i, 36]) + float(matrix[i, 43])
        except ValueError:
            print "error in line: ", i
        if gene_count[i] == 0.00:
            deleted_rows = deleted_rows +1

    low_count_rows = np.where(gene_count <= 0)[0]
    filtered_matrix = np.delete(matrix,low_count_rows,axis=0)
    filtered_genes = np.delete(gene_order,low_count_rows,axis=0)

    print "\nNo. of low count rows = ",np.shape(low_count_rows)
    print "\nDimensions of Filtered matrix = ",np.shape(filtered_matrix)
    print "\nDimensions of Filtered genes = ",np.shape(filtered_genes)

    return filtered_matrix, filtered_genes

