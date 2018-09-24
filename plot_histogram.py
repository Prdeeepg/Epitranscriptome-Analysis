import numpy as np
import matplotlib.pyplot as pyplot

"""
This function primarily works to create a histogram(s) from an input data.
"""

def histogram(matrix,gene_order,lib_order):


    #Print information about the dimensions of your counts array and the number of independent samples
    print np.shape(matrix), "size of counts array"
    print np.shape(gene_order), "size of gene order array"
    print lib_order, "Here are the replicates in this file"

    #Here, for every file, print descriptive statistics and plot histograms and boxplots
    for j in range(0,np.shape(matrix)[1]):
        print j, np.shape(matrix[:,j])
        print np.sum(matrix[:,j]), "here is the sum of counts in sample ", lib_order[j]

        pyplot.figure()
        pyplot.hist(matrix[:,j],20,normed=True)
        #pyplot.xlim(0,5000)
        pyplot.xlabel('Counts')
        pyplot.ylabel('Relative Frequency')
        pyplot.title('Probability distribution for counts in '+lib_order[j])
        pyplot.show()
