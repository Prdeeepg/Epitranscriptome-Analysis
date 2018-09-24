from scipy import stats
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as pyplot
from sklearn.cluster import KMeans


def kmeans(file_string,f,outfile):


    normed = np.loadtxt(file_string, usecols=(1, 2), dtype = 'float', skiprows=1)
    normed = normed.astype(float)
    gene_name = np.loadtxt(file_string, usecols=(0,), skiprows=1, dtype = 'string')
    print np.shape(normed)
    print normed[0:10,]
    np.random.seed(5)

    ###Transpose your data according to your question
    #X=np.transpose(X)
    #print np.shape(X)

   ###Here we need to normalize data to range from 0 to 1 or by zscore

    kmeans = KMeans(n_clusters=f, random_state=0).fit(normed)
    labels = kmeans.labels_
    centroid = kmeans.cluster_centers_
    colors = ['r', 'g', 'b', 'magenta', 'orange', 'grey', 'm', 'c']

   ###Above, depending on your question, determine how many clusters you should use
   ###Here we need to work with scikitlearn.cluster.KMeans to determine how to plot clusters on the scatterplot below

    kmeans.fit(normed)
    print np.shape(centroid)
    pyplot.figure()
    for i in range(0, f):
        values = np.where(kmeans.labels_ == i)[0]
        print np.shape(values)
        print "0 value = ",normed[values][:,0]
        print "1 value = ",normed[values][:,1]
        pyplot.scatter(normed[values][:, 0], normed[values][:, 1], color=colors[i])
    pyplot.title(('Kmeans clustering of the genes with ' + str(f) + ' clusters'))
    labels = ["Human Liver Carcinoma", "Heatshock", "UV", "Hepatocyte GF", "Interferons", "Brain"]
    for i in range(0,6):
        pyplot.annotate(labels[i],(normed[i,0], normed[i,1]))
    pyplot.xlabel("PC 1")
    pyplot.ylabel("PC 2")

    pyplot.savefig(outfile, bbox_inches='tight')

    #pyplot.figure()
    #pyplot.scatter(X[:,0],X[:,1])

    #pyplot.scatter(X[(np.where(kmeans.labels_=0)),0],X[(np.where(kmeans.labels_=0),1]),color='red')
    #pyplot.scatter(X[np.where(kmeans.labels_=1),0],X[np.where(kmeans.labels_=1),1],color='green')
    #pyplot.title('Kmeans clusters: ES vs. NPC')
    #pyplot.show()

    return









