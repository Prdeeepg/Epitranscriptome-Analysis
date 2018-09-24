import numpy as np
from compute_pca import compute_pca
from plot_pca import plot_pca
from plot_pca3D import plot_pca3D

def PCA( NormMatrix, file2D, file3D):

    transposed_matrix = NormMatrix.T

    proj, pve, pcs = compute_pca(transposed_matrix, scaled=True, logged=False, kernel=None,
                                 kernel_kwargs=None, variant='pca', pf=1)
    print '\nTransposed Matrix =', np.shape(transposed_matrix)
    print '\nProj =', np.shape(proj)
    print '\nPVE =', np.shape(pve)
    print '\nPCS =', np.shape(pcs)
    print proj
    experiment = [ "Human Liver Carcinoma", "Heatshock", "UV",
                  "Hepatocyte GF", "Interferons", "Brain"]
    plot_pca(file2D, proj, pcs=(0, 1), labels=experiment, label_points=True,
             levels=None, colors=None, legend=True, s=100)
    plot_pca3D(file3D, proj, pcs=(0, 1, 2), labels=experiment, label_points=True,
             levels=None, colors=None, legend=True, s=100)

