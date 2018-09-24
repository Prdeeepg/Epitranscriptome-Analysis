import numpy as np
from prettytable import from_csv
from prettytable import PrettyTable
x = PrettyTable()
import csv

def normtable(matrix, gene_order):

    lib_order = ["Experiments", "Human Liver Carcinoma", "Heatshock", "UV", "Hepatocyte GF", "Interferons", "Brain"]
    x = 1 + np.shape(matrix)[0]
    y = 1 + np.shape(matrix)[1]
    print x, y
    a = np.arange(x*y).reshape(x,y).astype('object')

    for k in range(0, y):
        a[0,k] = lib_order[k]

    for h in range(1, x):
        a[h,0] = gene_order[h-1]

    for h in range(1, x):
        for t in range(1,7):
            a[h,t] = matrix[(h-1),(t-1)]

    with open('TPMnorm_matrix.csv', 'wb') as hello:
        writer = csv.writer(hello, delimiter=',')
        data = a
        writer.writerows(data)

    fp = open("TPMnorm_matrix.csv", "r")
    file = from_csv(fp)
    file.align="l"
    fp.close()


