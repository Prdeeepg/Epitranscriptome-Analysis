import numpy as np
from prettytable import from_csv
from prettytable import PrettyTable
x = PrettyTable()
import csv

def Table(matrix, lib_order, gene_order):

    a = np.arange(250).reshape(50,5).astype('object')

    i=0
    a[0,0] = 'id'

    for k in range(1, 50):
        a[k,0] = lib_order[(k-1)]

    for j in range(0, 4):
        i=j+1
        a[0,i] = gene_order[j]
        for u in range(1, 50):
            a[u,i] = matrix[j,(u-1)]

    with open('SummaryStat.csv', 'wb') as hello:
        writer = csv.writer(hello, delimiter=',')
        data = a
        writer.writerows(data)

    fp = open("SummaryStat.csv", "r")
    file = from_csv(fp)
    file.align="l"
    fp.close()
    print(file)