import numpy as np
from prettytable import from_csv
from prettytable import PrettyTable
x = PrettyTable()
import csv



def PeaksAndReads(matrix):
    Summ = np.arange(70).reshape(7, 10).astype('object')
    Reads = np.arange(28).reshape(7, 4).astype('object')

    experiment = ["Experiments", "Human Liver Carcinoma", "Heatshock", "UV",
                  "Hepatocyte GF", "Interferons", "Brain"]
    Summary = ["Total Trans", "Trans with Peaks", "Total Peaks", "Coding trans",
               "Coding w/ peaks", "Total peaks in Coding", "Nnocoding trans",
               "Nnocoding w/ peaks", "Total peaks in noncoding"]
    read = ["Total Reads", "Noncoding Reads", "Coding Reads"]
    for i in range(0, 7):
        Summ[i,0] = experiment[i]
        Reads[i,0] = experiment[i]

    for i in range(0,9):
        Summ[0,(i+1)] = Summary[i]

    for i in range(0,3):
        Reads[0,(i+1)] = read[i]

    b = range(1,7)
    c = [8, 15, 22, 29, 36, 43]
    y = [10, 17, 24, 31, 38, 45]
    for h in range(0,6):
        noncoding =0
        coding = 0
        peakedgenenonc = 0
        noofpeaksnonc = 0.00
        peakedgenecod = 0
        noofpeakscod = 0.00
        codingreads = 0
        noncodingreads = 0
        for i in range(0, np.shape(matrix)[0]):
            if matrix[i,6] == 'noncoding':
                noncoding = noncoding + 1
                if float(matrix[i, c[h]]) > 0.00:
                    peakedgenenonc = peakedgenenonc + 1
                    noofpeaksnonc = noofpeaksnonc + float(matrix[i, c[h]])
                    noncodingreads = noncodingreads + float(matrix[i, y[h]].split( )[0])
            if matrix[i,6] == 'coding':
                coding = coding + 1
                if float(matrix[i, c[h]]) > 0.00:
                    peakedgenecod = peakedgenecod + 1
                    noofpeakscod = noofpeakscod + float(matrix[i, c[h]])
                    codingreads = codingreads + float((matrix[i, y[h]].split( )[0]))
        i = h
        Summ[b[i],1] = np.shape(matrix)[0]             # Total number of Transcripts       SAME
        Summ[b[i],2] = peakedgenenonc+peakedgenecod    # transcripts with peaks            DIFFERENT
        Summ[b[i],3] = noofpeaksnonc+noofpeakscod      # total no. of peaks in all trans.  DIFFERENT
        Summ[b[i],4] = coding                          # Coding transcripts                SAME
        Summ[b[i],5] = peakedgenecod                   # Coding trans. with peaks          DIFFERENT
        Summ[b[i],6] = noofpeakscod                    # Total Peaks in coding trans.      DIFFERENT
        Summ[b[i],7] = noncoding                       # Noncoding transcripts             SAME
        Summ[b[i],8] = peakedgenenonc                  # noncoding transcripts with peaks  DIFFERENT
        Summ[b[i],9] = noofpeaksnonc                   # Total peaks in noncoding trans.   DIFFERENT
        Reads[b[i],1] = codingreads+noncodingreads
        Reads[b[i],2] = codingreads
        Reads[b[i],3] = noncodingreads

    with open('SummaryStat2.csv', 'wb') as mello:
        writer = csv.writer(mello, delimiter=',')
        data2 = Summ
        writer.writerows(data2)

    fp = open("SummaryStat2.csv", "r")
    file = from_csv(fp)
    file.align = "l"
    fp.close()
    print(file)
    print "\n"

    with open('SummaryStat_Reads.csv', 'wb') as vello:
        writer = csv.writer(vello, delimiter=',')
        data3 = Reads
        writer.writerows(data3)

    fp = open("SummaryStat_Reads.csv", "r")
    file = from_csv(fp)
    file.align = "l"
    fp.close()
    print(file)