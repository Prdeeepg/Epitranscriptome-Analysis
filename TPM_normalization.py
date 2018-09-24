import numpy as np

def TPM_norm(filtered_matrix, Norm_Type = "transcript"):
    lengthnorm = np.arange((np.shape(filtered_matrix)[0])*6).reshape(np.shape(filtered_matrix)[0], 6).astype('object')
    rangeofgene = np.arange((np.shape(filtered_matrix)[0])*6).reshape(np.shape(filtered_matrix)[0], 6).astype('object')
    rangeofexon = np.arange((np.shape(filtered_matrix)[0])*6).reshape(np.shape(filtered_matrix)[0], 6).astype('object')

    c = [11, 18, 25, 32, 39, 46]    # Has all the exon length for each transcripts
    y = [10, 17, 24, 31, 38, 45]    # Has all the reads for each transcripts

    if Norm_Type == "transcript":
        for h in range(0,6):
            for i in range(0, np.shape(lengthnorm)[0]):

                # normalizing for transcript length
                try:
                    lengthnorm[i,h] = ((float(filtered_matrix[i,y[h]].split(" ")[0]))/((float(filtered_matrix[i,2])-float(filtered_matrix[i,1]))/1000))
                except AttributeError:
                    print "error"
                except ValueError:
                    lengthnorm[i, h] = float(0)

                # Collecting stats for range of transcript length
                rangeofgene[i,h] = (float(filtered_matrix[i,2]))-(float(filtered_matrix[i,1]))

                #Collecting stats for range of exon length
                try:
                    rangeofexon[i,h] = float(filtered_matrix[i, c[h]])
                except ValueError:
                    try:
                        rangeofexon[i,h] = float(filtered_matrix[i, c[h]].split(" ")[0])
                    except ValueError:
                        rangeofexon[i,h] = float(0)

    if Norm_Type == "exon":
        for h in range(0,6):
            for i in range(0, np.shape(lengthnorm)[0]):
                try:
                    rangeofexon[i,h] =(((float(filtered_matrix[i,y[h]].split(" ")[0])))/(float(filtered_matrix[i, c[h]])/1000))
                except ZeroDivisionError:
                    print "zero division error in line =", i, c[h],filtered_matrix[i, c[h]]
                except ValueError:
                    try:
                        rangeofexon[i,h] = (((float(filtered_matrix[i,y[h]].split(" ")[0])))/((float(filtered_matrix[i, c[h]].split(" ")[0]))/1000))
                    except ValueError:
                        rangeofexon[i,h] = float(0)




    rpm = np.nansum(lengthnorm, axis=0) / 1000000
    tpm_matrix = lengthnorm / rpm

    return tpm_matrix, rangeofexon, rangeofgene