import numpy as np

def load_m6Aseq_differential_counts(filename,delimiter='\t'):
    """
    Loads a RNA-seq counts file from disk.

    Parameters
    ----------
    filename : str
        The name of the RNA-seq counts file to parse
    delimiter : str
        The character by which the parser will separate the columns. 
        It is recommended to pass '\t' if your file is tab-delimited,
        or ' ' if your file is space-delimited.

    Returns
    -------
    (np.ndarray,List[str],List[str])
        Returns a tuple containing a 2 x 2 counts matrix, where
        each row contains gene counts for a given gene across all libraries and
        the columns are different libraries, and two lists containing gene and
        library names, corresponding to the row and column names of the counts
        matrix. 
    """
    lib_order = []
    with open(filename,'r') as handle:
        # the line.split(delimiter) is looking at each line of the txt file and then from these lines
        # it collects the columns in ().
        columns = zip(*[line.split(delimiter) for line in handle])

        # It saves the columns in[]
        columns = [list(column) for column in columns]

        # Now looking in each set of column in the collective columns in columns
        for column in columns:

            # The name gives the first column value which is the field name for each column.
            name = column[0]
            # this removes the first value from each column
            column.remove(name)
            for i in range(len(column)):
                column[i] = column[i].strip('\n').strip('\r')

            if 'id' in name.strip('#').strip('\n').lower():
                gene_order = column
                continue
            else:
                lib_order.append(name.strip('#').strip('\n').strip('\r'))

    return np.array([columns[i] for i in range(1,len(columns))]).T,gene_order,lib_order
