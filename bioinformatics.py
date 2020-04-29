def fasta_genome_to_pickle(fasta_filename, out_filename):
    '''
    Takes a FASTA file
    With one record per chromosome
    Returns a dict with each key as a chromosome
    Saves dict to a pickle file for future access
    '''

    genome_iter = SeqIO.parse(fasta_filename, 'fasta')
    genome_dict = SeqIO.to_dict(genome_iter)
    with open(out_filename, 'wb') as pickle_dump_file:
        pickle.dump(genome_dict, pickle_dump_file)

def fetch_seq(directory, chrom, start, end, UCSC=False, lineBreaks=True, header = True):
    '''
    Returns the DNA sequence from a region of the specified genome
    '''

    fn = directory + chrom + '.fa'
    fh = open(fn,'r')

    headerOffset = 0
    nStart = 0
    nEnd = 0

    if header:
        fh.seek(0)
        headerOffset = len(fh.readline())
    if lineBreaks:
        nStart = (start-1)/50
        nEnd = (end-1)/50
    if UCSC:
        fh.seek((start+nStart+headerOffset))
    else:
        fh.seek((start-1+nStart+headerOffset))
    span = ((end+nEnd-1)-(start+nStart-1))

    read = fh.read(span)
    if lineBreaks:
        read = read.replace('\n','')

    return read
    fh.close()
