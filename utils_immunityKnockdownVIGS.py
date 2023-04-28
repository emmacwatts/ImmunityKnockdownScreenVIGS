#Justif. - VIGS tool only available for Niben annotation sequences. 
#Papers use VIGs sequences targeting only one gene.

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def revComp(inputSeq):
    """
    This function takes an input sequence and returns the reverse complement.

    params: inputSeq in str format
    returns: revComp in str format

    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    
    revComp = ""
    for base in inputSeq[::-1]:
        revComp += complement[(base.upper())]

    return revComp

def localBLAST(query, outputFilename = 'unspecified', database = 'inputs/BLAST_LAB350cDNA/NbLAB350_Proteome.fasta'):
    """
    Run BLAST against LAB350 N benthamiana transcriptome for an input gene. 
    Return results for selection of sequences that correspond to input (potentially multiple).
    """
    import subprocess

    #If an output file name has not been provided, name it based on the query file
    if outputFilename == 'unspecified':
        queryFull = query.split("/") #Just use the filename (rather than the whole filepath)
        outputFilename = f"inputs/BLAST_LAB350cDNA/BLASToutputs/{queryFull[-1]}_blast"

    cmd = f"blastp -db {database} -query {query} -out {outputFilename}"

    subprocess.run(cmd, shell=True)

def findPrimerRegion(FWDprimer, REVprimer, geneSequenceFASTA):
    """
    Find region of a gene sequence that primers are designed to anneal to. 
    Will iterate and cut one off, re-searching to find a best-match region.
    Returns original full gene with highlight for best match region (should still manual check if this is a meaningful match).
    """
    from Bio import SeqIO

    for seq in SeqIO.parse(open(geneSequenceFASTA), 'fasta'):
        geneSequence = str(seq.seq)

    #Find start position of the forward primer
    FWDstart = geneSequence.find(FWDprimer)
    FWDprimerc = FWDprimer #Save original before manipulating
    while FWDstart == -1: #If substring not found, remove one letter and search again
        FWDprimerc = FWDprimerc[:-1]
        print(f"using truncated FWD primer to search {FWDprimerc}")
        FWDstart = geneSequence.find(FWDprimerc)

    #Find end position of the revComp reverse primer
    REVprimerrevcomp = revComp(REVprimer)
    REVstart = geneSequence.find(REVprimerrevcomp)
    while REVstart == -1:
        REVprimerrevcomp = REVprimerrevcomp[:-1]
        print(f"using truncated REV primer to search {REVprimerrevcomp}")
        REVstart = geneSequence.find(REVprimerrevcomp)

    REVend = REVstart + len(REVprimerrevcomp)

    highlightedSeq = geneSequence[:FWDstart] + f"{bcolors.BOLD}{geneSequence[FWDstart:REVend+1]}{bcolors.BOLD}" + geneSequence[REVend+1:]

    return highlightedSeq

def VIGSsequenceAlignment():
    """
    Align sequences to be targeted by VIGS. Define a common region to be used for VIGS, 
    and narrow this down further with VIGS target checks - fragment size, n-mers, mismatches.
    """

def TRV2ggVIGSprimers():
    """
    Designs primers for VIGs for a given target region. 
    """