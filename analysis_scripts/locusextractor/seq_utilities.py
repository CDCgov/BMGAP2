import Bio 
from Bio import SeqIO, Seq, SeqRecord
import re
import utilities
import os
from collections import defaultdict
unambig = "GATC"
unambig_set = set(unambig)
IUPAC_AmbiguousDNA = "GATCRYWSMKHBVDN"
ambig_all_set = set(IUPAC_AmbiguousDNA)
ambig_only_set = ambig_all_set.difference(unambig_set)

base_sets = {
    "DNA4": set([i.lower() for i in unambig_set] + [i.upper() for i in unambig_set]),
    "DNA4_upper": unambig_set,
    "DNA_ambig":set([i.lower() for i in ambig_all_set] + [i.upper() for i in ambig_all_set]),
    "Gaps":set('-')
}
script_version = 1.3 #27 Jan 2016

##Seq is a SeqRecord, Seq, or string
### Result is "None" or a location on the string (alt: should "none" be len(seq)?)
def convertToString(seq):
    s = seq
    if isinstance(s,Bio.SeqRecord.SeqRecord):
        s=s.seq
    if isinstance(s,Bio.Seq.Seq):
        s=str(s)
    return s
    
def find_first_stop(seq):
    s = convertToString(seq)
    if not isinstance(s, str):
        print("###Error###")
        print(seq)
        raise "Must have a string!"

    first = None
    for i in range(0,len(s),3):
        x=s[i:i+3]
        if x in ('TAA','TAG','TGA'): 
            first = i
            break
    return first

##Seq is a SeqRecord, Seq, or string        
def trim_at_first_stop(seq):
    stop = find_first_stop(seq)
    if stop is None:
        return seq
    else:
        return seq[:stop]
    
def internal_stop(my_seq,stop_codon=None):
    if stop_codon is None:
        stop_codon = find_first_stop(my_seq)
    return (stop_codon is not None) and (stop_codon < len(my_seq) - 3)    
        
##Returns a list of dicts containing results from analyzing each gene
##Defaults to bacterial.
from Bio.Data import CodonTable
def ORF_analysis(multifasta, translation_table = 11): 
    data = []
    with open(multifasta,'rt') as fin:
        for seq in SeqIO.parse(fin,'fasta'):
            seq_info = ORF_info(seq,translation_table)
            data.append(seq_info)
    return data
            
def ORF_info(seq,translation_table=11):
    trans_code = CodonTable.unambiguous_dna_by_id[translation_table]
    seq_info = {}
    try:
        seq_info['Allele'] =seq.name
    except AttributeError:
        seq_info['Allele'] = 'Unnamed'
    seq_info['Length'] = len(seq)
    first_stop = find_first_stop(seq)
    seq_info['FirstStop'] = first_stop
    seq_info['InternalStop'] = internal_stop(seq,first_stop)
    first_cdn = convertToString(seq)[0:3]
    has_start = first_cdn in trans_code.start_codons
    seq_info['HasStart'] = has_start
    seq_info['FullORF'] = has_start and isinstance(first_stop,int) and  (first_stop == len(seq) - 3)    
    return seq_info    
    
def unambiguous_sequence(seq):
    s = convertToString(seq)
    safe_char = base_sets['DNA4']
    seq_set = set(s)
    return seq_set <= safe_char
    
def describeSequences(sequenceFile):
    result = defaultdict(int)
    result['FileSize'] = os.path.getsize(sequenceFile)
    for s in seqs_guess_and_parse(sequenceFile):
        result['Sequences'] += 1
        result['Nucleotides'] += len(s)

##Trims all bases from edges of "FASTQ_seq" until reaching first base with Qual score <= "qual"
def trimFASTQtoFirstBase(FASTQ_seq, qual):
    start = 0
    quals =  FASTQ_seq.letter_annotations["phred_quality"]
    while quals[start] < qual:
        start += 1
        if start == len(FASTQ_seq):
            return None 
    stop = len(FASTQ_seq)
    while quals[stop-1] < qual:
        stop -= 1 #No need to test for over-run, since we already found a stop point when looking from left
    return FASTQ_seq[start:stop]
    
def standardize_contig_names(contig_list,base_ID):
    new_contigs = []
    name_set = set([None])
    c = 0 #counter for un-named contigs
    for contig in contig_list:
        name_match = re.search(r'(\d+)$',contig.id) ##Search for digits
        ctgID = int(name_match.group(1)) if name_match else None
        while ctgID in name_set:
            c += 1
            ctgID = c
        name_set.add(ctgID)
        contig.id = '{}_ctg_{:03d}'.format(base_ID,ctgID)
        new_contigs.append(contig)
    return new_contigs, c    

def seq_guess_and_read(filename):
    seq = None
    seq_format, compressed = utilities.guessFileFormat(filename)
    if seq_format is not None:
        with utilities.flexible_handle(filename, compressed, 'rt') as seq_in:
            seq = SeqIO.read(seq_in,seq_format)
    else:
        print("Cannot infer sequence format for file: " + filename)            
    return seq

def seqs_guess_and_parse(filename):
    seq = None
    seq_format, compressed = utilities.guessFileFormat(filename)
    if seq_format is not None:
        with utilities.flexible_handle(filename, compressed, 'rt') as seq_in:
            return SeqIO.parse(seq_in,seq_format)
    else:
        print("Cannot infer sequence format for file: " + filename)
    return seq

def seqs_guess_and_parse2list(filename):
    seq = None
    seq_format, compressed = utilities.guessFileFormat(filename)
    if seq_format is not None:
        with utilities.flexible_handle(filename, compressed, 'rt') as seq_in:
            seq = [x for x in SeqIO.parse(seq_in,seq_format)]
    else:
        print("Cannot infer sequence format for file: " + filename)
    return seq

##will throw an IOError if bad filename
def seqs_guess_and_parse2dict(filename):
    if not isinstance(filename,str):
        raise TypeError("Filename must be string, is {}".format(type(filename)))
    seq_dict = None
    seq_format, compressed = utilities.guessFileFormat(filename)
    if seq_format is not None:
        with utilities.flexible_handle(filename, compressed, 'rt') as seq_in:
            seq_dict = SeqIO.to_dict(SeqIO.parse(seq_in,seq_format))
    else:
        print("Cannot infer sequence format for file: " + filename)        
    return seq_dict

##Extract regions from seq, as defined by start and stop values in the tuples of region_list.
## If stop < start, extracted regions will be reverse complement of Seq
## Assume that values are indexed at 1, but allow them to be indexed at 0
def extractRegions(seq,region_list,basename='Region',index=1): 
    ##Validation
    if isinstance(seq,SeqRecord): ##Reduce sequence to string
        s = str(seq.seq)
    elif isinstance(seq,Seq):
        s = str(seq)
    elif isinstance(seq,str):
        s = seq
    else:
        raise TypeError("seq variable should be string, Seq, or SeqRecord, but is {}".format(type(seq)))
    if not isinstance(region_list,list):
        raise TypeError("Region list should be a list, but is {}".format(type(region_list)))
    for r in region_list: 
        r[0] -= index
        r[1] -= index
        if not isinstance(r,list):
            raise TypeError("Items in regions list should be list. Is {}".format(type(r)))##This should tolerate tuples too...
        if not (len(r) == 2):
            raise ValueError("Items in region list should be of length 2. Is {}".format(len(r)))
        if (r[0] >= len(s)) or (r[0] < 0):
            raise ValueError("Sequence indexes must be smaller than sequence length: {} and {}".format(r[0],len(s)))
        if (r[1] >= len(s)) or (r[1] < 1):
            raise ValueError("Sequence indexes must be smaller than sequence length: {} and {}".format(r[1],len(s)))
    if not isinstance(basename,str):
        raise TypeError('Basename must be string. is {}'.format(type(basename)))
    ##Extract
    seqs = []
    for r in region_list:
        beg = r[0]
        end = r[1]
        seq_name = '{}_{}_{}'.format(basename,beg,end)
        if beg < end:
            first = beg
            last = end
            Forward = True
        else:
            first = beg
            last = end
            Forward = False
        region = Seq.Seq(s[first,last+1])
        if not Forward:
            region = region.reverse_complement()
        sr = SeqRecord.SeqRecord(region,id=seq_name, name=seq_name, description='')
        seqs.append(sr) 
    return seqs
     