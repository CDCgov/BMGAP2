#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 14:48:32 2018

@author: ymw8
"""

from Bio.SubsMat import MatrixInfo
from collections import defaultdict
import os
import numpy as np
from math import log2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import pairwise2
import pandas as pd
import re

valid_mutations = ['Substitution','Insertion','Deletion']

Nm = 'Neisseria meningitidis'
Hi = 'Haemophilus influenzae'
AST_database_dict = {
    Hi:'Hi_db/',
    Nm:'Nm_db/'
}

##Loads reference sequences, translates, and validates; Returns copy of dict with sequences loaded
def load_locus_info(locus_info,variant_info,db_location):
    loaded_seqs = {}
    has_file = locus_info.ref_file.notnull()
    if sum(has_file) < len(locus_info):
        print("Only have sequence file for {}/{} genes".format(sum(has_file),len(locus_info)))
    locus_info = locus_info[has_file].copy()
    full_path = db_location + locus_info.ref_file
    locus_info['ref_file'] = full_path.apply(os.path.abspath)
    
    for locus in locus_info.index:##This is our name, not PMGA name
        locus_variants = variant_info[variant_info.Gene == locus].copy()#.set_index()
        deletions = locus_variants['Type of Mutation'] == 'Deletion'
        locus_variants.loc[deletions,'Resistant Variants'] = '-'
        if len(locus_variants) > 0: ##Nothing to validate if no variants are identified
            ref_file = locus_info.loc[locus,'ref_file']
            try:
                ref_raw = SeqIO.read(ref_file,'fasta')
            except:
                print(ref_file)
                raise
            else:
                print("{} is {}bp".format(locus,len(ref_raw)))
                if locus_info.loc[locus,'coding']:
                    if (len(ref_raw) % 3) != 0:
                        print("Warning: {} length is not multiple of three: {}bp".format(ref_raw.name,len(ref_raw)))
                    ref_seq = ref_raw.seq.translate() 
                    stop = ref_seq.find('*')
                    if stop > 0:
                        ref_seq = ref_seq[:stop]
                else:
                    ref_seq = ref_raw

                for i,r in locus_variants.iterrows():
                    start = int(r['First Position of Feature in Reference Protein Sequence']) - 1
                    assert (start > 0) and (start < len(ref_seq)),'{} Variant start position {} is not in sequence length {}'.format(locus,start+1,len(ref_seq))
                    stop  = int(r['Last Position of Feature in Reference Protein Sequence'])
                    assert (start > 0) and (start < len(ref_seq)),'{} Variant stop position {} is not in sequence length {}'.format(locus,stop+1,len(ref_seq))
                    if r['Feature Category'] == 'Mutation':
                        susc  = r['Susceptible Variants'].strip()
                        resist = r['Resistant Variants'].strip()
                        assert (stop-start) == len(susc),"{}-{} Variant must match length of locus. {} vs {}".format(locus,start+1,len(susc),(stop-start))
                        l_flank = start - 10
                        r_flank = stop + 10
                        if ref_seq[start:stop] != susc:
                            print("{}-{} Mutation reference must match sequence; {}_{}_{} vs {}".format(locus,start+1,ref_seq[l_flank:start],ref_seq[start:stop],ref_seq[stop:r_flank],susc))
                        
                    elif r['Feature Category'] != 'Region':
                        print("{}-{} The only valid feature categories are Mutation and Region".format(locus,start+1))
                    
                loaded_seqs[locus_info.loc[locus,'PMGA']]= ref_seq  
    return loaded_seqs
    
    
#loci_data is extract_loci_from_JSON_collection -- this must be performed with a locus list that is inclusive of locus_info keys
#locus_info is dataframe with variants
#locus_seqs is SeqRecords
# return table with each genome showing the known mutations plus a column for additional mutations at those loci
def identify_mutations_from_JSON_collection(loci_data,variant_info,locus_references, locus_rename=None):
    mutations_genomes = {}
    for filename,loci_dict in loci_data.items():
        mutation_dict = identify_mutations_in_single_JSON(loci_dict,variant_info,locus_references,locus_rename,filename.split("_")[0])
        mutations_genomes[filename] = mutation_dict
    return pd.DataFrame(mutations_genomes).transpose()
               
def get_seq_from_locus_list(locus_list):
    seq_str = ''
    if len(locus_list) > 0:
        seq_str = locus_list[0]['qseq']
        if len(locus_list) > 1: ##check if sequences are identical at all loci
            seqs_set = set([locus_item['qseq'] for locus_item in locus_list])
            try:
                locus = locus_list[0]['allele']
            except:
                locus = 'Unknown locus'  
            try:
                lid = locus_list[0]['contig'].split("_")[0]
            except:
                lid = "unknown genome"            
            if len(seqs_set) > 1:
                print("Warning: found {} sequences at {} loci for {} in {}".format(len(seqs_set),len(locus_list),locus,lid))            
                seq_str = ''
            else:
                print("Notice: found {} loci with {} sequence for {} in {}".format(len(locus_list),len(seqs_set),locus,lid))
    return seq_str

def validate_seq(locus_name,ref_alnseq,p,mutation_ref,true_p):
    if ref_alnseq[p] not in mutation_ref[true_p]['ref']+'-':
        print(ref_alnseq[p-10:p+10])
        print("Warning: {} Reference allele {} must match position in reference sequencese ({} at {}).".format(locus_name,mutation_ref[true_p]['ref'],
                                                                                                               ref_alnseq[p],true_p))
    assert len(mutation_ref[true_p]['ref']) == 1,"Multiple reference alleles not permitted (yet)."
    ##A gap can be tolerated, since that will prevent true_p from incrementing
    assert ref_alnseq[p] in (mutation_ref[true_p]['ref']+'-'),"{} Reference allele {} must match position in reference sequences ({} at {}).".format(locus_name,mutation_ref[true_p]['ref'],ref_alnseq[p],true_p)
    
##Used in identify mutations. Returns whether mutation is known. 
def document_mutation(mutation_string,locus_name,position,domain_mutations,mutation_dict,locus_regions,locus_variants,mutation_type,variant):
    known_mutation = False
    if not isinstance(mutation_type,str):
        raise ValueError("Mutation type not defined")
    if mutation_type not in valid_mutations:
        raise ValueError("Mutation type must be Ins Del or Subs")
    if mutation_string == '':
        print()
        print("ERROR! No mutation string defined {} {}".format(locus_name,position))
        print()                               
    ##Record mutations by regions
    domain_mutations['all'].append(mutation_string) ## all
    if len(locus_regions) > 0: ## each region
        for _, r in locus_regions.iterrows():
            domain = r['Name']
            if (position >= r['First Position of Feature in Reference Protein Sequence']) and (position <= r['Last Position of Feature in Reference Protein Sequence']):
                domain_mutations[domain].append(mutation_string)  
    ##Examine any variants defined for this position
    ref_p = position + 1 if mutation_type == 'Insertion' else position
    position_variants = locus_variants['First Position of Feature in Reference Protein Sequence'] == ref_p
    mutation_variants = locus_variants['Type of Mutation'] == mutation_type    
    this_p = locus_variants.loc[position_variants & mutation_variants]
    if len(this_p) > 0: #known mutation location
        ##Real calculations
        if variant not in this_p['Susceptible Variants'].values: ##Don't document known susceptibles in main columns
            if variant in this_p['Resistant Variants'].values:#known mutation
                mutation_name = '{}_{}'.format(locus_name,ref_p) ##Note: this key must match the summary at the end of the evaluation
                ##append a given name if possible
                exact = this_p.set_index('Resistant Variants')
                n = exact.loc[variant,'Name']
                if pd.notnull(n):
                    mutation_name += ' ({})'.format(n) ##Note: this must match the end of the evlautaion too...
                mutation_dict[mutation_name] = mutation_string
                known_mutation = True
            else:#novel mutation at existing location. 
                other_variant = "{}_{}".format(locus_name,mutation_string)
                mutation_dict['other'].append(other_variant)    
    return known_mutation

def identify_mutations_in_single_JSON(loci_dict,variant_info,locus_references,locus_rename,genome_id):
    mutation_dict = {}
    mutation_dict['other'] = []
    for locus,locus_list in loci_dict.items():##THis is PGMA nomenclature
        locus_name = locus_rename[locus] if isinstance(locus_rename,dict) else locus
        ##Only look at those with defined mutations
        locus_features = variant_info[variant_info.PMGA == locus]#.set_index()
        if len(locus_features) > 0: ##if this is true, there should be a sequence in locus_references
            locus_regions = locus_features[locus_features['Feature Category'] == 'Region']
            locus_variants = locus_features[locus_features['Feature Category'] == 'Mutation']
            seq_str = get_seq_from_locus_list(locus_list)
            if len(seq_str) == 0:
                print("Warning: could not identify single sequence for {} in {}".format(locus_name,genome_id))
                mutation_dict["{}_testLen".format(locus_name)] = len(seq_str)
            else:
                ##Start by performing alignment
                seq = Seq(seq_str)
                ref_seq = locus_references[locus]
                IUPAC_extendedProtein =  'ACDEFGHIKLMNPQRSTVWYBXZJUO*'
                IUPAC_extendedDNA = 'GATCBDSW'
                #checks if seq is protein or DNA sequence
                is_protein=all(char in IUPAC_extendedProtein for char in str(seq).upper())
                is_dna=all(char in IUPAC_extendedDNA for char in str(seq).upper())
                if is_protein:
                    stop = ref_seq.find('*') ##This should be redundant with load sequence
                    if stop > 0:
                        ref_seq = ref_seq[:stop]
                    seq = seq.translate()
                    stop = seq.find("*")
                    if stop > 0:
                        mutation_dict["{}_firstStop".format(locus_name)] = stop/len(seq)
                        seq = seq[:stop]   
                    else:
                        mutation_dict["{}_firstStop".format(locus_name)] = 'None'
                        # If the sequence starts with *, find any sequence between two stop codons.
                        match = re.search(r"\*([A-Z]+)\*?", str(seq))
                        if match:
                            seq = Seq(match.group(1))
                        # Otherwise if this is a string of just stop codons (really unexpected!),
                        # assign a dummy sequence to avoid biopython errors later
                        else:
                            seq = Seq("WWWWWWW")

                    alignments = pairwise2.align.globalds(str(seq),str(ref_seq),MatrixInfo.blosum95,-10,-0.1) ##default gap open and extend for culstal. 
                    
                elif is_dna:
                    print("DNA alignment has not been optimized")
                    alignments = pairwise2.align.globalxs(str(seq),str(ref_seq),-10,-0.1)  
                else:
                    raise ValueError("Must be protein or nucleotide!") # Not {}".format(ref_seq.alphabet))
                ##result format: list of (test_seq1,ref_seq2,score,begin,end)
                mutation_dict["{}_alnCount".format(locus_name)] = len(alignments)
                mutation_dict["{}_refLen".format(locus_name)] = len(ref_seq)
                mutation_dict["{}_testLen".format(locus_name)] = len(seq)
                max_aln = 0
                min_aln = len(seq) + len(ref_seq)
                for scan_aln in alignments:
                    aln_len = len(scan_aln[0])
                    min_aln = min(min_aln,aln_len)
                    max_aln = max(max_aln,aln_len)
                    if aln_len == min_aln:
                        aln = scan_aln
                mutation_dict["{}_maxLen".format(locus_name)] = max_aln
                mutation_dict["{}_minLen".format(locus_name)] = min_aln  
                ##Count mutations for this gene                
                domain_mutations = defaultdict(list)
                test_alnseq = aln[0]
                ref_alnseq = aln[1]
                true_p = deletion_start = 0 ##True_p is position on peptide sequence (not aligned with gaps)
                mutation_count = mismatches = known_mutations =0 
                insertion = deletion = ''
                for p in range(len(ref_alnseq)):
                    mutation_string = '' #make sure this is defined inside conditional statements
                    if ref_alnseq[p] == '-': ##insertion 
                        mismatches += 1
                        if insertion == '':
                            mutation_count += 1
                        insertion += test_alnseq[p] ##true_p reflects the position before the insertion
                        if len(deletion) > 0: ##end of deletion -- recorded below also
                            mutation_string = "{}del{}".format(deletion_start,deletion)
                            if document_mutation(mutation_string, locus_name, deletion_start,domain_mutations,mutation_dict,locus_regions,locus_variants,"Deletion",deletion):
                                known_mutations += 1
                            deletion = ''                            

                    else: ## match, deletion or substitution, move to next position in reference
                        if len(insertion) > 0: ##end of insertion
                            prior_base = ref_seq[true_p-1] if true_p > 0 else 'N-term'
                            subsequent_base = ref_seq[true_p] if true_p < len(ref_seq) else 'C-term'
                            if true_p == 0:
                                mutation_string = "{}_ins{}_{}{}".format(prior_base,insertion,true_p+1,subsequent_base)
                            elif true_p >= len(ref_seq):
                                mutation_string = "{}{}_ins{}_{}".format(true_p,prior_base,insertion,subsequent_base)
                            else:
                                mutation_string = "{}{}_ins{}_{}{}".format(true_p,prior_base,insertion,true_p+1,subsequent_base)
                            if document_mutation(mutation_string,locus_name,true_p,domain_mutations,mutation_dict,locus_regions,locus_variants,"Insertion",insertion):
                                known_mutations += 1                            
                            insertion = ''                        
                        true_p += 1
                        if test_alnseq[p] != ref_alnseq[p]: ##deletion or substitution                           
                            mismatches += 1
                            if test_alnseq[p] == '-': #deletion
                                if deletion == '': ##new deletion
                                    deletion_start = true_p
                                    mutation_count += 1 
                                deletion += ref_alnseq[p]
                            else: ## substitution 
                                mutation_count += 1 
                                if len(deletion) > 0: ##end of deletion
                                    mutation_string = "{}del{}".format(deletion_start,deletion)           
                                    if document_mutation(mutation_string,locus_name,deletion_start,domain_mutations,mutation_dict,locus_regions,locus_variants,"Deletion",deletion):
                                        known_mutations += 1
                                    deletion = ''   
                                ##                    
                                mutation_string = "{}{}{}".format(ref_alnseq[p],true_p,test_alnseq[p])
                                if document_mutation(mutation_string,locus_name,true_p,domain_mutations,mutation_dict,locus_regions,locus_variants,"Substitution",test_alnseq[p]):
                                    known_mutations += 1
                ###Wrap up any insertions or deletions
                if len(insertion) > 0: ##end of insertion
                    prior_base = ref_seq[true_p-1] if true_p > 0 else 'N-term'
                    subsequent_base = ref_seq[true_p] if true_p < len(ref_seq) else 'C-term'
                    mutation_string = "{}{}_{}{}ins{}".format(true_p,prior_base,true_p+1,subsequent_base,insertion)
                    if document_mutation(mutation_string,locus_name,true_p,domain_mutations,mutation_dict,locus_regions,locus_variants,"Insertion",insertion):
                        known_mutations += 1
                    insertion = ''
                if len(deletion) > 0: ##end of deletion -- recorded below also
                    mutation_string = "{}del{}".format(deletion_start,deletion)
                    document_mutation(mutation_string, locus_name, deletion_start,domain_mutations,mutation_dict,locus_regions,locus_variants,"Deletion",deletion)
                    deletion = ''              
                ##Fill in wild types
                for i in locus_variants.index:
                    ref_p = locus_variants.loc[i,'First Position of Feature in Reference Protein Sequence']
                    n = locus_variants.loc[i,'Name']
                    mutation_name = '{}_{}'.format(locus_name,ref_p) ##This should probably be a separate function.
                    if pd.notnull(n):
                        mutation_name += ' ({})'.format(n)                    
                    if mutation_name not in mutation_dict: ##This name must match document_mutations
                        mutation_dict[mutation_name] = '.'

                ##Summarize mutations for this gene
                if len(domain_mutations['all']) != mutation_count:
                    print("Error: {} {} did not keep all mutations {} vs {}".format(genome_id,locus_name,mutation_count,len(domain_mutations['all'])))
                    print(domain_mutations['all'])
                    print(aln)
                for domain, d_mutations in domain_mutations.items():
                    mutation_dict['{}_{}'.format(locus_name,domain)] = ', '.join(d_mutations) if len(d_mutations) > 0 else '.'
                mutation_dict["{}_mutation_count".format(locus_name)] = mutation_count
                mutation_dict["{}_mismatches".format(locus_name)] = mismatches
                mutation_dict["{}_known_mutations".format(locus_name)] = known_mutations

    mutation_dict['other'] = ', '.join(mutation_dict['other']) if len(mutation_dict['other']) > 0 else '.'
    return mutation_dict

#####Getting Json file from directory######
def get_json_file(directory):
    for file in os.listdir(directory):
        if file.endswith('final_results.json'):
            return file
##converting dictionary to fasta #######
def dict_to_fasta(seq_dict,outfile):
   # print(seq_dict)
    for id,seq in seq_dict.items():
        for key in seq:
            outfile.write('>'+id+'|'+seq[key])
            outfile.write("\n")
            outfile.write(key)
            outfile.write("\n")

###Making clean gene datatbase file########
def make_geneDB_file(gene_key):
    new_gene_key = gene_key.reset_index()
    new_gene_key2 = new_gene_key[['Gene','PMGA','Antimicrobic class']]
    new_gene_key_re = new_gene_key2[['PMGA','Gene','Antimicrobic class']]
    new_gene_key_sort = new_gene_key_re.iloc[new_gene_key_re.PMGA.str.lower().argsort()]
    return new_gene_key_sort
    
#summarize what genes were found######
def summarize_genes(locus_collection,gene_key):
    loci_set = set([k for loci in locus_collection.values() for k in loci.keys()])
    print('found')
    print(loci_set)
    print('did not find')
    print(set(gene_key.PMGA.tolist()).difference(loci_set))

#####fixing "internal stop codon call" ##########
def fix_internal_stop(locus_presence_table,locus_rename,gene_key):
    locus_presence_table.rename(columns=locus_rename,inplace=True)
    gene_key['coding']=gene_key['coding'].astype(bool) #change values to boolean to use ~ in next step
    not_coding = gene_key[~gene_key.coding].index.tolist()
    for c in locus_presence_table:
        if c in not_coding:
            isc = (locus_presence_table[c] == 'internal stop codon')
            locus_presence_table.loc[isc,c] = 'Present'

####Make Dictiorary of all the AMR genes found with their Antimicrobial class#####
def get_amr_genes(locus_table,allele_table,gene_key):
    amr_genes = {}
    locus_dict = locus_table.to_dict('list')
    allele_dict = allele_table.to_dict('list')
    #print(allele_dict)
    for index,row in gene_key.iterrows():
        for key,temp_value in locus_dict.items():
            if "(" in key:
                key_parse = re.search('\((.*)\)',key).group(1)
                key_parse2 = re.search('(.*) ',key).group(1)
                key_parse2.rstrip()
            else:
                key_parse = key
                key_parse2 = key
            value = temp_value[0]
            if index == key_parse:
                if key_parse2 not in amr_genes:
                    amr_genes[key_parse2] = {}
                amr_genes[key_parse2]["Gene_name"] = key_parse
                amr_genes[key_parse2]["status"] = value
                amr_genes[key_parse2]["Antimicrobial"] = row['Antimicrobic class']
                if len(allele_dict) != 0:
                    if key_parse2 in allele_table.keys():
                        amr_genes[key_parse2]["allele"] = allele_dict.get(key_parse2)[0]    
                        if amr_genes[key_parse2]["allele"] == 'temp_1':
                            amr_genes[key_parse2]["allele"] = 'New' #Since this is running on single genomes, there will only be a single allele for each genome. We want end-user to see 'New' rather than 'temp_1'...
                    else:amr_genes[key_parse2]["allele"] ="NA"
            else:
                if row['PMGA'] not in amr_genes:
                    amr_genes[row['PMGA']] = {}
                    amr_genes[row['PMGA']]["Gene_name"] = index
                    amr_genes[row['PMGA']]["status"] = "Absent"
                    amr_genes[row['PMGA']]["allele"] = "NA"
                    amr_genes[row['PMGA']]["Antimicrobial"] = row['Antimicrobic class']

    return amr_genes

######Making dictionary of all the mutations in the amr genes ################
def get_mutation(mutation_table,varinats):

    amr_mutations = {}
    mutation_dict = mutation_table.to_dict('list')
    for index,row in varinats.iterrows():
        pd.set_option('display.max_columns', None)
        var_gene = row['Gene']
        mutation = row['Susceptible Variants']+str(row['First Position of Feature in Reference Protein Sequence'])+row['Resistant Variants']
        amino_change = row['Susceptible Variants']+"->"+row['Resistant Variants']
        antibiotic = row['Associated Antibiotic Resistance(s)']
        mutation_type = row['Type of Mutation']
        for key,temp_value in mutation_dict.items():
            if "(" in key:
                key_parse = re.search('\((.*)\)',key).group(1)
                pmga_name = re.search('(.*) ',key).group(1)
            else:key_parse = key
            value = temp_value[0]
            if var_gene == key_parse:
                if "_all" in key:
                    if pmga_name not in amr_mutations:
                        amr_mutations[pmga_name] = {}
                    amr_mutations[pmga_name]['all_mutations'] = value
                    if mutation in value:
                        if mutation not in amr_mutations[pmga_name]:
                            amr_mutations[pmga_name][mutation] = {}
                        amr_mutations[pmga_name][mutation]["amino_acid_change"] = amino_change
                        amr_mutations[pmga_name][mutation]["type"] = mutation_type
                        amr_mutations[pmga_name][mutation]["resistance"] = antibiotic
    return amr_mutations

def combine_amr_mutations(amr,mutation):
    for key in amr:
        amr[key]['known_mutations'] = {}
        amr[key]['all_mutations'] = {}
        amr[key]['all_mutations'] = "."
        if key in mutation.keys():
            for key2 in  mutation[key].keys():
                if key2 == "all_mutations":
                    amr[key]['all_mutations'] = mutation[key]['all_mutations']
                else:
                    amr[key]['known_mutations'][key2] = mutation[key][key2]
    return amr

def get_interpretation(interpret,mutation,allele_table):
    final_marker = {};
    allele_dict = allele_table.to_dict('list')
    if "blaROB-1" in allele_dict.keys():final_marker["blaROB-1"] = "1"
    if "blaTEM-1" in allele_dict.keys():final_marker["blaTEM-1"] = "1"
    if "HAEM0118" in allele_dict.keys():final_marker["blaTEM-1"] = "1"
    if "NEIS2357" in allele_dict.keys():final_marker["blaTEM-1"] = "1"
    mutation_dict = mutation.to_dict('list')
    for sheet_name in interpret.sheet_names:
        gene_list=[]
        if "Nm_" in sheet_name:
            if "Fluoroquinolones" in sheet_name:gene_list.append("gyrA")
            if "Penicillins" in sheet_name:gene_list.append("penA")
        if "Hi_" in sheet_name:
            if "Fluoroquinolones" in sheet_name:
                gene_list.append("gyrA")
                gene_list.append("parC")
            if "Penicillins" in sheet_name:gene_list.append("ftsI")

        df = pd.read_excel(interpret,sheet_name = sheet_name).set_index('Name')
        if "Penicillins" in sheet_name:
            df.drop(columns=['blaROB-1','blaTEM-1'],inplace=True)
            df.drop(['blaROB-1','blaTEM-1'],inplace=True)
        df_drop = df.iloc[:,:-1]
        df_up = df_drop.fillna("TRUE")
        all_mutation_list = []
        for key,value in mutation_dict.items():
            for gene in gene_list:
                lookup = "("+gene+")"+"_all"
                if lookup in key:
                    all_mutation = value[0]
                    
                    all_mutation_list = all_mutation_list + all_mutation.split(", ")
        all_mutation_list.append('TRUE')
        df_res = df_up.isin(all_mutation_list)
        df_final = df_res.all(axis=1)
        for key2,value2 in df_final.items():
            if df_final[key2] == True:final_marker[key2] = "1"                
    return final_marker                 
def get_interpret_final(interpret_frame,interpret):
    antimicrobics = {}
    if not interpret_frame:
        return antimicrobics

    pp = "susceptible"
    for sheet_name in interpret.sheet_names:
        df = pd.read_excel(interpret,sheet_name = sheet_name).set_index('Name')
        for key in df.keys():
            if ":" in key:
                x = key.split(":")
                anti_agent = x[1]
                anti_class = x[0]
        for index,row in df.iterrows():
            if index in interpret_frame.keys():
                if anti_class not in antimicrobics:
                    antimicrobics[anti_class] = {}
                if anti_agent not in antimicrobics[anti_class]:
                    antimicrobics[anti_class][anti_agent] = {}
                if "markers" not in antimicrobics[anti_class][anti_agent]:
                    antimicrobics[anti_class][anti_agent]["markers"] = {}
                if index not in antimicrobics[anti_class][anti_agent]["markers"]:
                    antimicrobics[anti_class][anti_agent]["markers"][index] = {}
                antimicrobics[anti_class][anti_agent]["markers"][index] = row[-1]
              
    for key1,value1 in antimicrobics.items():
        for key2,value2 in value1.items():
            if "Resistant" in antimicrobics[key1][key2]["markers"].values():pp = "Resistant"
            elif "Intermediate" in antimicrobics[key1][key2]["markers"].values():pp = "Intermediate"
            elif "nonsusceptible" in antimicrobics[key1][key2]["markers"].values():pp = "nonsusceptible"
        antimicrobics[key1][key2]["predicted_phenotype"] = pp
    
    return antimicrobics
            
def get_summary(antimicrobics):
    from collections import defaultdict
    summ = defaultdict(list)
    pp = "testValue"
    if not antimicrobics:
        summ["predicted_resistance"] = "None"
        return summ

    for keyT,valueT in antimicrobics.items():
        for keyT2,valueT2 in valueT.items():
            if antimicrobics[keyT][keyT2]['predicted_phenotype'] == "Resistant":pp = "Resistance"
    
    if pp != "Resistance":
        summ["predicted_resistance"] = "None"
        return summ
    for key,value in antimicrobics.items():
        for key2,value2 in value.items():
            if antimicrobics[key][key2]['predicted_phenotype'] == "Resistant":
                summ['predicted_resistance'].append(key2)
                for key3,value3 in antimicrobics[key][key2]['markers'].items():
                    if value3 == "Resistant":
                        summ['resistance_genotype'].append(key3)
    for key,value in summ.items():
        if key == "predicted_resistance":
            resistance = ", ".join(value)
            summ["predicted_resistance"] = resistance
        if key == "resistance_genotype":
            markers = ", ".join(value)
            summ["resistance_genotype"] = markers
    return summ

def get_info(db_df,script_version,script_subversion):
    run_info = {};
    version_script = ("v{}.{}".format(script_version,script_subversion))
    db_df['date'] = pd.to_datetime(db_df['Date']).dt.date
    db_v = db_df['Version'].iloc[0]
    db_date = db_df['date'].iloc[0]
    date = db_date.strftime("%m-%d-%y")
    run_info["DB Version"] = db_v
    run_info["Updated"] = date
    run_info["Script_Version"] = version_script
    return run_info
