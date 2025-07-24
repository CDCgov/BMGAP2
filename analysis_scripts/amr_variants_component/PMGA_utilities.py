#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 14:48:54 2018
@author: ymw8

Edited on Fri Feb 07 2025
@author: ujz6
Made changes to accomodate RHEL8 migration
"""

from collections import defaultdict 
import json
import os
import utilities
import pandas as pd
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord

Nm = 'Neisseria meningitidis'
Hi = 'Haemophilus influenzae'

PMGA_species_nicknames = {
    Hi:'hinfluenzae',
    Nm:'neisseria'
}

##Return a dict with keys from this json file -- should simply be one  genome assembly.
## The value for the ruslt is a dict of loci with a list containing each instance of the locus (this is just the contents of the json, with contig key removed)
def extract_loci_from_JSON_file(json_file,locus_set=None,annotation_species='hinfluenzae'):
    parsed_json = json.load(open(json_file))
    if len(parsed_json.keys()) > 1:
        print("Warning: Parsing JSON with {} records".format(len(parsed_json.keys())))
        print(json_file)
    found_loci = defaultdict(list)
    result = {}
    for k in parsed_json.keys(): ##should be one key with one genome assembly
        json_contigs = parsed_json[k]['contigs']
        json_species = parsed_json[k]['species'] ##This is the scheme used for annotation
        if json_species == 'Haemophilus influenzae' : json_species = 'hinfluenzae'
        if json_species.upper() != annotation_species.upper():
            raise ValueError("Annotation file using the wrong species loci ({}). \n{}".format(json_species,json_file))
        for contig in json_contigs:
            locus_in_contig = set(json_contigs[contig].keys())
            locus_overlap = locus_in_contig if (locus_set is None) else locus_set.intersection(locus_in_contig)
            for locus in locus_overlap:
                if len(json_contigs[contig][locus]) > 0:
                    found_loci[locus] += json_contigs[contig][locus]
        result[k] = found_loci
    return result

##Default locus_set is every identified locus.
## Merges the dicts from extract_loci_from_JSON_file -- end result should be a dict with one record for each genome.
def extract_loci_from_JSON_collection(file_list,locus_set=None,annotation_species='hinfluenzae'):
    all_files = {}
    existing_keys = set()
    failures = []
    for filename in file_list:
        try:
            extracted = extract_loci_from_JSON_file(filename,locus_set,annotation_species) ##Dict with genome as key
        except Exception as e:
            failures.append(filename)
            print("Error")
            utilities.printExceptionDetails(e)
        else:
            extracted_keys = set(extracted.keys())
            if len(extracted_keys.intersection(existing_keys)) > 0:
                raise ValueError("Multiple files use the same assembly identifier")
            else:
                all_files.update(extracted)
                existing_keys = set(all_files.keys())
    for f in failures:
        basename = os.path.basename(f)
        print(basename)
        print(basename.split('_')[0])
    if isinstance(locus_set,set):
        found_set = set([k for locus in all_files.values() for k in locus.keys()])
        unfound_set = (locus_set.difference(found_set))
        if len(unfound_set) >0 :
            for unfound_locus in unfound_set:
                print("Unable to find any sequence for {}".format(unfound_locus))
    return all_files

###convert to presence/absence/flag table
def get_locus_presence_andAllele_JSON(extracted_collection,locus_rename=None):
    restructure = {}
    for a in extracted_collection.keys(): #file
        restructure[a] = {}
        for b in extracted_collection[a]: #locus
            locus_name = locus_rename[b] if isinstance(locus_rename,dict) else b
            output = {}
            output['alleles'] = ','.join([x['allele_id'] for x in extracted_collection[a][b]])
            if len(extracted_collection[a][b]) > 1:
                output['status'] = 'Multiple'
            elif len(extracted_collection[a][b]) == 0:
                output['status'] = 'Missing'
            else: # == 1
                if len(extracted_collection[a][b][0]['flags']) > 0:
                    output['status'] = ', '.join(extracted_collection[a][b][0]['flags'])
                else:
                    output['status'] = 'Present' #since there is not always a specific allele id
            for k,v in output.items():
                restructure[a][locus_name+'_'+k] = v
    result = pd.DataFrame(restructure).transpose().fillna('Missing')
    return result

###convert to presence/absence/flag table
##Updated to report maximum coverage for incomplete alleles; use "final=False" for 'raw' json files.
def make_locus_presence_table_JSON(extracted_collection,cov_filter = 0.9,final=True):
    restructure = {}
    for assembly in extracted_collection: #file
        restructure[assembly] = {}
        for locus_id, locus_list in extracted_collection[assembly].items(): #locus
            filtered_list = [x for x in locus_list if x['cov'] >= cov_filter]
            output = ''
            if len(filtered_list) > 1:
                output = 'Multiple'
            elif len(filtered_list) == 0:
                if (len(locus_list) == 0):
                   output = 'Missing'
                else: 
                   my_max = max([x['cov'] for x in locus_list]) 
                   output = 'Partial_only_{:0.3f}'.format(my_max)
            else: # == 1
                if final: ##"raw" json does not have flags.
                    if len(filtered_list[0]['flags']) > 0:
                        output = ', '.join(filtered_list[0]['flags'])
                    else:
                        output = 'Present' #since there is not always a specific allele id
                else:
                    output = 'Present'
            restructure[assembly][locus_id] = output
    result = pd.DataFrame(restructure).transpose().fillna('Missing')
    return result

##Updated to provide temp allele ID for new alleles
## Returns both a table with allele identifiers, and a dict defining new alleles
def get_allele_from_JSON(extracted_collection,allele_dict=None,locus_rename=None,cov_filter=0.9):
    if allele_dict is None:
        allele_dict = defaultdict(dict)
    if not isinstance(allele_dict,dict):##Should be a dict of dict: keys are locus, and sequences, value is ID
        raise ValueError("Illegal allele dictionary")
    restructure = {} ##We are collapsing multiple hits for the same gene
    for assembly in extracted_collection.keys(): #file
        restructure[assembly] = {}
        for locus_id, locus_list in extracted_collection[assembly].items(): #locus
            filtered_list = [x for x in locus_list if x['cov'] > cov_filter]
            locus_name = locus_rename[locus_id] if isinstance(locus_rename,dict) else locus_id
            allele_list = []
            for this_locus in filtered_list:
                if this_locus['allele_id'].startswith("new_allele"):
                    seq_str = ''
                    seq_str = this_locus['qseq']
                    if seq_str not in allele_dict[locus_id]:
                        incomplete = (this_locus['cov'] < 1.0) and (translate(seq_str[-3:]) != '*')
                        prefix = 'incomplete_' if incomplete else 'temp_'
                        allele_dict[locus_id][seq_str] = '{}{}'.format(prefix,len(allele_dict[locus_id])+1)
                        if incomplete:
                            print('{} {}'.format(locus_id,allele_dict[locus_id][seq_str]))
                        if (len(seq_str) % 3 != 0):
                            print("Warning: new sequence {} for {} is not multiple of three".format(allele_dict[locus_id][seq_str],locus_id))
                    allele_list.append(allele_dict[locus_id][seq_str]) 
                else:
                    allele_list.append(this_locus['allele_id'])
            if len(allele_list) > 0:
                restructure[assembly][locus_name] = ','.join(sorted(allele_list))
    result = pd.DataFrame(restructure).transpose().fillna('Missing')
    return result,allele_dict

###Troubleshoot: If you have trouble with the BMGAP JSON
script_header = [
'#!/bin/bash -l',
'#$ -cwd',
'#$ -q dbd.q,all.q ',
'#$ -pe smp 8-16',
'module load Python/3.4',
'module load ncbi-blast+/LATEST',
'module load Mash/1.1',
'echo Starting'    
]

##TODO: Doesn't look like used anymore, remove
##Note, the script and guide files must be relative to the same location. The script will contain 
##  a guide file location taht is relative to its own directory.--excel format
#def make_pmga_script(script_file,assembly_group,guide_file=None,sequence_directory=None):
#    with open(script_file,'wt') as script_out:
#        for h in script_header:
#            print(h,file=script_out)
#        pmga_command = 'python3 /scicomp/groups/OID/NCIRD/DBD/MVPDB/ML/tools/PMGA/blast_pubmlst.py -o pmga_{} -t $NSLOTS -sg '.format(assembly_group)
#        if isinstance(guide_file,str):
#            if os.path.isfile(guide_file):
#                relative_guide = os.path.relpath(guide_file,os.path.dirname(script_file))
#                pmga_command += '-x {}'.format(relative_guide)
#            else:
#                print("Error: guide file is not a file")
#        elif isinstance(sequence_directory,str):
#            if os.path.isdir(sequence_directory):
#                pmga_command += '-d {}'.format(sequence_directory)
#        print(pmga_command,file=script_out)
#    print('Wrote script to {}'.format(script_file))

##Methods for counting results in allele table
count_fields = ['incomplete','missing','multiple','known_allele','temp_allele']
def is_multiple(x):
    return ',' in x

def is_temp(x):
    return x.startswith('temp_')

def is_incomplete(x):
    return x.startswith("incomplete_")

def summarize_allele_table(allele_table):
    a = allele_table.copy()
    if "known_allele" in a.columns:
        a.drop('known_allele',axis=1,inplace=True)
    known_bool = (a.applymap(str.isnumeric))
    temp_bool = a.applymap(is_temp)
    multi_bool = (a.applymap(is_multiple))
    missing_bool = (a == "Missing")
    inc_bool = (a.applymap(is_incomplete))

    a['incomplete'] = inc_bool.sum(axis=1).astype(int)
    a['known_allele'] = known_bool.sum(axis=1).astype(int)
    a['multiple'] = multi_bool.sum(axis=1).astype(int)
    a['temp_allele'] = (temp_bool & ~multi_bool).sum(axis=1).astype(int)
    a['missing'] = missing_bool.sum(axis=1).astype(int)
    return a

def valid_alleles(x):
    result = x.isnumeric() or x.startswith("temp_") or (',' in x) ## I'm not really sure if I should be counting multipules. They seem
    return result

##Each column of df has values that are expected to match across rows. There should be no other columns or rows
## Pass a function for screening out valid values for matching
def proportion_mismatch(my_series,valid_func):
    def choose_two(x):
        return ((x**2)-x)/2

    valid_values = my_series.apply(valid_func)
    counts = my_series[valid_values].value_counts()
    total = choose_two(counts.sum())
    if total == 0:
        print("No comparisons for {}".format(my_series.name))
        return 'N/A'
    matching = sum([choose_two(x) for x in counts.values])
    mismatch = total-matching
    return mismatch/total

def valid_counts(my_series,valid_func):
    valid_values = my_series.apply(valid_func)
    counts = my_series[valid_values].value_counts()
    return counts.sum()

def extract_sequences(extracted_collection,locus_rename_dict):
    loci_sequences = defaultdict(list)
    for filename,loci_dict in extracted_collection.items():
        lid = filename.split("_")[0]
        for locus,locus_list in loci_dict.items():
            seq_str = ''
            if len(locus_list) > 0:
                seq_str = locus_list[0]['qseq']
                if len(locus_list) > 1:
                    seqs_set = set([locus_item['qseq'] for locus_item in locus_list])
                    if len(seqs_set) > 1:
                        print("Warning: found {} sequenes at {} loci for {} in {}".format(len(seqs_set),len(locus_list),locus,lid))            
                        seq_str = ''
                    else:
                        print("Notice: found {} loci with {} sequence for {} in {}".format(len(locus_list),len(seqs_set),locus,lid))
    
            if len(seq_str) > 0:
                seq = Seq(seq_str)
                seqR = SeqRecord(seq,id=lid,name='{}_{}'.format(lid,locus),description=locus_rename_dict[locus])
                loci_sequences[locus].append(seqR)
    return loci_sequences
