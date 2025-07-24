## requires python 3
## Type -h to see usage instructions

import os
import pandas as pd
import sys
import json
import argparse
import utilities
import AST_utilities
import PMGA_utilities
import openpyxl

SCRIPT_VERSION = 2 #Feb 2025
SCRIPT_SUBVERSION = 0

##getting relative path####
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

Hi_db_path = os.path.join(script_dir,"Hi_db/")
Nm_db_path = os.path.join(script_dir,"Nm_db/")

##Setup AST database:
AST_database_dict = {
    "Hi":Hi_db_path,
    "Nm":Nm_db_path
}

PMGA_species_nicknames = {
    "Hi":'hinfluenzae',
    "Nm":'neisseria'
}

def main():

    script_base = os.path.basename(__file__)
    _outputBase = '{}_v{}.{}'.format(os.path.splitext(script_base)[0],SCRIPT_VERSION,SCRIPT_SUBVERSION)
    default_logfile = os.path.splitext(script_base)[0]+'.log'

    parser = argparse.ArgumentParser(description='A program to report marker genes and AA substitutions')
    ### general info
    parser.add_argument('--version','-V',action='version',version='%(prog)s {}.{}'.format(SCRIPT_VERSION,SCRIPT_SUBVERSION))

    ### controls
    requiredArg = parser.add_argument_group('required arguments')

    requiredArg.add_argument('-o','--out_dir',help='Location to place results',required=True) 
    requiredArg.add_argument('-i','--input_dir',help='Location of PMGA json files',required=True) 
    requiredArg.add_argument('-s','--species', choices=['Nm','Hi'],help='Species to analyze',required=True)
    parser.add_argument('-x','--save_intermediate',action='store_true',help='Save intermediate excel files -- alleles, variants')
    
    args = parser.parse_args()
    try:
        if args.out_dir:
            if os.path.isdir(args.out_dir):
                outdir = args.out_dir
            else:
                os.system("mkdir {}".format(args.out_dir))
                print("Created output directory",args.out_dir)
                outdir = args.out_dir                   
        if args.input_dir: 
            if os.path.isdir(args.input_dir):
                indir = args.input_dir
            else:
                raise ValueError("Directory not present")          
        if args.species:
            specie = args.species            
        else:
            raise ValueError("Specie not known")

        ##Setup a log file.
        logFile = os.path.join(outdir,default_logfile)
        print("LogFile is : "+logFile)    
        sys.stdout = utilities.Logger(logFile)
        
        ##Record arguments
        print(_outputBase)
        print("Options are:")
        for arg in vars(args):
            print (arg, getattr(args, arg))        
        print()        
    
    except (ValueError, IOError) as e:
        print(e)
        parser.print_usage()
        
    ##Get the input data
    abs_path = os.path.abspath(indir)                            ##absolute path for input directory
    json_file = AST_utilities.get_json_file(indir)              ##getting the json file (note: will return only one json file, even if there are many)
    json_file_abs = os.path.join(abs_path,json_file)
    print("PMGA jason file: " + json_file)
    print("Absolute path: " + json_file_abs)
    temp_hold = json_file.split("_")                             ##extracting isolate id
    isolate_id = temp_hold[0]

    #Load genotyping database 
    gene_key = pd.read_excel(os.path.join(AST_database_dict[specie],'pmga_gene_key.xlsx'),engine='openpyxl').set_index('Gene')      
    print("Loaded gene key for {} with {} genes".format(specie,len(gene_key)))
    locus_rename = {r['PMGA']:'{} ({})'.format(r['PMGA'],gene) if gene != r['PMGA'] else r['PMGA'] for gene,r in gene_key.iterrows()}
    ### Load the variant list
    variants = pd.read_excel(os.path.join(AST_database_dict[specie],'AMR_variant_key.xlsx'),engine='openpyxl') ##Peptide variants
    db_version = pd.read_excel(os.path.join(AST_database_dict[specie],'AMR_variant_key.xlsx'),sheet_name= "version",engine='openpyxl') ## reading sheet with version info
    interpretation = pd.ExcelFile(os.path.join(AST_database_dict[specie],'AMR_interpretation.xlsx'),engine='openpyxl')
    variants = variants.dropna(how='all')
    #converting from float to integer
    variants['First Position of Feature in Reference Protein Sequence']= variants['First Position of Feature in Reference Protein Sequence'].astype(int)
    valid_features = ['Mutation','Region']
    valid_mutations = ['Substitution','Insertion','Deletion']
    print("Loaded {} valid features; {} valid mutations".format(sum(variants['Feature Category'].isin(valid_features)),sum(variants['Type of Mutation'].isin(valid_mutations))))
    print()
    print("Loading locus info")
   ###This will update the file path in the gene_key
    loaded_seqs = AST_utilities.load_locus_info(gene_key,variants,AST_database_dict[specie])
    
    ##Extract gene sequences for this isolate
    locus_collection = PMGA_utilities.extract_loci_from_JSON_file(json_file_abs,set(gene_key.PMGA.tolist()),annotation_species=PMGA_species_nicknames[specie])
    allele_table,new_seq_dict = PMGA_utilities.get_allele_from_JSON(locus_collection)
    if args.save_intermediate:
        PMGA_allele_file = "{}_PMGA_alleles.xlsx".format(isolate_id)
        allele_table.to_excel(os.path.join(outdir,PMGA_allele_file)) 
    ##Extracting new alleles to a fasta file 
    new_allele_fasta_file = "{}_new_allele.fasta".format(isolate_id)
    with open(os.path.join(outdir,new_allele_fasta_file),'w') as my_outfile: 
       AST_utilities.dict_to_fasta(new_seq_dict,my_outfile)
    ## Making clean gene key file for ouput#####
    clean_gene_key_file = AST_utilities.make_geneDB_file(gene_key)
    gene_database_file = "{}_gene_db.xlsx".format(isolate_id)
    clean_gene_key_file.to_excel(os.path.join(outdir,gene_database_file),index=False)

    print()
    AST_utilities.summarize_genes(locus_collection,gene_key) ###Summarize what genes were found - to sdtout 
    
    ## Get presence of gene in each genome
    locus_presence_table = PMGA_utilities.make_locus_presence_table_JSON(locus_collection, cov_filter = 1,final=True) ##ACR: we may need to discuss the cov_filter option.
    AST_utilities.fix_internal_stop(locus_presence_table,locus_rename,gene_key)     ##fix the "internal stop codon" call   SS: made this its own function
    ##Making amr genes dictionary
    amr_genes = AST_utilities.get_amr_genes(locus_presence_table,allele_table,gene_key) ##Note: new alleles are renamed as 'New' rather than 'temp_1' (in get_alleles_from_JSON)
    ## identify mutations. Note: this often generates a WARNING from Biopypthon
    mutation_frame = AST_utilities.identify_mutations_from_JSON_collection(locus_collection,variants,loaded_seqs,locus_rename)
    if args.save_intermediate:
        PMGA_mutation_file = "{}_PMGA_mutations.xlsx".format(isolate_id)
        mutation_frame.to_excel(os.path.join(outdir,PMGA_mutation_file))     
    #Making mutation dictionary for each gene
    amr_mutation = AST_utilities.get_mutation(mutation_frame,variants)
    ##getting relevant markers responsible for reduce suspetibility                
    interpret_frame = AST_utilities.get_interpretation(interpretation,mutation_frame,allele_table)
    ## making dictionary of antimicrobic class,drug and markers 
    antimicrobics = AST_utilities.get_interpret_final(interpret_frame,interpretation)
    ## making top level summary for resisatnt phenotype and responsible markers
    run_summary = AST_utilities.get_summary(antimicrobics)
    ##combining amr genes and mutation dictionaries
    combine_dict = AST_utilities.combine_amr_mutations(amr_genes,amr_mutation)
    ##Generating run summary
    run_info = AST_utilities.get_info(db_version,SCRIPT_VERSION,SCRIPT_SUBVERSION)
    ##Combining all the dictionaries  to one single Json file
    amr_data = {'amr_genes':combine_dict,'run_info':run_info,'antimicrobics':antimicrobics,'summary':run_summary}
    amr_data_file = "{}_amr_data.json".format(isolate_id)
    with open(os.path.join(outdir,amr_data_file),'w') as my_amrFile:
        json.dump(amr_data,my_amrFile)

if __name__ == "__main__":
    if not utilities.has_preferred_python():
        raise Exception("Upgrade your python version")
    main()