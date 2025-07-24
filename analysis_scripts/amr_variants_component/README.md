This project contains the functions and reference data for finding AMR genes and calling AMR-related genetic variants from genome
assemblies. There are databases for Hi and Nm. Configuration and usage are as follows:

###### Usage example:
module load Python/3.4\
python runAST.py -i 'inputDir containing PMGA json file' -s [Hi|Nm] -o 'outdir to store results'

###### Output Files:
1. Json  -  One Json file containing information on AMR genes found in the assemly and the variants present.&emsp; [Example Json output](docs/test_amr_data.json) 
2. Fasta -  fasta file of all the new alleles found in the assembly. &emsp; [Example fasta file](docs/test_new_allele.fasta)
3. Excel -  List of genes/markers evaluated by the database. &emsp; [Example Nm list](docs/test_Nm_gene_db.xlsx), [Example Hi list](docs/test_Hi_gene_db.xlsx)
4. Excel -  Excel tables with intermediate results (alleles, mutations). Exported with '-x' flag. &emsp; [Example allele file](docs/test_PMGA_alleles.xlsx)

###### Structure of the JSON output file:
The json file contains four sections/dictionaries (Fields marked asterisks will be in the BMGAP export table):
- summary: This dictionary contains the top-line, human readable results.\
  &emsp;         predicted_resistance*:            "a string with comma delimited list of antimicrobial agents predicted to be resistant"\
  &emsp;         resistance_genotype*:            "a string with comma delimited list of markers that are responsible for resistance to any antimicrobic."

- amr_genes: This dictionary contains the information on AMR genes and mutations/variants found in the genome assembly (keyed on 'gene id')\:   
  'gene id' :\
  &emsp;	 	Gene_name:         "name of the gene"\
  &emsp;        Antimicrobial: list of antimicrobic classes that may be affected by this gene (to guide the BMGAP display)\
  &emsp;        allele*:            "numeric identifier for the sequence of the gene ('NA' if gene is absent)"\
  &emsp;        status*:            "detection of full length gene"\
  &emsp;	 	all_mutations*:     "list of all amino acid replacements found in the gene ('.' if no mutation found)"\
  &emsp;	 	known_mutations:  "a dictionary of mutations that have been previously reported to affect resistance"\
  &emsp;&emsp;   	mutation_name:\
  &emsp;&emsp;&emsp;	type:              "type of mutation (substitution, insertion, etc)"\
  &emsp;&emsp;&emsp;	amino_acid_change: "amino acid replacement with respect to the reference"\
  &emsp;&emsp;&emsp;	resistance: "antimicrobic resistance potentially associated with this mutation"

- antimicrobics: This dictionary contains interpretations for individual antimicrobics (keyed on antimicrobic class).\
  'antimicrobic class' : (a class is only shown if the genome contains a marker associated with reduced susceptibility)\
  &emsp;    'antimicrobial agent' :\
  &emsp;&emsp;	 	predicted_phenotype*:         Dictionary with drug as key and phenotype as value ("susceptible/intermediate/resistant/nonsusceptible")\
  &emsp;&emsp;	 	markers*:         "relevant genetic markers resulting in reduced susceptibility for the drugs in this class"

- run_info : This dictionary contains the run information.\
  DB_version*:      "database version"\
  updated*:         "date when database last updated"\
  script_version*:  "version of the script" 


###### Database configuration files:
Each database consists of a subdirectory of the AMR_variants_component directory. It contains the following files.
1.  pmga_gene_key.xlsx: The critical fields here are: the gene name (used for reporting), the PMGA identifier (to locate gene in PMGA results), and the location of a reference allele file in the DB directory (optional: for aligning allels and identifying specific variants), Coding (for specifying whether sequences should be translated for alignment)
2.  AMR_variant_key.xlsx: The critical fields here are: gene name and PMGA identier (matching the pmga_gene_key), First/Last Position of Feature in Reference Protein Sequence, Feature Category/Type of Mutation (for specifying the specific mutation being examined); Susceptible/Resistant (for specifying which amino acids or nucleotides are potentially resistant)
3.  AMR_interpretation.xlsx: Data for final reporting of AMR genotype.
  
###### Database status (as of July 16 2020):

*Neisseria meningitidis* (Nm)
- 11 acquired resistance genes
- 5 modified target genes
- Focus: Penicillin, Ampicillin, Ciprofloxacin
- Unvalidated: Cefotaxime, Ceftriaxone, Meropenem, Rifampin


*Haemophilus influenzae* (Hi)
- 11 acquired resistance genes
- 9 modified target genes
- Focus: Ampicillin
- Unvalidated: Ampicillin/Sulbactam, Amoxicillin/Clavulanate, Cefotaxime, Ceftriaxone, Rifampin

###### Testing AMR_variant_component installation using known sequence data:
The folders "Nm_controls" and "Hi_controls" contain FASTA files with each of the genes that we are testing for, along with the analysis results (both PMGA and AMR modules). The subfolders are:
- sequences: FASTA files for each gene. Some genes have multiple alleles representing the important variation -- in this case, each allele is in a separate FASTA file.
- pmga: PMGA output for these sequences. Generate this folder by running 'pmga_test.sh' with your pmga installation.
- AMR: Output from the AMR testing module. These are the correct results for the sequences in the FASTA files. Generate this folder by running 'amr_test.sh'
