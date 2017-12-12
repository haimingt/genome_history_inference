# genome_history_inference
Scripts to infer the evolutionary history of gene gains and losses  by comparing phylogenetic gene tree with a reference species tree

PREREQUISITES:
-Perl 
-Perl module BioPerl http://bioperl.org/
-Add modules HAIMING/GENEHISTORY/CompareOne.pm and HAIMING/GENEHISTORY/Summarize.pm to PATH

Files in this project:

HAIMING/GENEHISTORY/PANTHER12.newick                        
The reference species tree for PANTHER12


PANTHER12_treetab/PTHR#####.newick                          
The phylogenetic gene tree for each PANTHER family

PANTHER12_treetab/PTHR#####.tab                             
Node info for corresponding phylogenetic gene tree	

PANTHER12_treeAsum/PTHR#####.sum.txt                        
Pre-computed results for gene gains and losses

Only one PANTHER12 family is listed here, for the other families, go to ftp://pantherdb.org/humangenomehistory 


HAIMING/GENEHISTORY/CompareOne.pm                           
Perl module for comparing one PANTHER gene tree with the reference species tree

HAIMING/GENEHISTORY/Summarize.pm                            
Perl module for running all PANTHER gene trees in PANTHER12_treetab directory using CompareOne.pm


runCompareOne.pl                                            
Sample Perl script for running CompareOne.pm

runSummarize.pl                                             
Sample Perl script for runnning Summarize.pm


For details, please check paper .....
If you find this useful, please cite.
