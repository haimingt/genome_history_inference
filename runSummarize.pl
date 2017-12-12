use HAIMING::GENEHISTORY::summarize;

my $s1tree = "./PANTHER12.newick"; #location of the PANTHER12 species tree
my $p12dir = "/dir/to/PANTHER12_treetap"; #location of the PANTHER12 gene trees and corresponding node files
my $outdir = "/dir/to/PANTHER12_treeAsum";

runEach($p12dir,$s1tree,$outdir); #run compareOne.pm on each gene tree of the directory $p12dir 


my $sumFile = "PANTHER12_evolutionary_history_summary.txt"; # the summary file for evolutionary history across all families
my $horizFile = "PANTHER12_horizontal_transfer_summary.txt" # the summary file for the horizontal transfers across all families

historySum($outdir,$sumFile,$horizFile);
