use HAIMING::GENEHISTORY::compareOne;

my $ptherlimit = shift; # perl runCompareOne.pl PTHR10000
unless ($pthrlimit){
  die "provide PTHR family id\n";
}

my $s1tree = "./PANTHER12.newick"; #location of the PANTHER12 species tree
my $p12dir = "/dir/to/PANTHER12_treetap"; 
my $treefile = "$p12dir/$pthrlimit.newick"; #location of the PANTHER12 gene trees
my $tabfile = "$p12dir/$pthrlimit.tab"; #location of the tab file for gene trees

my $outfile = "$pthrlimit.test"; # location of the output file

treeCompare($treefile,$tabfile,$s1tree,$outfile,$check);

