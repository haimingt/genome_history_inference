use HAIMING::GENEHISTORY::compareOne;

my $ptherlimit = shift; # perl runCompareOne.pl PTHR10000
unless ($pthrlimit){
  die "provide PTHR family id\n";
}

my $p12dir = "/dir/to/PANTHER12_treeAsum";
my $treefile = "$p12dir/$pthrlimit.newick";
my $tabfile = "$p12dir/$pthrlimit.tab";

my $outfile = "$pthrlimit.test";
my $outtree = "$pthrlimit.newick";

treeCompare($treefile,$tabfile,$s1tree,$outfile,$check);

