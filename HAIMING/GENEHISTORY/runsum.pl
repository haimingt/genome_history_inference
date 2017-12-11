use HAIMING::GENEHISTORY::summarize;
my $p12dir = "/home/pmd-02/pdt/pdthomas/users/haiming/PANTHER12_treeAsum";
my $treefile = "$s1dir/$pthrlimit.newick";
my $tabfile = "$s1dir/$pthrlimit.tab";

my $outfile = "$pthrlimit.test";
my $outtree = "$pthrlimit.newick";


my $oritree = "/home/pmd-02/pdt/pdthomas/panther/famlib/dev/UPL/SuperPANTHER1.0/lib_1.0/books/$pthrlimit/orig/tree/$pthrlimit.orig.tree";

unless (-e $outtree){
  tonewick($oritree,$outtree);
}
treeCompare($treefile,$tabfile,$s1tree,$outfile,$check);

