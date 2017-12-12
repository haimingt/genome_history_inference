#!/usr/bin/perl
package HAIMING::GENEHISTORY::compareOne;

# This is the new version of tree comparison module
# the old one seems to have some probelm with loss in an ancestral branch #use warnings;

#should input a tree file in newick format
# a plain file which links the tree's node ids to species name:
# first column node id, second column species name, third column node type
# node type: leaf, speciation, duplication, horizontal transfer
# it is good if correct ancestral species is found for the duplication node, speciation node or horizontal transfer, but it is OK if not inputted
# a species tree with which the gene tree structure is consistent
# besides the species name should in the plain file and species tree should match

use strict;
use warnings;
use Bio::TreeIO;
use Data::Dumper;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(treeCompare);
my %nodetype;
my %nodespecies;
my %infostore;
my $gtree;
my $stree;
my $groot;
my %momchild;
my %visited;
my $pthr;
my $check;
my $childrenaref = ();

sub treeCompare{

  %nodetype = ();
  %nodespecies = ();
  %infostore = ();
  %momchild = ();
  %visited = ();

  my $treefile = shift;
  my $desfile = shift;
  my $speciesfile = shift;
  my $outfile = shift;
  $check = shift;

  if ($treefile =~ /(PTHR[0-9]+)/){
    $pthr = $1;
  }
  open OUT, "> $outfile" or die "cannot output to $outfile\n";
  open DIN, "< $desfile" or die "cannot open the node description file $desfile\n";
  while(<DIN>){
    chomp;
    my @array = split(/\t/);
    unless ($array[2]){
      die "missing type for $array[0]\n";
    }
    $nodetype{$array[0]} = $array[2];
    if ($array[1]){
      $nodespecies{$array[0]} = $array[1];
    }
    else{
      if ($array[2] =~ /leaf/i){
	die "no species name for leaf node $array[0]\n";
      }
    }
  }
  close DIN;
  
  my $sinput = Bio::TreeIO->new(-file => $speciesfile, -format => "newick");
  $stree = $sinput->next_tree;
  my $ginput = Bio::TreeIO->new(-file => $treefile, -format => "newick");
  $gtree = $ginput->next_tree;
  $groot = $gtree->get_root_node();
  my $gid = $groot->id;
  unless ($gid){
    $groot = ($groot->each_Descendent)[0];
    #print STDERR "The $treefile 's root has no id, start with its child\n";
  }
  &getspecies;
  &updatedesfile($desfile);
  &compareUseNode($groot,1);
  &validate($treefile);

  foreach my $key (sort keys %infostore){
    print OUT "$key\t";
    foreach my $type (sort keys %{$infostore{$key}}){
      my $value = $infostore{$key}{$type};
      print OUT "$type:$value\t";
    }
    print OUT "\n";
  }
  close OUT;
}
sub updatedesfile{
  my $desfile = shift;
  my $target = $desfile.".up";
  open DOUT, "> $target" or die;
  open DIN, "< $desfile" or die "cannot open the node description file $desfile\n";
  while(<DIN>){
    chomp;
    my @array = split(/\t/);
    unless ($array[1]){
      my $s = $nodespecies{$array[0]};
      print DOUT "$array[0]\t$s\t$array[2]\n";
      print STDERR "findspecies: $array[0]\t$s\t$array[2]\n";
      unless (($array[2] =~ /ROOT/i) and ($array[2] =~ /horiz/i)){
	unless ($s){
	  die "no species found for $array[0]\n";
	}
      }
    }
    else{
      print DOUT "$_\n";
    }
  }
  close DIN;
  close DOUT;
}
sub getspecies{
  my @horiznodes;
  my @dupnodes;
  my @nodes = $groot->get_all_Descendents;
  @nodes =reverse @nodes;
  push(@nodes,$groot);
  foreach my $node (@nodes){
    my $nid = $node->id;
    my $ntype = $nodetype{$nid};
    unless ($node->id){
      die "existing node without id\n";
    }
    unless (exists $nodetype{$node->id}){
      die "existing node without type\n";
    }
    unless (exists $nodespecies{$node->id}){

      if ($ntype =~ /horiz/i){
	unless ($ntype =~ /root/i){
	  push(@horiznodes,$node);
	}
      }
      elsif ($ntype =~ /dup/i){
	push(@dupnodes,$node);
      }
      else{
	die "why $ntype $nid has no species?\n";
      }
    }
  }

  foreach my $node (@dupnodes){
    my $nid = $node->id;
    print STDERR "find species for dup node $nid\n";
    my $nspecies = &duplabel($node);
    unless ($nspecies){
      die "cannot find the species for $nid\n";
    }
    else{
      $nodespecies{$nid} = $nspecies;
    }
  }

  foreach my $node (@horiznodes){
    my $nid = $node->id;
    # print STDERR "find species for horiz node $nid\n";
    my $nspecies = &horizlabel($node);
    unless ($nspecies){
      die "cannot find the species for $nid\n";
    }
    else{
      $nodespecies{$nid} = $nspecies;
    }
  }
}

sub intricate{
  my $node = shift;
  my $mspecies = shift;
  my @children = $node->each_Descendent;
  my $mom = $node->ancestor || return; # so this is the root of the tree
  return unless ($mom->id); #so this is the root of the tree
  my @siblings;

  foreach my $child ($mom->each_Descendent){
    next if ($child eq $node);
    my $s = $nodespecies{$child->id};
    if ($s){
      push(@siblings,$s);
    }
  }
  my @sets;

  foreach my $child (@children){
    my $indicator =0;
    foreach my $sibling (@siblings){
      if (&isMom($sibling,$nodespecies{$child->id})){
	$indicator =1;
	last;
      }
    }
    if ($indicator ne 1){
      push(@sets,$nodespecies{$child->id});
    }
  }
  
  my $hspecies = &commonancestor(\@sets);
  if (($hspecies eq $mspecies) or isMom($hspecies,$mspecies)){
    return $sets[0];
  }
  return $hspecies;
}


sub horizlabel{
  my $node = shift;
  print "\n\n********find label for horiz".$node->id."\n";

  if ($node->id eq "AN77"){
    #    print "stop here\n";
  }
  my $current = $node;
  my $mspecies;
  my $n =0;
  my $mom;
  while($n<20){
    $mom = $current ->ancestor || last ;
    $mspecies = $nodespecies{$mom->id};
    if ($mspecies){
      last;
    }
    else{
      $current = $mom;
    }
    $n++;
  }

  unless ($mspecies){
    #it means the mom is perhaps a duplication node
    my $mom = $node->ancestor;

    foreach my $child ($mom->each_Descendent){
      next if ($child eq $node);
      if (exists $nodespecies{$child->id}){
	$mspecies = $nodespecies{$child->id};
	last;
      }
    }
  }
  if ($mspecies ne "LUCA"){
    $childrenaref = ();
    &getleavestopdown($node);

    my @sets; my$n = 0;
    my $size = @$childrenaref;

    foreach my $nodeid (@$childrenaref){
      if (isMom($mspecies,$nodeid) or ($mspecies eq $nodeid)){
	push(@sets,$nodeid);
	$n++;
      }
    }
    print "all horizlabel find:\t".join(";",@sets);
    print "\nmomspecies is $mspecies\t";

    unless ($n eq $size){
      my $hspecies = &commonancestor(\@sets);
      print "\nhorizspecies is $hspecies\n";
      return $hspecies;
    }
    else{
      my $hs = &intricate($node,$mspecies);
      if ($hs){
	return $hs;
      }
      die "really really need manual check\n";
    }
  }
  else{
    # we cannot get the species for the common ancestor,
    # we need to find the sibling species
    my $hs = &intricate($node,$mspecies);
    if ($hs){
      return $hs;
    }
    die "really really need manual check! HORIZ label\n";
  }
}


sub getleavestopdown{
  my $node = shift;

  foreach my $des ($node->each_Descendent){
    my $desid = $des->id;
    
    my $species = $nodespecies{$desid};
    my $nodetype = $nodetype{$desid};
    print "getleavestopdown: $desid\t$nodetype\t$species\n";
    if ($species){
      push(@$childrenaref,$species);
    }
    else{
      unless($nodetype =~ /horiz/i){
	&getleavestodown($des);
      }
    }
  }
}


sub duplabel{
  my $node = shift;

  my @sets;
  foreach my $child ($node->each_Descendent){
    my $cspecies = $nodespecies{$child->id};
    
    unless ($cspecies){
      if ($nodetype{$child->id} =~ /dup/i){
	$cspecies = &duplabel($child);
      }
      elsif ($nodetype{$child->id} =~ /spe/i){
	$childrenaref = ();
	&getleavestopdown($child);
	my @sset;
	foreach my $nodeid (@$childrenaref){
	  my $species = $nodespecies{$nodeid};
	  if ($species){
	    push(@sset,$species);
	  }
	}
	$cspecies = &commonancestor(\@sset);
      }
    }
    if ($cspecies){
      push(@sets,$cspecies);
    }
  }

  my $hspecies = &commonancestor(\@sets);
  return $hspecies;
}


sub commonancestor{
 
  # should start from the youngest!
  my $aref= shift;
  my $size = @$aref;
  my %diffs;
  foreach my $nodeid(@$aref){
    $diffs{$nodeid} =1;
  }
  my @keys = keys %diffs;
  my $ksize = @keys;
  if ($ksize ==1){
    return $keys[0];
  }
  elsif ($ksize <1){
    die "input $aref is empty!!Check\n";
  }
  my %store;

  foreach my $nodeid (@$aref){
    my $node = $stree->find_node(-id=>$nodeid);
    $store{$nodeid}++;
    my $n = 0; my $current = $node;
    unless ($node){
      die "no node found for id $nodeid\n";
    }
    while($n<50){
      my $mom = $current->ancestor || last;
      my $momid = $mom->id;
#      last unless $momid;
      $store{$momid}++;
      $current = $mom;
      $n++;
    }
  }
  print Dumper(\%diffs);
  print Dumper(\%store);
  foreach my $nodeid (@$aref){
    my $node = $stree->find_node(-id=>$nodeid);
    if ($store{$nodeid} >= $size){
      print "return $nodeid\n";
      return $nodeid;
    }
    my $n =0;
    while($n<50){
      my $mom = $node->ancestor ||last;
      my $momid = $mom->id;
      if ($store{$momid} >= $size){
	print "return $nodeid\n";
	return $momid;
      }
      $node = $mom;
      $n++;
    }
  }
}

sub validate{
  # just a validate, how many leaves with the sum file return and if the number is consistent with the tree file
  my $treefile = shift;
  my %genenum;
  my %sumnum;
  my %fromsum;

  my @nodes = $groot->get_all_Descendents;
  @nodes =reverse @nodes;
  push(@nodes,$groot);

  foreach my $node (@nodes){
    next unless $node ->is_Leaf;
    my $species = $nodespecies{$node->id};
    $genenum{$species}++;
  }
  foreach my $node (keys %infostore){
    next unless ($node);
    foreach my $type (keys %{$infostore{$node}}){
      my $value = $infostore{$node}{$type};
      my @des;
      my $treenode = $stree->find_node(-id=> $node);
      unless ($treenode->is_Leaf){
	foreach my $child ($treenode->get_all_Descendents){
	  next unless ($child->is_Leaf);
	  push(@des,$child->id);
	}
      }
      else{
	push(@des,$node);
      }
      if ($type eq "L"){
	foreach my $des (@des){
	  $sumnum{$des} -= $value;
	  $fromsum{$des}{$node} .= ";$type:".$value;
	}
      }
      else{
	foreach my $des(@des){
	  $sumnum{$des} += $value;
	  $fromsum{$des}{$node} .= ";$type:".$value;
	}
      }
    }    
  }

#  print Dumper(\%infostore);
 # print Dumper(\%genenum);
 # print Dumper(\%sumnum);
  #print Dumper(\%fromsum);

  my $indicator =0;
  foreach my $key (keys %genenum){
    my $gn = $genenum{$key};
    my $sn = $sumnum{$key};
    if ($gn eq $sn){
      #print STDERR "$pthr CORRECT numbers !in $key should be $gn, actual is $sn\n";
    }
    else{
      print STDERR "$pthr WRONG numbers !in $key should be $gn, actual is $sn\n";
      print STDERR Dumper(\$fromsum{$key});

      $indicator ++;
    }
  }

  if ($indicator > 4){
    die "validation didnot pass for $treefile";
  }
  print STDERR "validation passed\n";
}

sub compareUseNode{
  my $cnode = shift;
  my $isroot = shift;
  my $cid = $cnode->id; 
  #print "compareUseNode on $cid\n";
  my $cspecies = $nodespecies{$cid};
  my $ctype = $nodetype{$cid};
  if ($ctype =~ /dup/i){
    my @children = $cnode->each_Descendent;
    foreach my $child (@children){
      &lossCompare($child,$cspecies);
    }
    my $N = scalar @children;
    $infostore{$cspecies}{"G"} += $N-1;
    if ($isroot == 1){
      $infostore{$cspecies}{"R"} ++;
    }
  }
  elsif ($ctype =~ /spe/i){
    if ($isroot == 1){
      $infostore{$cspecies}{"R"} ++;
    }
    &lossCompare($cnode,$cspecies);
  }
  elsif ($ctype =~ /leaf/i){
    return;
  }
  elsif ($ctype =~ /horiz/i){
    my @children = $cnode->each_Descendent;
    if ($isroot){ # so if the root is horizontal transfer event
      # then consider all of the children come from horizontal transfer
      foreach my $child (@children){
	my $childid = $child->id;
	my $cspecies = $nodespecies{$childid};
	$infostore{$cspecies}{"HR"}++;
      }
    }
    else{
      # so the species for a horizontal transfer event is the one without the transfered event
      foreach my $child (@children){
	my $cid = $child->id;
	my $chspecies = $nodespecies{$cid};
	#my $momspecies = $nodespecies{$child->ancestor->id};
	
	if (isMom($cspecies,$chspecies) or ($cspecies eq $chspecies)){
	  &lossCompare($child,$cspecies);
	}
	else{ # this is gained from horizontal transfer
	  print OUT "HorizTO\t$chspecies\tHorizFrom\t$cspecies\n";
	  $infostore{$chspecies}{"H"} ++;
	  # as this horizontal transfer node is not root, it must has a mom!!
	  
	}
      }
    }
  }
  else{
    die "unknown type: $ctype for node $cid\n";
  }
  
  foreach my $childnode ($cnode->each_Descendent){
    &compareUseNode($childnode,0);
  }
}

sub lossCompare{
  my $node = shift;
  my $node_br = shift;
  my $specialcase = shift;
  if (exists $visited{$node->id}){
    unless ($specialcase){
      return;
    }
  }
  else{
    $visited{$node->id} =1;
  }

  unless ($node_br){
    die "error, no branch sepecies found for a node, check again!!! \n";
  }

  my %setB;
  my $snode = $stree->find_node(-id => $node_br);
  unless ($snode){
    die "cannot find the internal node for species $node_br\n";
    return;
  }

  my $ntype = $nodetype{$node->id};
  my $nspecies = $nodespecies{$node->id};
 
  print "lossCompare: ".$node->id." with $node_br\n";

  if (isMom($node_br,$nspecies)){
    print "losstrace: ".$node->id." $nspecies, with $node_br\n";
    &losstrace($node_br,$nspecies);
    &lossCompare($node,$nspecies,1);
    return;
  }
  if ($node_br ne $nspecies){
    die "BR $node_br and node species $nspecies should be the same\n";
  }

  if ($ntype =~ /dup|hori/i){
    return;
  }

  my %setA; # node 's direct descendants and their corresponding species
  foreach my $child ($node->each_Descendent){
    my $cid = $child->id;
    my $species = $nodespecies{$cid};
    $setA{$species} ++;
  }

  foreach my $child ($snode->each_Descendent){
    my $species = $child->id;
    $setB{$species}++;
  }

  my $match;
  my %match;


  if (%setA ~~ %setB){
    return;
  }

  my $indicator =1;
  foreach my $sb (keys %setB){ #direct descendent of B, if any of them is mom of node species
    last if $sb eq $nspecies;
    if (isMom($sb,$nspecies)){
      &losstrace($sb,$nspecies);
      $match{$sb} =1;
      $indicator =2;
      last;
    }
  }
  if ($indicator == 2){
    foreach my $sb (keys %setB){
      unless (exists $match{$sb}){
	$infostore{$sb}{"L"}++;
	print "For other: record 1 loss at $sb\n"; #so for other descendent, we will record loss
      }
    }
  }
  else{ # we should map setA to setB 
    # %match should be empty now
    foreach my $sb (keys %setB){
      my $sbi = 1;
      foreach my $sa (keys %setA){
	if (isMom($sb,$sa)){
	  &losstrace($sb,$sa);
	  $sbi = 2;
	  last;
	}
	elsif ($sb eq $sa){
	  $sbi =2 ;
	}
      }
      unless ($sbi == 2){
	$infostore{$sb}{"L"}++;
	print "record 1 loss at $sb\n";
      }
    }
  }
}
sub isMom{
  my $first = shift;
  my $second = shift;

  unless ($first and  $second){
    return 0;
  }
  unless (%momchild){
    foreach my $node ($stree->get_nodes){
      my $id = $node->id;
      next unless ($id);
      next if ($node->is_Leaf);
      foreach my $child ($node->get_all_Descendents){
	my $childid = $child->id;
	$momchild{$id}{$childid} =1;
      }
    }
  }
  return $momchild{$first}{$second};
}

sub isMomold{
  my $first = shift;
  my $second = shift;

  my $fnode = $stree->find_node(-id=>$first);
  my $snode = $stree->find_node(-id=>$second);

  unless ($fnode){ die "cannt find node for $first";}
  unless ($snode){die "cannt find node for $second";}
  my $n =0;
  while($n<40){
    my $mom= $snode->ancestor ||last;
    if ($mom->id eq $second){
      return 1;
    }
    $snode = $mom;
    $n++;
  }
  return 0;
}
sub losstrace{
  my $sb = shift; # the node in species tree
  my $sa = shift;
  if ($sb eq $sa){
    return;
  }
  my $node = $stree->find_node(-id=>$sa);
  my $nloop;
  while(1){
    my $mom = $node->ancestor;
    foreach my $child ($mom->each_Descendent){
      my $cid = $child->id;
      unless ($cid eq $node->id){
	$infostore{$cid}{"L"}++;
	print "losstrace record loss at $cid\n";
      }
    }
    if ($mom->id eq $sb){
      last;
    }
    else{
      $node = $mom;
      $nloop++;
      if ($nloop > 20){
	die "infinite loop in losstrace! using $sb and $sa \n";
      }
    }
  }
}
1;
