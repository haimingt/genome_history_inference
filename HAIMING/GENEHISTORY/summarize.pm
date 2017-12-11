#!/usr/bin/perl

package HAIMING::GENEHISTORY::summarize;

use HAIMING::GENEHISTORY::compareOne;
use Data::Dumper;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(runEach historySum toFiles);


sub runEach{
  my $directory = shift;  # so it should be a directory with the sum.txt files
  my $speciestree = shift;
  my $odir = shift;
  my $pthrlimit = shift;

  opendir (DH, $directory) or die "cannot open dir $directory\n";
  my @files = readdir(DH);
  closedir DH;
  my $n = 0;
  foreach my $file (sort @files){
    next unless ($file =~ /(PTHR[0-9]+)/);
    if ($pthrlimit){
      next unless ($file =~ /$pthrlimit/);
    }
    my $pthr= $1;
    last if $n> 10;
#    my $treefile = $directory."/$1/orig/tree/$1.orig.tree";
    my $treefile = $directory."/$pthr.newick";
    my $tab = $directory."/$pthr.tab";
    my $tmp = "$odir/$pthr.sum.txt";
    my $ntab = $tab.".up";
    if (-s $tmp){
      my $indicator =0;
      open TMP, "< $tmp" or die;
      while(<TMP>){
	if ($_ =~ /L:/){
	  $indicator =1;
	  last;
	}
      }
      close TMP;
      if ($indicator ==1){
	next;
      }
    }
    print STDERR "workin on $treefile\n";
    treeCompare($treefile,$tab,$speciestree,$tmp);
    $n++;
  }
}

sub historySum{
  my %infostore;
  my %horizsum;

  my $directory = shift;
  my $outfile = shift;
  my $horizfile = shift;
  my $hd = $horizfile.".de";

  open HD, "> $hd" or die;
  open OUT, "> $outfile" or die;
  open HO, "> $horizfile" or die;
  opendir (DH, $directory) or die;
  my @files = readdir(DH);
  closedir DH;
  foreach my $file (@files){
    next unless ($file =~ /(PTHR[0-9]+.sum.txt)/);
    my $branch;
    my $f = $directory."/".$file;
    open IN, "< $f" or die;
    print "working on $f\n";
    while(<IN>){
      chomp;
      my @a = split(/\t/);
      if ($a[0] =~ /Horiz/){
	print HD "$file:\t$_\n";
	$horizsum{$a[1]}{$a[3]}++;
	$horizsum{$a[1]}{"total"}++;

      }
      else{
	my $size =scalar @a;
	foreach my $i (1..$size-1){
	  my @tmp = split(":",$a[$i]);
	  $infostore{$a[0]}{$tmp[0]} += $tmp[1];
	}
      }
    }
    close IN;
  }

  foreach my $key (keys %horizsum){
    $infostore{$key}{"H"} =$horizsum{$key}{"total"};
  }

  print OUT Dumper(\%infostore);
  close OUT;
  print HO Dumper(\%horizsum);
  close HO;
  close HD;

}

sub toFiles{
  my $infofile = shift;
  my $piechart = shift;
  my $branchlabel =shift;

  my %infostore;
  my $branch;
  open IN, "< $infofile" or die;
  
  while(<IN>){
    if ($_ =~ /\'(.*)\' => \{/){
      $branch = $1;
    }
    elsif ($_ =~ /\'([A-Z])\' => ([0-9]+)/){
      my $type = $1;
      my $n = $2;
      $infostore{$branch}{$type} += $n;
    }
  }
  close IN;

  open PIE, "> $piechart" or die;
  print PIE "DATASET_PIECHART\nSEPARATOR TAB\nDATASET_LABEL\tGene evolutionary history (PANTHER12)\nCOLOR\t#ff0000\n";
  print PIE "FIELD_LABELS\tDup\tLoss\tHoriz\tRoot\nFIELD_COLORS\trgba(255,0,0,0.5)\trgba(0,255,0,0.5)\trgba(0,0,255,0.5)\trgba(255,0,255,0.5)\n";
  
  print PIE "LEGEND_TITLE\tG/L/R/H representations\n";
  print PIE "LEGEND_SHAPES\t1\t1\t1\t1\n";
  print PIE "LEGEND_COLORS\trgba(255,0,0,0.5)\trgba(0,255,0,0.5)\trgba(0,0,255,0.5)\trgba(255,0,255,0.5)\n";
  print PIE "LEGEND_LABELS\tDuplication\tLoss\tHorizontal Transfer\tDe novo creation\n";

  print PIE "MARGIN\t1.50\n";
  print PIE "HEIGHT_FACTOR\t1\nDATA\n";


  open BR, "> $branchlabel" or die;
  print BR "DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\tAncestral branches (PANTHER12)\nCOLOR\t#ff00ff\n";
  print BR "DATA\n";


  my $largest =0 ;
  my $smallest = 10000;

  foreach my $key (keys %infostore){
    next unless ($key);
    my ($G,$L,$H,$R) = (0,0,0,0);

    $G = $infostore{$key}{"G"};
    $L = $infostore{$key}{"L"};
    $H = $infostore{$key}{"H"};
    $R = $infostore{$key}{"R"};


    unless ($G){
      $G =0;
    }
    unless ($L){
      $L =0;
    }
    unless ($H){
      $H =0;
    }
    unless ($R){
      $R =0;
    }

    my $t = $G+$L+$H+R;
    if ($t < $smallest){
      $smallest = $t;
    }
    if ($t > $largest){
      $largest = $t;
    }
#    print "$smallest\t$largest\n";
  }
  foreach my $key (keys %infostore){
    next unless ($key);
    my ($G,$L,$H,$R) = (0,0,0,0);

    $G = $infostore{$key}{"G"};
    $L = $infostore{$key}{"L"};
    $H = $infostore{$key}{"H"};
    $R = $infostore{$key}{"R"};


    unless ($G){
      $G =0;
    }
    unless ($L){
      $L =0;
    }
    unless ($H){
      $H =0;
    }
    unless ($R){
      $R =0;
    }

    my $t = $G+$L+$H+R;
    my $s = &size($t,$smallest,$largest);
    my $l;
    if (($key =~ /[a-z]/) or ($key =~ /LUCA/)){
      $l = 0.8;
      my $length = length($key);
      my $sizefactor =1;
      if ($length > 10){
	$sizefactor = 10/$length;
      }

      print BR "$key\t$key\t0.3\t#000000\tnormal\t$sizefactor\t0\n";
    }
    else{
      $l = -1;
    }
    print PIE "$key\t$l\t$s\t$G\t$L\t$H\t$R\n";

    print "$key\t$t\t$s\n";

  }

  close PIE;  
  close BR;
}

sub size{
  my $t = shift;
  my $smallest = shift;
  my $largest = shift;

  my $s = 10;
  my $l = 200;
  
  my $size = ($t-$smallest)*($l-$s)/($largest-$smallest)+$s;
  return $size;
}
