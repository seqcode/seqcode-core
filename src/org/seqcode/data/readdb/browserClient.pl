#!/usr/bin/perl 
#BEGIN {push @INC, '/afs/bx.psu.edu/home/giardine/share/perl/5.20.2/'}


use strict;
use warnings;
use DBI;
use lib "./";
use ReadDBClient;

my $client = new ReadDBClient(); #Bio::DB::ReadDBClient->new();
my $cmd = shift(@ARGV);
my %props = read_mysqlconfig_file();
my ($username,$passwd,$conn) = @props{qw(user passwd dbiconnectstring)};
my $dbh = DBI->connect($conn, $username, $passwd, {RaiseError=>1, AutoCommit=>0});


if ($cmd eq 'getWeightHistogram') {
  #my ($self,$alignid,$chrom,$isType2,$isPaired,$extension,$binsize,$start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
  my ($align, $build, $chrom, $start, $stop) = @ARGV;
  my $binsize;
  if ($ARGV[5]) { $binsize = $ARGV[5]; }
  else { $binsize = 1; }
  if ($stop < $start) { die "Error stop is less than start $start\n"; }
  #no chr on chrom
  $chrom =~ s/chr//;
  #align may have some extra parameters 
  my $param = "$align";  #quote to make a copy
  $align =~ s/;.*//; #remove parameters from id
  my $isPlusStrand = undef;
  my $isLeft = undef;
  my $minWeight = 0;
  if ($param =~ /;/) { 
#print STDERR "TESTING have param $param and align $align\n";
     if ($param =~ /isPlusStrand=(\d)/) { $isPlusStrand=$1; }
     if ($param =~ /isLeft=(\d)/) { $isLeft=$1; } #boolean
     if ($param =~ /minWeight=(\d+)/) { $minWeight=$1; }
  }else { undef $param; }
  my $sth = $dbh->prepare("select chromosome.id from chromosome, genome where genome.id = chromosome.genome and genome.version = ? and chromosome.name = ?");
  $sth->execute($build, $chrom);
  my $chrid = $sth->fetchrow_array;
  my $sth2 = $dbh->prepare("select expttype.name, readlength from seqdata.seqexpt, expttype where expttype.id = seqdata.seqexpt.expttype and seqdata.seqexpt.id = (select expt from seqdata.seqalignment where id = ? limit 1)");
  $sth2->execute($align);
  my($expttype, $readlen) = $sth2->fetchrow_array;
  my $extension = 0;
  if ($expttype && $expttype eq 'CHIPSEQ') { $extension = 100; }
  #TODO use fragment length for RNA-seq
  #print STDERR "TESTING found type=$expttype, len=$readlen\n";
#print STDERR "TESTING sending command getWeightHistogram($align, $chrid, 0, 0, $extension, $binsize, $start, $stop, $minWeight, $isLeft, $isPlusStrand);\n";
  my $hits = $client->getWeightHistogram($align, $chrid, 0, 0, $extension, $binsize, $start, $stop, $minWeight, $isLeft, $isPlusStrand);
  #original value 0 or 1 based??? 1??
  foreach (sort { $a <=> $b } keys %$hits) {
    if ($binsize == 1) {
        if (defined $isPlusStrand && $isPlusStrand == 0) { #minus strand neg val
           print "chr$chrom\t", ($_-1), "\t$_\t-$hits->{$_}\n";
        }else {
	   print "chr$chrom\t", ($_-1), "\t$_\t$hits->{$_}\n";
        }
    }else { #if binned spread interval over full bin
        my $off = int($binsize/2);
        print "chr$chrom\t", ($_-$off), "\t", ($_+$off), "\t$hits->{$_}\n";
    }
  }
}elsif ($cmd eq 'getDetailPage') {
  my $align = shift @ARGV;
  #core tables have text for many of the fields
  my $cellsth = $dbh->prepare("select name from cellline where id = ?");
  my $labsth = $dbh->prepare("select name from lab where id = ?");
  my $typsth = $dbh->prepare("select name from expttype where id = ?");
  my $spsth = $dbh->prepare("select name from species where id = ?");
  my $consth = $dbh->prepare("select name from exptcondition where id = ?");
  my $tarsth = $dbh->prepare("select name from expttarget where id = ?");
  my $gnsth = $dbh->prepare("select version from genome where id = ?");
  my $sth = $dbh->prepare("select expttype, seqdata.seqexpt.name, replicate, species, readlength, numreads, cellline, expttarget from seqdata.seqexpt, seqdata.seqalignment where seqdata.seqalignment.expt = seqdata.seqexpt.id and seqdata.seqalignment.id = ?");
  $align =~ s/;.*//; #may have other parameters for the wiggles
  $sth->execute($align);
  my @row = $sth->fetchrow_array;
  if (!@row) {
     print "No data found for $align";
     return;
  }
  if ($row[0]) {
     $typsth->execute($row[0]);
     my @t = $typsth->fetchrow_array;
     $typsth->finish; #should be only one
     if (@t && $t[0]) { $row[0] = $t[0]; }
  }else { $row[0] = ''; }
  if (!$row[1]) { $row[1] = ''; }
  if (!$row[2]) { $row[2] = ''; }
  if ($row[3]) {
     $spsth->execute($row[3]);
     my @t = $spsth->fetchrow_array;
     $spsth->finish; #should be only one
     if (@t && $t[0]) { $row[3] = $t[0]; }
  }else { $row[3] = ''; }
  if (!$row[4]) { $row[4] = ''; }
  if (!$row[5]) { $row[5] = ''; }
  if ($row[6]) {
     $cellsth->execute($row[6]);
     my @t = $cellsth->fetchrow_array;
     $cellsth->finish; #should be only one
     if (@t && $t[0]) { $row[6] = $t[0]; }
  }else { $row[6] = ''; }
  if ($row[7]) {
     $tarsth->execute($row[7]);
     my @t = $tarsth->fetchrow_array;
     $tarsth->finish; #should be only one
     if (@t && $t[0]) { $row[7] = $t[0]; }
  }else { $row[7] = ''; }
  print "Experiment metadata for $align:<br>\n",
        "Experiment type: $row[0]<br>\n",
        "Experiment name: $row[1]<br>\n",
        "Cell line: $row[6]<br>\n",
        "Target: $row[7]<br>\n",
        "Replicate: $row[2]<br>\n",
        "Species: $row[3]<br>\n",
        "Read length: $row[4]<br>\n",
        "Number of reads: $row[5]<br>\n";
  $sth->finish;
}

$client->close();
$dbh->disconnect();

exit;

sub read_mysqlconfig_file {
  my ($self) = @_;
  my $homedir = './';
  my $basename = 'core_passwd';
  if ($ENV{'READDBROLE'}) {
    $basename = $ENV{'READDBROLE'} . $basename;
  }
    #open(PROPS,"${homedir}/.${basename}") or die "can't open config file .$basename in $homedir : $!";
  open(PROPS,"${basename}") or die "can't open config file $basename in current working directory";
  my %props = ();
  while(<PROPS>) {
    chomp;
    s/^\s*//;
    s/\s*$//;
    my ($k,$v) = ($_ =~ /([^=]*)\s*=\s*(.*)/);
    next unless ($k);
    $props{$k} = $v;
  }
  return %props;
}

