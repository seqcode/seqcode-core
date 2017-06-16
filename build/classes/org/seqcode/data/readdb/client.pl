#!/usr/bin/perl 
#BEGIN {push @INC, '/afs/bx.psu.edu/home/giardine/share/perl/5.20.2/'}


use strict;
use warnings;
use DBI;
#use lib "./";
use Bio::DB::ReadDBClient;

my $client = Bio::DB::ReadDBClient->new();
my $cmd = shift(@ARGV);
my %props = read_mysqlconfig_file();
my ($username,$passwd,$conn) = @props{qw(user passwd dbiconnectstring)};
my $dbh = DBI->connect($conn, $username, $passwd, {RaiseError=>1, AutoCommit=>0});


if ($cmd eq 'getChroms') {
  my $align = shift(@ARGV);
  my @c = $client->getChroms($align);
  print "chroms for $align are @c\n";
} elsif ($cmd eq 'getCount') {
  my $align = shift(@ARGV);
  my $c = $client->getCount($align);
  print "count for $align is $c\n";
} elsif ($cmd eq 'getHits') {
  my ($align, $chrom, $start, $stop) = @ARGV;
  my $hits;
  if ($start and $stop) {
    $hits = $client->getHitsRange($align,$chrom,$start,$stop);
  } else {
    $hits = $client->getHits($align,$chrom);
  }
  foreach (@$hits) {
    print "$_\n";
  }
} elsif ($cmd eq 'getWeights') {
  my ($align, $chrom, $start, $stop) = @ARGV;
  my $hits;
  if ($start and $stop) {
    #my ($self,$alignid,$chrom,$isType2,$isPaired,$start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
    #$hits = $client->getWeightRange($align,$chrom,$start,$stop);
    $hits = $client->getWeightRange($align, $chrom, 0, 0, $start, $stop, 0, 0, 1);
  } else {
    $hits = $client->getWeight($align,$chrom);
  }
    #foreach (@$hits) {
    #print "$_\n";
    #}
    print $hits;
} elsif ($cmd eq 'getWeightHistogram') {
  #my ($self,$alignid,$chrom,$isType2,$isPaired,$extension,$binsize,$start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
  my ($align, $build, $chrom, $start, $stop) = @ARGV;
  my $binsize;
  if ($ARGV[5]) { $binsize = $ARGV[5]; }
  else { $binsize = 1; }
  #no chr on chrom
  $chrom =~ s/chr//;
  my $sth = $dbh->prepare("select chromosome.id from chromosome, genome where genome.id = chromosome.genome and genome.version = ? and chromosome.name = ?");
  $sth->execute($build, $chrom);
  my $chrid = $sth->fetchrow_array;
  my $hits = $client->getWeightHistogram($align, $chrid, 0, 0, 0, $binsize, $start, $stop, 0, 0, 1);
  #TODO get offset to be added as inpud
  #original value 0 or 1 based???
  foreach (sort { $a <=> $b } keys %$hits) {
    if ($binsize == 1) {
    print "chr$chrom\t$_\t", ($_+1), "\t$hits->{$_}\n";
    }else {
        my $off = int($binsize/2);
        print "chr$chrom\t", ($_-$off), "\t", ($_+$off), "\t$hits->{$_}\n";
    }
  }
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

