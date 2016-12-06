#!/usr/bin/perl 
#BEGIN {push @INC, '/afs/bx.psu.edu/home/giardine/share/perl/5.20.2/'}


use strict;
use warnings;
use DBI;
use lib "./";
use ReadDBClient;

my $client = Bio::DB::ReadDBClient->new();
my $cmd = shift(@ARGV);
my %props = read_mysqlconfig_file();
my ($username,$passwd,$conn) = @props{qw(user passwd dbiconnectstring)};
my $dbh = DBI->connect($conn, $username, $passwd, {RaiseError=>1, AutoCommit=>0});


if ($cmd eq 'getWeightHistogram') {
  #my ($self,$alignid,$chrom,$isType2,$isPaired,$extension,$binsize,$start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
  #getHistogram, extension is used here
  #my ($self,$alignid,$chrom,$isType2,$isPaired,$extension,$binsize, $start,$stop,$minweight, $isLeft,$isPlusStrand) = @_;
  my ($align, $build, $chrom, $start, $stop) = @ARGV;
  my $binsize;
  if ($ARGV[5]) { $binsize = $ARGV[5]; }
  else { $binsize = 1; }
  #no chr on chrom
  $chrom =~ s/chr//;
  my $sth = $dbh->prepare("select chromosome.id from chromosome, genome where genome.id = chromosome.genome and genome.version = ? and chromosome.name = ?");
  $sth->execute($build, $chrom);
  my $chrid = $sth->fetchrow_array;
  #TODO set extension based on expttype
  my $sth2 = $dbh->prepare("select expttype.name, readlength from seqdata.seqexpt, expttype where expttype.id = seqdata.seqexpt.expttype and seqdata.seqexpt.id = (select expt from seqalignment where id = ? limit 1)");
  $sth2->execute($align);
  my($expttype, $readlen) = $sth2->fetchrow_array;
  my $extension = 0;
  #TODO do offset based on strand for chipseq? similar to macs
  if ($expttype && $expttype eq 'CHIPSEQ') { $extension = 100; }
  #print STDERR "TESTING found type=$expttype, len=$readlen\n";
  my $hits = $client->getWeightHistogram($align, $chrid, 0, 0, $extension, $binsize, $start, $stop, 0, 0, 1);
  #original value 0 or 1 based??? 1??
  foreach (sort { $a <=> $b } keys %$hits) {
    if ($binsize == 1) {
        print "chr$chrom\t", ($_-1), "\t$_\t$hits->{$_}\n";
    }else { #if binned spread interval over full bin
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

