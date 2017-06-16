#!/usr/bin/perl
#script to start a trackDb file from the mysql database for readdb
#NOTE: need to drop readdb* tables and reload through bbi when making changes

use strict;
use warnings;
use DBI;

my $tcnt = 0;
my $hub = 1; #flag to set whether script is for hub or resident tracks
my %props = read_mysqlconfig_file();
my ($username,$passwd,$conn) = @props{qw(user passwd dbiconnectstring)};
my $dbh = DBI->connect($conn, $username, $passwd, {RaiseError=>1, AutoCommit=>0});
my $build = shift @ARGV;

my $species = $dbh->selectcol_arrayref("select species from genome where id = ?", undef, $build);
print STDERR "Doing species $species->[0]\n";
#get values to go with foreign keys
my %cell = %{$dbh->selectall_hashref("select id, name from cellline", "id")};
my $cellSth = $dbh->prepare("select id, name from cellline where id in (select cellline from seqdata.seqexpt where expttype = ? and species = ?)");
my %target = %{$dbh->selectall_hashref("select id, name from expttarget", "id")};
my $targetSth = $dbh->prepare("select id, name from expttarget where id in (select expttarget from seqdata.seqexpt where expttype = ? and cellline = ? and species = ?)");
my %type = %{$dbh->selectall_hashref("select id, name from expttype", "id")};
my %lab = %{$dbh->selectall_hashref("select id, name from lab", "id")};
my $labSth = $dbh->prepare("select id, name from lab where id in (select lab from seqdata.seqexpt where cellline = ? and expttype = ? and species = ?)");
my %cond = %{$dbh->selectall_hashref("select id, name from exptcondition", "id")};
my %enclab;
my @bbi; #resident tracks need a table loaded
my $ucscBuild = $dbh->selectrow_arrayref("select version from genome where id = ?", undef, $build);

my $cell_cnt = scalar(keys %cell);
my $lab_cnt = scalar(keys %lab);
print STDERR "TESTING found $cell_cnt cells, and $lab_cnt labs\n";
#520 cells, 99 labs
#super track with cells as subTracks(composite)

#print super track foreach expttype
my $trCnt = 0;
my %super;  #only print super tracks that have wiggles
my %comp;   # "                                     "
foreach my $t (keys %type) {
   #Hubs can't handle parent tracks without children
   $super{"readdbType$t"} =
          "track readdbType$t\n" .
          "shortLabel ReadDb $type{$t}->{name}\n" .
          "longLabel ReadDb $type{$t}->{name}\n" .
          "superTrack on\n\n";
          #"group x\n\n"; #don't need group in hub
          #can have a color and an html page

   #foreach cell in this type do a composite track
   $cellSth->execute($type{$t}->{id}, $species->[0]);
   my $cellCnt = 0;
   while (my @cell = $cellSth->fetchrow_array) {
       #only print if used
       $comp{"readdbType${t}c$cell[0]"} =
             "track readdbType${t}c$cell[0]\n" .
             "shortLabel ReadDb $cell[1]\n" .
             "longLabel ReadDb $type{$t}->{name} $cell[1]\n" .
             "compositeTrack on\n" .
             "parent readdbType$t\n" .
             "type readdb\n";
       #print "group readdb\n";
       #this can have a color and an html page
       my $subgrp = 1;
       if ($type{$t}->{name} eq 'CHIPEXO' or $type{$t}->{name} eq 'RNASEQ') {
          $comp{"readdbType${t}c$cell[0]"} .=  "subGroup$subgrp view Views a=Plus_strand b=Minus_strand\n";
          $subgrp++;
       }
       $comp{"readdbType${t}c$cell[0]"} .=  "subGroup$subgrp lab Lab ";
       $subgrp++;
       $labSth->execute($cell[0], $t, $species->[0]);
       while (my @lab = $labSth->fetchrow_array) {
          $comp{"readdbType${t}c$cell[0]"} .=  " L$lab[0]=$lab[1]";
       }
       $comp{"readdbType${t}c$cell[0]"} .=  "\n";
       if ($type{$t}->{name} eq 'CHIPSEQ') {
          $comp{"readdbType${t}c$cell[0]"} .=  "subGroup$subgrp target Target ";
          $subgrp++;
          $targetSth->execute($t, $cell[0], $species->[0]);
          while (my @tar = $targetSth->fetchrow_array) {
             $tar[1] =~ s/ /_/g; 
             $comp{"readdbType${t}c$cell[0]"} .=  " t$tar[0]=$tar[1]"; 
          }
          $comp{"readdbType${t}c$cell[0]"} .=  "\n";
          $comp{"readdbType${t}c$cell[0]"} .=  "dimensions dimensionY=target dimensionX=lab\n";
       }

       #print "type readdb\n";
       $comp{"readdbType${t}c$cell[0]"} .=  "autoScale on\n";
       $comp{"readdbType${t}c$cell[0]"} .=  "\n";
       #print view subGroup if needed
       if ($type{$t}->{name} eq 'CHIPEXO' or $type{$t}->{name} eq 'RNASEQ') {
          $comp{"readdbType${t}c$cell[0]"} .=  
                "track readdbType${t}c$cell[0]a\n" .
                "shortLabel plus\n" .
                "longLabel plus strand\n" .
                "parent readdbType${t}c$cell[0]\n" .
                "visibility full\n" .
                "view a\n\n" .
                "track readdbType${t}c$cell[0]b\n" .
                "shortLabel minus\n" .
                "longLabel minus strand\n" .
                "visibility full\n" .
                "parent readdbType${t}c$cell[0]\n" .
                "view b\n\n";
       }
       $cellCnt++;
       $trCnt++;
   }
   print STDERR "Found $cellCnt cell tracks for $type{$t}{name}\n";
}
print STDERR "Found $trCnt tracks total\n";

print "#start of wiggle subtracks\n\n";
#fetch public alignments
$dbh->do("use seqdata"); #switch database
#use permissions not publicdbid for determining public
my $sth = $dbh->prepare("select e.name, e.expttype, e.lab, e.exptcondition, e.expttarget, e.cellline, e.publicdbid, a.id from seqalignment a, seqexpt e where a.expt = e.id and a.genome = ? and a.permissions like '\%public\%'");
$sth->execute($build);
my %used;
#TODO should we order the query as well?
my $priority = 0;  #needed to keep plus and minus strand tracks together
while (my @row = $sth->fetchrow_array) {
   $priority++;
   #print super and composite track if not done
   if (defined $super{"readdbType$row[1]"}) {
      print $super{"readdbType$row[1]"};
      undef $super{"readdbType$row[1]"}
   }
   if (defined $comp{"readdbType$row[1]c$row[5]"}) {
      print $comp{"readdbType$row[1]c$row[5]"};
      undef $comp{"readdbType$row[1]c$row[5]"};
   }
   #for exo tracks do plus and minus strand separately?
   my $bbiName = $row[7];
   if ($type{$row[1]}->{name} eq 'CHIPEXO' or $type{$row[1]}->{name} eq 'RNASEQ') {
      print "track readdb$row[7]a\n"; #a for plus strand so goes above minus
      print "parent readdbType$row[1]c$row[5]a\n"; #the parent is the view
      $bbiName .= ';isPlusStrand=1';
   }else {
      print "track readdb$row[7]\n"; #sequence alignment id
      print "parent readdbType$row[1]c$row[5]\n"; #the parent is the cell track
   }
   print "shortLabel $cell{$row[5]}->{name}";
   #if ($row[3]) { print " $cond{$row[3]}->{name}"; } #not helpful in most cases
   if ($row[4] && $target{$row[4]}->{name}) { print " $target{$row[4]}->{name}\n"; }
   else { print "\n"; }
   print "longLabel $row[0]\n";
   print "type readdb\n"; #required in all tracks for hubs
   $used{"readdbType$row[1]c$row[5]"}++;
   print "subGroups lab=L$row[2]";
   if ($type{$row[1]}->{name} eq 'CHIPSEQ') {
      print " target=t$row[4]";
   }elsif ($type{$row[1]}->{name} eq 'CHIPEXO' or $type{$row[1]}->{name} eq 'RNASEQ') {
      print " view=a";
   }
   print "\npriority $priority";
   #if this is for a hub need the bigDataUrl instead of bbi table
   if ($hub) { print "\nbigDataUrl $bbiName\n"; }
   print "\n\n";
   if ($type{$row[1]}->{name} eq 'CHIPEXO' or $type{$row[1]}->{name} eq 'RNASEQ') {
      $priority++;
      print "track readdb$row[7]b\n";
      print "shortLabel $cell{$row[5]}->{name}"; #ideally limited to 17 chars
      #if ($row[3]) { print " $cond{$row[3]}->{name}"; } #not helpful in most cases
      if ($row[4] && $target{$row[4]}->{name}) { print " $target{$row[4]}->{name} minus\n"; }
      else { print "minus\n"; }
      print "longLabel $row[0] minus strand\n";
      print "parent readdbType$row[1]c$row[5]b\n"; #the parent is the view
      $used{"readdbType$row[1]c$row[5]"}++;
      print "subGroups lab=L$row[2] view=b";
      print "\npriority $priority";
      print "\ntype readdb";
      if ($hub) { print "\nbigDataUrl $row[7];isPlusStrand=0"; }
      print "\n\n";
      if (!$hub) {
         push(@bbi, "\thgBbiDbLink $ucscBuild->[0] readdb$row[7]b $row[7];isPlusStrand=0\n");
         push(@bbi, "\thgBbiDbLink $ucscBuild->[0] readdb$row[7]a $row[7];isPlusStrand=1\n");
      }
   }else {
      if (!$hub) { push(@bbi, "\thgBbiDbLink $ucscBuild->[0] readdb$row[7] $row[7]\n"); }
   }
   $tcnt++;
}
print STDERR "have subtracks for ", scalar(keys %used), " tracks\n";


if (@bbi) {
   open(FH, ">", "bbiMakefile") or die "Couldn't open bbi file, $!\n";
   print FH "go:\n";
   foreach (@bbi) { print FH $_; }
   close FH or die "Couldn't finish bbi file, $!\n";
}
print STDERR "Found $tcnt wiggle tracks\n";
$dbh->disconnect;
exit(0);


#**************************************************
sub read_mysqlconfig_file {
  my ($self) = @_;
  my $homedir = './';
  my $basename = 'core_passwd';
  if ($ENV{'READDBROLE'}) {
    $basename = $ENV{'READDBROLE'} . $basename;
  }
  open(PROPS,"${homedir}/.${basename}") or die "can't open config file .$basename in $homedir : $!";
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

