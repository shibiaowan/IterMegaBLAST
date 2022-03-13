#!/usr/bin/perl 
#THis is to generate clusters where the sizes of all clusters are 2 (the last cluster may be 2 or three).
#The input should be the sequence files to be clustered. The output will be the ACs that have been clustered.

####################################################################
###You may need to adjust the E, G and W
####################################################################
# Usage:
#  ./get_size2_cluster.pl.pl <fasta_dir> <tmp dir> <results file>
# Example:
#./get_size2_cluster.pl ~/dnafasta_seq1 info_dnaalign1 dnaalign1_results.txt

####################################################################
#You may need to adjust the -W to adapt to different DNA datasets###
####################################################################

####################################################################
#For DNA dataset1: 
# The average distance is 10.6785714285714
# Wed Aug 17 17:57:11 2016  --  Wed Aug 17 17:57:15 2016
# Time in total: 4 seconds
#-E 5 -G 500 -W 200
####################################################################

####################################################################
# For DNA datasetII
# The average distance is 4.50537634408602
# Wed Aug 17 18:00:57 2016  --  Wed Aug 17 18:01:27 2016
# Time in total: 30 seconds
# for < 150, use -E 5 -G 500 -W 80
# for 
####################################################################

#by Shibiao WAN, 15-Aug-2016, Princeton, NJ, USA
#revised on 08-Mar-2022

use strict;
use warnings;
use Switch;
use Carp;
use File::Basename;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use IPC::System::Simple qw(system capture);
# use Math::NumberCruncher;

my $start = localtime();
my $t1 = time();
print "$start \n";

#################################################################
#the global variable
our %ncode;#code corresponding to different types of nucleotides
$ncode{'A'} = 'A';$ncode{'G'} = 'G';$ncode{'T'} = 'T';$ncode{'C'} = 'C';#Level 1 
$ncode{'R'} = 'AG';$ncode{'Y'} = 'CT';$ncode{'S'} = 'CG';$ncode{'W'} = 'AT';$ncode{'K'} = 'GT';$ncode{'M'} = 'AC';#Level 2
$ncode{'B'} = 'CGT';$ncode{'D'} = 'AGT';$ncode{'H'} = 'ACT';$ncode{'V'} = 'ACG';#Level 3
$ncode{'N'} = 'ACGT';#Level 4
#print Dumper(\%ncode);
# print "@{$ncode{'B'}}\n";#for array, we must use @
    
our %revncode = reverse %ncode;#the reverse ncodes
#print Dumper(\%revncode);

#################################################################


($#ARGV == 2) || die "Usage: $0 <fasta_dir> <tmp dir> <results file>\n"; 

# Find all .fasta files in fasta dir
my $dir = $ARGV[0];
opendir DIR, $dir or die "error: cannot open directory \"$dir\" $!";
my @seqs = sort grep (!/^\.$|^\.\.$/, readdir(DIR));#@seqs include the whole path
closedir (DIR);

@seqs = &sort_file_numerically(\@seqs);#sort numerically

# # Read .fasta file process sequences one by one
# my @seqs = read_all_sequences($fastafile,'fasta');

my $tmpdir = $ARGV[1];

my $resultsfile = $ARGV[2];
unless(open (RESULTSFILE, ">$resultsfile")){croak "cannot open $resultsfile\n";}

################################
#Saving final results
my @misnum;#number of mismatches between the pair sequences
my @tdist;#the distances of all the pair sequences
my @genpseq;#the generalization/parent sequences of all the pair sequences
my @misind;#the mismatch positions of all the pair sequences
my $avgdist = 0;#the average distance

my $num = @seqs;#number of sequences
# print "All seqs: @seqs\n";
my $t = $num;

#################################################
#This is to randomly select n sequences from @seq
# my @items;
# my $n = $num;
# my $n = 174;
#   for ( 1 .. $n )
#   {
#     push @items, splice @seqs, rand @seqs, 1;
#   }
# @seqs = @items;#already randomly sorted
#################################################

#################################################
#According to the paper, we selected subsets of the sequences in ascending order of Genbank IDs.
#################################################
my $n = $num;
my $n = 200;
# my $c = int(rand($num-$n));
my $c = 0;
@seqs = @seqs[$c .. ($n + $c-1)];#already randomly sorted


my @discardseq = ();
my @newdbseq;
my $q;#the query seq
my $j = 0;
my %idmap;
my %filemap;
while($t > 1){
    $j++;
    print "------------------------------------------------------------------------\n";
    print "The $j-th round sequence alignment:\n";
    if ($t == $num){
        $q = $seqs[0];#use the first sequence as the query sequence
        }else{
            $q = $newdbseq[0];
        }
   
    push @discardseq, $q;##############make a serious difference between . and ,
#      print "Discarded sequence: @discardseq\n";
        
    @newdbseq = &get_setdiff(\@seqs, \@discardseq);##################note when passing arrays, using /@seqs
#     print "The sequences for the database are: @newdbseq\n";
    my @newdb;
    foreach my $item(@newdbseq){
        my $seqn = "$dir/$item";#with the full path
        push @newdb, $seqn;#the new database for search
        
        my $in = Bio::SeqIO->new(-file => "$dir/$item", '-format' => 'Fasta');
        my $seq = $in->next_seq();
        $idmap{$item} = $seq->display_id;#create a map between the filename and the sequence ID
#         print "The seq id is $idmap{$item}\n";
        }
    %filemap = reverse %idmap;    
    
#     print "@newdb\n";
    my @cmdt;
    my $newdbname = "newdb-$j.fasta";#the new database for the $j-th time
    my $alignfile = "alignresults-$j.fasta";
    unless (-d $tmpdir){
        $cmdt[0] = "mkdir $tmpdir; cat @newdb > $tmpdir/$newdbname";
        }else{
            $cmdt[0] = "cat @newdb > $tmpdir/$newdbname";#rm -f $tmpdir/*;
        }
#     print "$cmdt[0]\n";
    system($cmdt[0]);
    $cmdt[1] = "cd $tmpdir; formatdb -i $newdbname -p F -o F";
    print "$cmdt[1]\n";
    system($cmdt[1]);
    
    my $initW = 220;
    $cmdt[2] = "nohup megablast -d $tmpdir/$newdbname  -i $dir/$q -E 5 -G 5 -W $initW > $tmpdir/$alignfile";
    print "$cmdt[2]\n";
    system($cmdt[2]);
    # Run a command, wait until it finishes, and make sure it works.
    # The output of this command is captured into $results.
#     my results = capture($^X, "get_dist_from_align.pl", "$tmpdir/$alignfile");
    my ($mislen, $newind, $totaldist, $parentseq, $queryid, $matchedid) = &get_dist_from_align("$tmpdir/$alignfile");
    
    my $r = 0;
    while ($matchedid eq "none"){#if cannot find the matched seqs, we use a lower criteria
        $r++;
        my $Wvalue = $initW - 20*$r; 
        $cmdt[2] = "nohup megablast -d $tmpdir/$newdbname  -i $dir/$q -E 5 -G 500 -W $Wvalue > $tmpdir/$alignfile";
        print "$cmdt[2]\n";
        system($cmdt[2]);
        
        ($mislen, $newind, $totaldist, $parentseq, $queryid, $matchedid) = &get_dist_from_align("$tmpdir/$alignfile");
        }
    
    my @newind = @$newind;#first deal with the returned array
    
    if (exists $filemap{$matchedid}){#the matched sequence
        push @discardseq, $filemap{$matchedid};#the matched sequence should be ignored in the next round
        @newdbseq = &get_setdiff(\@seqs, \@discardseq);##################this time, we should also remove the matched sequence
        }
    push @misnum, $mislen;
    push @tdist, $totaldist;
    $avgdist = $avgdist + $totaldist;
    push @genpseq, $parentseq;
    push @misind, @newind;
    print RESULTSFILE "------------------------------------------------------------------------\n";
    print RESULTSFILE "The $j-th pair sequence:\n\n";
    print RESULTSFILE "Query ID: $queryid; Matched ID: $matchedid\n";
    print RESULTSFILE "Mismatch number: $mislen\n";
    print RESULTSFILE "The distance: $totaldist\n";
    print RESULTSFILE "Mismatch positions: @newind\n";
    print RESULTSFILE "Generalized sequence:\n$parentseq\n";
    
#     print "The mismatch positions are: (in total $mislen mismatches)\n ";
#     print join(", ", @newind);
#     print "\n";
# 
#     print "The total distance between the parent seq and the two DNA sequences is: $totaldist\n";
#     print "$results\n";
    
#     $cmdt[3] = "rm -f $tmpdir/$newdbname.*";
#     print "$cmdt[3]\n";
#     system($cmdt[3]);

#     my $cmd = join(";", @cmdt);
#     print join("\n", @cmdt);
#     system($cmd);
    #my $cmd1 = "find . -maxdepth 1 -iname '*.fasta' -not -name '$q.txt' -exec cat {} +>database.txt";
    $t = @newdb;
    %filemap = ();
    %idmap = ();
}

my $sumdist = $avgdist;
$avgdist = $avgdist/(scalar @tdist);
print "The average distance is $avgdist (distance sum: $sumdist)\n";
my $finish = localtime();
my $t2 = time();
my $timedif = $t2 - $t1;
# my $costtime = $finish - $start;
print "$start  --  $finish\nTime in total: $timedif seconds\n";

print RESULTSFILE "------------------------------------------------------------------------\n";
print RESULTSFILE "The average distance: $avgdist\n";
print RESULTSFILE "Time cost: $timedif\n";
print RESULTSFILE "All distances: ($j distances in total)\n@tdist\n";

close RESULTSFILE or die "Cannot close $resultsfile\n";


###############################################################################################################################
########The following are subrountines
sub sort_file_numerically{
    my $aref = $_[0];#note the difference between $_[0] and @_
    
    my @a = @{$aref};#passing the arrays
    
    my %b;
    foreach my $item (@a) {#get the number-file hash
        if ($item =~ m/(\d+)/g){
            my $t = $1;
            $b{$t} = $item;
            }
        }
    
    my @sorta;
    foreach my $name (sort { $a <=> $b } keys %b) {#sort the file numerically
        push @sorta, $b{$name};
        }
    return @sorta;
}
    
sub get_setdiff{
    #This is to get the set difference between @a and @b, i.e., @a/@b 
    my ($aref, $bref) = @_;
    my @a = @{$aref};#passing the arrays
    my @b = @{$bref};
#     print "Array a: @a\n";
#     print "Array b: @b\n";
    my %b=map{$_=>1} @b;#the smaller subset
    my @restset = grep(!defined $b{$_}, @a);
    return @restset;
    }

################################################################################################################################
#######The following codes are derived from "get_dist_from_align.pl"###########
################################################################################################################################
sub get_dist_from_align{
#This is to calculate the distance between the query and subject sequences, meanwhile generate the generalization sequence.

# Define constant and input parameters to psi-blast
my $alignfile = $_[0];#the input file

my $totaldist = 0;#the distance between the generalized sequence and the two DNA sequences;initial value is 0
my $parentseq = ();#the generalized sequence


unless(open (ALIGNFILE, "$alignfile")){croak "Cannot open $alignfile\n";}

my $info = join("",<ALIGNFILE>);#all the contents

my $queryid = $1 if ($info =~ m/Query= ([\w\d\|]+)/); #the query sequence
# print "THe query ID is: $queryid\n";
my $matchedid = $1 if ($info =~ m/>([\w\d\|]+)/);#the matched sequence
# print "The matched ID is: $matchedid\n";

my $aligncontent = (split(/>/,$info))[1];#save the first sequence-alignment match (the second item of the content) 
# print "$aligncontent";

my @tmp;
if (defined $aligncontent){
  @tmp = split(/Query:/,$aligncontent);#the alignment sequence
}else{
  close ALIGNFILE or die "Cannot close $alignfile\n";

  my $matchedid = "none";
  my $mislen = "";
  my @newind = ();
  my $totaldist = ();
  my $parentseq = "";
  my $queryid = "";
  return ($mislen, \@newind, $totaldist, $parentseq, $queryid, $matchedid);
    }
my $l = @tmp;
chomp @tmp;

my @newind = ();
my @allgenseq;
for(my $i = 1; $i<$l; $i++){#starting from the second item (starting from the query sequence)
    my @pos =($tmp[$i]=~ m/(\d+)/g);#nucleotide position;note that perhaps there are more than 4 numbers, please ignore the rest numbers
#     print join(", ", @pos);
#     print "\n";
    my @tmp1 = split(/\n/,$tmp[$i]);#all info for a piece of sequence-segment alignment
#     print "Segment $i:\n";
#     print "$tmp1[0]\n";
#     print join(", ", @tmp1);
#     print "\n";    
   
    my $qs = $1 if ($tmp1[0] =~ m/\d+\s+([\w-]+)\s+\d+/);#the query sequence; note the difference between the matched positions and contents
#     print "$qs\n";
    my $len = length($qs);
    my $genseq = $qs;#the generalized sequence is based on the query sequence
    
   # print "\nThe length of the query segment is $len\n";
    my $ss = $1 if ($tmp1[2] =~ m/\d+\s+([\w-]+)\s+\d+/);#the subject sequence#########Note that gap (-) should not be forgottern
#     print "$ss\n";
    my $flag = substr($tmp1[1], -$len);#the match flag, only save the last $len characters
    #$flag =~ s/\n\w+://;#special word "Sbjct:"
    my $flaglen = length($flag);
    #print "The flag length is $flaglen\n";#the length is consistent with the alignment length, including whitespace for mismatch
    
    my @ind = ();
    while ($flag =~ m/\s/g){
        push @ind, $-[0];#find the starting indices of all mismatches
    }
    
#     print "Segment $i\n";
#     print join(", ", @ind);
#     print "\n";
    if (scalar @ind >0){#mismatch exists
        my $x = scalar @ind;
#         print "Mismatch number is: $x and are @ind\n";
      foreach my $item(@ind){
          my $qmis = substr($qs, $item, 1);#mismatch in the query sequence
          my $smis = substr($ss, $item, 1);#mismatch in the subject sequence
#         print "mismatch position at $item of the query is $qmis\n";
          my ($dist, $parent) = &get_dist_from_NNpair($qmis, $smis);
#         print "The distance is $dist\n";
          $totaldist = $totaldist + $dist; 
        
          substr($genseq, $item, 1) = $parent;
#           print "Pos: $item\n";
        if ($qmis ne 'n'){#########################################Note the new change here########################################
        my $k = $item + $pos[0];#the position in the whole DNA sequence
        push @newind, $k;
         
        }
      }
    
    #$parentseq = join('\n', $parentseq, $genseq);
    }
    push @allgenseq, $genseq; 
}

$parentseq = join("\n", values @allgenseq);

my $mislen = @newind;
# print "The mismatch positions are: (in total $mislen mismatches)\n ";
# print join(", ", @newind);
# print "\n";
# 
# print "The total distance between the parent seq and the two DNA sequences is: $totaldist\n";
# print "The generalized/parent sequence is:\n$parentseq\n";
            
close ALIGNFILE or die "Cannot close $alignfile\n";

return ($mislen, \@newind, $totaldist, $parentseq, $queryid, $matchedid);######################################returning an array is different from scalars
}

sub get_dist_from_NNpair {#nucleotide-nucleotide pairs
    my $n1 = $_[0];#the first N; the query one
    my $n2 =$_[1];#the second N; the matched one
    
    #convert to uppercases
    $n1 = uc $n1;
    $n2 = uc $n2;

    my $level1;
    my $level2;
    if ($n1 eq 'N'){#it is probably masked by the Megablast program
        my $parent = $n2; 
        my $dist = 0;
    
        return ($dist, $parent);
   }else{
            if (exists $ncode{$n1}){
       $level1 = length($ncode{$n1});#global varialbe for the levels of codes
    }else{
        $level1 = 3; #the gap
    }
    
    if (exists $ncode{$n2}){
       $level2 = length($ncode{$n2});
    }else{
        $level2 = 3;# the gap 
    }
    
    my ($parent, $parentlevel) = &get_parent_from_NNpair($n1, $n2);
#     print "The three distances are: $parentlevel, $level1, $level2\n";
    
    my $dist = 2*$parentlevel - $level1 - $level2;
    
    return ($dist, $parent);
}
}

sub get_parent_from_NNpair{#for generalization
    my $n1 = $_[0];
my $n2 = $_[1];
my $parent = ();
my $level = ();#the level of the code
    
#convert to uppercases
$n1 = uc $n1;
$n2 = uc $n2;
    
my $p = ();
    
if ((exists $ncode{$n1}) && (exists $ncode{$n2})){
#     push @p, @{$ncode{$n1}};#use @ to represent arrays
#     push @p, @{$ncode{$n2}};
    $p = join('',$ncode{$n1},$ncode{$n2});
#     print "$p\n";
    my $uniq = &get_uniq_array($p);#unique elements
#     print "The unique elements are: @uniq\n";
#     print "The unique sorted elements are: $uniq\n";
    if (exists $revncode{$uniq}){
        $parent = $revncode{$uniq};
        }
}elsif (($n1 eq '-') || ($n2 eq '-')){
    $parent = "N";#if any one of these two nucleotides is a gap, then their parent should be 'N''
    }
       
$level = length($ncode{$parent});       
# print "The parent code of $n1 and $n2 is $parent\nThe level the parent code is $level\n";
return ($parent, $level);
}
             
sub get_uniq_array{
    #This is to get the sorted unique elements of an array and return a scalar consisting of those sorted unique characters
    my $p = $_[0];#the original array
    my @tmp = split(//,$p);#characters split
    my %seen = ();
    my @uniq;
    my $uniqsort;
    foreach my $item(@tmp){
        push(@uniq, $item) unless $seen{$item}++;
    }
    %seen = ();#clear
    @uniq = sort @uniq;
    $uniqsort = join('',@uniq);
    
    return $uniqsort;
}
