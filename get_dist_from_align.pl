#!/usr/bin/perl 
# This is to calculate the distance between the query and subject sequences, meanwhile generate the generalization sequence.
# Usage:
#  ./get_dist_from_align.pl <ncalign_file> <dist>
# Example:
#  ./get_dist_from_align.pl 
#  ./get_dist_from_align.pl 

#by Shibiao WAN, 10-Aug-2016, Princeton, NJ, USA
#revised on 08-Mar-2022

use strict;
use warnings;
use Switch;
use Carp;
use Data::Dumper;
use File::Basename;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

($#ARGV == 0) || die "Usage: $0 <ncalign_file>\n"; 

# Define constant and input parameters to psi-blast
my $alignfile = $ARGV[0];
my $totaldist = 0;#the distance between the generalized sequence and the two DNA sequences;initial value is 0
my $parentseq = ();#the generalized sequence


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


unless(open (ALIGNFILE, "$alignfile")){croak "Cannot open $alignfile\n";}

my $info = join("",<ALIGNFILE>);#all the contents

my $aligncontent = (split(/>/,$info))[1];#save the first sequence-alignment match (the second item of the content) 
# print "$aligncontent";

my @tmp = split(/Query:/,$aligncontent);#the alignment sequence
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
#           print "Pos: $item\n";
        my $k = $item + $pos[0];#the position in the whole DNA sequence
        push @newind, $k;
         
        my $qmis = substr($qs, $item, 1);#mismatch in the query sequence
        my $smis = substr($ss, $item, 1);#mismatch in the subject sequence
#         print "mismatch position at $item of the query is $qmis\n";
        my ($dist, $parent) = &get_dist_from_NNpair($qmis, $smis);
#         print "The distance is $dist\n";
        $totaldist = $totaldist + $dist; 
        
        substr($genseq, $item, 1) = $parent; 
      }
    
    #$parentseq = join('\n', $parentseq, $genseq);
    }
    push @allgenseq, $genseq; 
}

$parentseq = join("\n", values @allgenseq);

my $mislen = @newind;
print "The mismatch positions are: (in total $mislen mismatches)\n ";
print join(", ", @newind);
print "\n";

print "The total distance between the parent seq and the two DNA sequences is: $totaldist\n";
print "The generalized/parent sequence is:\n$parentseq\n";
            
close ALIGNFILE or die "Cannot close $alignfile\n";

# return ($mislen, @newind, $totaldist, $parentseq);

sub get_dist_from_NNpair {#nucleotide-nucleotide pairs
    my $n1 = $_[0];#the first N
    my $n2 =$_[1];#the second N
    
    #convert to uppercases
    $n1 = uc $n1;
    $n2 = uc $n2;

    my $level1;
    my $level2;
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
}else{
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
    
