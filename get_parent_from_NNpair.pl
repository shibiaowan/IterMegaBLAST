#!/usr/bin/perl 
# This is to generate the generalization sequence (find the parent code of two nucleotides) and also the level of the parent code.
# Usage:
#  ./get_parent_from_NNpair.pl <ncode1> <ncode2>
# Example:
#  ./get_parent_from_NNpair.pl B T

#by Shibiao WAN, 13-Aug-2016, Princeton, NJ, USA

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
($#ARGV == 1) || die "Usage: $0 <ncode1> <ncode2>\n"; 

my $n1 = $ARGV[0];
my $n2 = $ARGV[1];
my $parent = ();
my $level = ();#the level of the code
    
#convert to uppercases
$n1 = uc $n1;
$n2 = uc $n2;
    
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
    
my $p = ();
    
if ((exists $ncode{$n1}) && (exists $ncode{$n2})){
#     push @p, @{$ncode{$n1}};#use @ to represent arrays
#     push @p, @{$ncode{$n2}};
    $p = join('',$ncode{$n1},$ncode{$n2});
    print "$p\n";
    my $uniq = &get_uniq_array($p);#unique elements
#     print "The unique elements are: @uniq\n";
    print "The unique sorted elements are: $uniq\n";
    if (exists $revncode{$uniq}){
        $parent = $revncode{$uniq};
        }
}else{
    $parent = "N";#if any one of these two nucleotides is a gap, then their parent should be 'N''
    }
       
$level = length($ncode{$parent});       
print "The parent code of $n1 and $n2 is $parent\nThe level the parent code is $level\n";
             
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
    