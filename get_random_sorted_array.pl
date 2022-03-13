#!/usr/bin/perl 
#This is to get the random order of an array of a folder

#by Shibiao WAN, 03-Sep-2016, Princeton, NJ, USA
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

($#ARGV == 0) || die "Usage: $0 <fasta_dir>\n"; 

# Find all .fasta files in fasta dir
my $dir = $ARGV[0];
opendir DIR, $dir or die "error: cannot open directory \"$dir\" $!";
my @seqs = sort grep (!/^\.$|^\.\.$/, readdir(DIR));#@seqs include the whole path

@seqs = &sort_file_numerically(\@seqs);
# my @seqs = sort { $a <=> $b } grep (!/^\.$|^\.\.$/, readdir(DIR));
# my @files = sort { $a <=> $b } readdir(DATA_DIR);
closedir (DIR);

print "Before randomization, the array is: @seqs\n";
print "------------------------------------------------------------------------\n";

my $num = @seqs;#number of sequences

#This is to randomly select n sequences from @seq
my @items;
# my $n = $num;
my $n = 10;
  for ( 1 .. $n )
  {
    push @items, splice @seqs, rand @seqs, 1;
  }
@seqs = @items;#already randomly sorted

print "Before randomization, the array is: @seqs\n";

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
    
    