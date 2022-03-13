#!/usr/bin/perl 
# This is to process the files in a directory.
# Usage:
#  scripts/process_files_in_dir.pl <fasta dir> 
# Example:
#  scripts/process_files_in_dir.pl ~/so/bioinfo/SeqAnalysis/SeqData/eu2423/train

use strict;
use warnings;
use Switch;
use Carp;
use Bio::SeqIO;

($#ARGV == 0) || die "Usage: $0 <fasta dir>\n"; 

# Find all .fasta files in fasta dir
my $dir = $ARGV[0];
opendir DIR, $dir or die "error: cannot open directory \"$dir\" $!";
my @files = sort grep (!/^\.$|^\.\.$/&& -f "$dir/$_", readdir(DIR));
closedir (DIR);

my $num = @files;

my @newseqs;
# Print the sequence id of all fasta files
for (my $i = 0; $i < $num; $i++) {
    if ($i<$num-2){
#     print "$files[$i]\n";
#     my $seqn = "$dir/$files[$i]";
    my $seqn = $files[$i];
    push @newseqs, $seqn;
    }
    
#     my $in = Bio::SeqIO->new(-file => "$dir/$file", '-format' => 'Fasta');
#     my $seq = $in->next_seq();
#     my @field = split(/_/,$file);
#     my $label = $field[0];
#     printf ANFILE "%s\n",$seq->display_id;
#     printf LBFILE "%s\n",$label;
}
    
# my $cmd = "cat @newseqs > newfile.fasta";
# system($cmd);

# the difference of @files and @newseqs (i.e., @files/@newseqs)
my %newseqs=map{$_=>1} @newseqs;
my @restseq=grep(!defined $newseqs{$_}, @files);
# print "The rest files are:\n";
print "@restseq\n";



