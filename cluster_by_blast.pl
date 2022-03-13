#!/usr/bin/perl 
# This is to use BLAST to do the clustering job. The input should be the sequence files to be clustered. The output will be the ACs that have been clustered.
# Usage:
#  ./cluster_by_blast.pl <fasta_file> <clustered AC>
# Example:
#  ./cluster_by_blast.pl ../SeqData/eu2423/train/1_1.fasta 1 100
#  ./cluster_by_blast.pl ~/so/java/Workspace/GoaSvmServer/WebContent/test_data/test.fasta 1 10

#by Shibiao WAN, 09-Aug-2016, Princeton, NJ, USA

use strict;
use warnings;
use Switch;
use Carp;
use File::Basename;
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

($#ARGV == 1) || die "Usage: $0 <fasta_file> <clustered AC>\n"; 

# Define constant and input parameters to psi-blast
my $fastafile = $ARGV[0];

my $evalue = $ARGV[2]; #For single-label ,Evalue =100; multi-label, Evalue = 10.
my $basename = fileparse($fastafile, ".fasta");
################Make sure these files are in the right directories for different deployments in different accounts############
my $blastdir = "/home4a/sbwan/blast";
#my $blastdir = "/home5a/mwmak/so/bioinfo/SeqAnalysis/blast";
my $logdir = "/home4a/sbwan/java/tmp3";
#my $logdir = "/home5a/mwmak/so/java/Workspace/mGoaSvmServer2/tmp";

################use a compact sequence database#############################
my $database = "/home4a/sbwan/Interpro/DBprocessing/proseq_DB/proseq.fasta";
#my $database = "/home4a/sbwan/databases/proseq_2013_07/proseq_sprot.fasta";

#my $acfile = "/home5a/mwmak/so/java/Workspace/mGoaSvmServer2/WebContent/new_ac.txt";
my $acfile = "/home4a/sbwan/java/workspace/mGoaSvmServer/WebContent/new_ac.txt";
my $maxitr = 1;
my $hvalue = 0.001;
my $maxlength = 1000;
##my $evalue = 100;

# Read .fasta file process sequences one by one
my @seqs = read_all_sequences($fastafile,'fasta');

my $num = @seqs;#number of sequences

while($num - 2*i > 3)

# Call PSI-BLAST and save homolog info in $logfile
my $j=1;
foreach my $seq (@seqs) {
    #print ">",$seq->display_id,"\n";
    #my $tempSeqFile = "/tmp/mwmak/${basename}-${j}.fasta";
    my $tempSeqFile = "${logdir}/${basename}-${j}.fasta";
    my $seqio_obj = Bio::SeqIO->new(-file => ">$tempSeqFile", -format => 'fasta' );
    $seq = $seq->trunc(1,min $maxlength,length($seq->seq()));
    $seqio_obj->write_seq($seq);
    my $logfile = "${logdir}/${basename}-${j}.log";
    my $cmd = "${blastdir}/bin/blastpgp -j $maxitr -h $hvalue -e $evalue -m 6 -d $database  -i $tempSeqFile -o $logfile > /dev/null 2>&1";
    #print "$cmd\n";
    system($cmd);
    my $ac;
    my $Evalue;
    ($ac,$Evalue) = &get_acEvalue_from_logfile($logfile);
    print "$ac   $Evalue\n";
    $j++;
}

sub get_acEvalue_from_logfile {
    my $logfile = $_[0];
    #my $acfile =$_[1];
    #my $mark = $_[2];
    my $ac = '--';
    my $Evalue = '--';
    unless(open (LOGFILE, "$logfile")){croak "Cannot open $logfile\n";}
    
    my $i=0;
    my $flag = 0;
    while (my $line=<LOGFILE>) {
        next if ($.<30);#skip the first 30 lines
        chomp $line;
	if ($line =~ m/\|[A-Z0-9]{6}\|/){
            $flag = 1;
            $ac = (split(/\|/,$line))[1]; 
	    $Evalue = (split(/\s{2,}/,$line))[2];
            #remove the whitespace on the right side
            $ac =~ s/\s*$//;#because the ac is surely valid, we do not need to do check
            $Evalue =~ s/\s*$//;
            last;
        }
    }

    if ($flag == 0){$ac = "XXXXXX"; $Evalue = -100;}
    
    my $j = $i + 1;
    #print "The $j -th homolog is used\n";
    close LOGFILE or die "Cannot close $logfile\n";
    return ($ac,$Evalue);
}

exit;

