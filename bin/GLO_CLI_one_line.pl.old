#!/usr/bin/perl

#####Documentation######
$num_args = $#ARGV + 1;
if ($num_args != 1) {
  print "\nUsage: GLO_CLI_one_line.pl DNA_sequence_fasta_file\n\n";
  print "\tSeq_fasta: Fasta file containing a single sequence to optimize\n";
  print "\tPlease note code was originally from Dan Dickinson but adapted for our needs. Also sequence_lib_scores.db file should be present\n";

  exit;
}
############
#Load libraries

use BerkeleyDB;
use Bio::Seq;
use Math::Random;
use IO::Handle;
use lib './bin/libs/';
use Seqscore;
use OptimizerTools;


##Reading and verifying input commands
my $Sqinps=$ARGV[0];

###If the following files are not present -> stop the program
if (!(-e $Sqinps."")){die $Sqinps." file not found!\n";}

##Initialize parameters
my $error = '';

# Get database
my $dbfilename = 'bin/sequence_lib_scores.db';
my $sequence_lib = new BerkeleyDB::Btree
	-Filename => $dbfilename
	or die "Cannot open $dbfilename: $! $BerkeleyDB::Error\n" ;
        
# Get input
my $dnaseq = '';
my $seqtype = 'DNA';
my $seqname = 'Input';

### If a nucleotide sequence was entered, calculate its score ###

my ( $input_sequence_score, $input_lowest_score, $input_n_w_lowest_score );

my $results = {}; #Get a pointer to an empty array that will hold the results

my @input_coding_sequence;
##Reading and hashing genomic sequence
open(op,$Sqinps) or die "cannot open Ref file\n";
while($line=<op>){
    chomp($line);
if($line =~ m/^>(.+)/){
    $seqname = $1;
    }else{
$dnaseq=$dnaseq.$line;
}}
close(op);

##Use BioPerl to get aminoacid sequence
my $seqobj = Bio::Seq->new( -seq => $dnaseq );

    my $trans = $seqobj->translate();
    $AAseq = $trans->seq();
    $DNAseq = $seqobj->seq();
    
    my $stopflag = 0;
    my $AAlength = length($AAseq);
    if ( $AAseq =~ /\*/ && $-[0] != $AAlength ) {
        #$warning = "You entered a nucleotide sequence that contains an internal stop codon.  Only the portion before the stop codon will be optimized.";
        $AAseq = substr($AAseq, 0, $-[0]);
        $DNAseq = substr($DNAseq, 0, $-[0] * 3 );
        $stopflag = 1;
    }


##Score sequence
#@input_coding_sequence = unpack("(A3)*", $dnaseq);
#( $results->{'input_sequence_score'}, $results->{'input_lowest_score'}, $results->{'input_n_w_lowest_score'} ) = Seqscore::score_sequence( \@input_coding_sequence, $sequence_lib );

#print "Input sequence:\n";        
#print $DNAseq;
#print "\nSequence score:\n";
#print "\t".$results->{'input_sequence_score'}."\n";

### Optimize the sequence ###
my $optimization_results = OptimizerTools::optimize($sequence_lib, $AAseq);
$results = {%$results, %$optimization_results};

#print "Output sequence:\n";                
if ($stopflag == 1){
print $optimization_results->{'Sequence'}."TAA\n";
}else{
  print $optimization_results->{'Sequence'}."\n";
}
#$results = {}; #Get a pointer to an empty array that will hold the results
#@input_coding_sequence = unpack("(A3)*", $optimization_results->{'Sequence'});
#( $results->{'input_sequence_score'}, $results->{'input_lowest_score'}, $results->{'input_n_w_lowest_score'} ) = Seqscore::score_sequence( \@input_coding_sequence, $sequence_lib );

#print "\nSequence score:\n";
#print "\t".$results->{'input_sequence_score'}."\n";

exit;
##########################################################################################
