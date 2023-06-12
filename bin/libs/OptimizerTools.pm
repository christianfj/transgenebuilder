#!/usr/bin/perl -w

# Package of functions to optimize an input sequence for germline expression, based on RNAseq expression data
# For now only amino acid sequence input is supported
#
#################################################################

package OptimizerTools;

use strict;
use warnings;
use BerkeleyDB;
use Math::Random;
use POSIX;
use Progress::Any;
use lib './libs/';
use Seqscore;
use Memoize;

#use Data::Dumper;

###################################################################

sub optimize {
    my $progress = Progress::Any -> get_indicator( task => 'optimize' );
    $progress -> update( message => 'Setting up...');
    
    # Get input
    my $sequence_lib = shift;
    my $AA = shift;
    my $AA_length = length($AA);
    my @AA = split //, $AA;
    
    # Set up data structure: Array of hashes, with each position in the array representing one amino acid
    my @seqdata; 
    B: for ( my $b = 0; $b < $AA_length - 3; $b++ ) {  #The last three residues are part of the last 4mer and so don't have explicit array indices
        $seqdata[$b] = { position => $b,                        # The amino acid position
                         AAword => join ('',@AA[$b..$b+3]),     # The amino acid 4mer beginning at this position
                         all_NTwords => {},                     # An empty hash that will eventually hold all possible coding 12mers and their scores
                         top_NTword => '',                      # The highest-scoring nucleotide sequence coding for this 4mer
                         top_NTword_score => 0,                 # The score for that nucleotide sequence
                         assigned_NTword => '',                 # The nucleotide sequence chosen for this 4mer in the final sequence
                         assigned_NTword_index => 0,            # The index (rank) of the assigned word
                         assigned_NTword_score => 0 }           # The score for that nucleotide sequence
    }
    
    
    # Set up codon table
    my %codons = (
        A => ['GCA', 'GCC', 'GCG', 'GCT'],
        C => ['TGC', 'TGT'],
        D => ['GAC', 'GAT'],
        E => ['GAA', 'GAG'],
        F => ['TTC', 'TTT'],
        G => ['GGA', 'GGC', 'GGG', 'GGT'],
        H => ['CAC', 'CAT'],
        I => ['ATA', 'ATC', 'ATT'],
        K => ['AAA', 'AAG'],
        L => ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
        M => ['ATG'],
        N => ['AAC', 'AAT'],
        P => ['CCA', 'CCC', 'CCG', 'CCT'],
        Q => ['CAA', 'CAG'],
        R => ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
        S => ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
        T => ['ACA', 'ACC', 'ACG', 'ACT'],
        V => ['GTA', 'GTC', 'GTG', 'GTT'],
        W => ['TGG'],
        Y => ['TAC','TAT']
    );
    
    
    #### Generate and score all possible 12mers coding for each 4mer in the amino acid sequence ####
    C: for ( my $c = 0; $c < $AA_length - 3; $c++ ) {   #Loops B and C could be combined, but it's ok like this and more intuitive
        my @AAword = split //, $seqdata[$c]->{'AAword'};
        
        # Generate and score all possible nucleotide sequences coding for this 4mer
        D1: foreach my $d1 ( @{$codons{$AAword[0]}} ) {
            my $coding_sequence1 = $d1;
            D2: foreach my $d2 ( @{$codons{$AAword[1]}} ) {
                my $coding_sequence2 = $coding_sequence1 . $d2;
                D3: foreach my $d3 ( @{$codons{$AAword[2]}} ) {
                    my $coding_sequence3 = $coding_sequence2 . $d3;
                    D4: foreach my $d4 ( @{$codons{$AAword[3]}} ) {
                        my $coding_sequence4 = $coding_sequence3 . $d4;
                        my $score = $sequence_lib ->{$coding_sequence4};
                        #$sequence_lib -> db_get($coding_sequence4, my $score);          # Get the score for the word
                        $seqdata[$c]->{'all_NTwords'}->{$coding_sequence4} = $score;    # Put the score in the big data structure
                    }
                }
            }
        }
        
        # Sort the words according to their scores and find the highest score
        $seqdata[$c]->{'all_NTwords_sorted'} = [ sort { $seqdata[$c]->{'all_NTwords'}->{$b} <=> $seqdata[$c]->{'all_NTwords'}->{$a} } keys(%{$seqdata[$c]->{'all_NTwords'}}) ];
        my $top_NTword = $seqdata[$c]->{'all_NTwords_sorted'}->[0];
        $seqdata[$c]->{'top_NTword'} = $top_NTword;
        $seqdata[$c]->{'assigned_NTword'} = $top_NTword;
        $seqdata[$c]->{'top_NTword_score'} = $seqdata[$c]->{'all_NTwords'}->{$top_NTword};
        $seqdata[$c]->{'assigned_NTword_score'} = $seqdata[$c]->{'all_NTwords'}->{$top_NTword};
        # print Dumper($seqdata[$c]);
    }
    
    
    #### Take a first guess at an optimized sequence by plugging in the highest-scoring words first, then filling in the gaps with the best available options ###
    $progress->update( message => 'Making a first guess at an optimized sequence...' );
    
    # Start by producing a list of indices with those pointing to the highest-scoring words at the top
    my @indices = 0 .. $AA_length-4;
    my @indices_sorted = sort { $seqdata[$b]->{'top_NTword_score'} <=> $seqdata[$a]->{'top_NTword_score'} } @indices;
    
    # Now assign words into the sequence, working from highest to lowest scores
    my @coding_sequence = ('') x $AA_length;
    E: while (@indices_sorted) {
        my $e = shift @indices_sorted;
        my @NT_word = unpack("(A3)*", $seqdata[$e]->{'assigned_NTword'});  # @NT_word will contain the four codons for the word we want to plug in
        
        if ( $coding_sequence[$e] || $coding_sequence[$e+1] || $coding_sequence[$e+2] || $coding_sequence[$e+3] ) {
            # If some of these positions are already filled, see if our top-scoring word matches what we already have
            # Figure out what we already have
            my @current_word = @coding_sequence[$e..$e+3];
            E1: for (my $e1 = 0; $e1 < @current_word; $e1++) {
                unless ($current_word[$e1]) {
                    $current_word[$e1] = '...';    
                }
            }
            my $currword_str = join '', @current_word;
            # Compare to the top word currently being considered
            my $NT_word_str = join '', @NT_word;
            if ( $NT_word_str =~ /$currword_str/) {
                # If it matches, it's easy - put the word into place
                @coding_sequence[$e..$e+3] = @NT_word;
            } else {
                # If it doesn't match, move down the list to the next highest-scoring possibility for this position
                my $index = $seqdata[$e]->{'assigned_NTword_index'};
                $index++;
                $seqdata[$e]->{'assigned_NTword_index'} = $index;
                my $new_assigned_word = $seqdata[$e]->{'all_NTwords_sorted'}->[$index];
                $seqdata[$e]->{'assigned_NTword'} = $new_assigned_word;
                $seqdata[$e]->{'assigned_NTword_score'} = $seqdata[$e]->{'all_NTwords'}->{$new_assigned_word};
                # Put this index back in the list and re-sort the list based on the new score
                my @new_indices = @indices_sorted;
                push @new_indices, $e;
                @indices_sorted = sort { $seqdata[$b]->{'assigned_NTword_score'} <=> $seqdata[$a]->{'assigned_NTword_score'} } @new_indices;
            }
        
        } else { # If all positions are empty, it's simpler
        @coding_sequence[$e..$e+3] = @NT_word;  # Put the word in place
        }
    }
    
    # Calculate the score for the first-pass sequence and save the results
    my ( $first_sequence_score, $first_lowest_score, $first_n_w_lowest_score ) = Seqscore::score_sequence( \@coding_sequence, $sequence_lib );
    my %first_results;
    $first_results{'Sequence'} = join('',@coding_sequence);
    $first_results{'Sequence_score'} = $first_sequence_score;
    $first_results{'Lowest_score'} = $first_lowest_score;
    $first_results{'Words_w_lowest_score'} = $first_n_w_lowest_score;
    
    memoize('Seqscore::score_word');
    ####  Interatively try to improve the sequence by changing one word at a time  ####
    my $refine_counter = 0;
    my $maxiter = -4.6052 / log( 1 - 1/($AA_length-3));  # Iterate enough times that every residue has a 99% chance of being hit at least once.  ln(0.01) = -4.6052
    my @optimized_sequences;
    
    H: for (my $h = 0; $h < 10; $h++ ) {  # Since the interative process is random, do it 10 times to make sure we don't get stuck in a local minimum
        $progress->update( message => join('', 'Refining sequence... optimization round ', $h+1, ' of 10'  ) );
        my $counter = 0;  # Keeps a count of the number of iterations done.
        
        my @new_seq = @coding_sequence; # Make a copy of @coding_sequence so we can reuse it later
        my $current_sequence_score = $first_sequence_score;
        F: while ($counter <= $maxiter) {
            # Pick a random position
            my $pos = random_uniform_integer(1, 0, $AA_length-4);
            
            # Calculate the context score for each possible word at that position
            my @context = @new_seq[ ($pos >= 3 ? $pos-3 : 0) .. ($pos <= $AA_length-7 ? $pos+6 : $AA_length-1) ];
            my $rel_pos = ( $pos >= 3 ? 3 : $pos );
            my %context_scores;
            G: foreach my $g (@{$seqdata[$pos]->{'all_NTwords_sorted'}}) {
                @context[$rel_pos..$rel_pos+3] = unpack("(A3)*", $g);
                $context_scores{$g} = Seqscore::score_word($rel_pos, \@context, $sequence_lib);
            }
            
            # Sort and find the highest scoring word given the context
            my @words_sorted = sort { $context_scores{$b} <=> $context_scores{$a} } keys(%context_scores);
            my $top_word = $words_sorted[0];
            
            if ($top_word eq join('', @new_seq[$pos..$pos+3])) {
                # If these are equal, we already have the best-scoring word; move on
                $counter++;
                next F;
            } else {
                # Otherwise, see if changing the word will improve the score for the whole sequence
                my @test_seq = @new_seq;
                @test_seq[$pos..$pos+3] = unpack("(A3)*",$top_word);
                my ( $test_score, $test_lowest_score, $test_n_w_lowest_score ) = Seqscore::score_sequence( \@test_seq, $sequence_lib );
                if ( $test_score > $current_sequence_score ) {
                    # If the change improves the score, keep it
                    $refine_counter++;
                    @new_seq = @test_seq;
                    $current_sequence_score = $test_score;
                    $counter = 0;
                } else {
                    # Otherwise, move on
                    $counter++;
                }
            }
        }
        
        # Calculate the score for the optimized sequence and save the results
        $progress->update(message => 'Finishing...');
        
        my ( $final_sequence_score, $final_lowest_score, $final_n_w_lowest_score ) = Seqscore::score_sequence( \@new_seq, $sequence_lib );
        my $new_seq = join('',@new_seq);
        I: foreach my $i (@optimized_sequences) {   # This loop ensures that we only save unique results.  Might not be the most efficient way to do it, but good enough
            next H if ( $i->{'Sequence'} eq $new_seq );
        }
        my %results;
        $results{'Sequence'} = $new_seq;
        $results{'Sequence_score'} = $final_sequence_score;
        $results{'Lowest_score'} = $final_lowest_score;
        $results{'Words_w_lowest_score'} = $final_n_w_lowest_score;
        push @optimized_sequences, \%results;
    }
    
    # Return the results
    my @optimization_results = sort { $b->{'Sequence_score'} <=> $a->{'Sequence_score'} } @optimized_sequences;
    my @holder;
    return ($optimization_results[0], @holder);
    
}

################################################################################################

sub addintrons {
    my $inseq = shift;
    my @CDS = split //, $inseq;
    
    # Figure out how many introns to add and their approximate positions
    my $seq_length = length($inseq);
    my $n_introns = ceil($seq_length / 300);
    if ($n_introns == 0 ) {  # Makes sure at least one intron is added, since we only call this function if we want to add introns
        $n_introns = 1;
    }
    my $exon_length = floor( $seq_length / ($n_introns+1) );
    
    # Split the sequence into exons, add an intron after each
    my $last_exon_end = 0;
    my $outseq = '';
    L: for (my $l = 0; $l < $n_introns; $l++) {
        
        # Find the beginning of the exon
        my $exon_start;
        if ($l == 0) {
            $exon_start = 0;
        } else {
            $exon_start = $last_exon_end + 1;
        } # If this is the first exon, it starts at the beginning of the sequence; otherwise, it starts at the next base after the end of the last exon.
        
        # Find the end of the exon
        my $exon_end;
        my $offset = 0;
        my $found_junction = 0;
        M: while (!$found_junction) { # This loop looks for the AGR motif
            my $start_point = $exon_length * ($l+1) + $offset;
            if ( $start_point <= $exon_start + 20 ) { # This means we can't find a place to put the intron without making the preceeding exon too short (<20 bp)
                $outseq = $outseq . $CDS[$last_exon_end];
                next L;
            }
            my $candidate_junction = join ('', @CDS[ $start_point .. $start_point+2 ] );
            if ( $candidate_junction =~ /AG[AG]/ ) {
                $exon_end = $start_point + 1;
                $found_junction = 1;
            } else {
                $offset = $offset - 1;
            }
        }
        $last_exon_end = $exon_end;
        
        # Add the exon to the output sequence
        $outseq = $outseq . join ('', @CDS[ $exon_start .. $exon_end]);
        
        # Generate the intron
        my $good_intron = 0;
        N: while (!$good_intron) {
            my $candidate_intron =  random_dna(0.3, 36); # Parameters: GC content (between 0 and 1), length
            next N if ($candidate_intron =~ /[ct]ag/ ); # Makes sure that the intron doesn't contain splice acceptors
            $good_intron = $candidate_intron;
        }
        my $intron = join ('', 'gtaagttt', $good_intron, 'ttttcag');
        
        # Add the intron to the output sequence
        $outseq = $outseq . $intron;
    }
    
    # Add the last exon to the sequence
    $outseq = $outseq . join ('', @CDS[ $last_exon_end+1 .. $#CDS] );
    
    # Return the results
    return $outseq;
}

################################################################################################

sub word_options {
    my $position = shift;
    my $seqdata_ref = shift;
    my $currword_ref = shift;
    my @currword_sub = @$currword_ref;
    Z: for (my $z = 0; $z < @currword_sub; $z++) {
        unless ($currword_sub[$z]) {
            $currword_sub[$z] = '...';    
        }
    }
    my $currword_str = join '', @currword_sub;
    
    my $options_sub = {};
    Y: foreach my $y ( @{$$seqdata_ref[$position]->{'all_NTwords_sorted'}} ) {
        if ( $y =~ /$currword_str/) {
            $options_sub->{$y} = 0;
        }
    }
    
    return $options_sub;
}

################################################################################################

sub random_dna {
    my $pctGC = shift;
    my $length = shift;
    my @sequence;
    W: for (my $w = 0; $w < $length; $w++) {
        my $letter;
        if (random_uniform() < $pctGC) {
            if (random_uniform() < 0.5) {
                $letter = 'g';
            } else {
                $letter = 'c';
            }
        } else {
            if (random_uniform() < 0.5) {
                $letter = 'a';
            } else {
                $letter = 't';
            }
        }
        push @sequence, $letter;
    }
    my $sequence = join('', @sequence);
    return $sequence;
}

################################################################################################

1;
