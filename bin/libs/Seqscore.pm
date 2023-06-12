#!/usr/bin/perl -w
#
# Package of functions for assigning a germline expression score to a sequence
#
#################################################################

package Seqscore;

use strict;
use warnings;
#use BerkeleyDB;

##########################################################################################

# Syntax: ( $sequence_score, $lowest_scoring_word, $numer_with_lowest_score ) = score_sequence( \@sequence, $pointer_to_BerkeleyDB_word_library )

sub score_sequence {
    my $CDSref = shift;
    my $CDS_length = @$CDSref;
    my $sequence_lib = shift;
    
    # Calculate the total score and the lowest word score for the sequence
    my $totalscore = 0;
    my $lowscore;
    my $words_w_lowscore;
    B: for ( my $b = 0; $b <= ($CDS_length); $b += 4 ) {
        my $wordscore = score_word($b, $CDSref, $sequence_lib);
        $totalscore += $wordscore;
        if ( $b == 0 || $wordscore < $lowscore) {
            $lowscore = $wordscore;
            $words_w_lowscore = 1;
        } elsif ($wordscore == $lowscore) {
            $words_w_lowscore++;
        }
    }
    
    # Normalize the score to the sequence length
    my $normscore = $totalscore / $CDS_length;
    
    # Return the results
    return ( $normscore, $lowscore, $words_w_lowscore );
}

###########################################################################################

#Syntax: $wordscore = score_word( $position, \@sequence, $pointer_to_BerkeleyDB_word_library )

sub score_word {
    my $position = shift;
    my $CDSref = shift;
    my $length = @$CDSref;
    my $sequence_lib = shift;
    local $" = '';

    #Score for word that overlaps by 3 nt on the 5' side
    my $L1score = 0;
    if ( $position >= 3 && $position <= $length - 1) {
        #my $L1word = join('', @$CDSref[$position-3..$position]);
        $L1score = $sequence_lib ->{join('', @$CDSref[$position-3..$position])};
        #local $" = '';
        #$L1score = $sequence_lib ->{@$CDSref[$position-3..$position]};
        #$L1score = $sequence_lib ->{ "@{$CDSref}[$position-3..$position]"  };
        #$sequence_lib -> db_get($L1word, $L1score);
    }
    
    #Score for word that overlaps by 6 nt on the 5' side
    my $L2score = 0;
    if ( $position >= 2 && $position <= $length - 2 ) {
        #my $L2word = join('', @$CDSref[$position-2..$position+1]);
        $L2score = $sequence_lib ->{join('', @$CDSref[$position-2..$position+1])};
        #$L2score = $sequence_lib ->{ "@{$CDSref}[$position-2..$position+1]"  };
        #$sequence_lib -> db_get($L2word, $L2score);
    }
    
    #Score for word that overlaps by 9 nt on the 5' side
    my $L3score = 0;
    if ( $position >= 1 && $position <= $length - 3 ) {
        #my $L3word = join('', @$CDSref[$position-1..$position+2]);
        $L3score = $sequence_lib ->{join('', @$CDSref[$position-1..$position+2])};
        #$L3score = $sequence_lib ->{ "@{$CDSref}[$position-1..$position+2]"  };
        #$sequence_lib -> db_get($L3word, $L3score);
    }
    
    #Score for central word
    my $mainscore = 0;
    if ( $position <= $length - 4 ) {
        #my $mainword = join('', @$CDSref[$position..$position+3]);
        $mainscore = $sequence_lib ->{join('', @$CDSref[$position..$position+3])};
        #$mainscore = $sequence_lib ->{ "@{$CDSref}[$position..$position+3]"  };
        #$sequence_lib -> db_get($mainword, $mainscore);
    }
        
    #Score for word that overlaps by 9 nt on the 3' side
    my $R3score = 0;
    if ( $position <= $length - 5 ) {
        #my $R3word = join('', @$CDSref[$position+1..$position+4]);
        $R3score = $sequence_lib ->{join('', @$CDSref[$position+1..$position+4])};
        #$R3score = $sequence_lib ->{ "@{$CDSref}[$position+1..$position+4]"  };
        #$sequence_lib -> db_get($R3word, $R3score);
    }
    
    #Score for word that overlaps by 6 nt on the 3' side
    my $R2score = 0;
    if ( $position <= $length - 6 ) {
        #my $R2word = join('', @$CDSref[$position+2..$position+5]);
        $R2score = $sequence_lib ->{join('', @$CDSref[$position+2..$position+5])};
        #$R3score = $sequence_lib ->{ "@{$CDSref}[$position+2..$position+5]"  };
        #$sequence_lib -> db_get($R2word, $R2score);
    }

    #Score for word that overlaps by 3 nt on the 3' side
    my $R1score = 0;
    if ( $position <= $length - 7 ) {
        #my $R1word = join('', @$CDSref[$position+3..$position+6]);
        $R1score = $sequence_lib ->{join('', @$CDSref[$position+3..$position+6])};
        #$R1score = $sequence_lib ->{ "@{$CDSref}[$position+3..$position+6]"  };
        #$sequence_lib -> db_get($R1word, $R1score);
    }
    
    #Add up scores
    my $wordscore = 0.25*$L1score + 0.5*$L2score + 0.75*$L3score + $mainscore + 0.75*$R3score + 0.5*$R2score + 0.25*$R1score;
    
    return $wordscore;
}

1;
