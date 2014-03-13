$|=1;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

$datalocation = "/panfs/sequencers/sequencerO/runs/140121_SN7001144_0462_AH8896ADXX/Data/Intensities/BaseCalls/Unaligned/Project_Na238/Sample_010114-4_012/010114-4_012_CTTGTA_L002_R1_001.fastq.gz";
@outputfilename = split("/",$datalocation);
$outputfilename = $outputfilename[$#outputfilename];
$outputfilename =~ s/.fastq.gz/.output.txt/;
open OUTFILE, ">" . $outputfilename;

$chr4barcode1 = "GAGCGA";
$chr4barcode2 = "GATGTG";
$CFTRbarcode1 = "CGGAAT";
$CFTRbarcode2 = "CGTGGC";
$MCFTRbarcode1 = "TGCGTA";
$MCFTRbarcode2 = "TTCTGG";

$CFTR_WT = "TGTTCTCAGTTTTCCTGGATTATGCCTGGCACCATTAAAGAAAATATCATTGGTGTTTCCTATGATGAATATAGATACAGAAGCGTCATCAAAGCATGCCA";
$CFTR_WT_interest = "AGAAAATATCATTGG";
$CFTR_mut = "TGTTCTCAGTTTTCCTGGATTATGCCTGGCACCATTAAGGAGAACATTATCTTTGGTGTTTCCTATGATGAATATAGATACAGAAGCGTCATCAAAGCATGCCA";
$CFTR_mut_interest = "GGAGAACATTATCTTTGG";
$MCFTR_WT = "TCTGCTCTCAATTTTCTTGGATTATGCCGGGTACTATCAAAGAAAATATCATTGGTGTTTCCTATGATGAGTACAGATATAAGAGTGTTGTCAAAGCTTGCC";
$MCFTR_WT_interest = "AGAAAATATCATTGG";
$MCFTR_mut = "TCTGCTCTCAATTTTCTTGGATTATGCCGGGTACTATCAAGGAGAACATTATCTTTGGTGTTTCCTATGATGAGTACAGATATAAGAGTGTTGTCAAAGCTTGCC";
$MCFTR_mut_interest = "GGAGAACATTATCTTTGG";
$chr4_WT = "TTCTACTAAAAGAAAACTTCTGTGTCCTAATTAAAGAAAATATCCTATTTGGTCATTTATCATGAAGATAATGATAAATGTTAAATAACGTCTTGTTTGCATTAAGTCTGTGGGATG";
$chr4_WT_interest = "ATTAAAGAAAATATCCTATTTGGTCATTTATCATGAAGATA";
$chr4_WT_interest_truncated = "TTTGGTCATTTATCATGAAGATA";

my $input = new IO::File "<$datalocation" or die "\nError opening $datalocation\n";
my $infile = new IO::Uncompress::Gunzip $input;

my $CFTR_WT_count = 0;
my $CFTR_mut_count = 0;
my $CFTR_nomatch_count = 0;
my $CFTR_WT_interest_mutation_count = 0;
my $CFTR_WT_flanking_mutation_count = 0;
my $CFTR_mut_interest_mutation_count = 0;
my $CFTR_mut_flanking_mutation_count = 0;
my $MCFTR_WT_count = 0;
my $MCFTR_mut_count = 0;
my $MCFTR_nomatch_count = 0;
my $MCFTR_WT_interest_mutation_count = 0;
my $MCFTR_WT_flanking_mutation_count = 0;
my $MCFTR_mut_interest_mutation_count = 0;
my $MCFTR_mut_flanking_mutation_count = 0;
my $chr4_WT_count = 0;
my $chr4_WT_interest_mutation_count = 0;
my $chr4_WT_flanking_mutation_count = 0;
my $chr4_nomatch_count = 0;

my $counter = 0;
while (my $header = <$infile>) { # This is for scanning all lines
#foreach (1..1000) { ######### This is for scanning 1000 lines
	if (($counter % 100) == 0) {
		print $counter, "\n";
	}
	$counter++;
#	my $header = <$infile>; ######## This is for scanning 1000 lines
	my $rawseq = <$infile>;
	my $garbage = <$infile>;
	my $garbage = <$infile>;
	chomp $rawseq;
	my $rawseq_barcode = substr($rawseq,0,6);
#	my $rawseq_right_barcode = substr($rawseq,(length($rawseq) - 6));
#	print length($rawseq) - 6, " ", $rawseq_barcode, " ", $rawseq_right_barcode, " ", reversecomplement($rawseq_right_barcode), "\n"; ########
#	if (countStringMatches($rawseq_barcode,reversecomplement($rawseq_right_barcode)) < 5) { # Chimera
#		next;
#	}
	my $seq_type;
	if (countStringMatches($rawseq_barcode,$chr4barcode1)  >= 5) {
		$seq_type = "chr4";
	} elsif (countStringMatches($rawseq_barcode,$chr4barcode2)  >= 5) {
		$seq_type = "chr4";
	} elsif (countStringMatches($rawseq_barcode,$CFTRbarcode1)  >= 5) {
		$seq_type = "CFTR";
	} elsif (countStringMatches($rawseq_barcode,$CFTRbarcode1)  >= 5) {
		$seq_type = "CFTR";
	} elsif (countStringMatches($rawseq_barcode,$MCFTRbarcode1)  >= 5) {
		$seq_type = "MCFTR";
	} elsif (countStringMatches($rawseq_barcode,$MCFTRbarcode2)  >= 5) {
		$seq_type = "MCFTR";
	}
	if ($seq_type eq "CFTR") {
		my ($WT_align1,$WT_align2) = alignSequences($rawseq,$CFTR_WT_interest);
		my $WT_align_score = countStringMatches($WT_align1,$WT_align2);
		my ($mut_align1,$mut_align2) = alignSequences($rawseq,$CFTR_mut_interest);
		my $mut_align_score = countStringMatches($mut_align1,$mut_align2);
		if ($WT_align_score >= (length($CFTR_WT_interest) - 2)) { # WT match
			$CFTR_WT_count++;
			$CFTR_WT_interest_mutation_count += (length($CFTR_WT_interest) - $WT_align_score);
			my ($full_align1,$full_align2) = alignSequences($rawseq,$CFTR_WT);
			my $WT_full_align_score = countStringMatches($full_align1,$full_align2);
			$CFTR_WT_flanking_mutation_count += ((length($CFTR_WT) - $WT_full_align_score) - (length($CFTR_WT_interest) - $WT_align_score));
		} elsif ($WT_align_score >= (length($CFTR_mut_interest) - 2)) { # mut match
			$CFTR_mut_count++;
			$CFTR_mut_interest_mutation_count += (length($CFTR_mut_interest) - $mut_align_score);
			my ($full_align1,$full_align2) = alignSequences($rawseq,$CFTR_mut);
			my $mut_full_align_score = countStringMatches($full_align1,$full_align2);
			$CFTR_mut_flanking_mutation_count += ((length($CFTR_mut) - $mut_full_align_score) - (length($CFTR_mut_interest) - $mut_align_score));
		} else { # No match
			$rawseq = reversecomplement($rawseq);
			my ($WT_align1,$WT_align2) = alignSequences($rawseq,$CFTR_WT_interest);
			my $WT_align_score = countStringMatches($WT_align1,$WT_align2);
			my ($mut_align1,$mut_align2) = alignSequences($rawseq,$CFTR_mut_interest);
			my $mut_align_score = countStringMatches($mut_align1,$mut_align2);
			if ($WT_align_score >= (length($CFTR_WT_interest) - 2)) { # WT match
				$CFTR_WT_count++;
				$CFTR_WT_interest_mutation_count += (length($CFTR_WT_interest) - $WT_align_score);
				my ($full_align1,$full_align2) = alignSequences($rawseq,$CFTR_WT);
				my $WT_full_align_score = countStringMatches($full_align1,$full_align2);
				$CFTR_WT_flanking_mutation_count += ((length($CFTR_WT) - $WT_full_align_score) - (length($CFTR_WT_interest) - $WT_align_score));
			} elsif ($WT_align_score >= (length($CFTR_mut_interest) - 2)) { # mut match
				$CFTR_mut_count++;
				$CFTR_mut_interest_mutation_count += (length($CFTR_mut_interest) - $mut_align_score);
				my ($full_align1,$full_align2) = alignSequences($rawseq,$CFTR_mut);
				my $mut_full_align_score = countStringMatches($full_align1,$full_align2);
				$CFTR_mut_flanking_mutation_count += ((length($CFTR_mut) - $mut_full_align_score) - (length($CFTR_mut_interest) - $mut_align_score));
			} else { # No match
				$CFTR_nomatch_count++;
			}
		}
	} elsif ($seq_type eq "MCFTR") {
#		if (substr($rawseq,46,15) eq $MCFTR_WT_interest) {
#			$WT_align_score = length($MCFTR_WT_interest);
#			$mut_align_score = 0;
#		} elsif (substr($rawseq,46,18) eq $MCFTR_mut_interest) {
#			$WT_align_score = 0;
#			$mut_align_score = length($MCFTR_mut_interest);
#		} else {
			my ($WT_align1,$WT_align2) = alignSequences($rawseq,$MCFTR_WT_interest);
			my $WT_align_score = countStringMatches($WT_align1,$WT_align2);
			my ($mut_align1,$mut_align2) = alignSequences($rawseq,$MCFTR_mut_interest);
			my $mut_align_score = countStringMatches($mut_align1,$mut_align2);
#		}
#		print $WT_align_score, " ", $mut_align_score, " ", $WT_align1, " ", $WT_align2, " ", $mut_align1, " ", $mut_align2, " ", $rawseq, "\n";
		if ($WT_align_score >= (length($MCFTR_WT_interest) - 2)) { # WT match
			$MCFTR_WT_count++;
			$MCFTR_WT_interest_mutation_count += (length($MCFTR_WT_interest) - $WT_align_score);
			my ($full_align1,$full_align2) = alignSequences($rawseq,$MCFTR_WT);
			my $WT_full_align_score = countStringMatches($full_align1,$full_align2);
			$MCFTR_WT_flanking_mutation_count += ((length($MCFTR_WT) - $WT_full_align_score) - (length($MCFTR_WT_interest) - $WT_align_score));
		} elsif ($WT_align_score >= (length($MCFTR_mut_interest) - 2)) { # mut match
			$MCFTR_mut_count++;
			$MCFTR_mut_interest_mutation_count += (length($MCFTR_mut_interest) - $mut_align_score);
			my ($full_align1,$full_align2) = alignSequences($rawseq,$MCFTR_mut);
			my $mut_full_align_score = countStringMatches($full_align1,$full_align2);
			$MCFTR_mut_flanking_mutation_count += ((length($MCFTR_mut) - $mut_full_align_score) - (length($MCFTR_mut_interest) - $mut_align_score));
		} else { # No match
			$rawseq = reversecomplement($rawseq);
			my ($WT_align1,$WT_align2) = alignSequences($rawseq,$MCFTR_WT_interest);
			my $WT_align_score = countStringMatches($WT_align1,$WT_align2);
			my ($mut_align1,$mut_align2) = alignSequences($rawseq,$MCFTR_mut_interest);
			my $mut_align_score = countStringMatches($mut_align1,$mut_align2);
			if ($WT_align_score >= (length($MCFTR_WT_interest) - 2)) { # WT match
				$MCFTR_WT_count++;
				$MCFTR_WT_interest_mutation_count += (length($MCFTR_WT_interest) - $WT_align_score);
				my ($full_align1,$full_align2) = alignSequences($rawseq,$MCFTR_WT);
				my $WT_full_align_score = countStringMatches($full_align1,$full_align2);
				$MCFTR_WT_flanking_mutation_count += ((length($MCFTR_WT) - $WT_full_align_score) - (length($MCFTR_WT_interest) - $WT_align_score));
			} elsif ($WT_align_score >= (length($MCFTR_mut_interest) - 2)) { # mut match
				$MCFTR_mut_count++;
				$MCFTR_mut_interest_mutation_count += (length($MCFTR_mut_interest) - $mut_align_score);
				my ($full_align1,$full_align2) = alignSequences($rawseq,$MCFTR_mut);
				my $mut_full_align_score = countStringMatches($full_align1,$full_align2);
				$MCFTR_mut_flanking_mutation_count += ((length($MCFTR_mut) - $mut_full_align_score) - (length($MCFTR_mut_interest) - $mut_align_score));
			} else { # No match
				$MCFTR_nomatch_count++;
			}
		}
	} elsif ($seq_type eq "chr4") {
#		if (substr($rawseq,35,41) eq $chr4_WT_interest) {
#			$WT_align_score = length($chr4_WT_interest);
#		} elsif (substr(reversecomplement($rawseq),0,23) eq $chr4_WT_interest_truncated) {
#			$WT_align_score = length($chr4_WT_interest);
#		} else {
			my ($WT_align1,$WT_align2) = alignSequences($rawseq,$chr4_WT_interest);
			my $WT_align_score = countStringMatches($WT_align1,$WT_align2);
#		}
#		print $WT_align_score, " ", $WT_align1, " ", $WT_align2, " ", $rawseq, "\n";
		if ($WT_align_score >= (length($chr4_WT_interest) - 2)) { # WT match
			$chr4_WT_count++;
			$chr4_WT_interest_mutation_count += (length($chr4_WT_interest) - $WT_align_score);
			my ($full_align1,$full_align2) = alignSequences($rawseq,$chr4_WT);
			my $WT_full_align_score = countStringMatches($full_align1,$full_align2);
			$chr4_WT_flanking_mutation_count += ((length($chr4_WT) - $WT_full_align_score) - (length($chr4_WT_interest) - $WT_align_score));
		} else {
			$rawseq = reversecomplement($rawseq);
			my ($WT_align1,$WT_align2) = alignSequences($rawseq,$chr4_WT_interest_truncated);
			my $WT_align_score = countStringMatches($WT_align1,$WT_align2);
			if ($WT_align_score >= (length($chr4_WT_interest_truncated) - 2)) { # WT match
				$chr4_WT_count++;
				$chr4_WT_interest_mutation_count += (length($chr4_WT_interest_truncated) - $WT_align_score);
				my ($full_align1,$full_align2) = alignSequences($rawseq,$chr4_WT);
				my $WT_full_align_score = countStringMatches($full_align1,$full_align2);
				$chr4_WT_flanking_mutation_count += ((length($chr4_WT) - $WT_full_align_score) - (length($chr4_WT_interest_truncated) - $WT_align_score));
			} else { # No match
				$chr4_nomatch_count++;
			}
		}
	}
}

print OUTFILE "CFTR WT seq\t", $CFTR_WT_count,"\n";
print OUTFILE "CFTR mut seq\t", $CFTR_mut_count,"\n";
print OUTFILE "CFTR no-match\t", $CFTR_nomatch_count,"\n";
print OUTFILE "CFTR ROI mutations in WT\t", $CFTR_WT_interest_mutation_count,"\n";
print OUTFILE "CFTR flanking mutations in WT\t", $CFTR_WT_flanking_mutation_count,"\n";
print OUTFILE "CFTR ROI mutations in WT\t", $CFTR_mut_interest_mutation_count,"\n";
print OUTFILE "CFTR flanking mutations in mut\t", $CFTR_mut_flanking_mutation_count,"\n";
print OUTFILE "MCFTR WT seq\t", $MCFTR_WT_count,"\n";
print OUTFILE "MCFTR mut seq\t", $MCFTR_mut_count,"\n";
print OUTFILE "MCFTR no-match\t", $MCFTR_nomatch_count,"\n";
print OUTFILE "MCFTR ROI mutations in WT\t", $MCFTR_WT_interest_mutation_count,"\n";
print OUTFILE "MCFTR flanking mutations in WT\t", $MCFTR_WT_flanking_mutation_count,"\n";
print OUTFILE "MCFTR ROI mutations in WT\t", $MCFTR_mut_interest_mutation_count,"\n";
print OUTFILE "MCFTR flanking mutations in mut\t", $MCFTR_mut_flanking_mutation_count,"\n";
print OUTFILE "Chr4 WT seq\t", $chr4_WT_count, "\n";
print OUTFILE "Chr4 ROI mutations in WT\t", $chr4_WT_interest_mutation_count,"\n";
print OUTFILE "Chr4 flanking mutations in WT\t", $chr4_WT_flanking_mutation_count,"\n";
print OUTFILE "Chr4 no-match\t", $chr4_nomatch_count ,"\n";

sub countStringMatches {
	my ($seq1, $seq2) = @_;
	if ($seq1 eq $seq2) {
		return length($seq1);
	}
	$n_matches = 0;
	foreach my $a (0..(length($seq1)-1)) {
		if (substr($seq1,$a,1) eq substr($seq2,$a,1)) {
			$n_matches++;
		}
	}
	return $n_matches;
}


sub alignSequences{
	my ($seq1, $seq2) = @_;

	# scoring scheme
	my $MATCH     =  1; # +1 for letters that match
	my $MISMATCH = -1; # -1 for letters that mismatch
	my $GAP       = -1; # -1 for any gap

	# initialization
	my @matrix;
	$matrix[0][0]{score}   = 0;
	$matrix[0][0]{pointer} = "none";
	for(my $j = 1; $j <= length($seq1); $j++) {
		 $matrix[0][$j]{score}   = 0;
		 $matrix[0][$j]{pointer} = "none";
	}
	for (my $i = 1; $i <= length($seq2); $i++) {
		 $matrix[$i][0]{score}   = 0;
		 $matrix[$i][0]{pointer} = "none";
	}

	# fill
	 my $max_i     = 0;
	 my $max_j     = 0;
	 my $max_score = 0;


	 for(my $i = 1; $i <= length($seq2); $i++) {
		 for(my $j = 1; $j <= length($seq1); $j++) {
			 my ($diagonal_score, $left_score, $up_score);
			 
			 # calculate match score
			 my $letter1 = substr($seq1, $j-1, 1);
			 my $letter2 = substr($seq2, $i-1, 1);      
			 if ($letter1 eq $letter2) {
				 $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
			  }
			 else {
				 $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
			  }
			 
			 # calculate gap scores
			 $up_score   = $matrix[$i-1][$j]{score} + $GAP;
			 $left_score = $matrix[$i][$j-1]{score} + $GAP;
			 
			 if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
				 $matrix[$i][$j]{score}   = 0;
				 $matrix[$i][$j]{pointer} = "none";
				 next; # terminate this iteration of the loop
			  }

			 
			 # choose best score
			 if ($diagonal_score >= $up_score) {
				 if ($diagonal_score >= $left_score) {
					 $matrix[$i][$j]{score}   = $diagonal_score;
					 $matrix[$i][$j]{pointer} = "diagonal";
				  }
				 else {
					 $matrix[$i][$j]{score}   = $left_score;
					 $matrix[$i][$j]{pointer} = "left";
				  }
			  } else {
				 if ($up_score >= $left_score) {
					 $matrix[$i][$j]{score}   = $up_score;
					 $matrix[$i][$j]{pointer} = "up";
				  }
				 else {
					 $matrix[$i][$j]{score}   = $left_score;
					 $matrix[$i][$j]{pointer} = "left";
				  }
			  }
			 
		   # set maximum score
			 if ($matrix[$i][$j]{score} > $max_score) {
				 $max_i     = $i;
				 $max_j     = $j;
				 $max_score = $matrix[$i][$j]{score};
			  }
		  }
	 }

	 # trace-back

	 my $align1 = "";
	 my $align2 = "";

	 my $j = $max_j;
	 my $i = $max_i;

	 while (1) {
		 last if $matrix[$i][$j]{pointer} eq "none";
		 
		 if ($matrix[$i][$j]{pointer} eq "diagonal") {
			 $align1 .= substr($seq1, $j-1, 1);
			 $align2 .= substr($seq2, $i-1, 1);
			 $i--; $j--;
		  }
		 elsif ($matrix[$i][$j]{pointer} eq "left") {
			 $align1 .= substr($seq1, $j-1, 1);
			 $align2 .= "-";
			 $j--;
		  }
		 elsif ($matrix[$i][$j]{pointer} eq "up") {
			 $align1 .= "-";
			 $align2 .= substr($seq2, $i-1, 1);
			 $i--;
		  }  
	 }

	 $align1 = reverse $align1;
	 $align2 = reverse $align2;
	 return ($align1,$align2);
}

sub reversecomplement {
  my $dna = $_[0];
  my $revcomp = reverse($dna);

  $revcomp =~ tr/ACGT/TGCA/;

  return $revcomp;
}