#!/usr/bin/env python

import time

# Bar codes for identifying relevant lines
chr4bar1	= "GAGCGA"
chr4bar2	= "GATGTG"
CFTRbar1	= "CGGAAT"
CFTRbar2	= "CGTGGC"
MCFTRbar1	= "TGCGTA"
MCFTRbar2	= "TTCTGG"

# ==== Begin sequence constants ====

# CFTR (what does WT stand for??)
interestCFTR = "AGAAAATATCATTGG"      # Portion of interest (same as mouse)

origCFTR  = "TGTTCTCAGTTTTCCTGGATTATGCCTGGCACCATTAA"
origCFTR += interestCFTR
origCFTR += "TGTTTCCTATGATGAATATAGATACAGAAGCGTCATCAAAGCATGCCA"

# Desired product post-mutation
mutInterestCFTR = "GGAGAACATTATCTTTGG"      # Mutation portion of interest (same as mouse)

mutCFTR  = "TGTTCTCAGTTTTCCTGGATTATGCCTGGCACCATTAA"
mutCFTR += mutInterestCFTR
mutCFTR += "TGTTTCCTATGATGAATATAGATACAGAAGCGTCATCAAAGCATGCCA"

# Mouse
interestMouse = "AGAAAATATCATTGG"     # Portion of interest (same as CFTR)
# Take 11 char slice.
origMouse  = "TCTGCTCTCAATTTTCTTGGATTATGCCGGGTACTATCAA" 
origMouse += interestMouse 
origMouse += "TGTTTCCTATGATGAGTACAGATATAAGAGTGTTGTCAAAGCTTGCC"

# Desired product post-mutation.
mutInterestMouse = "GGAGAACATTATCTTTGG"    # Mutation portion of interest (same as CFTR)

mutMouse  = "TCTGCTCTCAATTTTCTTGGATTATGCCGGGTACTATCAA" 
mutMouse += mutInterestMouse 
mutMouse += "TGTTTCCTATGATGAGTACAGATATAAGAGTGTTGTCAAAGCTTGCC"

# Chr4
interestChr4 = "ATTAAAGAAAATATCCTATTTGGTCATTTATCATGAAGATA"
interestChr4Trunc = "TTTGGTCATTTATCATGAAGATA"

origChr4  = "TTCTACTAAAAGAAAACTTCTGTGTCCTA" 
origChr4 += interestChr4 
origChr4 += "ATGATAAATGTTAAATAACGTCTTGTTTGCATTAAGTCTGTGGGATG"

# ==== End sequence constants ====

# ==== Utility functions ====

# This is a huge drain.
def strMatches(a, b):
	if a == b:
		return len(a)
	minLen = min(len(a), len(b))
	res = 0
	for i in range(0, minLen):
		if a[i] == b[i]:
			res += 1
	return res

# Return => match score of two sequences of same length
# If not same length, only process up to the length of the shorter sequence.
def matchScore(seq1, seq2):
	match = 1 		# match
	mismatch = -1 	# mismatch
	gap = -1 	# gap

	n = min(len(seq1), len(seq2))
	# Initialize dynamic programming array
	arr = [[0 for j in range(n)] for i in range(n)]

	maxScore = 0

	# Fill in array
	for i in range(1, n):
		# Only consider a slice of width 6 to save time.
		for j in range(max(1, i - 3), min(n, i + 3)):
			
			diag = arr[i - 1][j - 1]	# match/mismatch case

			if seq1[i] == seq2[j]:	# match
				diag += match
			else:					# mismatch
				diag += mismatch

			gap1 = arr[i - 1][j] + gap 	# deletion case
			gap2 = arr[i][j - 1] + gap 	# insertion case

			arr[i][j] = max(0, diag, gap1, gap2)

			maxScore = max(maxScore, arr[i][j])
	
	return maxScore
	
# ==== End Utility functions ====

infile = open("Sample_021714-1_005/021714-1_005_ACAGTG_L001_R1_001.fastq", "r")

startTime = time.time()

# Identify experiment by barcode.
chr4bar1	= "GAGCGA"
chr4bar2	= "GATGTG"

CFTRbar1	= "CGGAAT"
CFTRbar2	= "CGTGGC"

MCFTRbar1	= "TGCGTA"
MCFTRbar2	= "TTCTGG"

mouse = 0
count = 0

for line in infile:
	bar = line[0:6]
	# Chr4
	# (strMatches(bar, chr4bar1) > 4 or
	# 	strMatches(bar, chr4bar2) > 4):
	if (bar == chr4bar1 or
		bar == chr4bar2):
		continue

	# CFTR
	if (bar == CFTRbar1 or
		bar == CFTRbar2):
		continue
	
	# Mouse
	if (bar == MCFTRbar1 > 4 or
		bar == MCFTRbar2 > 4):
		mouse += 1
		if line[38:58].find("GGAGAACATTAT") > 0:
			print count
			print line
			count += 1




	
print "\n\n\n======================"
print "Total mouse: " + str(mouse)
print "Mutated: " + str(count)
print "Total time:"
print time.time() - startTime