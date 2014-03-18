#!/usr/bin/env python

import time
from copy import copy

# Bar codes for identifying relevant lines
chr4bar1	= "GAGCGA"
chr4bar2	= "GATGTG"
CFTRbar1	= "CGGAAT"
CFTRbar2	= "CGTGGC"
MCFTRbar1	= "TGCGTA"
MCFTRbar2	= "TTCTGG"

# === Begin sequence constants ===

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

# === End sequence constants ===

# === Begin output variables ===

# Enum to hold experiment types.
class E:
	CFTR 	= 0
	MCFTR 	= 1
	CHR4 	= 2

statList = {'mutInterest': "",	# Portion to be replaced.
			'total': 0,			# Total num of sequences from specified experiment.
			'unmut': 0,			# Unmutated sequences.
			'mut': 0,			# Mutated sequences.
			'noMatch': 0,		# No match sequences.
			'mutOrig': 0,		# Num mutations in original portion of interest.
			'flankOrig': 0,		# Num mutations in original flanking regions.
			'mutMut': 0,		# Num mutations in MUTATED portion of interest.
			'flankMut': 0 		# Num mutations in MUTATED flanking regions.
			}

# Object holding analysis variables.
var = {E.CFTR: copy(statList), E.MCFTR: copy(statList), E.CHR4: copy(statList)}

var[E.CFTR]['mutInterest'] 	= mutInterestCFTR
var[E.MCFTR]['mutInterest'] = mutInterestMouse
var[E.CHR4]['mutInterest'] 	= interestChr4

# === End output variables ===

# === Utility functions ===

# This is a huge drain.
def strMatches(a, b):
	if a == b:
		return len(a)
	minLen = min(len(a), len(b))
	res 
	for i in range(0, minLen):
		if a[i] == b[i]:
			res += 1
	return res

# Return > match score of two sequences of same length
# If not same length, only process up to the length of the shorter sequence.
def matchScore(seq1, seq2):
	match = 1 		# match
	mismatch = -1 	# mismatch
	gap = -1 	# gap

	n = min(len(seq1), len(seq2))
	# Initialize dynamic programming array
	arr = [[0 for j in range(n + 1)] for i in range(n + 1)]

	maxScore = 0

	# Fill in array
	for i in range(1, n + 1):
		# Only consider a slice of width 6 to save time.
		for j in range(max(1, i - 3), min(n + 1, i + 3)):
			
			diag = arr[i - 1][j - 1]	# match/mismatch case

			if seq1[i - 1] == seq2[j - 1]:	# match
				diag += match
			else:					# mismatch
				diag += mismatch

			gap1 = arr[i - 1][j] + gap 	# deletion case
			gap2 = arr[i][j - 1] + gap 	# insertion case

			arr[i][j] = max(0, diag, gap1, gap2)

			maxScore = max(maxScore, arr[i][j])
	
	return maxScore

# Return => void
# Primary processing function
def process(line, exp, start, end, scoreThresh):
	var[exp]['total'] += 1
	interest = line[start:end]

	score = matchScore(interest, var[exp]['mutInterest'])
	
	if score > scoreThresh:
		print "cool region: " + interest
		print "actual     : " + var[exp]['mutInterest']
		print score
		print exp
		var[exp]['mut'] += 1
	else:
		var[exp]['unmut'] += 1


# === End Utility functions ===

infile = open("Sample_021714-1_005/021714-1_005_ACAGTG_L001_R1_001.fastq", "r")

startTime = time.time()

# Identify experiment by barcode.
chr4bar1	= "GAGCGA"
chr4bar2	= "GATGTG"

CFTRbar1	= "CGGAAT"
CFTRbar2	= "CGTGGC"

MCFTRbar1	= "TGCGTA"
MCFTRbar2	= "TTCTGG"

mice = 0
mutMice = 0
unmutMice = 0

for line in infile:
	bar = line[0:6]
	# Chr4
	# (strMatches(bar, chr4bar1) > 4 or
	# 	strMatches(bar, chr4bar2) > 4):
	if (bar == chr4bar1 or
		bar == chr4bar2):
		process(line, E.CHR4, 29 + 6, 76, 10)

	# CFTR
	if (bar == CFTRbar1 or
		bar == CFTRbar2):
		process(line, E.CFTR, 38 + 6, 62, 8)	# errors
	
	# Mouse
	if (bar == MCFTRbar1 or
		bar == MCFTRbar2):
		process(line, E.MCFTR, 40 + 6, 64, 13)

	
print "\n\n\n======================"

print "CFTR total: " + str(var[E.CFTR]['total'])
print "CFTR mutated: " + str(var[E.CFTR]['mut'])
print "CFTR no-match: " + str(var[E.CFTR]['noMatch'])
print "CFTR ROI mutations (original): " + str(var[E.CFTR]['mutOrig'])
print "CFTR flanking mutations (original): " + str(var[E.CFTR]['flankOrig'])
print "CFTR ROI mutations (mutated): " + str(var[E.CFTR]['mutMut'])
print "CFTR flanking mutations (mutated): " + str(var[E.CFTR]['flankMut'])

print "MCFTR total: " + str(var[E.MCFTR]['total'])
print "MCFTR mutated: " + str(var[E.MCFTR]['mut'])
print "MCFTR no-match: " + str(var[E.MCFTR]['noMatch'])
print "MCFTR ROI mutations (original): " + str(var[E.MCFTR]['mutOrig'])
print "MCFTR flanking mutations (original): " + str(var[E.MCFTR]['flankOrig'])
print "MCFTR ROI mutations (mutated): " + str(var[E.MCFTR]['mutMut'])
print "MCFTR flanking mutations (mutated): " + str(var[E.MCFTR]['flankMut'])

print "Chr4 total: " + str(var[E.CHR4]['total'])
print "CHR4 mutated: " + str(var[E.CHR4]['mut'])
print "Chr4 ROI mutations: " + str(var[E.CHR4]['mutOrig'])
print "Chr4 flanking mutations: " + str(var[E.CHR4]['flankOrig'])
print "Chr4 no-match: " + str(var[E.CHR4]['noMatch'])


print "Total time: "
print time.time() - startTime